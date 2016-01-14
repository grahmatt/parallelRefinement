#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <mpi.h>
using namespace std;

struct grid {
    vector<vector<double> > rho;
    vector<vector<double> > Mx;
    vector<vector<double> > My;
    vector<vector<double> > rhoBCsLR;
    vector<vector<double> > MxBCsLR;
    vector<vector<double> > MyBCsLR;
    vector<vector<double> > rhoBCsBT;
    vector<vector<double> > MxBCsBT;
    vector<vector<double> > MyBCsBT;
    vector<vector<double> > rhoBCsLR_PlusMinus;
    vector<vector<double> > MxBCsLR_PlusMinus;
    vector<vector<double> > MyBCsLR_PlusMinus;
    vector<vector<double> > rhoBCsBT_PlusMinus;
    vector<vector<double> > MxBCsBT_PlusMinus;
    vector<vector<double> > MyBCsBT_PlusMinus;
    vector<vector<double> > rhoF_Plus;
    vector<vector<double> > rhoF_Minus;
    vector<vector<double> > rhoG_Plus;
    vector<vector<double> > rhoG_Minus;
    vector<vector<double> > MxF_Plus;
    vector<vector<double> > MxF_Minus;
    vector<vector<double> > MxG_Plus;
    vector<vector<double> > MxG_Minus;
    vector<vector<double> > MyF_Plus;
    vector<vector<double> > MyF_Minus;
    vector<vector<double> > MyG_Plus;
    vector<vector<double> > MyG_Minus;
    vector<double> oldRhoLeft;
    vector<double> oldMxLeft;
    vector<double> oldMyLeft;
    vector<double> oldRhoRight;
    vector<double> oldMxRight;
    vector<double> oldMyRight;
    double dt;
    double dx;
    double dy;
    bool connect [4];
    int xPoints;
    int yPoints;
};

#include "printData.H"
#include "fluxLimiter.H"
#include "reconstructFlux.H"
#include "roeSolver.H"
#include "updateBCs.H"
#include "updateAMR.H"
#include "MUSCL.H"

int main(int argc, char** argv) {

    int rank, size;
    double wallTime;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        wallTime = MPI_Wtime();
    }

    #include "timeControls.H"
    #include "createGrid.H"

    if (rank == 0) {
        cout << "Courant: " << shockSpeed*tunnel.dt/min(tunnel.dx,tunnel.dy) << "\n";
        // cin.get();
    }

    currentTime += dt;
    while (currentTime < endTime) {
        
        if (rank == 0) {
            cout << "t = " << currentTime << "\n";
        }

        updateAMR(refineTunnel, tunnel, refineFactor, refineColumns ,rhoShock, rhoStart, 
            rank, size, currentTime, sendLeft, recvLeft, sendRight, recvRight, vectorSize, vectorSizeSend);

        MUSCL(tunnel, soundSpeed, rank, sendLeft, recvLeft, sendRight, recvRight);

        if (refineColumns.size() > 0) {
            #include "AMR.H"
        }

        // printGridFiles(tunnel, refineTunnel, refineColumns, rank , size);
        // if (rank == 0) {cin.get();}

        currentTime += dt;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        wallTime = MPI_Wtime() - wallTime;
        cout << "ExecutionTime " << wallTime << "\n";
        ofstream output;
        output.precision(9);
        output.open ("output.txt");
        output << "Courant " << shockSpeed*tunnel.dt/min(tunnel.dx,tunnel.dy) 
            << "\nExecutionTime " << wallTime << "\nCells " << xPointsGlobal*yPointsGlobal
            << "\nNumberProcs " << size;
        output.close();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    printGridFiles(tunnel, refineTunnel, refineColumns, rank, size);

    MPI_Finalize();
    return 0;
}


