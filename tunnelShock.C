#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <time.h>
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
    string bounds [4];
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

    clock_t tStart;
    double wallTime;
    tStart = clock();

    #include "timeControls.H"
    #include "createGrid.H"

    cout << "Courant: " << shockSpeed*tunnel.dt/min(tunnel.dx,tunnel.dy) << "\n";
    // cin.get();

    currentTime += dt;
    while (currentTime < endTime) {

        cout << "t = " << currentTime << "\n";

        updateAMR(refineTunnel,tunnel,refineFactor,refineColumns,rhoShock,rhoStart);

        MUSCL(tunnel,soundSpeed);

        if (refineColumns.size() > 0) {
            #include "AMR.H"
        }

        // printGridFiles(tunnel,refineTunnel,refineColumns);
        // cin.get();

        currentTime += dt;
    }
    
    wallTime = (double)(clock() - tStart)/CLOCKS_PER_SEC;

    printGridFiles(tunnel,refineTunnel,refineColumns);

    cout << "ExecutionTime " << wallTime << "\n";
    ofstream output;
    output.precision(9);
    output.open ("output.txt");
    output << "Courant " << shockSpeed*tunnel.dt/min(tunnel.dx,tunnel.dy) 
        << "\nExecutionTime " << wallTime << "\nCells " << tunnel.xPoints*tunnel.yPoints
        << "\nNumberProcs " << 1;
    output.close();

    return 0;
}


