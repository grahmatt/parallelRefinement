#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
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
    // if (argc < 2){
    //     cout << "Missing matrix size: \"mpirun -np <#> hw3.bin <matrix size>\"" << endl;
    //     return 0;
    // }
    // else {
    //     n = atoi(argv[1]);
    // }

    #include "timeControls.H"
    #include "createGrid.H"

    cout << "Courant: " << shockSpeed*tunnel.dt/min(tunnel.dx,tunnel.dy) << "\n";
    cin.get();

    while (currentTime < endTime) {

        updateAMR(refineTunnel,tunnel,refineFactor,refineColumns,rhoShock,rhoStart);

        currentTime += dt;
        cout << "t = " << currentTime << "\n";

        MUSCL(tunnel,soundSpeed);

        cout << "Coarse Ran\n";

        if (refineColumns.size() > 0) {
            #include "AMR.H"
            cout << "Completed AMR\n";
        }

        // printGridFiles(tunnel);
        printGridFiles(tunnel,refineTunnel,refineColumns);
        cin.get();
    }
    
    // printGridFiles(tunnel);
    printGridFiles(tunnel,refineTunnel,refineColumns);

    return 0;
}


