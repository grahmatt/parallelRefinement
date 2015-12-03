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
    double dt;
    double dx;
    double dy;
    bool connect [4];
    string bounds [4];
    int xPoints;
    int yPoints;
};

// struct flux {
//     vector<vector<double> > rhoF_Plus;
//     vector<vector<double> > rhoF_Minus;
//     vector<vector<double> > rhoG_Plus;
//     vector<vector<double> > rhoG_Minus;
//     vector<vector<double> > MxF_Plus;
//     vector<vector<double> > MxF_Minus;
//     vector<vector<double> > MxG_Plus;
//     vector<vector<double> > MxG_Minus;
//     vector<vector<double> > MyF_Plus;
//     vector<vector<double> > MyF_Minus;
//     vector<vector<double> > MyG_Plus;
//     vector<vector<double> > MyG_Minus;
// };

#include "fluxLimiter.H"
#include "reconstructFlux.H"
#include "roeSolver.H"

int main(int argc, char** argv){

    // if (argc < 2){
    //     cout << "Missing matrix size: \"mpirun -np <#> hw3.bin <matrix size>\"" << endl;
    //     return 0;
    // }
    // else {
    //     n = atoi(argv[1]);
    // }

    int i, j, n;
    double rhoShock = 2.0, rhoStart = 1.0;

    grid tunnel;
    tunnel.dt = 0.01;
    tunnel.dx = 0.1;
    tunnel.dy = 0.1;
    tunnel.xPoints = 50;
    tunnel.yPoints = 10; 
    tunnel.rho[tunnel.xPoints][tunnel.yPoints];
    tunnel.Mx[tunnel.xPoints][tunnel.yPoints];
    tunnel.My[tunnel.xPoints][tunnel.yPoints];
    // BCs initialize
    for (i=0; i < 4; ++i) {
        if (i == 0) { // Left
            tunnel.connect[i] = false;
            tunnel.bounds[i] = "value";
            tunnel.rhoBCsLR[tunnel.yPoints][4];
            tunnel.MxBCsLR[tunnel.yPoints][4];
            tunnel.MyBCsLR[tunnel.yPoints][4];
            // For Roe Solver
            tunnel.rhoBCsLR_PlusMinus[tunnel.yPoints][2];
            tunnel.MxBCsLR_PlusMinus[tunnel.yPoints][2];
            tunnel.MyBCsLR_PlusMinus[tunnel.yPoints][2];
            for (n=0; n < tunnel.yPoints; ++n) {
                tunnel.rhoBCsLR[n][0] = rhoShock;
                tunnel.MxBCsLR[n][0] = 0.0;                // Should this be set to (rhoShock*tunnel.dx*tunnel.dy*soundSpeed)??
                tunnel.MyBCsLR[n][0] = 0.0;
                // For Roe Solver
                tunnel.rhoBCsLR[n][2] = rhoShock;
                tunnel.MxBCsLR[n][2] = 0.0;                // Should this be set to (rhoShock*tunnel.dx*tunnel.dy*soundSpeed)??
                tunnel.MyBCsLR[n][2] = 0.0;
            }
        }
        else if (i == 1) { // Right
            tunnel.connect[i] = true;
            tunnel.bounds[i] = "zeroGradient";
            for (n=0; n < tunnel.yPoints; ++n) {
                tunnel.rhoBCsLR[n][1] = rhoStart;
                tunnel.MxBCsLR[n][1] = 0.0;
                tunnel.MyBCsLR[n][1] = 0.0;
                // For Roe Solver
                tunnel.rhoBCsLR[n][3] = rhoStart;
                tunnel.MxBCsLR[n][3] = 0.0;
                tunnel.MyBCsLR[n][3] = 0.0;
            }
        }
        else if (i == 2) { // Bottom
            tunnel.connect[i] = false;
            tunnel.bounds[i] = "wall";
            tunnel.rhoBCsBT[4][tunnel.xPoints];
            tunnel.MxBCsBT[4][tunnel.xPoints];
            tunnel.MyBCsBT[4][tunnel.xPoints];
            // For Roe Solver
            tunnel.rhoBCsBT_PlusMinus[2][tunnel.xPoints];
            tunnel.MxBCsBT_PlusMinus[2][tunnel.xPoints];
            tunnel.MyBCsBT_PlusMinus[2][tunnel.xPoints];
            for (n=0; n < tunnel.xPoints; ++n) {
                tunnel.rhoBCsBT[0][n] = rhoStart;
                tunnel.MxBCsBT[0][n] = 0.0;
                tunnel.MyBCsBT[0][n] = 0.0;
                // For Roe Solver
                tunnel.rhoBCsBT[2][n] = rhoStart;
                tunnel.MxBCsBT[2][n] = 0.0;
                tunnel.MyBCsBT[2][n] = 0.0;
            }
        }
        else { // Top
            tunnel.connect[i] = false;
            tunnel.bounds[i] = "wall";
            for (n=0; n < tunnel.xPoints; ++n) {
                tunnel.rhoBCsBT[1][n] = rhoStart;
                tunnel.MxBCsBT[1][n] = 0.0;
                tunnel.MyBCsBT[1][n] = 0.0;
                // For Roe Solver
                tunnel.rhoBCsBT[3][n] = rhoStart;
                tunnel.MxBCsBT[3][n] = 0.0;
                tunnel.MyBCsBT[3][n] = 0.0;
            }
        }
    }
    // Grid initialize
    for (i=0; i < tunnel.yPoints; ++i) {
        for (j=0; j < tunnel.xPoints; ++j) {
            tunnel.rho[i][j] = 1.0;
            tunnel.Mx[i][j] = 0.0;
            tunnel.My[i][j] = 0.0;
        }
    }

    // Initialize flux variables
    // flux tunnelFlux;
    tunnel.rhoF_Plus[tunnel.xPoints][tunnel.yPoints];
    tunnel.rhoF_Minus[tunnel.xPoints][tunnel.yPoints];
    tunnel.rhoG_Plus[tunnel.xPoints][tunnel.yPoints];
    tunnel.rhoG_Minus[tunnel.xPoints][tunnel.yPoints];
    tunnel.MxF_Plus[tunnel.xPoints][tunnel.yPoints];
    tunnel.MxF_Minus[tunnel.xPoints][tunnel.yPoints];
    tunnel.MxG_Plus[tunnel.xPoints][tunnel.yPoints];
    tunnel.MxG_Minus[tunnel.xPoints][tunnel.yPoints];
    tunnel.MyF_Plus[tunnel.xPoints][tunnel.yPoints];
    tunnel.MyF_Minus[tunnel.xPoints][tunnel.yPoints];
    tunnel.MyG_Plus[tunnel.xPoints][tunnel.yPoints];
    tunnel.MyG_Minus[tunnel.xPoints][tunnel.yPoints];


    double soundSpeed = 1.0;
    double courant = 0.5*soundSpeed*tunnel.dt/tunnel.dx;  // Should this be sounds speed or u, v, or magnitude?? tunnel.dt/2

    #include "MUSCL.H"

}








