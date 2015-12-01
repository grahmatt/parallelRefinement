#include <iostream>
#include <cmath>
#include <vector>
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
    vector<vector<double> > rhoBCsTB;
    vector<vector<double> > MxBCsTB;
    vector<vector<double> > MyBCsTB;
    double dx;
    double dy;
    bool connect(4);
    char bounds(4);
    int xPoints;
    int yPoints;
};

struct flux {
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
};

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
    tunnel.rho[tunnelDim[0]][tunnelDim[1]];
    tunnel.Mx[tunnelDim[0]][tunnelDim[1]];
    tunnel.My[tunnelDim[0]][tunnelDim[1]];
    tunnel.dx = 0.1;
    tunnel.dy = 0.1;
    tunnel.xPoints = 50;
    tunnel.yPoints = 10; 
    // BCs initialize
    for (i=0; i < 4; ++i) {
        if (i == 0) { // Left
            tunnel.connect[i] = false;
            tunnel.bounds[i] = "value";
            tunnel.rhoBCsLR[tunnel.yPoints][2];
            tunnel.MxBCsLR[tunnel.yPoints][2];
            tunnel.MyBCsLR[tunnel.yPoints][2];
            for (n=0; n < tunnel.yPoints, ++n) {
                tunnel.rhoBCsLR[n][0] = rhoShock;
                tunnel.MxBCsLR[n][0] = 0.0;
                tunnel.MyBCsLR[n][0] = 0.0;
            }
        }
        else if (i == 1) { // Right
            tunnel.connect[i] = true;
            tunnel.bounds[i] = "outside";
            for (n=0; n < tunnel.yPoints, ++n) {
                tunnel.rhoBCsLR[n][1] = rhoStart;
                tunnel.MxBCsLR[n][1] = 0.0;
                tunnel.MyBCsLR[n][1] = 0.0;
            }
        }
        else if (i == 2) { // Bottom
            tunnel.connect[i] = false;
            tunnel.bounds[i] = "wall";
            tunnel.rhoBCsTB[2][tunnel.xPoints];
            tunnel.MxBCsTB[2][tunnel.xPoints];
            tunnel.MyBCsTB[2][tunnel.xPoints];
            for (n=0; n < tunnel.xPoints, ++n) {
                tunnel.rhoBCsTB[0][n] = rhoStart;
                tunnel.MxBCsTB[0][n] = 0.0;
                tunnel.MyBCsTB[0][n] = 0.0;
            }
        }
        else { // Top
            tunnel.connect[i] = false;
            tunnel.bounds[i] = "wall";
            for (n=0; n < tunnel.xPoints, ++n) {
                tunnel.rhoBCsTB[1][n] = rhoStart;
                tunnel.MxBCsTB[1][n] = 0.0;
                tunnel.MyBCsTB[1][n] = 0.0;
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
    flux tunnelFlux;
    tunnelFlux.rhoF_Plus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.rhoF_Minus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.rhoG_Plus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.rhoG_Minus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.MxF_Plus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.MxF_Minus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.MxG_Plus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.MxG_Minus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.MyF_Plus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.MyF_Minus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.MyG_Plus[tunnelDim[0]][tunnelDim[1]];
    tunnelFlux.MyG_Minus[tunnelDim[0]][tunnelDim[1]];

    // MUSCL Half Time
    for (i=0; i < tunnel.yPoints; ++i) {
        for (j=0; j < tunnel.xPoints; ++j) {
            xS_Plus = 
        }
    }


}








