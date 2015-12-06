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

void printGrid(grid& cells){
    cout.precision(4);
    for(int i=cells.yPoints-1; i >= 0; i--){
        for(int j=0; j < cells.xPoints; j++){
            cout << cells.rho[i][j] << "  ";
        }
    cout << "\n";
    }
}

#include "fluxLimiter.H"
#include "reconstructFlux.H"
#include "roeSolver.H"
#include "updateBCs.H"
// #include "gnuplot_i.hpp" //Gnuplot class handles POSIX-Pipe-communikation with Gnuplot

int main(int argc, char** argv){

    // if (argc < 2){
    //     cout << "Missing matrix size: \"mpirun -np <#> hw3.bin <matrix size>\"" << endl;
    //     return 0;
    // }
    // else {
    //     n = atoi(argv[1]);
    // }

    
    // Gnuplot g1("lines");
    // g1.savetops("test_output");
    // cout << "*** plotting slopes" << endl;
    // g1.set_title("Slopes\\nNew Line");
    // cout << "y = x" << endl;
    // g1.plot_slope(1.0,0.0,"y=x");
    // cout << "y = 2*x" << endl;
    // g1.plot_slope(2.0,0.0,"y=2x");
    // cout << "y = -x" << endl;
    // g1.plot_slope(-1.0,0.0,"y=-x");
    // g1.unset_title();

    #include "timeControls.H"
    #include "createGrid.H"

    cout << "Courant: " << shockSpeed*tunnel.dt/min(tunnel.dx,tunnel.dy) <<
        "\n Enter to continue...\n";
    cin.get();

    while (currentTime <= endTime) {

        // cout << "Printing Tunnel\n";
        // printGrid(tunnel);

        currentTime += dt;
        cout << "t = " << currentTime << "\n";

        #include "MUSCL.H"

        // cin.get();
    }
    
    ofstream rhoFile;
    ofstream xFile;
    ofstream yFile;
    rhoFile.precision(9);
    xFile.precision(9);
    yFile.precision(9);
    rhoFile.open ("RHO.txt");
    xFile.open ("X.txt");
    yFile.open ("Y.txt");
    xFile << "NumberOfColumns: " << tunnel.xPoints;
    xFile << " NumberOfRows: " << tunnel.yPoints << "\n";
    for(int i=0; i < tunnel.yPoints; i++){
        for(int j=0; j < tunnel.xPoints; j++){
            rhoFile << tunnel.rho[i][j] << " ";
            xFile << tunnel.dx*j + 0.5*tunnel.dx << " ";
            yFile << tunnel.dy*i + 0.5*tunnel.dy << " ";
        }
        if (i + 1 < tunnel.yPoints) {
            rhoFile << "\n";
            xFile << "\n";
            yFile << "\n";
        }
    }
    rhoFile.close();
    xFile.close();
    yFile.close();

    // cout << "Printing Tunnel\n";
    // printGrid(tunnel);

    // cout.precision(10);
    // cout << "Shock density: " << rhoShock << "\nDensity: " << tunnel.rho[0][0] << "\n";

}








