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

    while (currentTime <= endTime) {

        // cout << "Printing Tunnel\n";
        // printGrid(tunnel);

        currentTime += dt;
        cout << "t = " << currentTime << "\n";

        #include "MUSCL.H"

        // cin.get();
    }
    
    // ofstream myfile;
    // myfile.precision(15);
    // myfile.open ("output.txt");
    // myfile << "n = " << n << endl;
    // myfile << "Abs Sum = " << fixed << globalSum << endl;
    // myfile << "Minimum = " << fixed << globalMin << endl;
    // myfile << "wallTime = " << fixed << wallTime << endl;
    // myfile.close();
    // cout.precision(3);
    // cout << "n = " << n << endl;
    // cout << "globalMin = " << fixed << globalMin << endl;
    // cout << "globalAbsSum = " << fixed << globalSum << endl;
    // cout << "wallTime = " << fixed << wallTime << endl;

    cout << "Printing Tunnel\n";
    printGrid(tunnel);

    cout.precision(10);
    cout << "Shock density: " << rhoShock << "\nDensity: " << tunnel.rho[0][0] << "\n";

}








