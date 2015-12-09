int i, j, n, m, k;
bool updateRefCol, existingCol;
vector<int> refineColumns;
vector<double> rhoFixFPlus, rhoFixFMinus, MxFixFPlus, MxFixFMinus, 
    MyFixFPlus, MyFixFMinus;
grid tunnel, refineTunnel;

// Physical values
const double mach = 1.1;
const double pressAtm = 1.0;
const double soundSpeed = 1.0;
const double rhoStart = pressAtm;
const double rhoShock = pressAtm/soundSpeed*(pow(mach,2.0)/(pow(0.88128485,2.0)*
    pow((1.0 - 1.0/(7.0*pow(mach,2.0))),2.5)) - 1.0);
const double shockSpeed = soundSpeed*mach;

// Create coarse grid "tunnel"
tunnel.dt = dt;
tunnel.dx = 0.05;
tunnel.dy = 0.05;
tunnel.xPoints = 500;
tunnel.yPoints = 100; 
// Grid initialize
tunnel.rho.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,rhoStart)); //[tunnel.yPoints][tunnel.xPoints];
tunnel.Mx.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.My.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
// BCs initialize
for (i=0; i < 4; ++i) {
    if (i == 0) { // Left
        tunnel.connect[i] = false;
        tunnel.bounds[i] = "value";
        tunnel.rhoBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
        tunnel.MxBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
        tunnel.MyBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
        // For Roe Solver
        tunnel.rhoBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));
        tunnel.MxBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));
        tunnel.MyBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));
        for (n=0; n < tunnel.yPoints; ++n) {
            tunnel.rhoBCsLR[n][0] = rhoShock;
            tunnel.MxBCsLR[n][0] = 0.0; //rhoShock*soundSpeed;
            tunnel.MyBCsLR[n][0] = 0.0;
            // For Roe Solver
            tunnel.rhoBCsLR[n][2] = rhoShock;
            tunnel.MxBCsLR[n][2] = 0.0; //rhoShock*shockSpeed;
            tunnel.MyBCsLR[n][2] = 0.0;
        }
    }
    else if (i == 1) { // Right
        tunnel.connect[i] = false;
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
        tunnel.rhoBCsBT.resize(4, vector<double>(tunnel.xPoints,0));
        tunnel.MxBCsBT.resize(4, vector<double>(tunnel.xPoints,0));
        tunnel.MyBCsBT.resize(4, vector<double>(tunnel.xPoints,0));
        // For Roe Solver
        tunnel.rhoBCsBT_PlusMinus.resize(2, vector<double>(tunnel.xPoints,0));
        tunnel.MxBCsBT_PlusMinus.resize(2, vector<double>(tunnel.xPoints,0));
        tunnel.MyBCsBT_PlusMinus.resize(2, vector<double>(tunnel.xPoints,0));
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
// Initialize flux variables
tunnel.rhoF_Plus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.rhoF_Minus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.rhoG_Plus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.rhoG_Minus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.MxF_Plus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.MxF_Minus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.MxG_Plus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.MxG_Minus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.MyF_Plus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.MyF_Minus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.MyG_Plus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.MyG_Minus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));

// Create fine grid "refineTunnel"
const int refineFactor = 4;
refineColumns.push_back(0);
refineTunnel.dt = tunnel.dt/refineFactor;
refineTunnel.dx = tunnel.dx/refineFactor;
refineTunnel.dy = tunnel.dy/refineFactor;
// Grid initialize
refineTunnel.xPoints = 0;
refineTunnel.yPoints = tunnel.yPoints*refineFactor;
refineTunnel.connect[0] = false;
refineTunnel.connect[1] = true;
refineTunnel.connect[2] = tunnel.connect[2];
refineTunnel.connect[3] = tunnel.connect[3];
refineTunnel.rhoBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
refineTunnel.MxBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
refineTunnel.MyBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
refineTunnel.rhoBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));
refineTunnel.MxBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));
refineTunnel.MyBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));
// Initialize first column for first time step
refineTunnel.xPoints = refineColumns.size()*refineFactor;
refineTunnel.rho.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,rhoStart));
refineTunnel.Mx.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.My.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.rhoBCsBT.resize(4, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MxBCsBT.resize(4, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MyBCsBT.resize(4, vector<double>(refineTunnel.xPoints,0));
refineTunnel.rhoBCsBT_PlusMinus.resize(2, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MxBCsBT_PlusMinus.resize(2, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MyBCsBT_PlusMinus.resize(2, vector<double>(refineTunnel.xPoints,0));
refineTunnel.rhoF_Plus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.rhoF_Minus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.rhoG_Plus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.rhoG_Minus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MxF_Plus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MxF_Minus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MxG_Plus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MxG_Minus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MyF_Plus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MyF_Minus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MyG_Plus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
refineTunnel.MyG_Minus.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
cout << "HERE\n" << refineColumns.size() << " " << refineColumns.back() << "\n";
updateBCs(refineTunnel,tunnel,refineColumns,refineFactor);
cout << "HERE\n";

// Initialize flux fixes
rhoFixFPlus.resize(tunnel.yPoints,0);
rhoFixFMinus.resize(tunnel.yPoints,0);
MxFixFPlus.resize(tunnel.yPoints,0);
MxFixFMinus.resize(tunnel.yPoints,0);
MyFixFPlus.resize(tunnel.yPoints,0);
MyFixFMinus.resize(tunnel.yPoints,0);
