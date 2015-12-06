// double courant = 0.5*soundSpeed*tunnel.dt/tunnel.dx;  // Should this be sounds speed or u, v, or magnitude?? tunnel.dt/2
int i, j, n;
double pressAtm, mach, shockSpeed, rhoShock, rhoStart, soundSpeed;
grid tunnel, tunnel_halfTime;

mach = 1.1;
pressAtm = 1.0;
soundSpeed = 1.0;
rhoStart = pressAtm;
rhoShock = pressAtm/soundSpeed*(pow(mach,2.0)/(pow(0.88128485,2.0)*pow((1.0 - 1.0/(7.0*pow(mach,2.0))),2.5)) - 1.0);
shockSpeed = soundSpeed*mach;

tunnel.dt = dt;
tunnel.dx = 0.005;
tunnel.dy = 0.05;
tunnel.xPoints = 500;
tunnel.yPoints = 100; 
tunnel.rho.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0)); //[tunnel.yPoints][tunnel.xPoints];
tunnel.Mx.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
tunnel.My.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0));
// BCs initialize
for (i=0; i < 4; ++i) {
    if (i == 0) { // Left
        tunnel.connect[i] = false;
        tunnel.bounds[i] = "value";
        tunnel.rhoBCsLR.resize(tunnel.yPoints, vector<double>(4,0)); //[tunnel.yPoints][4];
        tunnel.MxBCsLR.resize(tunnel.yPoints, vector<double>(4,0)); //[tunnel.yPoints][4];
        tunnel.MyBCsLR.resize(tunnel.yPoints, vector<double>(4,0)); //[tunnel.yPoints][4];
        // For Roe Solver
        tunnel.rhoBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0)); //[tunnel.yPoints][2];
        tunnel.MxBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0)); //[tunnel.yPoints][2];
        tunnel.MyBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0)); //[tunnel.yPoints][2];
        for (n=0; n < tunnel.yPoints; ++n) {
            tunnel.rhoBCsLR[n][0] = rhoShock;
            tunnel.MxBCsLR[n][0] = rhoShock*soundSpeed; //*tunnel.dy*shockSpeed;                // Should this be set to (rhoShock*tunnel.dx*tunnel.dy*soundSpeed)??
            tunnel.MyBCsLR[n][0] = 0.0;
            // For Roe Solver
            tunnel.rhoBCsLR[n][2] = rhoShock;
            tunnel.MxBCsLR[n][2] = rhoShock*soundSpeed; //*tunnel.dy*shockSpeed;                // Should this be set to (rhoShock*tunnel.dx*tunnel.dy*soundSpeed)??
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
        tunnel.rhoBCsBT.resize(4, vector<double>(tunnel.xPoints,0)); //[4][tunnel.xPoints];
        tunnel.MxBCsBT.resize(4, vector<double>(tunnel.xPoints,0)); //[4][tunnel.xPoints];
        tunnel.MyBCsBT.resize(4, vector<double>(tunnel.xPoints,0)); //[4][tunnel.xPoints];
        // For Roe Solver
        tunnel.rhoBCsBT_PlusMinus.resize(2, vector<double>(tunnel.xPoints,0)); //[2][tunnel.xPoints];
        tunnel.MxBCsBT_PlusMinus.resize(2, vector<double>(tunnel.xPoints,0)); //[2][tunnel.xPoints];
        tunnel.MyBCsBT_PlusMinus.resize(2, vector<double>(tunnel.xPoints,0)); //[2][tunnel.xPoints];
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
        tunnel.rho[i][j] = rhoStart;
        tunnel.Mx[i][j] = 0.0;
        tunnel.My[i][j] = 0.0;
    }
}
// Initialize flux variables
tunnel.rhoF_Plus.resize(tunnel.yPoints, vector<double>(tunnel.xPoints,0)); //[tunnel.yPoints][tunnel.xPoints];
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