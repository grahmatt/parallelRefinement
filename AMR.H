// Active Mesh Refinement (AMR) only for x direction

int refineFactor = 4;
bool updateRefCol;
vector<int> refineColumns, oldRefCol;
grid refineTunnel, refineTunnel_halfTime;
refineTunnel.dt = tunnel.dt/refineFactor;
refineTunnel.dx = tunnel.dx/refineFactor;
refineTunnel.dy = tunnel.dy/refineFactor;
refineTunnel.xPoints = 0;
refineTunnel.yPoints = tunnel.yPoints*refineFactor;

i = 0;
n = 0;
updateRefCol = false;
for (j=0; j < tunnel.xPoints; ++j) {
    if (tunnel.rho[i][j] < 0.99*rhoShock && tunnel.rho[i][j] > 1.01*rhoStart) {
        if (!updateRefCol && oldRefCol[n] != j) {
            updateRefCol = true;
        }
        if (refinement.size() == 0 && j > 0) {
            refinement.push_back (j - 1);
        }
        refineColumns.push_back (j);
        if (refineColumns.size() > 1 && refineColumns.back() - 1 != refineColumns[refineColumns.end - 2]) {
            cout << "Error: Only one refinement zone available\n";
        }
        ++ n;
    }
}
if (refineColumns.back() != tunnel.xPoints - 1) {
    refineColumns.push_back (refineColumns.back() + 1);
}

if updateRefCol {
    refineTunnel_halfTime = refineTunnel;

    refineTunnel.xPoints = refineColumns.size()*refineFactor;
    refineTunnel.rho.resize(refineTunnel.yPoints, vector<double>(refineTunnel.xPoints,0));
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

    for (j=0; j < refineColumns.size(); ++j) {
        for (n=0; n < refineFactor; ++n) {
            for (i=0; i < refineTunnel.yPoints, ++i) {
                kdjflkdajf;
            }
        }
    }
}
