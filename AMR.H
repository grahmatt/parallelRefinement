// Active Mesh Refinement (AMR) only for x direction

int m, k, refineFactor = 4;
bool updateRefCol, existingCol;
vector<int> refineColumns, oldRefCol;
grid refineTunnel, refineTunnel_halfTime;
refineTunnel.dt = tunnel.dt/refineFactor;
refineTunnel.dx = tunnel.dx/refineFactor;
refineTunnel.dy = tunnel.dy/refineFactor;
refineTunnel.xPoints = 0;
refineTunnel.yPoints = tunnel.yPoints*refineFactor;
refineTunnel.bounds[2] = tunnel.bounds[2];
refineTunnel.bounds[3] = tunnel.bounds[3];
refineTunnel.rhoBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
refineTunnel.MxBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
refineTunnel.MyBCsLR.resize(tunnel.yPoints, vector<double>(4,0));
refineTunnel.rhoBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));
refineTunnel.MxBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));
refineTunnel.MyBCsLR_PlusMinus.resize(tunnel.yPoints, vector<double>(2,0));

vector<double> rhoFixFPlus, rhoFixFMinus, MxFixFPlus, MxFixFMinus, MyFixFPlus, MyFixFMinus;
rhoFixFPlus.resize(tunnel.yPoints,0);
rhoFixFMinus.resize(tunnel.yPoints,0);
MxFixFPlus.resize(tunnel.yPoints,0);
MxFixFMinus.resize(tunnel.yPoints,0);
MyFixFPlus.resize(tunnel.yPoints,0);
MyFixFMinus.resize(tunnel.yPoints,0);


i = 0;
n = -1;
updateRefCol = false;
oldRefCol = refineColumns;
refineColumns.resize(0);
for (j=0; j < tunnel.xPoints; ++j) {
    if (tunnel.rho[i][j] < 0.99*rhoShock && tunnel.rho[i][j] > 1.01*rhoStart) {
        if (!updateRefCol) {
            if (n + 1 == oldRefCol.size()) {
                updateRefCol = true;
            }
            else {
                ++n;
                if (oldRefCol[n] != j) {
                    updateRefCol = true;
                }
            }
        }
        if (refinement.size() == 0 && j > 0) {
            refinement.push_back (j - 1);
            if (j - 1 == 0) {
                refineTunnel.bounds[0] = tunnel.bounds[0];
            }
            else {
                refineTunnel.bounds[0] = true;
            }
        }
        refineColumns.push_back (j);
        if (refineColumns.size() > 1 && refineColumns.back() - 1 != refineColumns[refineColumns.end - 2]) {
            cout << "Error: Only one refinement zone available\n";
            return 0;
        }
    }
}
if (refineColumns.back() != tunnel.xPoints - 1) {
    refineColumns.push_back (refineColumns.back() + 1);
    if (refineColumns.back() == tunnel.xPoints - 1) {
        refineTunnel.bounds[1] = tunnel.bounds[1];
    }
    else {
        refineTunnel.bounds[1] = true;
    }
}
else {
    refineTunnel.bounds[1] = true;
}
if (!updateRefCol && refineColumns.size() != oldRefCol.size()) {
    updateRefCol = true;
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
        existingCol = false;
        for (m = 0; m < oldRefCol.size(); ++m) {
            if (refineColumns[j] == oldRefCol[m]) {
                existingCol = true;
                break;
            }
        }
        if existingCol {
            for (n=0; n < refineFactor; ++n) {
                for (i=0; i < refineTunnel.yPoints, ++i) {
                    refineTunnel.rho[i][j*refineFactor + n] = refineTunnel_halfTime.rho[i][m*refineFactor + n];
                    refineTunnel.Mx[i][j*refineFactor + n] = refineTunnel_halfTime.Mx[i][m*refineFactor + n];
                    refineTunnel.My[i][j*refineFactor + n] = refineTunnel_halfTime.My[i][m*refineFactor + n];
                }
            }
        }
        else {
            for (n=0; n < refineFactor; ++n) {
                for (i=0; i < tunnel.yPoints, ++i) {
                    for (k=0; k < refineFactor; ++k) {
                        refineTunnel.rho[i*refineFactor + k][j*refineFactor + n] = tunnel.rho[i][refineColumns[j]];
                        refineTunnel.Mx[i*refineFactor + k][j*refineFactor + n] = tunnel.Mx[i][refineColumns[j]];
                        refineTunnel.My[i*refineFactor + k][j*refineFactor + n] = tunnel.My[i][refineColumns[j]];
                    }
                }
            }
        }
    }

    updateBCs(refineTunnel,tunnel,refineColumns,refineFactor);


    
    // Start here
    for (i=0; i < tunnel.yPoints; ++i) {
        if (refineColumns[0] != 0) {
            rhoFixFPlus[i] = tunnel.rhoF_Plus[i][refineColumns[0] - 1];
            MxFixFPlus[i] = tunnel.MxF_Plus[i][refineColumns[0] - 1];
            MyFixFPlus[i] = tunnel.MyF_Plus[i][refineColumns[0] - 1];
        }
        if (refineColumns.back() != tunnel.xPoints - 1) {
            rhoFixFMinus[i] = - tunnel.rhoF_Minus[i][refineColumns.back() + 1];
            MxFixFMinus[i] = - tunnel.MxF_Minus[i][refineColumns.back() + 1];
            MyFixFMinus[i] = - tunnel.MyF_Minus[i][refineColumns.back() + 1];
        }
    }
    for (n=0; n < refineFactor; ++n) {
        MUSCL(refineTunnel,soundSpeed,tunnel,refineColumns,refineFactor);
        for (i=0; i < tunnel.yPoints, ++i) {
            for (k=0; k < refineFactor; ++k) {
                if (refineColumns[0] != 0) {
                    rhoFixFPlus[i] -= refineTunnel.rhoF_Minus[i*refineFactor + k][0]/pow(refineFactor,2);
                    MxFixFPlus[i] -= refineTunnel.MxF_Minus[i*refineFactor + k][0]/pow(refineFactor,2);
                    MyFixFPlus[i] -= refineTunnel.MyF_Minus[[i*refineFactor + k][0]/pow(refineFactor,2);
                }
                if (refineColumns.back() != tunnel.xPoints - 1) {
                    rhoFixFMinus[i] += refineTunnel.rhoF_Plus[i*refineFactor + k][refineTunnel.xPoints - 1]/pow(refineFactor,2);
                    MxFixFMinus[i] += refineTunnel.MxF_Plus[i*refineFactor + k][refineTunnel.xPoints - 1]/pow(refineFactor,2);
                    MyFixFMinus[i] += refineTunnel.MyF_Plus[i*refineFactor + k][refineTunnel.xPoints - 1]/pow(refineFactor,2);
                }
            }
        }
    }
    // Coarse border cell flux updates
    for (i=0; i < tunnel.yPoints; ++i) {
        if (refineColumns[0] != 0) {
            tunnel.rho[i][refineColumns[0] - 1] += tunnel.dt/tunnel.dx*rhoFixFPlus[i];
            tunnel.Mx[i][refineColumns[0] - 1] += tunnel.dt/tunnel.dx*MxFixFPlus[i];
            tunnel.My[i][refineColumns[0] - 1] += tunnel.dt/tunnel.dx*MyFixFPlus[i];
        }
        if (refineColumns.back() != tunnel.xPoints - 1) {
            tunnel.rho[i][refineColumns.back() + 1] += tunnel.dt/tunnel.dx*rhoFixFMinus[i];
            tunnel.Mx[i][refineColumns.back() + 1] += tunnel.dt/tunnel.dx*MxFixFMinus[i];
            tunnel.My[i][refineColumns.back() + 1] += tunnel.dt/tunnel.dx*MyFixFMinus[i];
        }
    }
    // Coarse inner cell center updates
    
}
