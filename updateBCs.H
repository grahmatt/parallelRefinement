// Coarse cell boundary condition updates
void updateBCs
(
    grid& cells
) 
{
    int i, n;
    for (i=0; i < 4; ++i) {
        if (!cells.connect[i]) {
            if (i == 0) { // Left
                for (n=0; n < cells.yPoints; ++n) {
                    cells.MxBCsLR[n][0] = cells.Mx[n][0];
                    // For Roe Solver
                    cells.MxBCsLR[n][2] = cells.Mx[n][0];
                }
            }
            else if (i == 1) { // Right
                for (n=0; n < cells.yPoints; ++n) {
                    cells.rhoBCsLR[n][1] = cells.rho[n][cells.xPoints-1];
                    cells.MxBCsLR[n][1] = cells.Mx[n][cells.xPoints-1];
                    cells.MyBCsLR[n][1] = cells.My[n][cells.xPoints-1];
                    // For Roe Solver
                    cells.rhoBCsLR[n][3] = cells.rho[n][cells.xPoints-2];
                    cells.MxBCsLR[n][3] = cells.Mx[n][cells.xPoints-2];
                    cells.MyBCsLR[n][3] = cells.My[n][cells.xPoints-2];
                }
            }
            else if (i == 2) { // Bottom
                for (n=0; n < cells.xPoints; ++n) {
                    cells.rhoBCsBT[0][n] = cells.rho[0][n];
                    cells.MxBCsBT[0][n] = cells.Mx[0][n];
                    cells.MyBCsBT[0][n] = -cells.My[0][n];
                    // For Roe Solver
                    cells.rhoBCsBT[2][n] = cells.rho[1][n];
                    cells.MxBCsBT[2][n] = cells.Mx[1][n];
                    cells.MyBCsBT[2][n] = -cells.My[1][n];
                }
            }
            else { // Top
                for (n=0; n < cells.xPoints; ++n) {
                    cells.rhoBCsBT[1][n] = cells.rho[cells.yPoints-1][n];
                    cells.MxBCsBT[1][n] = cells.Mx[cells.yPoints-1][n];
                    cells.MyBCsBT[1][n] = -cells.My[cells.yPoints-1][n];
                    // For Roe Solver
                    cells.rhoBCsBT[3][n] = cells.rho[cells.yPoints-2][n];
                    cells.MxBCsBT[3][n] = cells.Mx[cells.yPoints-2][n];
                    cells.MyBCsBT[3][n] = -cells.My[cells.yPoints-2][n];
                }
            }
        }
    }
}
// Refine cell boundary values storing coarse cell old time values
void updateBCs
(
    grid& cells,
    grid& coarseCells,
    vector<int> refineColumns,
    const int refineFactor
) 
{
    int i, j, m, n, k;
    for (m=0; m < 4; ++m) {
        if (!cells.connect[m]) {
            if (m == 0) { // Left
                for (i=0; i < coarseCells.yPoints; ++i) {
                    for (k=0; k < refineFactor; ++k) {
                        cells.rhoBCsLR[i*refineFactor + k][0] = coarseCells.rhoBCsLR[i][0];
                        cells.MxBCsLR[i*refineFactor + k][0] = cells.Mx[i*refineFactor + k][0];
                        cells.MyBCsLR[i*refineFactor + k][0] = 0.0;
                        // For Roe Solver
                        cells.rhoBCsLR[i*refineFactor + k][2] = coarseCells.rhoBCsLR[i][2];
                        cells.MxBCsLR[i*refineFactor + k][2] = cells.Mx[i*refineFactor + k][0];
                        cells.MyBCsLR[i*refineFactor + k][2] = 0.0;
                    }
                }
            }
            else if (m == 1) { // Right
                for (n=0; n < cells.yPoints; ++n) {
                    cells.rhoBCsLR[n][1] = cells.rho[n][cells.xPoints-1];
                    cells.MxBCsLR[n][1] = cells.Mx[n][cells.xPoints-1];
                    cells.MyBCsLR[n][1] = cells.My[n][cells.xPoints-1];
                    // For Roe Solver
                    cells.rhoBCsLR[n][3] = cells.rho[n][cells.xPoints-2];
                    cells.MxBCsLR[n][3] = cells.Mx[n][cells.xPoints-2];
                    cells.MyBCsLR[n][3] = cells.My[n][cells.xPoints-2];
                }
            }
            else if (m == 2) { // Bottom
                for (n=0; n < cells.xPoints; ++n) {
                    cells.rhoBCsBT[0][n] = cells.rho[0][n];
                    cells.MxBCsBT[0][n] = cells.Mx[0][n];
                    cells.MyBCsBT[0][n] = -cells.My[0][n];
                    // For Roe Solver
                    cells.rhoBCsBT[2][n] = cells.rho[1][n];
                    cells.MxBCsBT[2][n] = cells.Mx[1][n];
                    cells.MyBCsBT[2][n] = -cells.My[1][n];
                }
            }
            else { // Top
                for (n=0; n < cells.xPoints; ++n) {
                    cells.rhoBCsBT[1][n] = cells.rho[cells.yPoints-1][n];
                    cells.MxBCsBT[1][n] = cells.Mx[cells.yPoints-1][n];
                    cells.MyBCsBT[1][n] = -cells.My[cells.yPoints-1][n];
                    // For Roe Solver
                    cells.rhoBCsBT[3][n] = cells.rho[cells.yPoints-2][n];
                    cells.MxBCsBT[3][n] = cells.Mx[cells.yPoints-2][n];
                    cells.MyBCsBT[3][n] = -cells.My[cells.yPoints-2][n];
                }
            }
        }
        else {
            if (m == 0) { // Left
                for (i=0; i < coarseCells.yPoints; ++i) {
                    for (k=0; k < refineFactor; ++k) {
                        cells.rhoBCsLR[i*refineFactor + k][0] = coarseCells.rho[i][refineColumns[0] - 1];
                        cells.MxBCsLR[i*refineFactor + k][0] = coarseCells.rho[i][refineColumns[0] - 1];
                        cells.MyBCsLR[i*refineFactor + k][0] = coarseCells.rho[i][refineColumns[0] - 1];
                        // For Roe Solver
                        cells.rhoBCsLR[i*refineFactor + k][2] = coarseCells.rho[i][refineColumns[0] - 1];
                        cells.MxBCsLR[i*refineFactor + k][2] = coarseCells.rho[i][refineColumns[0] - 1];
                        cells.MyBCsLR[i*refineFactor + k][2] = coarseCells.rho[i][refineColumns[0] - 1];
                    }
                }
            }
            else if (m == 1) { // Right
                for (i=0; i < coarseCells.yPoints; ++i) {
                    for (k=0; k < refineFactor; ++k) {
                        cells.rhoBCsLR[i*refineFactor + k][1] = coarseCells.rho[i][refineColumns.back() + 1];
                        cells.MxBCsLR[i*refineFactor + k][1] = cells.Mx[i][refineColumns.back() + 1];
                        cells.MyBCsLR[i*refineFactor + k][1] = cells.My[i][refineColumns.back() + 1];
                        // For Roe Solver
                        cells.rhoBCsLR[i*refineFactor + k][3] = cells.rho[i][refineColumns.back() + 1];
                        cells.MxBCsLR[i*refineFactor + k][3] = cells.Mx[i][refineColumns.back() + 1];
                        cells.MyBCsLR[i*refineFactor + k][3] = cells.My[i][refineColumns.back() + 1];
                    }
                }
            }
            else {
                cout << "Error in updating refinement cells in updateBCs:line 132\n";
            }
        }
    }
}
// Refine cell interpolated boundary updates
void updateBCs
(
    grid& cells,
    grid& coarseCells,
    vector<int> refineColumns,
    vector<int> oldRhoValues,
    vector<int> oldMxValues,
    vector<int> oldMyValues,
    const int refineFactor,
    int timeSpot
) 
{
    int i, j, m, n, k;
    for (m=0; m < 4; ++m) {
        if (!cells.connect[m]) {
            if (m == 0) { // Left
                for (i=0; i < coarseCells.yPoints; ++i) {
                    for (k=0; k < refineFactor; ++k) {
                        cells.rhoBCsLR[i*refineFactor + k][0] = coarseCells.rhoBCsLR[i][0];
                        cells.MxBCsLR[i*refineFactor + k][0] = cells.Mx[i*refineFactor + k][0];
                        cells.MyBCsLR[i*refineFactor + k][0] = 0.0;
                        // For Roe Solver
                        cells.rhoBCsLR[i*refineFactor + k][2] = coarseCells.rhoBCsLR[i][2];
                        cells.MxBCsLR[i*refineFactor + k][2] = cells.Mx[i*refineFactor + k][0];
                        cells.MyBCsLR[i*refineFactor + k][2] = 0.0;
                    }
                }
            }
            else if (m == 1) { // Right
                for (n=0; n < cells.yPoints; ++n) {
                    cells.rhoBCsLR[n][1] = cells.rho[n][cells.xPoints-1];
                    cells.MxBCsLR[n][1] = cells.Mx[n][cells.xPoints-1];
                    cells.MyBCsLR[n][1] = cells.My[n][cells.xPoints-1];
                    // For Roe Solver
                    cells.rhoBCsLR[n][3] = cells.rho[n][cells.xPoints-2];
                    cells.MxBCsLR[n][3] = cells.Mx[n][cells.xPoints-2];
                    cells.MyBCsLR[n][3] = cells.My[n][cells.xPoints-2];
                }
            }
            else if (m == 2) { // Bottom
                for (n=0; n < cells.xPoints; ++n) {
                    cells.rhoBCsBT[0][n] = cells.rho[0][n];
                    cells.MxBCsBT[0][n] = cells.Mx[0][n];
                    cells.MyBCsBT[0][n] = -cells.My[0][n];
                    // For Roe Solver
                    cells.rhoBCsBT[2][n] = cells.rho[1][n];
                    cells.MxBCsBT[2][n] = cells.Mx[1][n];
                    cells.MyBCsBT[2][n] = -cells.My[1][n];
                }
            }
            else { // Top
                for (n=0; n < cells.xPoints; ++n) {
                    cells.rhoBCsBT[1][n] = cells.rho[cells.yPoints-1][n];
                    cells.MxBCsBT[1][n] = cells.Mx[cells.yPoints-1][n];
                    cells.MyBCsBT[1][n] = -cells.My[cells.yPoints-1][n];
                    // For Roe Solver
                    cells.rhoBCsBT[3][n] = cells.rho[cells.yPoints-2][n];
                    cells.MxBCsBT[3][n] = cells.Mx[cells.yPoints-2][n];
                    cells.MyBCsBT[3][n] = -cells.My[cells.yPoints-2][n];
                }
            }
        }
        else {
            double rhoNow, MxNow, MyNow, rhoNowPlus, MxNowPlus, MyNowPlus, 
                rhoY, MxY, MyY;
            if (m == 0) { // Left
                // for (i=-1; i < coarseCells.yPoints; ++i) {
                //     if (i == -1) {
                //     }
                //     else if (i == coarseCells.yPoints) {

                //     }
                //     else {
                //         rhoNow = (coarseCells.rho[i][refineColumns[0]-1] - oldRhoValues[i])/refineFactor*timeSpot + oldRhoValues[i];
                //         MxNow = (coarseCells.Mx[i][refineColumns[0]-1] - oldMxValues[i])/refineFactor*timeSpot + oldMxValues[i];
                //         MyNow = (coarseCells.My[i][refineColumns[0]-1] - oldMyValues[i])/refineFactor*timeSpot + oldMyValues[i];
                //         rhoNowPlus = (coarseCells.rho[i+1][refineColumns[0]-1] - oldRhoValues[i+1])/refineFactor*timeSpot + oldRhoValues[i+1];
                //         MxNowPlus = (coarseCells.Mx[i+1][refineColumns[0]-1] - oldMxValues[i+1])/refineFactor*timeSpot + oldMxValues[i+1];
                //         MyNowPlus = (coarseCells.My[i+1][refineColumns[0]-1] - oldMyValues[i+1])/refineFactor*timeSpot + oldMyValues[i+1];
                //         for (k=1; k < 2*refineFactor; k = 2 + k) {
                //             rhoY = (rhoNowPlus - rhoNow)/(2*refineFactor)*k + rhoNow;
                //             MxY = (MxNowPlus - MxNow)/(2*refineFactor)*k + MxNow;
                //             MyY = (MyNowPlus - MyNow)/(2*refineFactor)*k + MyNow;
                //             cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                //             cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                //             cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                //             // For Roe Solver
                //             cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                //             cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                //             cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                //         }
                //     }
                    
                // }
            }
            else if (m == 1) { // Right
                for (i=-1; i < coarseCells.yPoints; ++i) {
                    if (i == -1) {
                    }
                    else if (i == coarseCells.yPoints) {

                    }
                    else {
                        rhoNow = (coarseCells.rho[i][refineColumns[0]-1] - oldRhoValues[i])/refineFactor*timeSpot + oldRhoValues[i];
                        MxNow = (coarseCells.Mx[i][refineColumns[0]-1] - oldMxValues[i])/refineFactor*timeSpot + oldMxValues[i];
                        MyNow = (coarseCells.My[i][refineColumns[0]-1] - oldMyValues[i])/refineFactor*timeSpot + oldMyValues[i];
                        rhoNowPlus = (coarseCells.rho[i+1][refineColumns[0]-1] - oldRhoValues[i+1])/refineFactor*timeSpot + oldRhoValues[i+1];
                        MxNowPlus = (coarseCells.Mx[i+1][refineColumns[0]-1] - oldMxValues[i+1])/refineFactor*timeSpot + oldMxValues[i+1];
                        MyNowPlus = (coarseCells.My[i+1][refineColumns[0]-1] - oldMyValues[i+1])/refineFactor*timeSpot + oldMyValues[i+1];
                        for (k=1; k < 2*refineFactor; k = 2 + k) {
                            rhoY = (rhoNowPlus - rhoNow)/(2*refineFactor)*k + rhoNow;
                            MxY = (MxNowPlus - MxNow)/(2*refineFactor)*k + MxNow;
                            MyY = (MyNowPlus - MyNow)/(2*refineFactor)*k + MyNow;
                            cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                            cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                            cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                            // For Roe Solver
                            cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                            cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                            cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                        }
                    }
                    
                }
            }
            else {
                cout << "Error in updating refinement cells in updateBCs:line 132\n";
            }
        }
    }
}