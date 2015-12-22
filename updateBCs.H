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
        else {
            // transfer left right bcs

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
            if (m == 0 && refineColumns[0] != 0) { // Left
                if (refineColumns[0] != 0) {
                    for (i=0; i < coarseCells.yPoints+2; ++i) {
                        if (i == 0) {
                            coarseCells.oldRhoLeft[i] = coarseCells.rhoBCsBT[0][refineColumns[0]-1];
                            coarseCells.oldMxLeft[i] = coarseCells.MxBCsBT[0][refineColumns[0]-1];
                            coarseCells.oldMyLeft[i] = coarseCells.MyBCsBT[0][refineColumns[0]-1];
                        }
                        else if (i == coarseCells.yPoints+1) {
                            coarseCells.oldRhoLeft[i] = coarseCells.rhoBCsBT[1][refineColumns[0]-1];
                            coarseCells.oldMxLeft[i] = coarseCells.MxBCsBT[1][refineColumns[0]-1];
                            coarseCells.oldMyLeft[i] = coarseCells.MyBCsBT[1][refineColumns[0]-1];
                        }
                        else {
                            coarseCells.oldRhoLeft[i] = coarseCells.rho[i-1][refineColumns[0]-1];
                            coarseCells.oldMxLeft[i] = coarseCells.Mx[i-1][refineColumns[0]-1];
                            coarseCells.oldMyLeft[i] = coarseCells.My[i-1][refineColumns[0]-1];
                        }
                    }
                }
                else {
                    // Exchange refineColumns

                }
            }
            else if (m == 1) { // Right
                if (refineColumns.back == coarseCells.xPoints-1) {
                    for (i=0; i < coarseCells.yPoints+2; ++i) {
                        if (i == 0) {
                            coarseCells.oldRhoRight[i] = coarseCells.rhoBCsBT[0][refineColumns.back()+1];
                            coarseCells.oldMxRight[i] = coarseCells.MxBCsBT[0][refineColumns.back()+1];
                            coarseCells.oldMyRight[i] = coarseCells.MyBCsBT[0][refineColumns.back()+1];
                        }
                        else if (i == coarseCells.yPoints+1) {
                            coarseCells.oldRhoRight[i] = coarseCells.rhoBCsBT[1][refineColumns.back()+1];
                            coarseCells.oldMxRight[i] = coarseCells.MxBCsBT[1][refineColumns.back()+1];
                            coarseCells.oldMyRight[i] = coarseCells.MyBCsBT[1][refineColumns.back()+1];
                        }
                        else {
                            coarseCells.oldRhoRight[i] = coarseCells.rho[i-1][refineColumns.back()+1];
                            coarseCells.oldMxRight[i] = coarseCells.Mx[i-1][refineColumns.back()+1];
                            coarseCells.oldMyRight[i] = coarseCells.My[i-1][refineColumns.back()+1];
                        }
                    }
                }
                else {
                    //Exchange refineColumns

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
    const int refineFactor,
    double timeSpot
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
                if (refineColumns[0] != 0) {
                    for (i=-1; i < coarseCells.yPoints; ++i) {
                        if (i == -1) {
                            rhoNow = (coarseCells.rhoBCsBT[0][refineColumns[0]-1] - coarseCells.oldRhoLeft[i+1])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+1];
                            MxNow = (coarseCells.MxBCsBT[0][refineColumns[0]-1] - coarseCells.oldMxLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+1];
                            MyNow = (coarseCells.MyBCsBT[0][refineColumns[0]-1] - coarseCells.oldMyLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+1];
                            rhoNowPlus = (coarseCells.rho[i+1][refineColumns[0]-1] - coarseCells.oldRhoLeft[i+2])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+2];
                            MxNowPlus = (coarseCells.Mx[i+1][refineColumns[0]-1] - coarseCells.oldMxLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+2];
                            MyNowPlus = (coarseCells.My[i+1][refineColumns[0]-1] - coarseCells.oldMyLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+2];
                            for (k=1+refineFactor; k < 2*refineFactor; k = 2 + k) {
                                rhoY = (rhoNowPlus - rhoNow)/(2*refineFactor)*k + rhoNow;
                                MxY = (MxNowPlus - MxNow)/(2*refineFactor)*k + MxNow;
                                MyY = (MyNowPlus - MyNow)/(2*refineFactor)*k + MyNow;
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                // For Roe Solver
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*4 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*4 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0];
                            }
                        }
                        else if (i == coarseCells.yPoints-1) {
                            rhoNow = (coarseCells.rho[i][refineColumns[0]-1] - coarseCells.oldRhoLeft[i+1])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+1];
                            MxNow = (coarseCells.Mx[i][refineColumns[0]-1] - coarseCells.oldMxLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+1];
                            MyNow = (coarseCells.My[i][refineColumns[0]-1] - coarseCells.oldMyLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+1];
                            rhoNowPlus = (coarseCells.rhoBCsBT[1][refineColumns[0]-1] - coarseCells.oldRhoLeft[i+2])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+2];
                            MxNowPlus = (coarseCells.MxBCsBT[1][refineColumns[0]-1] - coarseCells.oldMxLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+2];
                            MyNowPlus = (coarseCells.MyBCsBT[1][refineColumns[0]-1] - coarseCells.oldMyLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+2];
                            for (k=1; k < refineFactor; k = 2 + k) {
                                rhoY = (rhoNowPlus - rhoNow)/(2*refineFactor)*k + rhoNow;
                                MxY = (MxNowPlus - MxNow)/(2*refineFactor)*k + MxNow;
                                MyY = (MyNowPlus - MyNow)/(2*refineFactor)*k + MyNow;
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                // For Roe Solver
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*4 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*4 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0];
                            }
                        }
                        else {
                            rhoNow = (coarseCells.rho[i][refineColumns[0]-1] - coarseCells.oldRhoLeft[i+1])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+1];
                            MxNow = (coarseCells.Mx[i][refineColumns[0]-1] - coarseCells.oldMxLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+1];
                            MyNow = (coarseCells.My[i][refineColumns[0]-1] - coarseCells.oldMyLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+1];
                            rhoNowPlus = (coarseCells.rho[i+1][refineColumns[0]-1] - coarseCells.oldRhoLeft[i+2])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+2];
                            MxNowPlus = (coarseCells.Mx[i+1][refineColumns[0]-1] - coarseCells.oldMxLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+2];
                            MyNowPlus = (coarseCells.My[i+1][refineColumns[0]-1] - coarseCells.oldMyLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+2];
                            for (k=1; k < 2*refineFactor; k = 2 + k) {
                                rhoY = (rhoNowPlus - rhoNow)/(2*refineFactor)*k + rhoNow;
                                MxY = (MxNowPlus - MxNow)/(2*refineFactor)*k + MxNow;
                                MyY = (MyNowPlus - MyNow)/(2*refineFactor)*k + MyNow;
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][0] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                // For Roe Solver
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*4 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*4 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][0];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][2] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][0];
                            }
                        }
                    }
                }
                else {
                    // Exchange refineCol cells

                }
            }
            else if (m == 1) { // Right
                if (refineColumns.back == coarseCells.xPoints-1) {
                    for (i=-1; i < coarseCells.yPoints; ++i) {
                        if (i == -1) {
                            rhoNow = (coarseCells.rhoBCsBT[0][refineColumns.back()+1] - coarseCells.oldRhoRight[i+1])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+1];
                            MxNow = (coarseCells.MxBCsBT[0][refineColumns.back()+1] - coarseCells.oldMxRight[i+1])/refineFactor*timeSpot + coarseCells.oldMxRight[i+1];
                            MyNow = (coarseCells.MyBCsBT[0][refineColumns.back()+1] - coarseCells.oldMyRight[i+1])/refineFactor*timeSpot + coarseCells.oldMyRight[i+1];
                            rhoNowPlus = (coarseCells.rho[i+1][refineColumns.back()+1] - coarseCells.oldRhoRight[i+2])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+2];
                            MxNowPlus = (coarseCells.Mx[i+1][refineColumns.back()+1] - coarseCells.oldMxRight[i+2])/refineFactor*timeSpot + coarseCells.oldMxRight[i+2];
                            MyNowPlus = (coarseCells.My[i+1][refineColumns.back()+1] - coarseCells.oldMyRight[i+2])/refineFactor*timeSpot + coarseCells.oldMyRight[i+2];
                            for (k=1+refineFactor; k < 2*refineFactor; k = 2+k) {
                                rhoY = (rhoNowPlus - rhoNow)/(2*refineFactor)*k + rhoNow;
                                MxY = (MxNowPlus - MxNow)/(2*refineFactor)*k + MxNow;
                                MyY = (MyNowPlus - MyNow)/(2*refineFactor)*k + MyNow;
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                // For Roe Solver
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                            }
                        }
                        else if (i == coarseCells.yPoints-1) {
                            rhoNow = (coarseCells.rho[i][refineColumns.back()+1] - coarseCells.oldRhoRight[i+1])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+1];
                            MxNow = (coarseCells.Mx[i][refineColumns.back()+1] - coarseCells.oldMxRight[i+1])/refineFactor*timeSpot + coarseCells.oldMxRight[i+1];
                            MyNow = (coarseCells.My[i][refineColumns.back()+1] - coarseCells.oldMyRight[i+1])/refineFactor*timeSpot + coarseCells.oldMyRight[i+1];
                            rhoNowPlus = (coarseCells.rhoBCsBT[1][refineColumns.back()+1] - coarseCells.oldRhoRight[i+2])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+2];
                            MxNowPlus = (coarseCells.MxBCsBT[1][refineColumns.back()+1] - coarseCells.oldMxRight[i+2])/refineFactor*timeSpot + coarseCells.oldMxRight[i+2];
                            MyNowPlus = (coarseCells.MyBCsBT[1][refineColumns.back()+1] - coarseCells.oldMyRight[i+2])/refineFactor*timeSpot + coarseCells.oldMyRight[i+2];
                            for (k=1; k < refineFactor; k = 2 + k) {
                                rhoY = (rhoNowPlus - rhoNow)/(2*refineFactor)*k + rhoNow;
                                MxY = (MxNowPlus - MxNow)/(2*refineFactor)*k + MxNow;
                                MyY = (MyNowPlus - MyNow)/(2*refineFactor)*k + MyNow;
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                // For Roe Solver
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                            }
                        }
                        else {
                            rhoNow = (coarseCells.rho[i][refineColumns.back()+1] - coarseCells.oldRhoRight[i+1])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+1];
                            MxNow = (coarseCells.Mx[i][refineColumns.back()+1] - coarseCells.oldMxRight[i+1])/refineFactor*timeSpot + coarseCells.oldMxRight[i+1];
                            MyNow = (coarseCells.My[i][refineColumns.back()+1] - coarseCells.oldMyRight[i+1])/refineFactor*timeSpot + coarseCells.oldMyRight[i+1];
                            rhoNowPlus = (coarseCells.rho[i+1][refineColumns.back()+1] - coarseCells.oldRhoRight[i+2])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+2];
                            MxNowPlus = (coarseCells.Mx[i+1][refineColumns.back()+1] - coarseCells.oldMxRight[i+2])/refineFactor*timeSpot + coarseCells.oldMxRight[i+2];
                            MyNowPlus = (coarseCells.My[i+1][refineColumns.back()+1] - coarseCells.oldMyRight[i+2])/refineFactor*timeSpot + coarseCells.oldMyRight[i+2];
                            for (k=1; k < 2*refineFactor; k = 2 + k) {
                                rhoY = (rhoNowPlus - rhoNow)/(2*refineFactor)*k + rhoNow;
                                MxY = (MxNowPlus - MxNow)/(2*refineFactor)*k + MxNow;
                                MyY = (MyNowPlus - MyNow)/(2*refineFactor)*k + MyNow;
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][1] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                // For Roe Solver
                                cells.rhoBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (rhoY - cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.rho[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MxBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (MxY - cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*4 + cells.Mx[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                                cells.MyBCsLR[i*refineFactor+refineFactor/2+(k-1)/2][3] = (MyY - cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1])/(refineFactor + 1)*2 + cells.My[i*refineFactor+refineFactor/2+(k-1)/2][cells.xPoints-1];
                            }
                        }   
                    }
                }
                else {
                    // Exchange refineCol cells
                    
                }
            }
            else {
                cout << "Error in updating refinement cells in updateBCs:line 132\n";
            }
        }
    }
}