// Coarse cell boundary condition updates
void updateBCs
(
    grid& cells,
    int rank,
    double sendLeft[],
    double recvLeft[],
    double sendRight[],
    double recvRight[]
) 
{
    int i, n, ierr, phi, count = 0;
    MPI_Status  status[2];
    MPI_Request send_request[2], recv_request[2];
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
            if (i == 0) { // Left
                for (phi=0; phi < 6; ++phi) {
                    for (n=0; n < cells.yPoints; ++n) {
                        switch (phi) {
                            case 0: sendLeft[n] = cells.rho[n][0];
                                        break;
                            case 1: sendLeft[n + cells.yPoints] = cells.Mx[n][0];
                                        break;
                            case 2: sendLeft[n + phi*cells.yPoints] = cells.My[n][0];
                                        break;
                            case 3: sendLeft[n + phi*cells.yPoints] = cells.rho[n][1];
                                        break;
                            case 4: sendLeft[n + phi*cells.yPoints] = cells.Mx[n][1];
                                        break;
                            case 5: sendLeft[n + phi*cells.yPoints] = cells.My[n][1];
                                        break;
                        }
                    }
                }
                ierr = MPI_Isend(&sendLeft[0], 6*cells.yPoints, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &send_request[0]);
                ierr = MPI_Irecv(&recvLeft[0], 6*cells.yPoints, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &recv_request[0]);
                ++count;
            }
            else if (i == 1) { // Right
                for (phi=0; phi < 6; ++phi) {
                    for (n=0; n < cells.yPoints; ++n) {
                        switch (phi) {
                            case 0: sendRight[n] = cells.rho[n][cells.xPoints-1];
                                        break; 
                            case 1: sendRight[n + cells.yPoints] = cells.Mx[n][cells.xPoints-1];
                                        break;
                            case 2: sendRight[n + phi*cells.yPoints] = cells.My[n][cells.xPoints-1];
                                        break;
                            case 3: sendRight[n + phi*cells.yPoints] = cells.rho[n][cells.xPoints-2];
                                        break;
                            case 4: sendRight[n + phi*cells.yPoints] = cells.Mx[n][cells.xPoints-2];
                                        break;
                            case 5: sendRight[n + phi*cells.yPoints] = cells.My[n][cells.xPoints-2];
                                        break;
                        }
                    }
                }
                if (count == 0) {
                    ierr = MPI_Isend(&sendRight[0], 6*cells.yPoints, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &send_request[0]);
                    ierr = MPI_Irecv(&recvRight[0], 6*cells.yPoints, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &recv_request[0]);
                }
                else {
                    ierr = MPI_Isend(&sendRight[0], 6*cells.yPoints, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &send_request[1]);
                    ierr = MPI_Irecv(&recvRight[0], 6*cells.yPoints, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &recv_request[1]);
                }
                ++count;
            }
            else {
                cout << "Error: Trying to transfer top and bottom cells!\n"
                    << "updateBCs.H line 112 \n";
            }
        }
    }
    if (count > 0) {
        ierr = MPI_Waitall(count, recv_request, status);
        if (cells.connect[0]) {
            for (phi=0; phi < 6; ++phi) {
                for (n=0; n < cells.yPoints; ++n) {
                    switch (phi) {
                        case 0: cells.rhoBCsLR[n][0] = recvLeft[n];
                                    break;
                        case 1: cells.MxBCsLR[n][0] = recvLeft[n + cells.yPoints];
                                    break;
                        case 2: cells.MyBCsLR[n][0] = recvLeft[n + phi*cells.yPoints];
                                    break;
                        case 3: cells.rhoBCsLR[n][2] = recvLeft[n + phi*cells.yPoints];
                                    break;
                        case 4: cells.MxBCsLR[n][2] = recvLeft[n + phi*cells.yPoints];
                                    break;
                        case 5: cells.MyBCsLR[n][2] = recvLeft[n + phi*cells.yPoints];
                                    break;
                    }
                }
            }
        }
        if (cells.connect[1]) {
            for (phi=0; phi < 6; ++phi) {
                for (n=0; n < cells.yPoints; ++n) {
                    switch (phi) {
                        case 0: cells.rhoBCsLR[n][1] = recvRight[n];
                                    break;
                        case 1: cells.MxBCsLR[n][1] = recvRight[n + cells.yPoints];
                                    break;
                        case 2: cells.MyBCsLR[n][1] = recvRight[n + phi*cells.yPoints];
                                    break;
                        case 3: cells.rhoBCsLR[n][3] = recvRight[n + phi*cells.yPoints];
                                    break;
                        case 4: cells.MxBCsLR[n][3] = recvRight[n + phi*cells.yPoints];
                                    break;
                        case 5: cells.MyBCsLR[n][3] = recvRight[n + phi*cells.yPoints];
                                    break;
                    }
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
    const int refineFactor,
    int rank,
    double sendLeft[],
    double recvLeft[],
    double sendRight[],
    double recvRight[],
    int vectorSize[],
    int vectorSizeSend[]
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
                else if (vectorSize[0] == 0 && vectorSizeSend[0] == 0) {
                    for (i=0; i < coarseCells.yPoints+2; ++i) {
                        if (i == 0) {
                            coarseCells.oldRhoLeft[i] = coarseCells.rhoBCsLR[0][0];
                            coarseCells.oldMxLeft[i] = coarseCells.MxBCsLR[0][0];
                            coarseCells.oldMyLeft[i] = -coarseCells.MyBCsLR[0][0];
                        }
                        else if (i == coarseCells.yPoints+1) {
                            coarseCells.oldRhoLeft[i] = coarseCells.rhoBCsLR[coarseCells.yPoints-1][0];
                            coarseCells.oldMxLeft[i] = coarseCells.MxBCsLR[coarseCells.yPoints-1][0];
                            coarseCells.oldMyLeft[i] = -coarseCells.MyBCsLR[coarseCells.yPoints-1][0];
                        }
                        else {
                            coarseCells.oldRhoLeft[i] = coarseCells.rhoBCsLR[i-1][0];
                            coarseCells.oldMxLeft[i] = coarseCells.MxBCsLR[i-1][0];
                            coarseCells.oldMyLeft[i] = coarseCells.MyBCsLR[i-1][0];
                        }
                    }
                }
            }
            else if (m == 1) { // Right
                if (refineColumns.back() != coarseCells.xPoints-1) {
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
                else if (vectorSize[1] == 0 && vectorSizeSend[1] == 0) {
                    for (i=0; i < coarseCells.yPoints+2; ++i) {
                        if (i == 0) {
                            coarseCells.oldRhoLeft[i] = coarseCells.rhoBCsLR[0][1];
                            coarseCells.oldMxLeft[i] = coarseCells.MxBCsLR[0][1];
                            coarseCells.oldMyLeft[i] = -coarseCells.MyBCsLR[0][1];
                        }
                        else if (i == coarseCells.yPoints+1) {
                            coarseCells.oldRhoLeft[i] = coarseCells.rhoBCsLR[coarseCells.yPoints-1][1];
                            coarseCells.oldMxLeft[i] = coarseCells.MxBCsLR[coarseCells.yPoints-1][1];
                            coarseCells.oldMyLeft[i] = -coarseCells.MyBCsLR[coarseCells.yPoints-1][1];
                        }
                        else {
                            coarseCells.oldRhoLeft[i] = coarseCells.rhoBCsLR[i-1][1];
                            coarseCells.oldMxLeft[i] = coarseCells.MxBCsLR[i-1][1];
                            coarseCells.oldMyLeft[i] = coarseCells.MyBCsLR[i-1][1];
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

// Refine cell interpolated boundary updates
void updateBCs
(
    grid& cells,
    grid& coarseCells,
    vector<int> refineColumns,
    const int refineFactor,
    double timeSpot,
    int rank,
    double sendLeft[],
    double recvLeft[],
    double sendRight[],
    double recvRight[],
    int vectorSize[],
    int vectorSizeSend[]
) 
{
    int i, j, m, n, k, ierr, phi, count = 0;
    MPI_Status  status[2];
    MPI_Request send_request[2], recv_request[2];
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
                else if (vectorSize[0] == 0 && vectorSizeSend[0] == 0) {
                    for (i=-1; i < coarseCells.yPoints; ++i) {
                        if (i == -1) {
                            rhoNow = (coarseCells.rhoBCsLR[0][0] - coarseCells.oldRhoLeft[i+1])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+1];
                            MxNow = (coarseCells.MxBCsLR[0][0] - coarseCells.oldMxLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+1];
                            MyNow = (-coarseCells.MyBCsLR[0][0] - coarseCells.oldMyLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+1];
                            rhoNowPlus = (coarseCells.rhoBCsLR[i+1][0] - coarseCells.oldRhoLeft[i+2])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+2];
                            MxNowPlus = (coarseCells.MxBCsLR[i+1][0] - coarseCells.oldMxLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+2];
                            MyNowPlus = (coarseCells.MyBCsLR[i+1][0] - coarseCells.oldMyLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+2];
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
                            rhoNow = (coarseCells.rhoBCsLR[i][0] - coarseCells.oldRhoLeft[i+1])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+1];
                            MxNow = (coarseCells.MxBCsLR[i][0] - coarseCells.oldMxLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+1];
                            MyNow = (coarseCells.MyBCsLR[i][0] - coarseCells.oldMyLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+1];
                            rhoNowPlus = (coarseCells.rhoBCsLR[i][0] - coarseCells.oldRhoLeft[i+2])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+2];
                            MxNowPlus = (coarseCells.MxBCsLR[i][0] - coarseCells.oldMxLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+2];
                            MyNowPlus = (-coarseCells.MyBCsLR[i][0] - coarseCells.oldMyLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+2];
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
                            rhoNow = (coarseCells.rhoBCsLR[i][0] - coarseCells.oldRhoLeft[i+1])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+1];
                            MxNow = (coarseCells.MxBCsLR[i][0] - coarseCells.oldMxLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+1];
                            MyNow = (coarseCells.MyBCsLR[i][0] - coarseCells.oldMyLeft[i+1])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+1];
                            rhoNowPlus = (coarseCells.rhoBCsLR[i+1][0] - coarseCells.oldRhoLeft[i+2])/refineFactor*timeSpot + coarseCells.oldRhoLeft[i+2];
                            MxNowPlus = (coarseCells.MxBCsLR[i+1][0] - coarseCells.oldMxLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMxLeft[i+2];
                            MyNowPlus = (coarseCells.MyBCsLR[i+1][0] - coarseCells.oldMyLeft[i+2])/refineFactor*timeSpot + coarseCells.oldMyLeft[i+2];
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
                else { // Exchange refine cells across processors
                    for (phi=0; phi < 6; ++phi) {
                        for (n=0; n < cells.yPoints; ++n) {
                            switch (phi) {
                                case 0: sendLeft[n] = cells.rho[n][0];
                                            break;
                                case 1: sendLeft[n + cells.yPoints] = cells.Mx[n][0];
                                            break;
                                case 2: sendLeft[n + phi*cells.yPoints] = cells.My[n][0];
                                            break;
                                case 3: sendLeft[n + phi*cells.yPoints] = cells.rho[n][1];
                                            break;
                                case 4: sendLeft[n + phi*cells.yPoints] = cells.Mx[n][1];
                                            break;
                                case 5: sendLeft[n + phi*cells.yPoints] = cells.My[n][1];
                                            break;
                            }
                        }
                    }
                    ierr = MPI_Isend(&sendLeft[0], 6*cells.yPoints, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &send_request[0]);
                    ierr = MPI_Irecv(&recvLeft[0], 6*cells.yPoints, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &recv_request[0]);
                    ++count;
                }
            }
            else if (m == 1) { // Right
                if (refineColumns.back() != coarseCells.xPoints-1) {
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
                else if (vectorSize[1] == 0 && vectorSizeSend[1] == 0) {
                    for (i=-1; i < coarseCells.yPoints; ++i) {
                        if (i == -1) {
                            rhoNow = (coarseCells.rhoBCsLR[0][1] - coarseCells.oldRhoRight[i+1])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+1];
                            MxNow = (coarseCells.MxBCsLR[0][1] - coarseCells.oldMxRight[i+1])/refineFactor*timeSpot + coarseCells.oldMxRight[i+1];
                            MyNow = (-coarseCells.MyBCsLR[0][1] - coarseCells.oldMyRight[i+1])/refineFactor*timeSpot + coarseCells.oldMyRight[i+1];
                            rhoNowPlus = (coarseCells.rhoBCsLR[i+1][1] - coarseCells.oldRhoRight[i+2])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+2];
                            MxNowPlus = (coarseCells.MxBCsLR[i+1][1] - coarseCells.oldMxRight[i+2])/refineFactor*timeSpot + coarseCells.oldMxRight[i+2];
                            MyNowPlus = (coarseCells.MyBCsLR[i+1][1] - coarseCells.oldMyRight[i+2])/refineFactor*timeSpot + coarseCells.oldMyRight[i+2];
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
                            rhoNow = (coarseCells.rhoBCsLR[i][1] - coarseCells.oldRhoRight[i+1])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+1];
                            MxNow = (coarseCells.MxBCsLR[i][1] - coarseCells.oldMxRight[i+1])/refineFactor*timeSpot + coarseCells.oldMxRight[i+1];
                            MyNow = (coarseCells.MyBCsLR[i][1] - coarseCells.oldMyRight[i+1])/refineFactor*timeSpot + coarseCells.oldMyRight[i+1];
                            rhoNowPlus = (coarseCells.rhoBCsLR[i][1] - coarseCells.oldRhoRight[i+2])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+2];
                            MxNowPlus = (coarseCells.MxBCsLR[i][1] - coarseCells.oldMxRight[i+2])/refineFactor*timeSpot + coarseCells.oldMxRight[i+2];
                            MyNowPlus = (-coarseCells.MyBCsLR[i][1] - coarseCells.oldMyRight[i+2])/refineFactor*timeSpot + coarseCells.oldMyRight[i+2];
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
                            rhoNow = (coarseCells.rhoBCsLR[i][1] - coarseCells.oldRhoRight[i+1])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+1];
                            MxNow = (coarseCells.MxBCsLR[i][1] - coarseCells.oldMxRight[i+1])/refineFactor*timeSpot + coarseCells.oldMxRight[i+1];
                            MyNow = (coarseCells.MyBCsLR[i][1] - coarseCells.oldMyRight[i+1])/refineFactor*timeSpot + coarseCells.oldMyRight[i+1];
                            rhoNowPlus = (coarseCells.rhoBCsLR[i+1][1] - coarseCells.oldRhoRight[i+2])/refineFactor*timeSpot + coarseCells.oldRhoRight[i+2];
                            MxNowPlus = (coarseCells.MxBCsLR[i+1][1] - coarseCells.oldMxRight[i+2])/refineFactor*timeSpot + coarseCells.oldMxRight[i+2];
                            MyNowPlus = (coarseCells.MyBCsLR[i+1][1] - coarseCells.oldMyRight[i+2])/refineFactor*timeSpot + coarseCells.oldMyRight[i+2];
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
                else { // Exchange refine cells across processors
                    for (phi=0; phi < 6; ++phi) {
                        for (n=0; n < cells.yPoints; ++n) {
                            switch (phi) {
                                case 0: sendRight[n] = cells.rho[n][cells.xPoints-1];
                                            break;
                                case 1: sendRight[n + cells.yPoints] = cells.Mx[n][cells.xPoints-1];
                                            break;
                                case 2: sendRight[n + phi*cells.yPoints] = cells.My[n][cells.xPoints-1];
                                            break;
                                case 3: sendRight[n + phi*cells.yPoints] = cells.rho[n][cells.xPoints-2];
                                            break;
                                case 4: sendRight[n + phi*cells.yPoints] = cells.Mx[n][cells.xPoints-2];
                                            break;
                                case 5: sendRight[n + phi*cells.yPoints] = cells.My[n][cells.xPoints-2];
                                            break;
                            }
                        }
                    }
                    if (count == 0) {
                        ierr = MPI_Isend(&sendRight[0], 6*cells.yPoints, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &send_request[0]);
                        ierr = MPI_Irecv(&recvRight[0], 6*cells.yPoints, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &recv_request[0]);
                    }
                    else {
                        ierr = MPI_Isend(&sendRight[0], 6*cells.yPoints, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &send_request[1]);
                        ierr = MPI_Irecv(&recvRight[0], 6*cells.yPoints, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &recv_request[1]);
                    }
                    ++count;
                }
            }
            else {
                cout << "Error in updating refinement cells in updateBCs:line 132\n";
            }
        }
    }
    ierr = MPI_Waitall(count, recv_request, status);

    if (count > 0) {
        if (vectorSize[0] != 0 || vectorSizeSend[0] > 0) {
            for (phi=0; phi < 6; ++phi) {
                for (n=0; n < cells.yPoints; ++n) {
                    switch (phi) {
                        case 0: cells.rhoBCsLR[n][0] = recvLeft[n];
                                    break;
                        case 1: cells.MxBCsLR[n][0] = recvLeft[n + cells.yPoints];
                                    break;
                        case 2: cells.MyBCsLR[n][0] = recvLeft[n + phi*cells.yPoints];
                                    break;
                        case 3: cells.rhoBCsLR[n][2] = recvLeft[n + phi*cells.yPoints];
                                    break;
                        case 4: cells.MxBCsLR[n][2] = recvLeft[n + phi*cells.yPoints];
                                    break;
                        case 5: cells.MyBCsLR[n][2] = recvLeft[n + phi*cells.yPoints];
                                    break;
                    }
                }
            }
        }
        if (vectorSize[1] != 0 || vectorSizeSend[1] > 0) {
            for (phi=0; phi < 6; ++phi) {
                for (n=0; n < cells.yPoints; ++n) {
                    switch (phi) {
                        case 0: cells.rhoBCsLR[n][1] = recvRight[n];
                                    break;
                        case 1: cells.MxBCsLR[n][1] = recvRight[n + cells.yPoints];
                                    break;
                        case 2: cells.MyBCsLR[n][1] = recvRight[n + phi*cells.yPoints];
                                    break;
                        case 3: cells.rhoBCsLR[n][3] = recvRight[n + phi*cells.yPoints];
                                    break;
                        case 4: cells.MxBCsLR[n][3] = recvRight[n + phi*cells.yPoints];
                                    break;
                        case 5: cells.MyBCsLR[n][3] = recvRight[n + phi*cells.yPoints];
                                    break;
                    }
                }
            }
        }
    }
}