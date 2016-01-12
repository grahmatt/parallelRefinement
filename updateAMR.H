void updateAMR
(
    grid& refineCells,
    grid& coarseCells,
    const int refineFactor,
    vector<int>& refineColumns,
    const double rhoShock,
    const double rhoStart,
    int rank,
    int size,
    double currentTime,
    double sendLeft[],
    double recvLeft[],
    double sendRight[],
    double recvRight[],
    int vectorSize[],
    int vectorSizeSend[]
) 
{
    int a, b, i, j, k, n, m, rowHalf, ierr, count;
    bool updateRefCol, existingCol;
    vector<int> oldRefCol;
    MPI_Status  status[2];
    MPI_Request send_request[2], recv_request[2];
    vector<int> refineColumnsSendLeft, refineColumnsSendRight;
    int refineColumnsRecvLeft[coarseCells.yPoints], refineColumnsRecvRight[coarseCells.yPoints];

    i = 0;
    n = -1;
    rowHalf = round(coarseCells.yPoints/2);
    updateRefCol = false;
    oldRefCol = refineColumns;
    refineColumns.resize(0);
    for (j=0; j < coarseCells.xPoints; ++j) {
        if (coarseCells.rho[rowHalf][j] < 0.99*rhoShock && 
            coarseCells.rho[rowHalf][j] > 1.01*rhoStart) {
            if (!updateRefCol) {
                if (n + 1 == oldRefCol.size()) {
                    updateRefCol = true;
                }
                else {
                    ++n;
                    if (oldRefCol[n] < j - 5) {
                        updateRefCol = true;
                    }
                }
            }
            if (refineColumns.size() == 0) {
                for (a=5; a >= 0; a--) {
                    if (j - a >= 0) {
                        refineColumns.push_back (j - a);
                    }
                    else if (coarseCells.connect[0]) {    // For parallel
                        refineColumnsSendLeft.push_back(j - a);
                    }
                }     
            }
            else {
                refineColumns.push_back (j);
            }
            if (refineColumns.size() > 1 && 
                refineColumns.back() - 1 != refineColumns[refineColumns.size() - 2]) {
                cout << "Error: Only one refinement zone available\n";
                bool clearVector = true;
                for (a=refineColumns[0]; a < j+1; a++) {
                    if (clearVector) {
                        refineColumns.clear();
                        clearVector = false; 
                    }
                    refineColumns.push_back(a);
                }
            }
        }
    }

    if (rank == 0 && currentTime < 5*coarseCells.dt && refineColumns.empty()) {
        refineColumns.push_back(0);
    }

    if (!refineColumns.empty()) {
        for (a=0; a < 6; a++) {
            if (refineColumns.back() != coarseCells.xPoints - 1 && !refineColumns.empty()) {
                refineColumns.push_back (refineColumns.back() + 1);
            }
            else if (coarseCells.connect[1]) {              // For parallel
                for (b=0; b + a < 6; b++) {
                    refineColumnsSendRight.push_back(b);
                }
                break;
            }
            else {
                break;
            }
        }

        if (refineColumns.back() == coarseCells.xPoints - 1) {
            refineCells.connect[1] = coarseCells.connect[1];
        }
        else {
            refineCells.connect[1] = true;
        }

        if (refineColumns[0] == 0) {
            refineCells.connect[0] = coarseCells.connect[0];
        }
        else {
            refineCells.connect[0] = true;
        }
    }

    if (!updateRefCol && refineColumns.size() != oldRefCol.size()) {
        updateRefCol = true;
    }

    // Send refine cells across processors
    if (rank == 0) {
        vectorSizeSend[1] = refineColumnsSendRight.size();
        ierr = MPI_Isend(&vectorSizeSend[1], 1, MPI_INT, rank+1, 0 , MPI_COMM_WORLD, &send_request[0]);
        ierr = MPI_Irecv(&vectorSize[1], 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &recv_request[0]);
        ierr = MPI_Wait(recv_request, status);

        if (vectorSizeSend[1] > 0) {
            ierr = MPI_Isend(&refineColumnsSendRight[0], refineColumnsSendRight.size(), MPI_INT, rank+1, 0, MPI_COMM_WORLD, &send_request[0]);
        }

        if (vectorSize[1] > 0) {
            // refineColumnsRecvRight.resize(vectorSize[0]);
            ierr = MPI_Irecv(&refineColumnsRecvRight[0], vectorSize[0], MPI_INT, rank+1, 0, MPI_COMM_WORLD, &recv_request[0]);
            count = 1;
        }
        else {
            count = 0;
        }
    }
    else if (rank == size - 1) {
        vectorSizeSend[0] = refineColumnsSendLeft.size();
        ierr = MPI_Isend(&vectorSizeSend[0], 1, MPI_INT, rank-1, 0 , MPI_COMM_WORLD, &send_request[0]);
        ierr = MPI_Irecv(&vectorSize[0], 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &recv_request[0]);
        ierr = MPI_Wait(recv_request, status);

        if (vectorSizeSend[0] > 0) {
            ierr = MPI_Isend(&refineColumnsSendLeft[0], refineColumnsSendLeft.size(), MPI_INT, rank-1, 0, MPI_COMM_WORLD, &send_request[0]);
        }

        if (vectorSize[0] > 0) {
            // refineColumnsRecvLeft.resize(vectorSize[0]);
            ierr = MPI_Irecv(&refineColumnsRecvLeft[0], vectorSize[0], MPI_INT, rank-1, 0, MPI_COMM_WORLD, &recv_request[0]);
            count = 1;
        }
        else {
            count = 0;
        }
    }
    else {
        vectorSizeSend[0] = refineColumnsSendLeft.size();
        ierr = MPI_Isend(&vectorSizeSend[0], 1, MPI_INT, rank-1, 0 , MPI_COMM_WORLD, &send_request[0]);
        ierr = MPI_Irecv(&vectorSize[0], 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &recv_request[0]);
        vectorSizeSend[1] = refineColumnsSendRight.size();
        ierr = MPI_Isend(&vectorSizeSend[1], 1, MPI_INT, rank+1, 0 , MPI_COMM_WORLD, &send_request[1]);
        ierr = MPI_Irecv(&vectorSize[1], 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &recv_request[1]);
        ierr = MPI_Waitall(2, recv_request, status);

        if (vectorSizeSend[0] > 0) {
            ierr = MPI_Isend(&refineColumnsSendLeft[0], refineColumnsSendLeft.size(), MPI_INT, rank-1, 0, MPI_COMM_WORLD, &send_request[0]);
        }
        if (vectorSizeSend[1] > 0) {
            ierr = MPI_Isend(&refineColumnsSendRight[0], refineColumnsSendRight.size(), MPI_INT, rank+1, 0, MPI_COMM_WORLD, &send_request[1]);
        }

        count = 0;
        if (vectorSize[0] > 0) {
            // refineColumnsRecvLeft.resize(vectorSize[0]);
            ierr = MPI_Irecv(&refineColumnsRecvLeft[0], vectorSize[0], MPI_INT, rank-1, 0, MPI_COMM_WORLD, &recv_request[0]);
            ++count;
        }
        if (vectorSize[1] > 0) {
            if (count == 0) {
                // refineColumnsRecvRight.resize(vectorSize[1]);
                ierr = MPI_Irecv(&refineColumnsRecvRight[0], vectorSize[1], MPI_INT, rank+1, 0, MPI_COMM_WORLD, &recv_request[0]);
            }
            else {
                // refineColumnsRecvRight.resize(vectorSize[1]);
                ierr = MPI_Irecv(&refineColumnsRecvRight[0], vectorSize[1], MPI_INT, rank+1, 0, MPI_COMM_WORLD, &recv_request[1]);
            }
            ++count;
        }
    }
    ierr = MPI_Waitall(count, recv_request, status);

    if (refineColumns.empty()) {
        if (vectorSize[0] > 0) {
            for (a=0; a < vectorSize[0]; ++a) {
                refineColumns.push_back(refineColumnsRecvLeft[a]);
            }
            updateRefCol = true;
        }
        else if (vectorSize[1] > 0) {
            cout << "rank: " << rank << " recvRight\n";
            // TEST THIS
            for (a=0; a < vectorSize[1]; a++) {
                refineColumns.push_back(coarseCells.xPoints + refineColumnsRecvRight[a]);
            }
            updateRefCol = true;
        }
    }


    // Grid update
    if (updateRefCol || refineColumns.size() != refineCells.xPoints/refineFactor) {
        grid refineCells_oldTime = refineCells;

        refineCells.xPoints = refineColumns.size()*refineFactor;
        refineCells.rho.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.Mx.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.My.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.rhoBCsBT.resize(4, vector<double>(refineCells.xPoints,0));
        refineCells.MxBCsBT.resize(4, vector<double>(refineCells.xPoints,0));
        refineCells.MyBCsBT.resize(4, vector<double>(refineCells.xPoints,0));
        refineCells.rhoBCsBT_PlusMinus.resize(2, vector<double>(refineCells.xPoints,0));
        refineCells.MxBCsBT_PlusMinus.resize(2, vector<double>(refineCells.xPoints,0));
        refineCells.MyBCsBT_PlusMinus.resize(2, vector<double>(refineCells.xPoints,0));
        refineCells.rhoF_Plus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.rhoF_Minus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.rhoG_Plus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.rhoG_Minus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.MxF_Plus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.MxF_Minus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.MxG_Plus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.MxG_Minus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.MyF_Plus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.MyF_Minus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.MyG_Plus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        refineCells.MyG_Minus.resize(refineCells.yPoints, vector<double>(refineCells.xPoints,0));
        for (j=0; j < refineColumns.size(); ++j) {
            existingCol = false;
            for (m = 0; m < oldRefCol.size(); ++m) {
                if (refineColumns[j] == oldRefCol[m]) {
                    existingCol = true;
                    break;
                }
            }
            if (existingCol) {
                for (n=0; n < refineFactor; ++n) {
                    for (i=0; i < refineCells.yPoints; ++i) {
                        refineCells.rho[i][j*refineFactor + n] = 
                            refineCells_oldTime.rho[i][m*refineFactor + n];
                        refineCells.Mx[i][j*refineFactor + n] = 
                            refineCells_oldTime.Mx[i][m*refineFactor + n];
                        refineCells.My[i][j*refineFactor + n] = 
                            refineCells_oldTime.My[i][m*refineFactor + n];
                    }
                }
            }
            else {
                for (n=0; n < refineFactor; ++n) {
                    for (i=0; i < coarseCells.yPoints; ++i) {
                        for (k=0; k < refineFactor; ++k) {
                            refineCells.rho[i*refineFactor + k][j*refineFactor + n] 
                                = coarseCells.rho[i][refineColumns[j]];
                            refineCells.Mx[i*refineFactor + k][j*refineFactor + n] 
                                = coarseCells.Mx[i][refineColumns[j]];
                            refineCells.My[i*refineFactor + k][j*refineFactor + n] 
                                = coarseCells.My[i][refineColumns[j]];
                        }
                    }
                }
            }
        }
    }

    if (!refineColumns.empty()) {
        updateBCs(refineCells, coarseCells, refineColumns, refineFactor,
            rank, sendLeft, recvLeft, sendRight, recvRight,vectorSize);
    }
}