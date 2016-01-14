void printGrid
(
    grid& cells
)
{
    cout.precision(4);
    for(int i=cells.yPoints-1; i >= 0; i--){
        for(int j=0; j < cells.xPoints; j++){
            cout << cells.rho[i][j] << "  ";
        }
    cout << "\n";
    }
}


void printGridFiles
(
    grid& cells
)
{
    ofstream rhoFile, xFile, yFile;
    rhoFile.precision(9);
    xFile.precision(9);
    yFile.precision(9);
    rhoFile.open ("RHO.txt");
    xFile.open ("X.txt");
    yFile.open ("Y.txt");
    xFile << "NumberOfColumns: " << cells.xPoints;
    xFile << " NumberOfRows: " << cells.yPoints << "\n";
    for(int i=0; i < cells.yPoints; i++) {
        for(int j=0; j < cells.xPoints; j++) {
            rhoFile << cells.rho[i][j] << " ";
            xFile << cells.dx*j + 0.5*cells.dx << " ";
            yFile << cells.dy*i + 0.5*cells.dy << " ";
        }
        if (i + 1 < cells.yPoints) {
            rhoFile << "\n";
            xFile << "\n";
            yFile << "\n";
        }
    }
    rhoFile.close();
    xFile.close();
    yFile.close();
}


void printGridFiles
(
    grid& coarseCells,
    grid& refinedCells,
    vector<int> refineColumns
)
{
    ofstream rhoFile, xFile, yFile, rhoFileRefine, xFileRefine, yFileRefine;
    // Coarse cells
    rhoFile.precision(9);
    xFile.precision(9);
    yFile.precision(9);
    rhoFile.open ("RHO.txt");
    xFile.open ("X.txt");
    yFile.open ("Y.txt");
    xFile << "NumberOfColumns: " << coarseCells.xPoints;
    xFile << " NumberOfRows: " << coarseCells.yPoints << "\n";
    for(int i=0; i < coarseCells.yPoints; i++) {
        for(int j=0; j < coarseCells.xPoints; j++) {
            rhoFile << coarseCells.rho[i][j] << " ";
            xFile << coarseCells.dx*j + 0.5*coarseCells.dx << " ";
            yFile << coarseCells.dy*i + 0.5*coarseCells.dy << " ";
        }
        if (i + 1 < coarseCells.yPoints) {
            rhoFile << "\n";
            xFile << "\n";
            yFile << "\n";
        }
    }
    rhoFile.close();
    xFile.close();
    yFile.close();
    // Refined cells
    rhoFileRefine.precision(9);
    xFileRefine.precision(9);
    yFileRefine.precision(9);
    rhoFileRefine.open ("refineRHO.txt");
    xFileRefine.open ("refineX.txt");
    yFileRefine.open ("refineY.txt");
    xFileRefine << "NumberOfColumns: " << refinedCells.xPoints;
    xFileRefine << " NumberOfRows: " << refinedCells.yPoints << "\n";
    for(int i=0; i < refinedCells.yPoints; i++) {
        for(int j=0; j < refinedCells.xPoints; j++) {
            rhoFileRefine << refinedCells.rho[i][j] << " ";
            xFileRefine << refineColumns[0]*coarseCells.dx + refinedCells.dx*j + 0.5*refinedCells.dx << " ";
            yFileRefine << refinedCells.dy*i + 0.5*refinedCells.dy << " ";
        }
        if (i + 1 < refinedCells.yPoints) {
            rhoFileRefine << "\n";
            xFileRefine << "\n";
            yFileRefine << "\n";
        }
    }
    rhoFileRefine.close();
    xFileRefine.close();
    yFileRefine.close();
}


void printGridFiles
(
    grid& coarseCells,
    grid& refinedCells,
    vector<int> refineColumns,
    int rank,
    int size
)
{
    int xPointsGlobal, ierr, xPointsProc, xPointsCurrent;
    ofstream rhoFile, xFile, yFile, rhoFileRefine, xFileRefine, yFileRefine;
    MPI_Status  status[1];
    double recv_row[coarseCells.xPoints];

    // Coarse cells
    ierr = MPI_Reduce(&coarseCells.xPoints, &xPointsGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        rhoFile.precision(9);
        xFile.precision(9);
        yFile.precision(9);
        rhoFile.open ("RHO.txt");
        xFile.open ("X.txt");
        yFile.open ("Y.txt");
        xFile << "NumberOfColumns: " << xPointsGlobal;
        xFile << " NumberOfRows: " << coarseCells.yPoints << "\n";
    }
    for (int i=0; i < coarseCells.yPoints; i++) {
        xPointsCurrent = 0;
        if (rank == 0) {
            for (int p=0; p < size; ++p) {
                if (p == 0) {
                    for (int j=0; j < coarseCells.xPoints; j++) { 
                        rhoFile << coarseCells.rho[i][j] << " ";
                        xFile << coarseCells.dx*j + 0.5*coarseCells.dx << " ";
                        yFile << coarseCells.dy*i + 0.5*coarseCells.dy << " ";
                    }
                    xPointsCurrent += coarseCells.xPoints;
                }
                else {
                    ierr = MPI_Recv(&xPointsProc, 1, MPI_INT, p, 0, MPI_COMM_WORLD, status);
                    ierr = MPI_Recv(&recv_row, xPointsProc, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, status);
                    for (int j=0; j < xPointsProc; j++) {
                        rhoFile << recv_row[j] << " ";
                        xFile << coarseCells.dx*j + 0.5*coarseCells.dx + xPointsCurrent*coarseCells.dx << " ";
                        yFile << coarseCells.dy*i + 0.5*coarseCells.dy << " ";
                    }
                    xPointsCurrent += xPointsProc;
                }
            }
            if (i + 1 < coarseCells.yPoints) {
                rhoFile << "\n";
                xFile << "\n";
                yFile << "\n";
            }
        }
        else {
            ierr = MPI_Ssend(&coarseCells.xPoints, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            ierr = MPI_Ssend(&coarseCells.rho[i][0], coarseCells.xPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    if (rank == 0) {
        rhoFile.close();
        xFile.close();
        yFile.close();
    }

    // Refined cells
    int procRefCellsSize[size], refineColStart[size], maxFineCells;
    ierr = MPI_Reduce(&refinedCells.xPoints, &maxFineCells, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    double recv_rowFine[maxFineCells];

    ierr = MPI_Reduce(&refinedCells.xPoints, &xPointsGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        rhoFileRefine.precision(9);
        xFileRefine.precision(9);
        yFileRefine.precision(9);
        rhoFileRefine.open ("refineRHO.txt");
        xFileRefine.open ("refineX.txt");
        yFileRefine.open ("refineY.txt");
        xFileRefine << "NumberOfColumns: " << xPointsGlobal;
        xFileRefine << " NumberOfRows: " << refinedCells.yPoints << "\n";
    }
    if (rank == 0) {
        for (int p=0; p < size; ++p) {
            if (p == 0) {
                procRefCellsSize[p] = refinedCells.xPoints;
                if (refinedCells.xPoints > 0) {
                    refineColStart[p] = refineColumns[0];
                }
                else {
                    refineColStart[p] = 0;
                }
            }
            else {
                ierr = MPI_Recv(&procRefCellsSize[p], 1, MPI_INT, p, 0, MPI_COMM_WORLD, status);
                if (procRefCellsSize[p] > 0) {
                    ierr = MPI_Recv(&refineColStart[p], 1, MPI_INT, p, 0, MPI_COMM_WORLD, status);
                }
                else {
                    refineColStart[p] = 0;
                }
            }
        }
    }
    else {
        ierr = MPI_Ssend(&refinedCells.xPoints, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        if (refinedCells.xPoints > 0) {
            ierr = MPI_Ssend(&refineColumns[0], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    for(int i=0; i < refinedCells.yPoints; i++) {
        if (rank == 0) {
             for (int p=0; p < size; ++p) {
                if (p == 0 && refinedCells.xPoints > 0) {
                    for(int j=0; j < refinedCells.xPoints; j++) {
                        rhoFileRefine << refinedCells.rho[i][j] << " ";
                        xFileRefine << refineColumns[0]*coarseCells.dx + refinedCells.dx*j + 0.5*refinedCells.dx << " ";
                        yFileRefine << refinedCells.dy*i + 0.5*refinedCells.dy << " ";
                    }
                }
                else if (p != 0 && procRefCellsSize[p] > 0) {
                    ierr = MPI_Recv(&recv_rowFine, procRefCellsSize[p], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, status);
                    for (int j=0; j < procRefCellsSize[p]; j++) {
                        rhoFileRefine << recv_rowFine[j] << " ";
                        xFileRefine << p*coarseCells.xPoints*coarseCells.dx + refineColStart[p]*coarseCells.dx + refinedCells.dx*j + 0.5*refinedCells.dx << " ";
                        yFileRefine << refinedCells.dy*i + 0.5*refinedCells.dy << " ";
                    }
                }
            }
            if (i + 1 < refinedCells.yPoints) {
                rhoFileRefine << "\n";
                xFileRefine << "\n";
                yFileRefine << "\n";
            }
        }
        else if (refinedCells.xPoints > 0) {
            ierr = MPI_Ssend(&refinedCells.rho[i][0], refinedCells.xPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    if (rank == 0) {
        rhoFileRefine.close();
        xFileRefine.close();
        yFileRefine.close();
    }
}