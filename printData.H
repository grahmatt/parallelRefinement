void printGrid(grid& cells){
    cout.precision(4);
    for(int i=cells.yPoints-1; i >= 0; i--){
        for(int j=0; j < cells.xPoints; j++){
            cout << cells.rho[i][j] << "  ";
        }
    cout << "\n";
    }
}

void printGridFiles(grid& cells) {
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
void printGridFiles(grid& coarseCells, grid& refinedCells, vector<int> refineColumns) {
    int rows, cols;
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
    // // Combined cells
    // rhoFileCombine.precision(9);
    // xFileCombine.precision(9);
    // yFileCombine.precision(9);
    // rhoFileCombine.open ("combineRHO.txt");
    // xFileCombine.open ("combineX.txt");
    // yFileCombine.open ("combineY.txt");
    // cols =  coarseCells.xPoints - refineColumns.size() + refineColumns.size()*refinedCells.xPoints;
    // xFileCombine << "NumberOfColumns: " << cols;
    // rows = refinedCells.yPoints + coarseCells.yPoints;
    // xFileCombine << " NumberOfRows: " << rows << "\n";
    // for(int i=0; i < rows; i++) {
    //     for(int j=0; j < cols; j++) {
    //         rhoFileCombine << refinedCells.rho[i][j] << " ";
    //         xFileCombine << refinedCells.dx*j + 0.5*refinedCells.dx << " ";
    //         yFileCombine << refinedCells.dy*i + 0.5*refinedCells.dy << " ";
    //     }
    //     if (i + 1 < rows) {
    //         rhoFileCombine << "\n";
    //         xFileCombine << "\n";
    //         yFileCombine << "\n";
    //     }
    // }
    // rhoFileCombine.close();
    // xFileCombine.close();
    // yFileCombine.close();
}