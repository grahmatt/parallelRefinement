void updateAMR
(
    grid& refineCells,
    grid& coarseCells,
    const int refineFactor,
    vector<int>& refineColumns,
    const double rhoShock,
    const double rhoStart
) 
{
    int i, j, k, n, m;
    bool updateRefCol, existingCol;
    vector<int> oldRefCol;

    i = 0;
    n = -1;
    updateRefCol = false;
    oldRefCol = refineColumns;
    refineColumns.resize(0);
    for (j=0; j < coarseCells.xPoints; ++j) {
        if (coarseCells.rho[i][j] < 0.99*rhoShock && 
            coarseCells.rho[i][j] > 1.01*rhoStart) {
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
            if (refineColumns.size() == 0 && j > 0) {
                refineColumns.push_back (j - 1);
                if (j - 1 == 0) {
                    refineCells.connect[0] = coarseCells.connect[0];
                }
                else {
                    refineCells.connect[0] = true;
                }
            }
            refineColumns.push_back (j);
            if (refineColumns.size() > 1 && 
                refineColumns.back() - 1 != refineColumns[refineColumns.size() - 2]) {
                cout << "Error: Only one refinement zone available\n";
            }
        }
    }
    if (refineColumns.back() != coarseCells.xPoints - 1) {
        refineColumns.push_back (refineColumns.back() + 1);
        if (refineColumns.back() == coarseCells.xPoints - 1) {
            refineCells.connect[1] = coarseCells.connect[1];
        }
        else {
            refineCells.connect[1] = true;
        }
    }
    else {
        refineCells.connect[1] = true;
    }
    if (!updateRefCol && refineColumns.size() != oldRefCol.size()) {
        updateRefCol = true;
    }
    if (updateRefCol) {
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
    updateBCs(refineCells,coarseCells,refineColumns,refineFactor);
}
