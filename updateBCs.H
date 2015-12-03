void updateBCs
(
    grid& cells
) 
{
    int i, n;
    for (i=0; i < 4; ++i) {
        if (!cells.connect[i]) {
            if (i == 0) { // Left
                // for (n=0; n < cells.yPoints; ++n) {
                //     cells.rhoBCsLR[n][0] = rhoShock;
                //     cells.MxBCsLR[n][0] = 0.0;                // Should this be set to (rhoShock*cells.dx*cells.dy*soundSpeed)??
                //     cells.MyBCsLR[n][0] = 0.0;
                //     // For Roe Solver
                //     cells.rhoBCsLR[n][2] = rhoShock;
                //     cells.MxBCsLR[n][2] = 0.0;                // Should this be set to (rhoShock*cells.dx*cells.dy*soundSpeed)??
                //     cells.MyBCsLR[n][2] = 0.0;
                // }
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
