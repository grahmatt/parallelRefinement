void MUSCL
(
    grid& cells,
    const double soundSpeed
) 
{
    int i, j;
    grid cells_halfTime = cells;

    // MUSCL Half Time
    for (i=0; i < cells.yPoints; ++i) {
        for (j=0; j < cells.xPoints; ++j) {
            reconstructFlux(i,j,cells_halfTime,soundSpeed);
            cells_halfTime.rho[i][j] = cells_halfTime.rho[i][j] - 0.5*cells_halfTime.dt*
            (
                (cells_halfTime.rhoF_Plus[i][j] - cells_halfTime.rhoF_Minus[i][j])/cells_halfTime.dx +
                (cells_halfTime.rhoG_Plus[i][j] - cells_halfTime.rhoG_Minus[i][j])/cells_halfTime.dy
            );
            cells_halfTime.Mx[i][j] = cells_halfTime.Mx[i][j] - 0.5*cells_halfTime.dt*
            (
                (cells_halfTime.MxF_Plus[i][j] - cells_halfTime.MxF_Minus[i][j])/cells_halfTime.dx +
                (cells_halfTime.MxG_Plus[i][j] - cells_halfTime.MxG_Minus[i][j])/cells_halfTime.dy
            );
            cells_halfTime.My[i][j] = cells_halfTime.My[i][j] - 0.5*cells_halfTime.dt*
            (
                (cells_halfTime.MyF_Plus[i][j] - cells_halfTime.MyF_Minus[i][j])/cells_halfTime.dx +
                (cells_halfTime.MyG_Plus[i][j] - cells_halfTime.MyG_Minus[i][j])/cells_halfTime.dy
            );
        }
    }
    updateBCs(cells_halfTime);

    // MUSCL Full Time
    for (i=0; i < cells.yPoints; ++i) {
        for (j=0; j < cells.xPoints; ++j) {
            reconstructFlux(i,j,cells_halfTime,cells,soundSpeed);   
        }
    }
    for (i=0; i <= cells.yPoints; ++i) {
        for (j=0; j <= cells.xPoints; ++j) {
            roeSolver(i,j,cells_halfTime,cells,soundSpeed);   
        }
    }
    for (i=0; i < cells.yPoints; ++i) {
        for (j=0; j < cells.xPoints; ++j) {
            if (i+1 == cells.yPoints && j+1 == cells.xPoints) {
                cells.rho[i][j] = cells.rho[i][j] - cells.dt*
                (
                    (cells.rhoF_Plus[i][j] - cells.rhoF_Minus[i][j])/cells.dx +
                    (cells.rhoG_Plus[i][j] - cells.rhoG_Minus[i][j])/cells.dy
                );
                cells.Mx[i][j] = cells.Mx[i][j] - cells.dt*
                (
                    (cells.MxF_Plus[i][j] - cells.MxF_Minus[i][j])/cells.dx +
                    (cells.MxG_Plus[i][j] - cells.MxG_Minus[i][j])/cells.dy
                );
                cells.My[i][j] = cells.My[i][j] - cells.dt*
                (
                    (cells.MyF_Plus[i][j] - cells.MyF_Minus[i][j])/cells.dx +
                    (cells.MyG_Plus[i][j] - cells.MyG_Minus[i][j])/cells.dy
                );
            }
            else if (j+1 == cells.xPoints) {
                cells.rho[i][j] = cells.rho[i][j] - cells.dt*
                (
                    (cells.rhoF_Plus[i][j] - cells.rhoF_Minus[i][j])/cells.dx +
                    (cells.rhoG_Minus[i+1][j] - cells.rhoG_Minus[i][j])/cells.dy
                );
                cells.Mx[i][j] = cells.Mx[i][j] - cells.dt*
                (
                    (cells.MxF_Plus[i][j] - cells.MxF_Minus[i][j])/cells.dx +
                    (cells.MxG_Minus[i+1][j] - cells.MxG_Minus[i][j])/cells.dy
                );
                cells.My[i][j] = cells.My[i][j] - cells.dt*
                (
                    (cells.MyF_Plus[i][j] - cells.MyF_Minus[i][j])/cells.dx +
                    (cells.MyG_Minus[i+1][j] - cells.MyG_Minus[i][j])/cells.dy
                );
            }
            else if (i+1 == cells.yPoints) {
                cells.rho[i][j] = cells.rho[i][j] - cells.dt*
                (
                    (cells.rhoF_Minus[i][j+1] - cells.rhoF_Minus[i][j])/cells.dx +
                    (cells.rhoG_Plus[i][j] - cells.rhoG_Minus[i][j])/cells.dy
                );
                cells.Mx[i][j] = cells.Mx[i][j] - cells.dt*
                (
                    (cells.MxF_Minus[i][j+1] - cells.MxF_Minus[i][j])/cells.dx +
                    (cells.MxG_Plus[i][j] - cells.MxG_Minus[i][j])/cells.dy
                );
                cells.My[i][j] = cells.My[i][j] - cells.dt*
                (
                    (cells.MyF_Minus[i][j+1] - cells.MyF_Minus[i][j])/cells.dx +
                    (cells.MyG_Plus[i][j] - cells.MyG_Minus[i][j])/cells.dy
                );
            }
            else {
                cells.rho[i][j] = cells.rho[i][j] - cells.dt*
                (
                    (cells.rhoF_Minus[i][j+1] - cells.rhoF_Minus[i][j])/cells.dx +
                    (cells.rhoG_Minus[i+1][j] - cells.rhoG_Minus[i][j])/cells.dy
                );
                cells.Mx[i][j] = cells.Mx[i][j] - cells.dt*
                (
                    (cells.MxF_Minus[i][j+1] - cells.MxF_Minus[i][j])/cells.dx +
                    (cells.MxG_Minus[i+1][j] - cells.MxG_Minus[i][j])/cells.dy
                );
                cells.My[i][j] = cells.My[i][j] - cells.dt*
                (
                    (cells.MyF_Minus[i][j+1] - cells.MyF_Minus[i][j])/cells.dx +
                    (cells.MyG_Minus[i+1][j] - cells.MyG_Minus[i][j])/cells.dy
                );
            }
        }
    }
    updateBCs(cells);
}

void MUSCL
(
    grid& refineCells,
    const double soundSpeed,
    grid& coarseCells,
    vector<int> refineColumns,
    const int refineFactor
) 
{
    int i, j;
    grid refineCells_halfTime = refineCells;
    // MUSCL Half Time
    for (i=0; i < refineCells.yPoints; ++i) {
        for (j=0; j < refineCells.xPoints; ++j) {
            // cout << "HERE 1: i = " << i << " j = " << j << "\n";
            // cout << refineCells.xPoints << " " << refineCells.yPoints << "\n";
            reconstructFlux(i,j,refineCells_halfTime,soundSpeed);
            refineCells_halfTime.rho[i][j] = refineCells_halfTime.rho[i][j] - 0.5*refineCells_halfTime.dt*
            (
                (refineCells_halfTime.rhoF_Plus[i][j] - refineCells_halfTime.rhoF_Minus[i][j])/refineCells_halfTime.dx +
                (refineCells_halfTime.rhoG_Plus[i][j] - refineCells_halfTime.rhoG_Minus[i][j])/refineCells_halfTime.dy
            );
            refineCells_halfTime.Mx[i][j] = refineCells_halfTime.Mx[i][j] - 0.5*refineCells_halfTime.dt*
            (
                (refineCells_halfTime.MxF_Plus[i][j] - refineCells_halfTime.MxF_Minus[i][j])/refineCells_halfTime.dx +
                (refineCells_halfTime.MxG_Plus[i][j] - refineCells_halfTime.MxG_Minus[i][j])/refineCells_halfTime.dy
            );
            refineCells_halfTime.My[i][j] = refineCells_halfTime.My[i][j] - 0.5*refineCells_halfTime.dt*
            (
                (refineCells_halfTime.MyF_Plus[i][j] - refineCells_halfTime.MyF_Minus[i][j])/refineCells_halfTime.dx +
                (refineCells_halfTime.MyG_Plus[i][j] - refineCells_halfTime.MyG_Minus[i][j])/refineCells_halfTime.dy
            );
        }
    }
    updateBCs(refineCells_halfTime,coarseCells,refineColumns,refineFactor);

    // MUSCL Full Time
    for (i=0; i < refineCells.yPoints; ++i) {
        for (j=0; j < refineCells.xPoints; ++j) {
            reconstructFlux(i,j,refineCells_halfTime,refineCells,soundSpeed);   
        }
    }
    for (i=0; i <= refineCells.yPoints; ++i) {
        for (j=0; j <= refineCells.xPoints; ++j) {
            roeSolver(i,j,refineCells_halfTime,refineCells,soundSpeed);   
        }
    }
    for (i=0; i < refineCells.yPoints; ++i) {
        for (j=0; j < refineCells.xPoints; ++j) {
            if (i+1 == refineCells.yPoints && j+1 == refineCells.xPoints) {
                refineCells.rho[i][j] = refineCells.rho[i][j] - refineCells.dt*
                (
                    (refineCells.rhoF_Plus[i][j] - refineCells.rhoF_Minus[i][j])/refineCells.dx +
                    (refineCells.rhoG_Plus[i][j] - refineCells.rhoG_Minus[i][j])/refineCells.dy
                );
                refineCells.Mx[i][j] = refineCells.Mx[i][j] - refineCells.dt*
                (
                    (refineCells.MxF_Plus[i][j] - refineCells.MxF_Minus[i][j])/refineCells.dx +
                    (refineCells.MxG_Plus[i][j] - refineCells.MxG_Minus[i][j])/refineCells.dy
                );
                refineCells.My[i][j] = refineCells.My[i][j] - refineCells.dt*
                (
                    (refineCells.MyF_Plus[i][j] - refineCells.MyF_Minus[i][j])/refineCells.dx +
                    (refineCells.MyG_Plus[i][j] - refineCells.MyG_Minus[i][j])/refineCells.dy
                );
            }
            else if (j+1 == refineCells.xPoints) {
                refineCells.rho[i][j] = refineCells.rho[i][j] - refineCells.dt*
                (
                    (refineCells.rhoF_Plus[i][j] - refineCells.rhoF_Minus[i][j])/refineCells.dx +
                    (refineCells.rhoG_Minus[i+1][j] - refineCells.rhoG_Minus[i][j])/refineCells.dy
                );
                refineCells.Mx[i][j] = refineCells.Mx[i][j] - refineCells.dt*
                (
                    (refineCells.MxF_Plus[i][j] - refineCells.MxF_Minus[i][j])/refineCells.dx +
                    (refineCells.MxG_Minus[i+1][j] - refineCells.MxG_Minus[i][j])/refineCells.dy
                );
                refineCells.My[i][j] = refineCells.My[i][j] - refineCells.dt*
                (
                    (refineCells.MyF_Plus[i][j] - refineCells.MyF_Minus[i][j])/refineCells.dx +
                    (refineCells.MyG_Minus[i+1][j] - refineCells.MyG_Minus[i][j])/refineCells.dy
                );
            }
            else if (i+1 == refineCells.yPoints) {
                refineCells.rho[i][j] = refineCells.rho[i][j] - refineCells.dt*
                (
                    (refineCells.rhoF_Minus[i][j+1] - refineCells.rhoF_Minus[i][j])/refineCells.dx +
                    (refineCells.rhoG_Plus[i][j] - refineCells.rhoG_Minus[i][j])/refineCells.dy
                );
                refineCells.Mx[i][j] = refineCells.Mx[i][j] - refineCells.dt*
                (
                    (refineCells.MxF_Minus[i][j+1] - refineCells.MxF_Minus[i][j])/refineCells.dx +
                    (refineCells.MxG_Plus[i][j] - refineCells.MxG_Minus[i][j])/refineCells.dy
                );
                refineCells.My[i][j] = refineCells.My[i][j] - refineCells.dt*
                (
                    (refineCells.MyF_Minus[i][j+1] - refineCells.MyF_Minus[i][j])/refineCells.dx +
                    (refineCells.MyG_Plus[i][j] - refineCells.MyG_Minus[i][j])/refineCells.dy
                );
            }
            else {
                refineCells.rho[i][j] = refineCells.rho[i][j] - refineCells.dt*
                (
                    (refineCells.rhoF_Minus[i][j+1] - refineCells.rhoF_Minus[i][j])/refineCells.dx +
                    (refineCells.rhoG_Minus[i+1][j] - refineCells.rhoG_Minus[i][j])/refineCells.dy
                );
                refineCells.Mx[i][j] = refineCells.Mx[i][j] - refineCells.dt*
                (
                    (refineCells.MxF_Minus[i][j+1] - refineCells.MxF_Minus[i][j])/refineCells.dx +
                    (refineCells.MxG_Minus[i+1][j] - refineCells.MxG_Minus[i][j])/refineCells.dy
                );
                refineCells.My[i][j] = refineCells.My[i][j] - refineCells.dt*
                (
                    (refineCells.MyF_Minus[i][j+1] - refineCells.MyF_Minus[i][j])/refineCells.dx +
                    (refineCells.MyG_Minus[i+1][j] - refineCells.MyG_Minus[i][j])/refineCells.dy
                );
            }
        }
    }    
    updateBCs(refineCells,coarseCells,refineColumns,refineFactor);
}