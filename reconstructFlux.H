void reconstructFlux
(
    int i,
    int j,
    grid& cellHalf,
    double soundSpeed
) 
{
    double rho_S, Mx_S, My_S, rho_Plus, rho_Minus, Mx_Plus, Mx_Minus, My_Plus, 
        My_Minus, courant;
    
    // F - Reconstuction and Fluxes
    if (j == 0) {
        rho_S = fluxLimiter
        (
            cellHalf.rhoBCsLR[i][0],
            cellHalf.rho[i][j],
            cellHalf.rho[i][j+1],
            cellHalf.dx
        );
        Mx_S = fluxLimiter
        (
            cellHalf.MxBCsLR[i][0],
            cellHalf.Mx[i][j],
            cellHalf.Mx[i][j+1],
            cellHalf.dx
        );
        My_S = fluxLimiter
        (
            cellHalf.MyBCsLR[i][0],
            cellHalf.My[i][j],
            cellHalf.My[i][j+1],
            cellHalf.dx
        );
    }
    else if (j+1 == cellHalf.xPoints) {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i][j-1],
            cellHalf.rho[i][j],
            cellHalf.rhoBCsLR[i][1],
            cellHalf.dx
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i][j-1],
            cellHalf.Mx[i][j],
            cellHalf.MxBCsLR[i][1],
            cellHalf.dx
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i][j-1],
            cellHalf.My[i][j],
            cellHalf.MyBCsLR[i][1],
            cellHalf.dx
        );
    }
    else {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i][j-1],
            cellHalf.rho[i][j],
            cellHalf.rho[i][j+1],
            cellHalf.dx
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i][j-1],
            cellHalf.Mx[i][j],
            cellHalf.Mx[i][j+1],
            cellHalf.dx
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i][j-1],
            cellHalf.My[i][j],
            cellHalf.My[i][j+1],
            cellHalf.dx
        );
    }
    courant = 0.5*cellHalf.Mx[i][j]/cellHalf.rho[i][j]*cellHalf.dt/cellHalf.dx;  // Should this be sounds speed or u, v, or magnitude?? tunnel.dt/2
    rho_Plus = soundSpeed*(cellHalf.rho[i][j] + 0.5*cellHalf.dx*(1.0 - courant*rho_S));
    rho_Minus = soundSpeed*(cellHalf.rho[i][j] - 0.5*cellHalf.dx*(1.0 - courant*rho_S));
    Mx_Plus = soundSpeed*(cellHalf.Mx[i][j] + 0.5*cellHalf.dx*(1.0 - courant*Mx_S));
    Mx_Minus = soundSpeed*(cellHalf.Mx[i][j] - 0.5*cellHalf.dx*(1.0 - courant*Mx_S));
    My_Plus = soundSpeed*(cellHalf.My[i][j] + 0.5*cellHalf.dx*(1.0 - courant*My_S));
    My_Minus = soundSpeed*(cellHalf.My[i][j] - 0.5*cellHalf.dx*(1.0 - courant*My_S));

    cellHalf.rhoF_Plus[i][j] = Mx_Plus;
    cellHalf.rhoF_Minus[i][j] = Mx_Minus;
    cellHalf.MxF_Plus[i][j] = pow(Mx_Plus,2.0)/rho_Plus + 
        pow(soundSpeed, 2.0)*rho_Plus;
    cellHalf.MxF_Minus[i][j] = pow(Mx_Minus,2.0)/rho_Minus + 
        pow(soundSpeed, 2.0)*rho_Minus;
    cellHalf.MyF_Plus[i][j] = Mx_Plus*My_Plus/rho_Plus;
    cellHalf.MyF_Minus[i][j] = Mx_Minus*My_Minus/rho_Minus;


    // G - Reconstruction and Fluxes
    if (i == 0) {
        rho_S = fluxLimiter
        (
            cellHalf.rhoBCsBT[0][j],
            cellHalf.rho[i][j],
            cellHalf.rho[i+1][j],
            cellHalf.dy
        );
        Mx_S = fluxLimiter
        (
            cellHalf.MxBCsBT[0][j],
            cellHalf.Mx[i][j],
            cellHalf.Mx[i+1][j],
            cellHalf.dy
        );
        My_S = fluxLimiter
        (
            cellHalf.MyBCsBT[0][j],
            cellHalf.My[i][j],
            cellHalf.My[i+1][j],
            cellHalf.dy
        );
    }
    else if (i+1 == cellHalf.yPoints) {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i-1][j],
            cellHalf.rho[i][j],
            cellHalf.rhoBCsBT[1][j],
            cellHalf.dy
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i-1][j],
            cellHalf.Mx[i][j],
            cellHalf.MxBCsBT[1][j],
            cellHalf.dy
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i-1][j],
            cellHalf.My[i][j],
            cellHalf.MyBCsBT[1][j],
            cellHalf.dy
        );
    }
    else {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i-1][j],
            cellHalf.rho[i][j],
            cellHalf.rho[i+1][j],
            cellHalf.dy
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i-1][j],
            cellHalf.Mx[i][j],
            cellHalf.Mx[i+1][j],
            cellHalf.dy
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i-1][j],
            cellHalf.My[i][j],
            cellHalf.My[i+1][j],
            cellHalf.dy
        );
    }
    courant = 0.5*cellHalf.My[i][j]/cellHalf.rho[i][j]*cellHalf.dt/cellHalf.dy;
    rho_Plus = soundSpeed*(cellHalf.rho[i][j] + 0.5*cellHalf.dy*(1.0 - courant*rho_S));
    rho_Minus = soundSpeed*(cellHalf.rho[i][j] - 0.5*cellHalf.dy*(1.0 - courant*rho_S));
    Mx_Plus = soundSpeed*(cellHalf.Mx[i][j] + 0.5*cellHalf.dy*(1.0 - courant*Mx_S));
    Mx_Minus = soundSpeed*(cellHalf.Mx[i][j] - 0.5*cellHalf.dy*(1.0 - courant*Mx_S));
    My_Plus = soundSpeed*(cellHalf.My[i][j] + 0.5*cellHalf.dy*(1.0 - courant*My_S));
    My_Minus = soundSpeed*(cellHalf.My[i][j] - 0.5*cellHalf.dy*(1.0 - courant*My_S));

    cellHalf.rhoG_Plus[i][j] = My_Plus;
    cellHalf.rhoG_Minus[i][j] = My_Minus;
    cellHalf.MxG_Plus[i][j] = Mx_Plus*My_Plus/rho_Plus;
    cellHalf.MxG_Minus[i][j] = Mx_Minus*My_Minus/rho_Minus;
    cellHalf.MyG_Plus[i][j] = pow(My_Plus,2.0)/rho_Plus + 
        pow(soundSpeed, 2.0)*rho_Plus;
    cellHalf.MyG_Minus[i][j] = pow(My_Minus,2.0)/rho_Minus + 
        pow(soundSpeed, 2.0)*rho_Minus;
}

void reconstructFlux
(
    int i,
    int j,
    grid& cellHalf,
    grid& cellFull,
    double soundSpeed
) 
{
    double rho_S, Mx_S, My_S, courant;
    
    // Reconstuction for left/right values into F half time to use in Roe Solver
    if (j == 0) {
        rho_S = fluxLimiter
        (
            cellHalf.rhoBCsLR[i][0],
            cellHalf.rho[i][j],
            cellHalf.rho[i][j+1],
            cellHalf.dx
        );
        Mx_S = fluxLimiter
        (
            cellHalf.MxBCsLR[i][0],
            cellHalf.Mx[i][j],
            cellHalf.Mx[i][j+1],
            cellHalf.dx
        );
        My_S = fluxLimiter
        (
            cellHalf.MyBCsLR[i][0],
            cellHalf.My[i][j],
            cellHalf.My[i][j+1],
            cellHalf.dx
        );
    }
    else if (j+1 == cellHalf.xPoints) {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i][j-1],
            cellHalf.rho[i][j],
            cellHalf.rhoBCsLR[i][1],
            cellHalf.dx
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i][j-1],
            cellHalf.Mx[i][j],
            cellHalf.MxBCsLR[i][1],
            cellHalf.dx
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i][j-1],
            cellHalf.My[i][j],
            cellHalf.MyBCsLR[i][1],
            cellHalf.dx
        );
    }
    else {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i][j-1],
            cellHalf.rho[i][j],
            cellHalf.rho[i][j+1],
            cellHalf.dx
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i][j-1],
            cellHalf.Mx[i][j],
            cellHalf.Mx[i][j+1],
            cellHalf.dx
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i][j-1],
            cellHalf.My[i][j],
            cellHalf.My[i][j+1],
            cellHalf.dx
        );
    }
    courant = 0.5*cellHalf.Mx[i][j]/cellHalf.rho[i][j]*cellHalf.dt/cellHalf.dx;
    cellHalf.rhoF_Plus[i][j] = soundSpeed*(cellHalf.rho[i][j] + 0.5*cellHalf.dx*(1.0 - courant*rho_S));
    cellHalf.rhoF_Minus[i][j] = soundSpeed*(cellHalf.rho[i][j] - 0.5*cellHalf.dx*(1.0 - courant*rho_S));
    cellHalf.MxF_Plus[i][j] = soundSpeed*(cellHalf.Mx[i][j] + 0.5*cellHalf.dx*(1.0 - courant*Mx_S));
    cellHalf.MxF_Minus[i][j] = soundSpeed*(cellHalf.Mx[i][j] - 0.5*cellHalf.dx*(1.0 - courant*Mx_S));
    cellHalf.MyF_Plus[i][j] = soundSpeed*(cellHalf.My[i][j] + 0.5*cellHalf.dx*(1.0 - courant*My_S));
    cellHalf.MyF_Minus[i][j] = soundSpeed*(cellHalf.My[i][j] - 0.5*cellHalf.dx*(1.0 - courant*My_S));
    // Boundary condition left/right values
    if (j == 0) {
        rho_S = fluxLimiter
        (
            cellHalf.rhoBCsLR[i][2],
            cellHalf.rhoBCsLR[i][0],
            cellHalf.rho[i][j],
            cellHalf.dx
        );
        Mx_S = fluxLimiter
        (
            cellHalf.MxBCsLR[i][2],
            cellHalf.MxBCsLR[i][0],
            cellHalf.Mx[i][j],
            cellHalf.dx
        );
        My_S = fluxLimiter
        (
            cellHalf.MyBCsLR[i][2],
            cellHalf.MyBCsLR[i][0],
            cellHalf.My[i][j],
            cellHalf.dx
        );
        courant = 0.5*cellHalf.MxBCsLR[i][0]/cellHalf.rhoBCsLR[i][0]*cellHalf.dt/cellHalf.dx;
        cellHalf.rhoBCsLR_PlusMinus[i][0] = soundSpeed*(cellHalf.rhoBCsLR[i][0] + 0.5*cellHalf.dx*(1.0 - courant*rho_S));
        cellHalf.MxBCsLR_PlusMinus[i][0] = soundSpeed*(cellHalf.MxBCsLR[i][0] + 0.5*cellHalf.dx*(1.0 - courant*Mx_S));
        cellHalf.MyBCsLR_PlusMinus[i][0] = soundSpeed*(cellHalf.MyBCsLR[i][0] + 0.5*cellHalf.dx*(1.0 - courant*My_S));
    }
    else if (j+1 == cellHalf.xPoints) {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i][j],
            cellHalf.rhoBCsLR[i][1],
            cellHalf.rhoBCsLR[i][3],
            cellHalf.dx
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i][j],
            cellHalf.MxBCsLR[i][1],
            cellHalf.MxBCsLR[i][3],
            cellHalf.dx
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i][j],
            cellHalf.MyBCsLR[i][1],
            cellHalf.MyBCsLR[i][3],
            cellHalf.dx
        );
        courant = 0.5*cellHalf.MxBCsLR[i][1]/cellHalf.rhoBCsLR[i][1]*cellHalf.dt/cellHalf.dx;
        cellHalf.rhoBCsLR_PlusMinus[i][1] = soundSpeed*(cellHalf.rhoBCsLR[i][1] - 0.5*cellHalf.dx*(1.0 - courant*rho_S));
        cellHalf.MxBCsLR_PlusMinus[i][1] = soundSpeed*(cellHalf.MxBCsLR[i][1] - 0.5*cellHalf.dx*(1.0 - courant*Mx_S));
        cellHalf.MyBCsLR_PlusMinus[i][1] = soundSpeed*(cellHalf.MyBCsLR[i][1] - 0.5*cellHalf.dx*(1.0 - courant*My_S));
    }

    // Reconstuction for bottom/top values into G half time to use in Roe Solver
    if (i == 0) {
        rho_S = fluxLimiter
        (
            cellHalf.rhoBCsBT[0][j],
            cellHalf.rho[i][j],
            cellHalf.rho[i+1][j],
            cellHalf.dy
        );
        Mx_S = fluxLimiter
        (
            cellHalf.MxBCsBT[0][j],
            cellHalf.Mx[i][j],
            cellHalf.Mx[i+1][j],
            cellHalf.dy
        );
        My_S = fluxLimiter
        (
            cellHalf.MyBCsBT[0][j],
            cellHalf.My[i][j],
            cellHalf.My[i+1][j],
            cellHalf.dy
        );
    }
    else if (i+1 == cellHalf.yPoints) {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i-1][j],
            cellHalf.rho[i][j],
            cellHalf.rhoBCsBT[1][j],
            cellHalf.dy
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i-1][j],
            cellHalf.Mx[i][j],
            cellHalf.MxBCsBT[1][j],
            cellHalf.dy
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i-1][j],
            cellHalf.My[i][j],
            cellHalf.MyBCsBT[1][j],
            cellHalf.dy
        );
    }
    else {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i-1][j],
            cellHalf.rho[i][j],
            cellHalf.rho[i+1][j],
            cellHalf.dy
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i-1][j],
            cellHalf.Mx[i][j],
            cellHalf.Mx[i+1][j],
            cellHalf.dy
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i-1][j],
            cellHalf.My[i][j],
            cellHalf.My[i+1][j],
            cellHalf.dy
        );
    }
    courant = 0.5*cellHalf.My[i][j]/cellHalf.rho[i][j]*cellHalf.dt/cellHalf.dy;
    cellHalf.rhoG_Plus[i][j] = soundSpeed*(cellHalf.rho[i][j] + 0.5*cellHalf.dy*(1.0 - courant*rho_S));
    cellHalf.rhoG_Minus[i][j] = soundSpeed*(cellHalf.rho[i][j] - 0.5*cellHalf.dy*(1.0 - courant*rho_S));
    cellHalf.MxG_Plus[i][j] = soundSpeed*(cellHalf.Mx[i][j] + 0.5*cellHalf.dy*(1.0 - courant*Mx_S));
    cellHalf.MxG_Minus[i][j] = soundSpeed*(cellHalf.Mx[i][j] - 0.5*cellHalf.dy*(1.0 - courant*Mx_S));
    cellHalf.MyG_Plus[i][j] = soundSpeed*(cellHalf.My[i][j] + 0.5*cellHalf.dy*(1.0 - courant*My_S));
    cellHalf.MyG_Minus[i][j] = soundSpeed*(cellHalf.My[i][j] - 0.5*cellHalf.dy*(1.0 - courant*My_S));
    // Boundary condition left/right values
    if (i == 0) {
        rho_S = fluxLimiter
        (
            cellHalf.rhoBCsBT[2][j],
            cellHalf.rhoBCsBT[0][j],
            cellHalf.rho[i][j],
            cellHalf.dy
        );
        Mx_S = fluxLimiter
        (
            cellHalf.MxBCsBT[2][j],
            cellHalf.MxBCsBT[0][j],
            cellHalf.Mx[i][j],
            cellHalf.dy
        );
        My_S = fluxLimiter
        (
            cellHalf.MyBCsBT[2][j],
            cellHalf.MyBCsBT[0][j],
            cellHalf.My[i][j],
            cellHalf.dy
        );
        courant = 0.5*cellHalf.MyBCsBT[0][j]/cellHalf.rhoBCsBT[0][j]*cellHalf.dt/cellHalf.dy;
        cellHalf.rhoBCsBT_PlusMinus[0][j] = soundSpeed*(cellHalf.rhoBCsBT[0][j] + 0.5*cellHalf.dy*(1.0 - courant*rho_S));
        cellHalf.MxBCsBT_PlusMinus[0][j] = soundSpeed*(cellHalf.MxBCsBT[0][j] + 0.5*cellHalf.dy*(1.0 - courant*Mx_S));
        cellHalf.MyBCsBT_PlusMinus[0][j] = soundSpeed*(cellHalf.MyBCsBT[0][j] + 0.5*cellHalf.dy*(1.0 - courant*My_S));
    }
    else if (i+1 == cellHalf.yPoints) {
        rho_S = fluxLimiter
        (
            cellHalf.rho[i][j],
            cellHalf.rhoBCsBT[1][j],
            cellHalf.rhoBCsBT[3][j],
            cellHalf.dy
        );
        Mx_S = fluxLimiter
        (
            cellHalf.Mx[i][j],
            cellHalf.MxBCsBT[1][j],
            cellHalf.MxBCsBT[3][j],
            cellHalf.dy
        );
        My_S = fluxLimiter
        (
            cellHalf.My[i][j],
            cellHalf.MyBCsBT[1][j],
            cellHalf.MyBCsBT[3][j],
            cellHalf.dy
        );
        courant = 0.5*cellHalf.MyBCsBT[1][j]/cellHalf.rhoBCsBT[1][j]*cellHalf.dt/cellHalf.dy;
        cellHalf.rhoBCsBT_PlusMinus[1][j] = soundSpeed*(cellHalf.rhoBCsBT[1][j] - 0.5*cellHalf.dy*(1.0 - courant*rho_S));
        cellHalf.MxBCsBT_PlusMinus[1][j] = soundSpeed*(cellHalf.MxBCsBT[1][j] - 0.5*cellHalf.dy*(1.0 - courant*Mx_S));
        cellHalf.MyBCsBT_PlusMinus[1][j] = soundSpeed*(cellHalf.MyBCsBT[1][j] - 0.5*cellHalf.dy*(1.0 - courant*My_S));
    }
}
