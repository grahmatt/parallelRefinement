void reconstructFlux
(
    int i,
    int j,
    grid& cell,
    double courant,
    double soundSpeed
) 
{
    double rho_S, Mx_S, My_S, rho_Plus, rho_Minus, Mx_Plus, Mx_Minus, My_Plus, 
        My_Minus;
    
    // F - Reconstuction and Fluxes
    if (j == 0) {
        rho_S = fluxLimiter
        (
            cell.rhoBCsLR[i][0],
            cell.rho[i][j],
            cell.rho[i][j+1],
            cell.dx
        );
        Mx_S = fluxLimiter
        (
            cell.MxBCsLR[i][0],
            cell.Mx[i][j],
            cell.Mx[i][j+1],
            cell.dx
        );
        My_S = fluxLimiter
        (
            cell.MyBCsLR[i][0],
            cell.My[i][j],
            cell.My[i][j+1],
            cell.dx
        );
    }
    else if (j+1 == cell.xPoints) {
        rho_S = fluxLimiter
        (
            cell.rho[i][j-1],
            cell.rho[i][j],
            cell.rhoBCsLR[i][1],
            cell.dx
        );
        Mx_S = fluxLimiter
        (
            cell.Mx[i][j-1],
            cell.Mx[i][j],
            cell.MxBCsLR[i][1],
            cell.dx
        );
        My_S = fluxLimiter
        (
            cell.My[i][j-1],
            cell.My[i][j],
            cell.MyBCsLR[i][1],
            cell.dx
        );
    }
    else {
        rho_S = fluxLimiter
        (
            cell.rho[i][j-1],
            cell.rho[i][j],
            cell.rho[i][j+1],
            cell.dx
        );
        Mx_S = fluxLimiter
        (
            cell.Mx[i][j-1],
            cell.Mx[i][j],
            cell.Mx[i][j+1],
            cell.dx
        );
        My_S = fluxLimiter
        (
            cell.My[i][j-1],
            cell.My[i][j],
            cell.My[i][j+1],
            cell.dx
        );
    }
    rho_Plus = soundSpeed*(cell.rho[i][j] + 0.5*cell.dx*(1.0 - courant*rho_S));
    rho_Minus = soundSpeed*(cell.rho[i][j] - 0.5*cell.dx*(1.0 - courant*rho_S));
    Mx_Plus = soundSpeed*(cell.Mx[i][j] + 0.5*cell.dx*(1.0 - courant*Mx_S));
    Mx_Minus = soundSpeed*(cell.Mx[i][j] - 0.5*cell.dx*(1.0 - courant*Mx_S));
    My_Plus = soundSpeed*(cell.My[i][j] + 0.5*cell.dx*(1.0 - courant*My_S));
    My_Minus = soundSpeed*(cell.My[i][j] - 0.5*cell.dx*(1.0 - courant*My_S));

    cell.rhoF_Plus[i][j] = Mx_Plus;
    cell.rhoF_Minus[i][j] = Mx_Minus;
    cell.MxF_Plus[i][j] = pow(Mx_Plus,2.0)/rho_Plus + 
        pow(soundSpeed, 2.0)*rho_Plus;
    cell.MxF_Minus[i][j] = pow(Mx_Minus,2.0)/rho_Minus + 
        pow(soundSpeed, 2.0)*rho_Minus;
    cell.MyF_Plus[i][j] = Mx_Plus*My_Plus/rho_Plus;
    cell.MyF_Minus[i][j] = Mx_Minus*My_Minus/rho_Minus;


    // G - Reconstruction and Fluxes
    if (i == 0) {
        rho_S = fluxLimiter
        (
            cell.rhoBCsBT[0][j],
            cell.rho[i][j],
            cell.rho[i+1][j],
            cell.dy
        );
        Mx_S = fluxLimiter
        (
            cell.MxBCsBT[0][j],
            cell.Mx[i][j],
            cell.Mx[i+1][j],
            cell.dy
        );
        My_S = fluxLimiter
        (
            cell.MyBCsBT[0][j],
            cell.My[i][j],
            cell.My[i+1][j],
            cell.dy
        );
    }
    else if (i+1 == cell.yPoints) {
        rho_S = fluxLimiter
        (
            cell.rho[i-1][j],
            cell.rho[i][j],
            cell.rhoBCsBT[1][j],
            cell.dy
        );
        Mx_S = fluxLimiter
        (
            cell.Mx[i-1][j],
            cell.Mx[i][j],
            cell.MxBCsBT[1][j],
            cell.dy
        );
        My_S = fluxLimiter
        (
            cell.My[i-1][j],
            cell.My[i][j],
            cell.MyBCsBT[1][j],
            cell.dy
        );
    }
    else {
        rho_S = fluxLimiter
        (
            cell.rho[i-1][j],
            cell.rho[i][j],
            cell.rho[i+1][j],
            cell.dy
        );
        Mx_S = fluxLimiter
        (
            cell.Mx[i-1][j],
            cell.Mx[i][j],
            cell.Mx[i+1][j],
            cell.dy
        );
        My_S = fluxLimiter
        (
            cell.My[i-1][j],
            cell.My[i][j],
            cell.My[i+1][j],
            cell.dy
        );
    }
    rho_Plus = soundSpeed*(cell.rho[i][j] + 0.5*cell.dx*(1.0 - courant*rho_S));
    rho_Minus = soundSpeed*(cell.rho[i][j] - 0.5*cell.dx*(1.0 - courant*rho_S));
    Mx_Plus = soundSpeed*(cell.Mx[i][j] + 0.5*cell.dx*(1.0 - courant*Mx_S));
    Mx_Minus = soundSpeed*(cell.Mx[i][j] - 0.5*cell.dx*(1.0 - courant*Mx_S));
    My_Plus = soundSpeed*(cell.My[i][j] + 0.5*cell.dx*(1.0 - courant*My_S));
    My_Minus = soundSpeed*(cell.My[i][j] - 0.5*cell.dx*(1.0 - courant*My_S));

    cell.rhoG_Plus[i][j] = My_Plus;
    cell.rhoG_Minus[i][j] = My_Minus;
    cell.MxG_Plus[i][j] = Mx_Plus*My_Plus/rho_Plus;
    cell.MxG_Minus[i][j] = Mx_Minus*My_Minus/rho_Minus;
    cell.MyG_Plus[i][j] = pow(My_Plus,2.0)/rho_Plus + 
        pow(soundSpeed, 2.0)*rho_Plus;
    cell.MyG_Minus[i][j] = pow(My_Minus,2.0)/rho_Minus + 
        pow(soundSpeed, 2.0)*rho_Minus;
}
void reconstructFlux
(
    int i,
    int j,
    grid& cell,
    double courant,
    double soundSpeed,
    bool roe
) 
{
    if (!roe) {
        cout << "roe = FALSE, no other option solving using roe solver\n";
    }
    double rho_S, Mx_S, My_S, rho_Plus, rho_Minus, Mx_Plus, Mx_Minus, My_Plus, 
        My_Minus;
    
    // F* - Reconstuction and Fluxes using Roe Solver
    if (j == 0) {
        rho_S = fluxLimiter
        (
            cell.rhoBCsLR[i][0],
            cell.rho[i][j],
            cell.rho[i][j+1],
            cell.dx
        );
        Mx_S = fluxLimiter
        (
            cell.MxBCsLR[i][0],
            cell.Mx[i][j],
            cell.Mx[i][j+1],
            cell.dx
        );
        My_S = fluxLimiter
        (
            cell.MyBCsLR[i][0],
            cell.My[i][j],
            cell.My[i][j+1],
            cell.dx
        );
    }
    else if (j+1 == cell.xPoints) {
        rho_S = fluxLimiter
        (
            cell.rho[i][j-1],
            cell.rho[i][j],
            cell.rhoBCsLR[i][1],
            cell.dx
        );
        Mx_S = fluxLimiter
        (
            cell.Mx[i][j-1],
            cell.Mx[i][j],
            cell.MxBCsLR[i][1],
            cell.dx
        );
        My_S = fluxLimiter
        (
            cell.My[i][j-1],
            cell.My[i][j],
            cell.MyBCsLR[i][1],
            cell.dx
        );
    }
    else {
        rho_S = fluxLimiter
        (
            cell.rho[i][j-1],
            cell.rho[i][j],
            cell.rho[i][j+1],
            cell.dx
        );
        Mx_S = fluxLimiter
        (
            cell.Mx[i][j-1],
            cell.Mx[i][j],
            cell.Mx[i][j+1],
            cell.dx
        );
        My_S = fluxLimiter
        (
            cell.My[i][j-1],
            cell.My[i][j],
            cell.My[i][j+1],
            cell.dx
        );
    }
    rho_Plus = soundSpeed*(cell.rho[i][j] + 0.5*cell.dx*(1.0 - courant*rho_S));
    rho_Minus = soundSpeed*(cell.rho[i][j] - 0.5*cell.dx*(1.0 - courant*rho_S));
    Mx_Plus = soundSpeed*(cell.Mx[i][j] + 0.5*cell.dx*(1.0 - courant*Mx_S));
    Mx_Minus = soundSpeed*(cell.Mx[i][j] - 0.5*cell.dx*(1.0 - courant*Mx_S));
    My_Plus = soundSpeed*(cell.My[i][j] + 0.5*cell.dx*(1.0 - courant*My_S));
    My_Minus = soundSpeed*(cell.My[i][j] - 0.5*cell.dx*(1.0 - courant*My_S));

    cell.rhoF_Plus[i][j] = Mx_Plus;
    cell.rhoF_Minus[i][j] = Mx_Minus;
    cell.MxF_Plus[i][j] = pow(Mx_Plus,2.0)/rho_Plus + 
        pow(soundSpeed, 2.0)*rho_Plus;
    cell.MxF_Minus[i][j] = pow(Mx_Minus,2.0)/rho_Minus + 
        pow(soundSpeed, 2.0)*rho_Minus;
    cell.MyF_Plus[i][j] = Mx_Plus*My_Plus/rho_Plus;
    cell.MyF_Minus[i][j] = Mx_Minus*My_Minus/rho_Minus;


    // G* - Reconstruction and Fluxes using Roe Solver
    if (i == 0) {
        rho_S = fluxLimiter
        (
            cell.rhoBCsBT[0][j],
            cell.rho[i][j],
            cell.rho[i+1][j],
            cell.dy
        );
        Mx_S = fluxLimiter
        (
            cell.MxBCsBT[0][j],
            cell.Mx[i][j],
            cell.Mx[i+1][j],
            cell.dy
        );
        My_S = fluxLimiter
        (
            cell.MyBCsBT[0][j],
            cell.My[i][j],
            cell.My[i+1][j],
            cell.dy
        );
    }
    else if (i+1 == cell.yPoints) {
        rho_S = fluxLimiter
        (
            cell.rho[i-1][j],
            cell.rho[i][j],
            cell.rhoBCsBT[1][j],
            cell.dy
        );
        Mx_S = fluxLimiter
        (
            cell.Mx[i-1][j],
            cell.Mx[i][j],
            cell.MxBCsBT[1][j],
            cell.dy
        );
        My_S = fluxLimiter
        (
            cell.My[i-1][j],
            cell.My[i][j],
            cell.MyBCsBT[1][j],
            cell.dy
        );
    }
    else {
        rho_S = fluxLimiter
        (
            cell.rho[i-1][j],
            cell.rho[i][j],
            cell.rho[i+1][j],
            cell.dy
        );
        Mx_S = fluxLimiter
        (
            cell.Mx[i-1][j],
            cell.Mx[i][j],
            cell.Mx[i+1][j],
            cell.dy
        );
        My_S = fluxLimiter
        (
            cell.My[i-1][j],
            cell.My[i][j],
            cell.My[i+1][j],
            cell.dy
        );
    }
    rho_Plus = soundSpeed*(cell.rho[i][j] + 0.5*cell.dx*(1.0 - courant*rho_S));
    rho_Minus = soundSpeed*(cell.rho[i][j] - 0.5*cell.dx*(1.0 - courant*rho_S));
    Mx_Plus = soundSpeed*(cell.Mx[i][j] + 0.5*cell.dx*(1.0 - courant*Mx_S));
    Mx_Minus = soundSpeed*(cell.Mx[i][j] - 0.5*cell.dx*(1.0 - courant*Mx_S));
    My_Plus = soundSpeed*(cell.My[i][j] + 0.5*cell.dx*(1.0 - courant*My_S));
    My_Minus = soundSpeed*(cell.My[i][j] - 0.5*cell.dx*(1.0 - courant*My_S));

    cell.rhoG_Plus[i][j] = My_Plus;
    cell.rhoG_Minus[i][j] = My_Minus;
    cell.MxG_Plus[i][j] = Mx_Plus*My_Plus/rho_Plus;
    cell.MxG_Minus[i][j] = Mx_Minus*My_Minus/rho_Minus;
    cell.MyG_Plus[i][j] = pow(My_Plus,2.0)/rho_Plus + 
        pow(soundSpeed, 2.0)*rho_Plus;
    cell.MyG_Minus[i][j] = pow(My_Minus,2.0)/rho_Minus + 
        pow(soundSpeed, 2.0)*rho_Minus;
}