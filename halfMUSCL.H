double rho_S, Mx_S, My_S, rho_Plus, rho_Minus, Mx_Plus, Mx_Minus, My_Plus, My_Minus, courant, soundSpeed;

soundSpeed = 1.0;
courant = soundSpeed*tunnel.dt/tunnel.dx;  // Should this be sounds speed or u, v, or magnitude??

// F Reconstuction and Fluxes
if (j == 0) {
    rho_S = fluxLimiter
    (
        tunnel.rhoBCsLR[i][0],
        tunnel.rho[i][j],
        tunnel.rho[i][j+1],
        tunnel.dx
    );
    Mx_S = fluxLimiter
    (
        tunnel.MxBCsLR[i][0],
        tunnel.Mx[i][j],
        tunnel.Mx[i][j+1],
        tunnel.dx
    );
    My_S = fluxLimiter
    (
        tunnel.MyBCsLR[i][0],
        tunnel.My[i][j],
        tunnel.My[i][j+1],
        tunnel.dx
    );
}
else if (j+1 == tunnel.xPoints) {
    rho_S = fluxLimiter
    (
        tunnel.rho[i][j-1],
        tunnel.rho[i][j],
        tunnel.rhoBCsLR[i][1],
        tunnel.dx
    );
    Mx_S = fluxLimiter
    (
        tunnel.Mx[i][j-1],
        tunnel.Mx[i][j],
        tunnel.MxBCsLR[i][1],
        tunnel.dx
    );
    My_S = fluxLimiter
    (
        tunnel.My[i][j-1],
        tunnel.My[i][j],
        tunnel.MyBCsLR[i][1],
        tunnel.dx
    );
}
else {
    rho_S = fluxLimiter
    (
        tunnel.rho[i][j-1],
        tunnel.rho[i][j],
        tunnel.rho[i][j+1],
        tunnel.dx
    );
    Mx_S = fluxLimiter
    (
        tunnel.Mx[i][j-1],
        tunnel.Mx[i][j],
        tunnel.Mx[i][j+1],
        tunnel.dx
    );
    My_S = fluxLimiter
    (
        tunnel.My[i][j-1],
        tunnel.My[i][j],
        tunnel.My[i][j+1],
        tunnel.dx
    );
}
rho_Plus = soundSpeed*(tunnel.rho[i][j] + 0.5*tunnel.dx*(1.0 - courant*rho_S));
rho_Minus = soundSpeed*(tunnel.rho[i][j] - 0.5*tunnel.dx*(1.0 - courant*rho_S));
Mx_Plus = soundSpeed*(tunnel.Mx[i][j] + 0.5*tunnel.dx*(1.0 - courant*Mx_S));
Mx_Minus = soundSpeed*(tunnel.Mx[i][j] - 0.5*tunnel.dx*(1.0 - courant*Mx_S));
My_Plus = soundSpeed*(tunnel.My[i][j] + 0.5*tunnel.dx*(1.0 - courant*My_S));
My_Minus = soundSpeed*(tunnel.My[i][j] - 0.5*tunnel.dx*(1.0 - courant*My_S));

tunnel.rhoF_Plus[i][j] = Mx_Plus;
tunnel.rhoF_Minus[i][j] = Mx_Minus;
tunnel.MxF_Plus[i][j] = pow(Mx_Plus,2.0)/rho_Plus + pow(soundSpeed, 2.0)*rho_Plus;
tunnel.MxF_Minus[i][j] = pow(Mx_Minus,2.0)/rho_Minus + pow(soundSpeed, 2.0)*rho_Minus;
tunnel.MyF_Plus[i][j] = Mx_Plus*My_Plus/rho_Plus;
tunnel.MyF_Minus[i][j] = Mx_Minus*My_Minus/rho_Minus;


// G Reconstruction and Fluxes
if (i == 0) {
    rho_S = fluxLimiter
    (
        tunnel.rhoBCsBT[0][j],
        tunnel.rho[i][j],
        tunnel.rho[i+1][j],
        tunnel.dy
    );
    Mx_S = fluxLimiter
    (
        tunnel.MxBCsBT[0][j],
        tunnel.Mx[i][j],
        tunnel.Mx[i+1][j],
        tunnel.dy
    );
    My_S = fluxLimiter
    (
        tunnel.MyBCsBT[0][j],
        tunnel.My[i][j],
        tunnel.My[i+1][j],
        tunnel.dy
    );
}
else if (i+1 == tunnel.yPoints) {
    rho_S = fluxLimiter
    (
        tunnel.rho[i-1][j],
        tunnel.rho[i][j],
        tunnel.rhoBCsBT[1][j],
        tunnel.dy
    );
    Mx_S = fluxLimiter
    (
        tunnel.Mx[i-1][j],
        tunnel.Mx[i][j],
        tunnel.MxBCsBT[1][j],
        tunnel.dy
    );
    My_S = fluxLimiter
    (
        tunnel.My[i-1][j],
        tunnel.My[i][j],
        tunnel.MyBCsBT[1][j],
        tunnel.dy
    );
}
else {
    rho_S = fluxLimiter
    (
        tunnel.rho[i-1][j],
        tunnel.rho[i][j],
        tunnel.rho[i+1][j],
        tunnel.dy
    );
    Mx_S = fluxLimiter
    (
        tunnel.Mx[i-1][j],
        tunnel.Mx[i][j],
        tunnel.Mx[i+1][j],
        tunnel.dy
    );
    My_S = fluxLimiter
    (
        tunnel.My[i-1][j],
        tunnel.My[i][j],
        tunnel.My[i+1][j],
        tunnel.dy
    );
}
rho_Plus = soundSpeed*(tunnel.rho[i][j] + 0.5*tunnel.dx*(1.0 - courant*rho_S));
rho_Minus = soundSpeed*(tunnel.rho[i][j] - 0.5*tunnel.dx*(1.0 - courant*rho_S));
Mx_Plus = soundSpeed*(tunnel.Mx[i][j] + 0.5*tunnel.dx*(1.0 - courant*Mx_S));
Mx_Minus = soundSpeed*(tunnel.Mx[i][j] - 0.5*tunnel.dx*(1.0 - courant*Mx_S));
My_Plus = soundSpeed*(tunnel.My[i][j] + 0.5*tunnel.dx*(1.0 - courant*My_S));
My_Minus = soundSpeed*(tunnel.My[i][j] - 0.5*tunnel.dx*(1.0 - courant*My_S));

tunnel.rhoG_Plus[i][j] = My_Plus;
tunnel.rhoG_Minus[i][j] = My_Minus;
tunnel.MxG_Plus[i][j] = Mx_Plus*My_Plus/rho_Plus;
tunnel.MxG_Minus[i][j] = Mx_Minus*My_Minus/rho_Minus;
tunnel.MyG_Plus[i][j] = pow(My_Plus,2.0)/rho_Plus + pow(soundSpeed, 2.0)*rho_Plus;
tunnel.MyG_Minus[i][j] = pow(My_Minus,2.0)/rho_Minus + pow(soundSpeed, 2.0)*rho_Minus;

// Half time step
tunnel_halfTime.rho[i][j] = tunnel.row[i][j] - 0.5*tunnel.dt*( \
    (tunnel.rhoF_Plus - tunnel.rhoF_Minus)/tunnel.dx + \
    (tunnel.rhoG_Plus - tunnel.rhoG_Minus)/tunnel.dy);
tunnel_halfTime.Mx[i][j] = ;
tunnel_halfTime.My[i][j] = ;
