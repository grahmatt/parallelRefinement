// Update BCs for tunnel

// MUSCL Half Time
grid tunnel_halfTime = tunnel;
for (i=0; i < tunnel.yPoints; ++i) {
    for (j=0; j < tunnel.xPoints; ++j) {
        reconstructFlux(i,j,tunnel_halfTime,courant,soundSpeed);
        tunnel_halfTime.rho[i][j] = tunnel_halfTime.rho[i][j] - 0.5*tunnel_halfTime.dt*
        (
            (tunnel_halfTime.rhoF_Plus[i][j] - tunnel_halfTime.rhoF_Minus[i][j])/tunnel_halfTime.dx +
            (tunnel_halfTime.rhoG_Plus[i][j] - tunnel_halfTime.rhoG_Minus[i][j])/tunnel_halfTime.dy
        );
        tunnel_halfTime.Mx[i][j] = tunnel_halfTime.Mx[i][j] - 0.5*tunnel_halfTime.dt*
        (
            (tunnel_halfTime.MxF_Plus[i][j] - tunnel_halfTime.MxF_Minus[i][j])/tunnel_halfTime.dx +
            (tunnel_halfTime.MxG_Plus[i][j] - tunnel_halfTime.MxG_Minus[i][j])/tunnel_halfTime.dy
        );
        tunnel_halfTime.My[i][j] = tunnel_halfTime.My[i][j] - 0.5*tunnel_halfTime.dt*
        (
            (tunnel_halfTime.MyF_Plus[i][j] - tunnel_halfTime.MyF_Minus[i][j])/tunnel_halfTime.dx +
            (tunnel_halfTime.MyG_Plus[i][j] - tunnel_halfTime.MyG_Minus[i][j])/tunnel_halfTime.dy
        );
    }
}

// Update BCs for tunnel_halfTime

// MUSCL Full Time
for (i=0; i < tunnel.yPoints; ++i) {
    for (j=0; j < tunnel.xPoints; ++j) {
        reconstructFlux(i,j,tunnel_halfTime,tunnel,courant,soundSpeed);
        
    }
}
