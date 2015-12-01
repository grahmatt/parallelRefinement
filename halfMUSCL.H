// F
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

rho_Plus = soundSpeed*(cell.rho + 0.5*cell.dx*(1.0 - courant*rho_S));
rho_Minus = soundSpeed*(cell.rho - 0.5*cell.dx*(1.0 - courant*rho_S));
Mx_Plus = soundSpeed*(cell.Mx + 0.5*cell.dx*(1.0 - courant*Mx_S));
Mx_Minus =