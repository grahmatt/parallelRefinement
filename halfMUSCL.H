// F
rhoS = fluxLimiter
(
    cell.rho[i][j-1],
    cell.rho[i][j],
    cell.rho[i][j+1],
    cell.dx
);
xMomentS = fluxLimiter
(
    cell.xMomentum[i][j-1],
    cell.xMomentum[i][j],
    cell.xMomentum[i][j+1],
    cell.dx
);
yMomentS = fluxLimiter
(
    cell.yMomentum[i][j-1],
    cell.yMomentum[i][j],
    cell.yMomentum[i][j+1],
    cell.dx
);

rho_Plus = soundSpeed*(cell.rho + 0.5*cell.dx*(1.0 - courant*rhoS));
rho_Minus = soundSpeed*(cell.rho - 0.5*cell.dx*(1.0 - courant*rhoS));
xMoment_Plus = 