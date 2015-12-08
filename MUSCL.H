void MUSCL
(
    grid& tunnel,
    double soundSpeed
) 
{
    int i, j;
    grid tunnel_halfTime = tunnel;

    // MUSCL Half Time
    for (i=0; i < tunnel.yPoints; ++i) {
        for (j=0; j < tunnel.xPoints; ++j) {
            reconstructFlux(i,j,tunnel_halfTime,soundSpeed);
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
    updateBCs(tunnel_halfTime);

    // MUSCL Full Time
    for (i=0; i < tunnel.yPoints; ++i) {
        for (j=0; j < tunnel.xPoints; ++j) {
            reconstructFlux(i,j,tunnel_halfTime,tunnel,soundSpeed);   
        }
    }
    for (i=0; i <= tunnel.yPoints; ++i) {
        for (j=0; j <= tunnel.xPoints; ++j) {
            roeSolver(i,j,tunnel_halfTime,tunnel,soundSpeed);   
        }
    }
    for (i=0; i < tunnel.yPoints; ++i) {
        for (j=0; j < tunnel.xPoints; ++j) {
            if (i+1 == tunnel.yPoints && j+1 == tunnel.xPoints) {
                tunnel.rho[i][j] = tunnel.rho[i][j] - tunnel.dt*
                (
                    (tunnel.rhoF_Plus[i][j] - tunnel.rhoF_Minus[i][j])/tunnel.dx +
                    (tunnel.rhoG_Plus[i][j] - tunnel.rhoG_Minus[i][j])/tunnel.dy
                );
                tunnel.Mx[i][j] = tunnel.Mx[i][j] - tunnel.dt*
                (
                    (tunnel.MxF_Plus[i][j] - tunnel.MxF_Minus[i][j])/tunnel.dx +
                    (tunnel.MxG_Plus[i][j] - tunnel.MxG_Minus[i][j])/tunnel.dy
                );
                tunnel.My[i][j] = tunnel.My[i][j] - tunnel.dt*
                (
                    (tunnel.MyF_Plus[i][j] - tunnel.MyF_Minus[i][j])/tunnel.dx +
                    (tunnel.MyG_Plus[i][j] - tunnel.MyG_Minus[i][j])/tunnel.dy
                );
            }
            else if (j+1 == tunnel.xPoints) {
                tunnel.rho[i][j] = tunnel.rho[i][j] - tunnel.dt*
                (
                    (tunnel.rhoF_Plus[i][j] - tunnel.rhoF_Minus[i][j])/tunnel.dx +
                    (tunnel.rhoG_Minus[i+1][j] - tunnel.rhoG_Minus[i][j])/tunnel.dy
                );
                tunnel.Mx[i][j] = tunnel.Mx[i][j] - tunnel.dt*
                (
                    (tunnel.MxF_Plus[i][j] - tunnel.MxF_Minus[i][j])/tunnel.dx +
                    (tunnel.MxG_Minus[i+1][j] - tunnel.MxG_Minus[i][j])/tunnel.dy
                );
                tunnel.My[i][j] = tunnel.My[i][j] - tunnel.dt*
                (
                    (tunnel.MyF_Plus[i][j] - tunnel.MyF_Minus[i][j])/tunnel.dx +
                    (tunnel.MyG_Minus[i+1][j] - tunnel.MyG_Minus[i][j])/tunnel.dy
                );
            }
            else if (i+1 == tunnel.yPoints) {
                tunnel.rho[i][j] = tunnel.rho[i][j] - tunnel.dt*
                (
                    (tunnel.rhoF_Minus[i][j+1] - tunnel.rhoF_Minus[i][j])/tunnel.dx +
                    (tunnel.rhoG_Plus[i][j] - tunnel.rhoG_Minus[i][j])/tunnel.dy
                );
                tunnel.Mx[i][j] = tunnel.Mx[i][j] - tunnel.dt*
                (
                    (tunnel.MxF_Minus[i][j+1] - tunnel.MxF_Minus[i][j])/tunnel.dx +
                    (tunnel.MxG_Plus[i][j] - tunnel.MxG_Minus[i][j])/tunnel.dy
                );
                tunnel.My[i][j] = tunnel.My[i][j] - tunnel.dt*
                (
                    (tunnel.MyF_Minus[i][j+1] - tunnel.MyF_Minus[i][j])/tunnel.dx +
                    (tunnel.MyG_Plus[i][j] - tunnel.MyG_Minus[i][j])/tunnel.dy
                );
            }
            else {
                tunnel.rho[i][j] = tunnel.rho[i][j] - tunnel.dt*
                (
                    (tunnel.rhoF_Minus[i][j+1] - tunnel.rhoF_Minus[i][j])/tunnel.dx +
                    (tunnel.rhoG_Minus[i+1][j] - tunnel.rhoG_Minus[i][j])/tunnel.dy
                );
                tunnel.Mx[i][j] = tunnel.Mx[i][j] - tunnel.dt*
                (
                    (tunnel.MxF_Minus[i][j+1] - tunnel.MxF_Minus[i][j])/tunnel.dx +
                    (tunnel.MxG_Minus[i+1][j] - tunnel.MxG_Minus[i][j])/tunnel.dy
                );
                tunnel.My[i][j] = tunnel.My[i][j] - tunnel.dt*
                (
                    (tunnel.MyF_Minus[i][j+1] - tunnel.MyF_Minus[i][j])/tunnel.dx +
                    (tunnel.MyG_Minus[i+1][j] - tunnel.MyG_Minus[i][j])/tunnel.dy
                );
            }
        }
    }    
    updateBCs(tunnel);
}

void MUSCL
(
    grid& tunnel,
    double soundSpeed,
    grid& coarseTunnel,
    vector<int> refineColumns,
    int refineFactor
) 
{
    int i, j;
    grid tunnel_halfTime = tunnel;

    // MUSCL Half Time
    for (i=0; i < tunnel.yPoints; ++i) {
        for (j=0; j < tunnel.xPoints; ++j) {
            reconstructFlux(i,j,tunnel_halfTime,soundSpeed);
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
    updateBCs(tunnel_halfTime,coarseTunnel,refineColumns,refineFactor);

    // MUSCL Full Time
    for (i=0; i < tunnel.yPoints; ++i) {
        for (j=0; j < tunnel.xPoints; ++j) {
            reconstructFlux(i,j,tunnel_halfTime,tunnel,soundSpeed);   
        }
    }
    for (i=0; i <= tunnel.yPoints; ++i) {
        for (j=0; j <= tunnel.xPoints; ++j) {
            roeSolver(i,j,tunnel_halfTime,tunnel,soundSpeed);   
        }
    }
    for (i=0; i < tunnel.yPoints; ++i) {
        for (j=0; j < tunnel.xPoints; ++j) {
            if (i+1 == tunnel.yPoints && j+1 == tunnel.xPoints) {
                tunnel.rho[i][j] = tunnel.rho[i][j] - tunnel.dt*
                (
                    (tunnel.rhoF_Plus[i][j] - tunnel.rhoF_Minus[i][j])/tunnel.dx +
                    (tunnel.rhoG_Plus[i][j] - tunnel.rhoG_Minus[i][j])/tunnel.dy
                );
                tunnel.Mx[i][j] = tunnel.Mx[i][j] - tunnel.dt*
                (
                    (tunnel.MxF_Plus[i][j] - tunnel.MxF_Minus[i][j])/tunnel.dx +
                    (tunnel.MxG_Plus[i][j] - tunnel.MxG_Minus[i][j])/tunnel.dy
                );
                tunnel.My[i][j] = tunnel.My[i][j] - tunnel.dt*
                (
                    (tunnel.MyF_Plus[i][j] - tunnel.MyF_Minus[i][j])/tunnel.dx +
                    (tunnel.MyG_Plus[i][j] - tunnel.MyG_Minus[i][j])/tunnel.dy
                );
            }
            else if (j+1 == tunnel.xPoints) {
                tunnel.rho[i][j] = tunnel.rho[i][j] - tunnel.dt*
                (
                    (tunnel.rhoF_Plus[i][j] - tunnel.rhoF_Minus[i][j])/tunnel.dx +
                    (tunnel.rhoG_Minus[i+1][j] - tunnel.rhoG_Minus[i][j])/tunnel.dy
                );
                tunnel.Mx[i][j] = tunnel.Mx[i][j] - tunnel.dt*
                (
                    (tunnel.MxF_Plus[i][j] - tunnel.MxF_Minus[i][j])/tunnel.dx +
                    (tunnel.MxG_Minus[i+1][j] - tunnel.MxG_Minus[i][j])/tunnel.dy
                );
                tunnel.My[i][j] = tunnel.My[i][j] - tunnel.dt*
                (
                    (tunnel.MyF_Plus[i][j] - tunnel.MyF_Minus[i][j])/tunnel.dx +
                    (tunnel.MyG_Minus[i+1][j] - tunnel.MyG_Minus[i][j])/tunnel.dy
                );
            }
            else if (i+1 == tunnel.yPoints) {
                tunnel.rho[i][j] = tunnel.rho[i][j] - tunnel.dt*
                (
                    (tunnel.rhoF_Minus[i][j+1] - tunnel.rhoF_Minus[i][j])/tunnel.dx +
                    (tunnel.rhoG_Plus[i][j] - tunnel.rhoG_Minus[i][j])/tunnel.dy
                );
                tunnel.Mx[i][j] = tunnel.Mx[i][j] - tunnel.dt*
                (
                    (tunnel.MxF_Minus[i][j+1] - tunnel.MxF_Minus[i][j])/tunnel.dx +
                    (tunnel.MxG_Plus[i][j] - tunnel.MxG_Minus[i][j])/tunnel.dy
                );
                tunnel.My[i][j] = tunnel.My[i][j] - tunnel.dt*
                (
                    (tunnel.MyF_Minus[i][j+1] - tunnel.MyF_Minus[i][j])/tunnel.dx +
                    (tunnel.MyG_Plus[i][j] - tunnel.MyG_Minus[i][j])/tunnel.dy
                );
            }
            else {
                tunnel.rho[i][j] = tunnel.rho[i][j] - tunnel.dt*
                (
                    (tunnel.rhoF_Minus[i][j+1] - tunnel.rhoF_Minus[i][j])/tunnel.dx +
                    (tunnel.rhoG_Minus[i+1][j] - tunnel.rhoG_Minus[i][j])/tunnel.dy
                );
                tunnel.Mx[i][j] = tunnel.Mx[i][j] - tunnel.dt*
                (
                    (tunnel.MxF_Minus[i][j+1] - tunnel.MxF_Minus[i][j])/tunnel.dx +
                    (tunnel.MxG_Minus[i+1][j] - tunnel.MxG_Minus[i][j])/tunnel.dy
                );
                tunnel.My[i][j] = tunnel.My[i][j] - tunnel.dt*
                (
                    (tunnel.MyF_Minus[i][j+1] - tunnel.MyF_Minus[i][j])/tunnel.dx +
                    (tunnel.MyG_Minus[i+1][j] - tunnel.MyG_Minus[i][j])/tunnel.dy
                );
            }
        }
    }    
    updateBCs(tunnel,coarseTunnel,refineColumns,refineFactor);
}