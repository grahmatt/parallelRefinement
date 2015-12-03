void roeSolver
(
    int i,
    int j,
    grid& cellHalf,
    grid& cellFull,
    double soundSpeed
) 
{
    double rho_Plus, rho_Minus, Mx_Plus, Mx_Minus, My_Plus, My_Minus,
        rhoF_Plus, rhoF_Minus, MxF_Plus, MxF_Minus, MyF_Plus, MyF_Minus,
        uAveRoe, vAveRoe, lambda_1, lambda_2, lambda_3, alpha_1, alpha_2, alpha_3;

    // Roe Solver x-dir
    if (i < cellHalf.yPoints) {
        if (j == 0) {
            rho_Plus = cellHalf.rhoBCsLR_PlusMinus[i][0];
            rho_Minus = cellHalf.rhoF_Minus[i][j];
            Mx_Plus = cellHalf.MxBCsLR_PlusMinus[i][0];
            Mx_Minus = cellHalf.MxF_Minus[i][j];
            My_Plus = cellHalf.MyBCsLR_PlusMinus[i][0];
            My_Minus = cellHalf.MyF_Minus[i][j];
        }
        else if (j == cellHalf.xPoints) {
            rho_Plus = cellHalf.rhoF_Plus[i][j-1];
            rho_Minus = cellHalf.rhoBCsLR_PlusMinus[i][1];
            Mx_Plus = cellHalf.MxF_Plus[i][j-1];
            Mx_Minus = cellHalf.MxBCsLR_PlusMinus[i][1];
            My_Plus = cellHalf.MyF_Plus[i][j-1];
            My_Minus = cellHalf.MyBCsLR_PlusMinus[i][1];
        }
        else {
            rho_Plus = cellHalf.rhoF_Plus[i][j-1];
            rho_Minus = cellHalf.rhoF_Minus[i][j];
            Mx_Plus = cellHalf.MxF_Plus[i][j-1];
            Mx_Minus = cellHalf.MxF_Minus[i][j];
            My_Plus = cellHalf.MyF_Plus[i][j-1];
            My_Minus = cellHalf.MyF_Minus[i][j];
        }
        rhoF_Plus = Mx_Plus;
        rhoF_Minus = Mx_Minus;
        MxF_Plus = pow(Mx_Plus,2.0)/rho_Plus + 
            pow(soundSpeed, 2.0)*rho_Plus;
        MxF_Minus = pow(Mx_Minus,2.0)/rho_Minus + 
            pow(soundSpeed, 2.0)*rho_Minus;
        MyF_Plus = Mx_Plus*My_Plus/rho_Plus;
        MyF_Minus = Mx_Minus*My_Minus/rho_Minus;
        // F* Calculations
        uAveRoe = (Mx_Plus/sqrt(rho_Plus) + Mx_Minus/sqrt(rho_Minus))/(sqrt(rho_Plus) + sqrt(rho_Minus));
        vAveRoe = (My_Plus/sqrt(rho_Plus) + My_Minus/sqrt(rho_Minus))/(sqrt(rho_Plus) + sqrt(rho_Minus));
        lambda_1 = abs(uAveRoe);
        lambda_2 = abs(uAveRoe + soundSpeed);
        lambda_3 = abs(uAveRoe - soundSpeed);
        alpha_1 = vAveRoe*(rho_Plus - rho_Minus) - (My_Plus - My_Minus);
        alpha_2 = 0.5*vAveRoe/soundSpeed*(soundSpeed - uAveRoe)*(rho_Minus - rho_Plus) + 
            0.5*vAveRoe/soundSpeed*(Mx_Minus - Mx_Plus);
        alpha_3 = 0.5*vAveRoe/soundSpeed*(uAveRoe + soundSpeed)*(rho_Minus - rho_Plus) -
            0.5*vAveRoe/soundSpeed*(Mx_Minus - Mx_Plus);
        if (j != cellHalf.xPoints) {
            cellFull.rhoF_Minus[i][j] = 0.5*(rhoF_Plus + rhoF_Minus) -
                0.5*(lambda_2*(0.5/soundSpeed*(soundSpeed - uAveRoe)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*(Mx_Minus - Mx_Plus)) + 
                lambda_3*(0.5/soundSpeed*(uAveRoe + soundSpeed)*(rho_Minus - rho_Plus) - 
                0.5/soundSpeed*(Mx_Minus - Mx_Plus)));
            cellFull.MxF_Minus[i][j] = 0.5*(MxF_Plus + MxF_Minus) - 
                0.5*(lambda_2*(0.5/soundSpeed*(soundSpeed - uAveRoe)*(uAveRoe + soundSpeed)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*(uAveRoe + soundSpeed)*(Mx_Minus - Mx_Plus)) + 
                lambda_3*(0.5/soundSpeed*(uAveRoe + soundSpeed)*(uAveRoe - soundSpeed)*(rho_Minus - rho_Plus) - 
                0.5/soundSpeed*(uAveRoe - soundSpeed)*(Mx_Minus - Mx_Plus)));
            cellFull.MyF_Minus[i][j] = 0.5*(MyF_Plus + MyF_Minus) - 
                0.5*(lambda_1*alpha_1 + lambda_2*alpha_2 + lambda_3*alpha_3);
        }
        else {
            cellFull.rhoF_Plus[i][j-1] = 0.5*(rhoF_Plus + rhoF_Minus) -
                0.5*(lambda_2*(0.5/soundSpeed*(soundSpeed - uAveRoe)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*(Mx_Minus - Mx_Plus)) + 
                lambda_3*(0.5/soundSpeed*(uAveRoe + soundSpeed)*(rho_Minus - rho_Plus) - 
                0.5/soundSpeed*(Mx_Minus - Mx_Plus)));
            cellFull.MxF_Plus[i][j-1] = 0.5*(MxF_Plus + MxF_Minus) - 
                0.5*(lambda_2*(0.5/soundSpeed*(soundSpeed - uAveRoe)*(uAveRoe + soundSpeed)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*(uAveRoe + soundSpeed)*(Mx_Minus - Mx_Plus)) + 
                lambda_3*(0.5/soundSpeed*(uAveRoe + soundSpeed)*(uAveRoe - soundSpeed)*(rho_Minus - rho_Plus) - 
                0.5/soundSpeed*(uAveRoe - soundSpeed)*(Mx_Minus - Mx_Plus)));
            cellFull.MyF_Plus[i][j-1] = 0.5*(MyF_Plus + MyF_Minus) - 
                0.5*(lambda_1*alpha_1 + lambda_2*alpha_2 + lambda_3*alpha_3);
        }
    }

    // Roe Solver y-dir
    if (j < cellHalf.xPoints) {
        if (i == 0) {
            rho_Plus = cellHalf.rhoBCsBT_PlusMinus[0][j];
            rho_Minus = cellHalf.rhoG_Minus[i][j];
            Mx_Plus = cellHalf.MxBCsBT_PlusMinus[0][j];
            Mx_Minus = cellHalf.MxG_Minus[i][j];
            My_Plus = cellHalf.MyBCsBT_PlusMinus[0][j];
            My_Minus = cellHalf.MyG_Minus[i][j];
        }
        else if (i == cellHalf.yPoints) {
            rho_Plus = cellHalf.rhoG_Plus[i-1][j];
            rho_Minus = cellHalf.rhoBCsBT_PlusMinus[1][j];
            Mx_Plus = cellHalf.MxG_Plus[i-1][j];
            Mx_Minus = cellHalf.MxBCsBT_PlusMinus[1][j];
            My_Plus = cellHalf.MyG_Plus[i-1][j];
            My_Minus = cellHalf.MyBCsBT_PlusMinus[1][j];
        }
        else {
            rho_Plus = cellHalf.rhoG_Plus[i-1][j];
            rho_Minus = cellHalf.rhoG_Minus[i][j];
            Mx_Plus = cellHalf.MxG_Plus[i-1][j];
            Mx_Minus = cellHalf.MxG_Minus[i][j];
            My_Plus = cellHalf.MyG_Plus[i-1][j];
            My_Minus = cellHalf.MyG_Minus[i][j];
        }
        rhoF_Plus = Mx_Plus;
        rhoF_Minus = Mx_Minus;
        MxF_Plus = pow(Mx_Plus,2.0)/rho_Plus + 
            pow(soundSpeed, 2.0)*rho_Plus;
        MxF_Minus = pow(Mx_Minus,2.0)/rho_Minus + 
            pow(soundSpeed, 2.0)*rho_Minus;
        MyF_Plus = Mx_Plus*My_Plus/rho_Plus;
        MyF_Minus = Mx_Minus*My_Minus/rho_Minus;
        // G* Calculations
        uAveRoe = (Mx_Plus/sqrt(rho_Plus) + Mx_Minus/sqrt(rho_Minus))/(sqrt(rho_Plus) + sqrt(rho_Minus));
        vAveRoe = (My_Plus/sqrt(rho_Plus) + My_Minus/sqrt(rho_Minus))/(sqrt(rho_Plus) + sqrt(rho_Minus));
        lambda_1 = abs(vAveRoe);
        lambda_2 = abs(vAveRoe + soundSpeed);
        lambda_3 = abs(vAveRoe - soundSpeed);
        alpha_1 = uAveRoe*(rho_Plus - rho_Minus) - (Mx_Plus - Mx_Minus);
        alpha_2 = 0.5/soundSpeed*(vAveRoe + soundSpeed)*(soundSpeed - vAveRoe)*(rho_Minus - rho_Plus) + 
            0.5/soundSpeed*(vAveRoe + soundSpeed)*(My_Minus - My_Plus);
        alpha_3 = 0.5/soundSpeed*(vAveRoe + soundSpeed)*(vAveRoe - soundSpeed)*(rho_Minus - rho_Plus) + 
            0.5/soundSpeed*(soundSpeed - vAveRoe)*(My_Minus - My_Plus);
        if (i != cellHalf.yPoints) {
            cellFull.rhoG_Minus[i][j] = 0.5*(rhoF_Plus + rhoF_Minus) -
                0.5*(lambda_2*(0.5/soundSpeed*(soundSpeed - vAveRoe)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*(My_Minus - My_Plus)) + 
                lambda_3*(0.5/soundSpeed*(vAveRoe + soundSpeed)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*(My_Minus - My_Plus)));
            cellFull.MxG_Minus[i][j] = 0.5*(MxF_Plus + MxF_Minus) - 
                0.5*(lambda_1*alpha_1 + 
                lambda_2*(0.5/soundSpeed*uAveRoe*(soundSpeed - vAveRoe)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*uAveRoe*(My_Minus - My_Plus)) + 
                lambda_3*(0.5/soundSpeed*uAveRoe*(vAveRoe + soundSpeed)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*uAveRoe*(My_Minus - My_Plus)));
            cellFull.MyG_Minus[i][j] = 0.5*(MyF_Plus + MyF_Minus) - 
                0.5*(lambda_2*alpha_2 + lambda_3*alpha_3);
        }
        else {
            cellFull.rhoG_Plus[i-1][j] = 0.5*(rhoF_Plus + rhoF_Minus) -
                0.5*(lambda_2*(0.5/soundSpeed*(soundSpeed - vAveRoe)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*(My_Minus - My_Plus)) + 
                lambda_3*(0.5/soundSpeed*(vAveRoe + soundSpeed)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*(My_Minus - My_Plus)));
            cellFull.MxG_Plus[i-1][j] = 0.5*(MxF_Plus + MxF_Minus) - 
                0.5*(lambda_1*alpha_1 + 
                lambda_2*(0.5/soundSpeed*uAveRoe*(soundSpeed - vAveRoe)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*uAveRoe*(My_Minus - My_Plus)) + 
                lambda_3*(0.5/soundSpeed*uAveRoe*(vAveRoe + soundSpeed)*(rho_Minus - rho_Plus) + 
                0.5/soundSpeed*uAveRoe*(My_Minus - My_Plus)));
            cellFull.MyG_Plus[i-1][j] = 0.5*(MyF_Plus + MyF_Minus) - 
                0.5*(lambda_2*alpha_2 + lambda_3*alpha_3);
        }
    }
}
