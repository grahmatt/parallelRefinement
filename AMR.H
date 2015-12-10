// Active Mesh Refinement (AMR) only for x direction
for (i=0; i < tunnel.yPoints; ++i) {
    if (refineColumns[0] != 0) {
        rhoFixFPlus[i] = tunnel.rhoF_Plus[i][refineColumns[0] - 1];
        MxFixFPlus[i] = tunnel.MxF_Plus[i][refineColumns[0] - 1];
        MyFixFPlus[i] = tunnel.MyF_Plus[i][refineColumns[0] - 1];
    }
    if (refineColumns.back() != tunnel.xPoints - 1) {
        rhoFixFMinus[i] = - tunnel.rhoF_Minus[i][refineColumns.back() + 1];
        MxFixFMinus[i] = - tunnel.MxF_Minus[i][refineColumns.back() + 1];
        MyFixFMinus[i] = - tunnel.MyF_Minus[i][refineColumns.back() + 1];
    }
}
for (n=0; n < refineFactor; ++n) {
    MUSCL(refineTunnel,soundSpeed,tunnel,refineColumns,refineFactor);
    for (i=0; i < tunnel.yPoints; ++i) {
        for (k=0; k < refineFactor; ++k) {
            if (refineColumns[0] != 0) {
                rhoFixFPlus[i] -= 
                    refineTunnel.rhoF_Minus[i*refineFactor + k][0]/
                    pow(refineFactor,2);
                MxFixFPlus[i] -= 
                    refineTunnel.MxF_Minus[i*refineFactor + k][0]/
                    pow(refineFactor,2);
                MyFixFPlus[i] -= 
                    refineTunnel.MyF_Minus[i*refineFactor + k][0]/
                    pow(refineFactor,2);
            }
            if (refineColumns.back() != tunnel.xPoints - 1) {
                rhoFixFMinus[i] += 
                    refineTunnel.rhoF_Plus[i*refineFactor + k][refineTunnel.xPoints - 1]/
                    pow(refineFactor,2);
                MxFixFMinus[i] += 
                    refineTunnel.MxF_Plus[i*refineFactor + k][refineTunnel.xPoints - 1]/
                    pow(refineFactor,2);
                MyFixFMinus[i] += 
                    refineTunnel.MyF_Plus[i*refineFactor + k][refineTunnel.xPoints - 1]/
                    pow(refineFactor,2);
            }
        }
    }
}
// Coarse border cell flux updates
for (i=0; i < tunnel.yPoints; ++i) {
    if (refineColumns[0] != 0) {
        tunnel.rho[i][refineColumns[0] - 1] += tunnel.dt/tunnel.dx*rhoFixFPlus[i];
        tunnel.Mx[i][refineColumns[0] - 1] += tunnel.dt/tunnel.dx*MxFixFPlus[i];
        tunnel.My[i][refineColumns[0] - 1] += tunnel.dt/tunnel.dx*MyFixFPlus[i];
    }
    if (refineColumns.back() != tunnel.xPoints - 1) {
        tunnel.rho[i][refineColumns.back() + 1] -= tunnel.dt/tunnel.dx*rhoFixFMinus[i];
        tunnel.Mx[i][refineColumns.back() + 1] -= tunnel.dt/tunnel.dx*MxFixFMinus[i];
        tunnel.My[i][refineColumns.back() + 1] -= tunnel.dt/tunnel.dx*MyFixFMinus[i];
    }
}
// Coarse inner cell center updates
double sumRho, sumMx, sumMy;
for (i=0; i < tunnel.yPoints; ++i) {
    for (j=0; j < refineColumns.size(); ++j) {
        sumRho = 0;
        sumMx = 0;
        sumMy = 0;
        for (k=0; k < refineFactor; ++k) {
            for (n=0; n < refineFactor; ++n) {
                sumRho += refineTunnel.rho[i*refineFactor + k][j*refineFactor + n];
                sumMx += refineTunnel.Mx[i*refineFactor + k][j*refineFactor + n];
                sumMy += refineTunnel.My[i*refineFactor + k][j*refineFactor + n];
            }
        }
        tunnel.rho[i][j] = sumRho/pow(refineFactor,2);
        tunnel.Mx[i][j] = sumMx/pow(refineFactor,2);
        tunnel.My[i][j] = sumMy/pow(refineFactor,2);
    }
}