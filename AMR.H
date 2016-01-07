// Active Mesh Refinement (AMR) only for x direction and no flux corrections on coarse boundary with refine cells
for (n=0; n < refineFactor; ++n) {
    MUSCL(refineTunnel, soundSpeed, tunnel, refineColumns, refineFactor, n, 
        rank, sendLeft, recvLeft, sendRight, recvRight);
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
        tunnel.rho[i][refineColumns[j]] = sumRho/pow(refineFactor,2.0);
        tunnel.Mx[i][refineColumns[j]] = sumMx/pow(refineFactor,2.0);
        tunnel.My[i][refineColumns[j]] = sumMy/pow(refineFactor,2.0);
    }
}