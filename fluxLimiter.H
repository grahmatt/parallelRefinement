double fluxLimiter
(
    double& u_Minus,
    double& u,
    double& u_Plus,
    double& distance // char& method
) 
{
    double sjPlus;
    double sjMinus;

    // Superbee
    sjPlus = (u_Plus - u)/distance;
    sjMinus = (u - u_Minus)/distance;

    return copysign(1.0,(sjMinus+sjPlus))*max
        (
            0.0,max
            (
                min(2.0*abs(sjMinus),abs(sjPlus)),
                min(abs(sjMinus),2.0*abs(sjPlus))
            )
        );
}