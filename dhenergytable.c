//********************************************************************************************/
// this script generates array of debye-huckel potential in an intervals of r^2
//
//********************************************************************************************/
#include <stdio.h>
#include <math.h>

void dhenergytable_(float *screeningFactor,
                     float *saltCoefficient,
                     float *deConstant,
                     float *esCutoffDistance,
                     float DebyeHuckelPotentials[],
                     float DebyeHuckelForces_over_r[])

{
    float DistanceSquared=0.01,Dist,interval=0.01,cutofflength,K=332.0;
    int i,intcutoff;
    
    cutofflength=*esCutoffDistance/interval;
    intcutoff=(int) cutofflength;     //this is the number of segments (the length of the array)
    for (i=0 ; i<=intcutoff ; i++)
    {
        Dist=sqrt(DistanceSquared);
        DebyeHuckelPotentials[i]= K*(*saltCoefficient)*exp(-(*screeningFactor)*(Dist))/((*deConstant)*(Dist));
        DebyeHuckelForces_over_r[i]=DebyeHuckelPotentials[i]*(1/DistanceSquared + *screeningFactor/Dist);
        DistanceSquared+= interval;
    }
    fflush(stdout);
}
