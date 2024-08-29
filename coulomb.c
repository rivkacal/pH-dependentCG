//********************************************************************************************/
// electrostatics energy term
//
//********************************************************************************************/
#include <stdio.h>
#include <math.h>
void coulomb_(int IC[],int JC[],float Q1Q2[], float X[], float Y[],float Z[],float Fx[],float Fy[],float Fz[] ,int *NE,float *Eelec, float *deConstant)
{
  int i;
  float K = 332.0;//?? not sure this is the correct number.
  float epsilon = *deConstant;
  for(i =0;i < *NE;i++)
  {
      int C1 = IC[i]-1;
      int C2 = JC[i]-1;
      float dx = X[C1] - X[C2];
      float dy = Y[C1] - Y[C2];
      float dz = Z[C1] - Z[C2];

      float r2 = pow(dx,2) + pow(dy,2) + pow(dz,2);
      float r = sqrt(r2);
      float K2 = K*Q1Q2[i]/epsilon;
      *Eelec+= K2/r;

      // force in the direction C1 to C2 ,devided by r (is used to compute forces in x,y,z directions)
      float F_over_r = K2/(r2*r);
		
      Fx[C1]+=  F_over_r*dx;
      Fx[C2]-=  F_over_r*dx;
      Fy[C1]+=  F_over_r*dy;
      Fy[C2]-=  F_over_r*dy;
      Fz[C1]+=  F_over_r*dz;
      Fz[C2]-=  F_over_r*dz;
     }
}
void coulombfactor_(float *sigma, float *deConstant, float *esEnergy)
{
  float K = 332.0;//?? not sure this is the correct number.
  
  *esEnergy = K/((*deConstant)*(sqrt(*sigma)));
}

