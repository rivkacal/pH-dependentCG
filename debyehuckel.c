//********************************************************************************************/
// electrostatics energy term
//
//********************************************************************************************/
#include <stdio.h>
#include <math.h>
void debyehuckel_(int *useESCutoff,
		  int *useDHTable,
		  int IC[],
		  int JC[],
		  float Q1Q2[],
		  float X[],
		  float Y[],
		  float Z[],
		  float Fx[],
		  float Fy[],
		  float Fz[] ,
		  int *NE,
                  float *cutoffDist,
		  float *Eelec,
		  float *deConstant,
		  float *screeningFactor,
		  float *saltCoefficient,
		  float *DebyeHuckelPotentials,
		  float *DebyeHuckelForces) 
{

  float K = 332.0; //?? not sure this is the correct number.
  float epsilon = *deConstant;
  float kappa = *screeningFactor;
  float B_kappa = *saltCoefficient;
  int i,intDistSquared;
  float helper;

  for(i =0;i < *NE;i++)
  {
      int C1 = IC[i]-1;
      int C2 = JC[i]-1;
      float dx = X[C1] - X[C2];
      float dy = Y[C1] - Y[C2];
      float dz = Z[C1] - Z[C2];

      float r2 = pow(dx,2) + pow(dy,2) + pow(dz,2);
      if ((r2< *cutoffDist) || (*useESCutoff == 0))
      {
        float F_over_r;
        if (*useDHTable == 0)
        {
          float r = sqrt(r2); 
          float K2 = K*B_kappa*Q1Q2[i]*exp(-kappa*r)/epsilon;
          *Eelec+= K2/r;
          // force in the direction C1 to C2 ,devided by r (is used to compute forces in x,y,z directions)
          F_over_r = K2*(1/(r2*r) + kappa/r2);
        }
	else
	{
          helper=100*r2;   //make sure to multiple at the inverse of interval in dhEnergyTable.c
          intDistSquared=(int) helper;
          *Eelec+= DebyeHuckelPotentials[intDistSquared]*Q1Q2[i];
          // force in the direction C1 to C2 ,devided by r (is used to compute forces in x,y,z directions)
          F_over_r = DebyeHuckelForces[intDistSquared]*Q1Q2[i];
	}	
        Fx[C1]+=  F_over_r*dx;
        Fx[C2]-=  F_over_r*dx;
        Fy[C1]+=  F_over_r*dy;
        Fy[C2]-=  F_over_r*dy;
        Fz[C1]+=  F_over_r*dz;
        Fz[C2]-=  F_over_r*dz;
      }
  }
}

void debyehuckelfactor_(float *sigma,
			float *deConstant,
			float *screeningFactor,
			float *saltCoefficient,
			float *esEnergy)
{
  float K = 332.0; //?? not sure this is the correct number.
  *esEnergy = K*(*saltCoefficient)*exp(-(*screeningFactor)*(sqrt(*sigma)))/((*deConstant)*(sqrt(*sigma)));
}

