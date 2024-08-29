//********************************************************************************************/
// creates a square box with dimentions Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
//
//********************************************************************************************/
#include <stdio.h>
#include <math.h>

void box_(float X[], float Y[],float Z[],float Fx[],float Fy[],float Fz[] ,int DynAtomRange[], int *DynAtomlen,float *Ebox,float boxMin[],float boxMax[],float *Kforce) {
       float K  = *Kforce;
       float K2 = K*2;
       float Xmin  = boxMin[0];
       float Xmax = boxMax[0];
       float Ymin  = boxMin[1];
       float Ymax = boxMax[1];
       float Zmin  = boxMin[2];
       float Zmax = boxMax[2];
       int CCCi = 0;
       int i = 0;
       //energy in every direction
       float ExBox = 0;
       float EyBox = 0;
       float EzBox = 0;

//       F = K*dR
//       E = (K*(dR**2))/2
// R is the distance of the penetration of the atom to the wall of the box
      for (i = 0;i+1 < *DynAtomlen;i=i+2){
       for (CCCi =DynAtomRange[i];CCCi< DynAtomRange[i+1];CCCi++){
		if (X[CCCi] < Xmin ){
			Fx[CCCi]+= -(K2* (X[CCCi] - Xmin));
			ExBox = K*pow((X[CCCi] - Xmin),2);
		}
		if (X[CCCi] > Xmax){
			Fx[CCCi]+= -(K2* (X[CCCi] - Xmax));
			ExBox = K*pow((X[CCCi] - Xmax),2);				
		}
		
		if (Y[CCCi] < Ymin ){
			Fy[CCCi]+= -(K2* (Y[CCCi] - Ymin));
			EyBox = K*pow((Y[CCCi] - Ymin),2);
		}
		if (Y[CCCi] > Ymax){
			Fy[CCCi]+= -(K2* (Y[CCCi] - Ymax));
			EyBox = K*pow((Y[CCCi] - Ymax),2);
		}	
		if (Z[CCCi] < Zmin ){
			Fz[CCCi]+= -(K2* (Z[CCCi] - Zmin));
			EzBox = K*pow((Z[CCCi] - Zmin),2);	
		}
		if (Z[CCCi] > Zmax){
			Fz[CCCi]+= -(K2* (Z[CCCi] - Zmax));
			EzBox = K*pow((Z[CCCi] - Zmax),2);
		}
		*Ebox += ExBox + EyBox + EzBox;	
       }
       }
}
