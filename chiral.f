!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* chiral computes the relative distance of the Cbeta bead from chirality 
!* and the resulting force to applied on the bead 
!**********************************************************************

        SUBROUTINE chiral(E)

        include 'MD.com'

      double precision lfac
      integer I3, J3, K3, L3
      real XI_J,YI_J,ZI_J,XI_K,YI_K,ZI_K,XI_L,YI_L,ZI_L,
     +     XIJxIK,YIJxIK,ZIJxIK,TP,deltaTP,dE_dTP,
     +     dYIJxIK_dXI, dZIJxIK_dXI,dXIJxIK_dYI,
     +     dZIJxIK_dYI,dXIJxIK_dZI,dYIJxIK_dZI,
     +     dYIJxIK_dXJ,dZIJxIK_dXJ,dXIJxIK_dYJ,
     +     dZIJxIK_dYJ,dXIJxIK_dZJ,dYIJxIK_dZJ,
     +     dYIJxIK_dXK,dZIJxIK_dXK,dZIJxIK_dYK,
     +     dXIJxIK_dYK,dXIJxIK_dZK,dYIJxIK_dZK,
     +     dTP_dXI,dTP_dXJ,dTP_dXK,dTP_dXL,
     +     dTP_dYI,dTP_dYJ,dTP_dYK,dTP_dYL,
     +     dTP_dZI,dTP_dZJ,dTP_dZK,dTP_dZL

      DIMENSION GMUL(10)
      DATA GMUL/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00,
     +          0.0d+00,8.0d+00,0.0d+00,10.0d+00/
      DATA TM24,TM06,tenm3/1.0d-18,1.0d-06,1.0d-03/
      data zero,one,two,four,six,twelve/0.d0,1.d0,2.d0,4.d0,6.d0,12.d0/
      E = 0.0

C
C     ----- GRAND LOOP FOR THE DIHEDRAL STUFF -----
C
          DO JN = 1,nChirals
            I3 = Ichiral(JN)
            J3 = Jchiral(JN)
            K3 = Kchiral(JN)
            L3 = Lchiral(JN)
C
C           ----- CALCULATION OF ij, ik, il VECTORS -----
C
            XI_J = X(J3)-X(I3)
            YI_J = Y(J3)-Y(I3)
            ZI_J = Z(J3)-Z(I3)
            XI_K = X(K3)-X(I3)
            YI_K = Y(K3)-Y(I3)
            ZI_K = Z(K3)-Z(I3)
            XI_L = X(L3)-X(I3)
            YI_L = Y(L3)-Y(I3)
            ZI_L = Z(L3)-Z(I3)
C
C           ----- CALCULATION OF normal VECTOR -----
C
            XIJxIK = YI_J*ZI_K-ZI_J*YI_K
            YIJxIK = ZI_J*XI_K-XI_J*ZI_K
            ZIJxIK = XI_J*YI_K-YI_J*XI_K
C
C           ----- CALCULATION OF Triple product -----
C
	    TP = XI_L*XIJxIK + YI_L*YIJxIK + ZI_L*ZIJxIK
	    deltaTP = TP-chiralValue(JN)
	    dYIJxIK_dXI = ZI_K-ZI_J
	    dZIJxIK_dXI = YI_J-YI_K
	    dTP_dXI = -XIJxIK + YI_L*dYIJxIK_dXI + ZI_L*dZIJxIK_dXI

	    dXIJxIK_dYI = ZI_J-ZI_K
	    dZIJxIK_dYI = XI_K-XI_J
	    dTP_dYI = -YIJxIK + XI_L*dXIJxIK_dYI + ZI_L*dZIJxIK_dYI

	    dXIJxIK_dZI = YI_K-YI_J
	    dYIJxIK_dZI = XI_J-XI_K
	    dTP_dZI = -ZIJxIK + XI_L*dXIJxIK_dZI + YI_L*dYIJxIK_dZI

	    dYIJxIK_dXJ = -ZI_K
	    dZIJxIK_dXJ =  YI_K
            dTP_dXJ = YI_L*dYIJxIK_dXJ+ZI_L*dZIJxIK_dXJ

	    dXIJxIK_dYJ =  ZI_K
	    dZIJxIK_dYJ = -XI_K
            dTP_dYJ = XI_L*dXIJxIK_dYJ + ZI_L*dZIJxIK_dYJ

	    dXIJxIK_dZJ = -YI_K
	    dYIJxIK_dZJ = XI_K
            dTP_dZJ = XI_L*dXIJxIK_dZJ + YI_L*dYIJxIK_dZJ


	    dYIJxIK_dXK = ZI_J
	    dZIJxIK_dXK = -YI_J
            dTP_dXK = YI_L*dYIJxIK_dXK + ZI_L*dZIJxIK_dXK

	    dXIJxIK_dYK = -ZI_J
	    dZIJxIK_dYK = XI_J
            dTP_dYK = XI_L*dXIJxIK_dYK + ZI_L*dZIJxIK_dYK

	    dXIJxIK_dZK = YI_J
	    dYIJxIK_dZK = -XI_J
            dTP_dZK = XI_L*dXIJxIK_dZK + YI_L*dYIJxIK_dZK

	    dTP_dXL = XIJxIK
	    dTP_dYL = YIJxIK
	    dTP_dZL = ZIJxIK


            dE_dTP = chiralCoeff(JN)*(deltaTP)

            FX(I3) = FX(I3)- dE_dTP*dTP_dXI
            FY(I3) = FY(I3)- dE_dTP*dTP_dYI
            FZ(I3) = FZ(I3)- dE_dTP*dTP_dZI
            FX(J3) = FX(J3)- dE_dTP*dTP_dXJ
            FY(J3) = FY(J3)- dE_dTP*dTP_dYJ
            FZ(J3) = FZ(J3)- dE_dTP*dTP_dZJ
            FX(K3) = FX(K3)- dE_dTP*dTP_dXK
            FY(K3) = FY(K3)- dE_dTP*dTP_dYK
            FZ(K3) = FZ(K3)- dE_dTP*dTP_dZK
            FX(L3) = FX(L3)- dE_dTP*dTP_dXL
            FY(L3) = FY(L3)- dE_dTP*dTP_dYL
            FZ(L3) = FZ(L3)- dE_dTP*dTP_dZL

             E =  E + (1.0/2.0)*chiralCoeff(JN)*(deltaTP**2)
	
	  enddo

          END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END of Chiral ^^^^^^^^^^^^^^^^^^^^^^^^^^^
