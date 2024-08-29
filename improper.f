!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* computes the improper dihedral (using triple product)   *
!**********************************************************************

        SUBROUTINE improperDihedral(E)

        include 'MD.com'

      double precision lfac
      integer I3, J3, K3, L3
      real XI_J,YI_J,ZI_J,XI_K,YI_K,ZI_K,XI_L,YI_L,ZI_L,
     +     XIJxIK,YIJxIK,ZIJxIK,TP,
     +     size_IJxIK,size_IL,one_over_size_IL,one_over_size_IL_3,
     +     one_over_size_IJxIK,one_over_size_IJxIK_3,factor_IJxIK_IL,
     +     factor_TP_IL,factor_TP_IJxIK,cos_phi,phi,deltaPhi,dE_dphi,
     +     dphi_dcosphi,const_factor,dYIJxIK_dXI, dZIJxIK_dXI,
     +     dXIJxIK_dYI,dZIJxIK_dYI,dXIJxIK_dZI,dYIJxIK_dZI,dYIJxIK_dXJ,
     +     dZIJxIK_dXJ,dXIJxIK_dYJ,dZIJxIK_dYJ,dXIJxIK_dZJ,dYIJxIK_dZJ,
     +     dYIJxIK_dXK,dZIJxIK_dXK,dXIJxIK_dYK,dZIJxIK_dYK,dXIJxIK_dZK,
     +     dYIJxIK_dZK,dTP_dXI,doosIJxIK_dXI,doosIL_dXI,
     +     dTP_dYI,doosIJxIK_dYI,doosIL_dYI,dTP_dZI,doosIJxIK_dZI,
     +     doosIL_dZI,dTP_dXJ,doosIJxIK_dXJ,dTP_dYJ,doosIJxIK_dYJ,
     +     dTP_dZJ,doosIJxIK_dZJ,dTP_dXK,doosIJxIK_dXK,dTP_dYK,
     +     doosIJxIK_dYK,dTP_dZK,doosIJxIK_dZK,dTP_dXL,doosIL_dXL,
     +     dTP_dYL,doosIL_dYL,dTP_dZL,doosIL_dZL,
     +     dcosphi_dXI,dcosphi_dYI,dcosphi_dZI,dcosphi_dXJ,dcosphi_dYJ,
     +     dcosphi_dZJ,dcosphi_dXK,dcosphi_dYK,dcosphi_dZK,dcosphi_dXL,
     +     dcosphi_dYL,dcosphi_dZL

      DIMENSION GMUL(10)
      DATA GMUL/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00,
     +          0.0d+00,8.0d+00,0.0d+00,10.0d+00/
      DATA TM24,TM06,tenm3/1.0d-18,1.0d-06,1.0d-03/
      data zero,one,two,four,six,twelve/0.d0,1.d0,2.d0,4.d0,6.d0,12.d0/
      E = 0.0

C
C     ----- GRAND LOOP FOR THE DIHEDRAL STUFF -----
C
          DO JN = 1,nImpropers
            I3 = Iimproper(JN)
            J3 = Jimproper(JN)
            K3 = Kimproper(JN)
            L3 = Limproper(JN)
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
C
C           ----- CALCULATION OF VECTOR sizes -----
C
            size_IJxIK = SQRT(XIJxIK**2 + 
     Q                            YIJxIK**2 +
     Q                            ZIJxIK**2+TM24)
            size_IL = SQRT(XI_L**2+YI_L**2+ZI_L**2+TM24)

            one_over_size_IJxIK = one/size_IJxIK
            one_over_size_IJxIK_3 = one_over_size_IJxIK**3		
            one_over_size_IL = one/size_IL
            one_over_size_IL_3 = one_over_size_IL**3

            factor_IJxIK_IL = one_over_size_IJxIK * one_over_size_IL
	    factor_TP_IL = TP * one_over_size_IL
	    factor_TP_IJxIK = TP * one_over_size_IJxIK
C
C           ----- CALCULATION OF angle -----
C
            cos_phi = MIN(one,TP*factor_IJxIK_IL)
            cos_phi = MAX(-one,cos_phi)
            phi = ACOS(cos_phi)
	    phi = SIGN(phi,ASIN(cos_phi))
C
C           ----- CALCULATION OF energy and forces -----
C
	    deltaPhi = improperAngle(JN)-phi
            if (minimalDeltaAngle .lt. deltaPhi) then

	      dYIJxIK_dXI = ZI_K-ZI_J
	      dZIJxIK_dXI = YI_J-YI_K

	      dTP_dXI = -XIJxIK + YI_L*dYIJxIK_dXI + ZI_L*dZIJxIK_dXI
	      doosIJxIK_dXI = -1*one_over_size_IJxIK_3*
     Q                         (YIJxIK*dYIJxIK_dXI +
     Q                          ZIJxIK*dZIJxIK_dXI)
	      doosIL_dXI = one_over_size_IL_3*XI_L

	      dXIJxIK_dYI = ZI_J-ZI_K
	      dZIJxIK_dYI = XI_K-XI_J

	      dTP_dYI = -YIJxIK + XI_L*dXIJxIK_dYI + ZI_L*dZIJxIK_dYI
	      doosIJxIK_dYI = -1*one_over_size_IJxIK_3*
     Q                         (XIJxIK*dXIJxIK_dYI +
     Q                          ZIJxIK*dZIJxIK_dYI)
	      doosIL_dYI = one_over_size_IL_3 * YI_L

	      dXIJxIK_dZI = YI_K-YI_J
	      dYIJxIK_dZI = XI_J-XI_K

	      dTP_dZI = -ZIJxIK + XI_L*dXIJxIK_dZI + YI_L*dYIJxIK_dZI
	      doosIJxIK_dZI = -1 * one_over_size_IJxIK_3 *
     Q                         (XIJxIK*dXIJxIK_dZI +
     Q                          YIJxIK*dYIJxIK_dZI)
	      doosIL_dZI = one_over_size_IL_3 * ZI_L

	      dYIJxIK_dXJ = -ZI_K
	      dZIJxIK_dXJ =  YI_K

              dTP_dXJ = YI_L*dYIJxIK_dXJ+ZI_L*dZIJxIK_dXJ
	      doosIJxIK_dXJ = -1 * one_over_size_IJxIK_3 *
     Q                        (YIJxIK*dYIJxIK_dXJ + ZIJxIK*dZIJxIK_dXJ)

	      dXIJxIK_dYJ =  ZI_K
	      dZIJxIK_dYJ = -XI_K

              dTP_dYJ = XI_L*dXIJxIK_dYJ + ZI_L*dZIJxIK_dYJ
              doosIJxIK_dYJ = -1 * one_over_size_IJxIK_3 *
     Q                        (XIJxIK*dXIJxIK_dYJ + ZIJxIK*dZIJxIK_dYJ)

	      dXIJxIK_dZJ = -YI_K
	      dYIJxIK_dZJ = XI_K

              dTP_dZJ = XI_L*dXIJxIK_dZJ + YI_L*dYIJxIK_dZJ
              doosIJxIK_dZJ = -1 * one_over_size_IJxIK_3 *
     Q                        (XIJxIK*dXIJxIK_dZJ + YIJxIK*dYIJxIK_dZJ)


	      dYIJxIK_dXK = ZI_J
	      dZIJxIK_dXK = -YI_J

              dTP_dXK = YI_L*dYIJxIK_dXK + ZI_L*dZIJxIK_dXK
	      doosIJxIK_dXK = -1 * one_over_size_IJxIK_3 * 
     Q                        (YIJxIK*dYIJxIK_dXK + ZIJxIK*dZIJxIK_dXK)

	      dXIJxIK_dYK = -ZI_J
	      dZIJxIK_dYK = XI_J

              dTP_dYK = XI_L*dXIJxIK_dYK + ZI_L*dZIJxIK_dYK
              doosIJxIK_dYK = -1 * one_over_size_IJxIK_3 * 
     Q                        (XIJxIK*dXIJxIK_dYK + ZIJxIK*dZIJxIK_dYK)

	      dXIJxIK_dZK = YI_J
	      dYIJxIK_dZK = -XI_J

              dTP_dZK = XI_L*dXIJxIK_dZK + YI_L*dYIJxIK_dZK
              doosIJxIK_dZK = -1 * one_over_size_IJxIK_3 * 
     Q                         (XIJxIK*dXIJxIK_dZK + YIJxIK*dYIJxIK_dZK)

	      dTP_dXL = XIJxIK
	      doosIL_dXL = -1*one_over_size_IL_3*XI_L 

	      dTP_dYL = YIJxIK
	      doosIL_dYL = -1*one_over_size_IL_3*YI_L 

	      dTP_dZL = ZIJxIK
	      doosIL_dZL = -1*one_over_size_IL_3*ZI_L 


              dcosphi_dXI = dTP_dXI * factor_IJxIK_IL + 
     Q                      factor_TP_IL * doosIJxIK_dXI + 
     Q                      factor_TP_IJxIK * doosIL_dXI 
              dcosphi_dYI = dTP_dYI * factor_IJxIK_IL + 
     Q                      factor_TP_IL * doosIJxIK_dYI + 
     Q                      factor_TP_IJxIK * doosIL_dYI 
              dcosphi_dZI = dTP_dZI * factor_IJxIK_IL + 
     Q                      factor_TP_IL * doosIJxIK_dZI + 
     Q                      factor_TP_IJxIK * doosIL_dZI 

              dcosphi_dXJ = dTP_dXJ*factor_IJxIK_IL + 
     Q                      factor_TP_IL*doosIJxIK_dXJ 
              dcosphi_dYJ = dTP_dYJ*factor_IJxIK_IL + 
     Q                      factor_TP_IL*doosIJxIK_dYJ 
              dcosphi_dZJ = dTP_dZJ*factor_IJxIK_IL + 
     Q                      factor_TP_IL*doosIJxIK_dZJ 

              dcosphi_dXK = dTP_dXK*factor_IJxIK_IL + 
     Q                      factor_TP_IL*doosIJxIK_dXK 
              dcosphi_dYK = dTP_dYK*factor_IJxIK_IL + 
     Q                      factor_TP_IL*doosIJxIK_dYK 
              dcosphi_dZK = dTP_dZK*factor_IJxIK_IL + 
     Q                      factor_TP_IL*doosIJxIK_dZK 

              dcosphi_dXL = dTP_dXL*factor_IJxIK_IL + 
     Q                      factor_TP_IJxIK*doosIL_dXL 
              dcosphi_dYL = dTP_dYL*factor_IJxIK_IL + 
     Q                      factor_TP_IJxIK*doosIL_dYL 
              dcosphi_dZL = dTP_dZL*factor_IJxIK_IL + 
     Q                      factor_TP_IJxIK*doosIL_dZL 


	      dphi_dcosphi = -1.0/(sqrt(1-cos_phi**2)) 
              dE_dPhi = improperCoeff(JN)*(deltaPhi)
 	      const_factor = dE_dPhi*dphi_dcosphi

              FX(I3) = FX(I3)- const_factor*dcosphi_dXI
              FY(I3) = FY(I3)- const_factor*dcosphi_dYI
              FZ(I3) = FZ(I3)- const_factor*dcosphi_dZI
              FX(J3) = FX(J3)- const_factor*dcosphi_dXJ
              FY(J3) = FY(J3)- const_factor*dcosphi_dYJ
              FZ(J3) = FZ(J3)- const_factor*dcosphi_dZJ
              FX(K3) = FX(K3)- const_factor*dcosphi_dXK
              FY(K3) = FY(K3)- const_factor*dcosphi_dYK
              FZ(K3) = FZ(K3)- const_factor*dcosphi_dZK
              FX(L3) = FX(L3)- const_factor*dcosphi_dXL
              FY(L3) = FY(L3)- const_factor*dcosphi_dYL
              FZ(L3) = FZ(L3)- const_factor*dcosphi_dZL

              E =  E + (1.0/2.0)*improperCoeff(JN)*(deltaPhi**2)
	
	    endif
	  enddo

          END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END of Improper Dihedral^^^^^^^^^^^^^^^^^^^^^^^^^^^
