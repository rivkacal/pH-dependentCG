!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* ANGL  computes the Force due to the bond angles                     *
!* This code is taken from AMBER, and modified                         *
!***********************************************************************

      SUBROUTINE ANGL(E)
      include 'MD.com'

! "XX, IX," was right after "X," but has been taken out for the time being 
!      IMPLICIT _REAL_ (A-H,O-Z)
      LOGICAL SKIP,NOCRST
C
C     ----- ROUTINE TO GET THE ANGLE ENERGIES AND FORCES FOR THE
C           POTENTIAL OF THE TYPE CT*(T-T0)**2
C
      integer I3, J3, K3
      real  CST,EAW,RIJ,RKJ,RIK,DFW,ANT, ebw
      dimension ebw(190)
      DIMENSION XIJ(nta),YIJ(nta),ZIJ(nta),XKJ(nta),YKJ(nta),
     + ZKJ(nta),CST(nta),EAW(nta),RIJ(nta),RKJ(nta),RIK(nta),
     + DFW(nta),ANT(nta)
      real RIJ0, RKJ0, RIK0, ANT0, DA, ST, STH,
     + CIK, CII, CKK, DT1, DT2, DT3, DT4, DT5, DT6, DT7, DT8, DT9, pt999
     Q , ebal

! These are all replaced with global arrays that don't need to declared
! I think AMBER uses this method since it has such a large memory
! If I need to, I will use this method later.
!      DIMENSION IT(*),JT(*),KT(*),ICT(*),X(*),F(*)

! ",XX(*),IX(*)" was removed from Dimension above

      data pt999 /0.9990d0/
      ebal= 0.0d0

C
C     ----- GRAND LOOP FOR THE angle STUFF -----
C

          DO JN = 1, nTA
            I3 = IT(JN)
            J3 = JT(JN)
            K3 = KT(JN)
C
C           ----- CALCULATION OF THE angle -----
C
            XIJ(JN) = X(I3)-X(J3)
            YIJ(JN) = Y(I3)-Y(J3)
            ZIJ(JN) = Z(I3)-Z(J3)
            XKJ(JN) = X(K3)-X(J3)
            YKJ(JN) = Y(K3)-Y(J3)
            ZKJ(JN) = Z(K3)-Z(J3)
          END DO
C
          DO JN = 1,nTA
            RIJ0 = XIJ(JN)*XIJ(JN)+YIJ(JN)*YIJ(JN)+ZIJ(JN)*ZIJ(JN)
            RKJ0 = XKJ(JN)*XKJ(JN)+YKJ(JN)*YKJ(JN)+ZKJ(JN)*ZKJ(JN)
            RIK0 = SQRT(RIJ0*RKJ0)
            CT0 = (XIJ(JN)*XKJ(JN)+YIJ(JN)*YKJ(JN)+ZIJ(JN)*ZKJ(JN))/RIK0
            CT1 = MAX(-pt999,CT0)
            CT2 = MIN(pt999,CT1)
            CST(JN) = CT2
            ANT(JN) = ACOS(CT2)
            RIJ(JN) = RIJ0
            RKJ(JN) = RKJ0
            RIK(JN) = RIK0
          END DO

! end of insertion


C
C         ----- CALCULATION OF THE ENERGY AND DER -----
C
          
          DO JN = 1,nTA
            ANT0 = ANT(JN)
            DA = ANT0 - ANTC(JN)
            DF = TK(JN)*DA

! These lines were in AMBER, but I don't need them... yet...
!            if(idecomp.eq.1 .or. idecomp.eq.2) then
!             II = (IT(JN) + 3)/3
!             JJ = (JT(JN) + 3)/3
!             KK = (KT(JN) + 3)/3
!              call decangle(XX,IX,II,JJ,KK,EAW(JN))
!            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            DFW(JN) = -(DF+DF)/SIN(ANT0)
          END DO
C
C         ----- CALCULATION OF THE FORCE -----
C
          DO JN = 1,nTA
            I3 = IT(JN)
            J3 = JT(JN)
            K3 = KT(JN)
C
            ST = DFW(JN)
            STH = ST*CST(JN)
            CIK = ST/RIK(JN)
            CII = STH/RIJ(JN)
            CKK = STH/RKJ(JN)
            DT1 = CIK*XKJ(JN)-CII*XIJ(JN)
            DT2 = CIK*YKJ(JN)-CII*YIJ(JN)
            DT3 = CIK*ZKJ(JN)-CII*ZIJ(JN)
            DT7 = CIK*XIJ(JN)-CKK*XKJ(JN)
            DT8 = CIK*YIJ(JN)-CKK*YKJ(JN)
            DT9 = CIK*ZIJ(JN)-CKK*ZKJ(JN)
            DT4 = -DT1-DT7
            DT5 = -DT2-DT8
            DT6 = -DT3-DT9
C

            Fx(I3) = Fx(I3)-DT1
            Fy(I3) = Fy(I3)-DT2
            Fz(I3) = Fz(I3)-DT3
            Fx(J3) = Fx(J3)-DT4
            Fy(J3) = Fy(J3)-DT5
            Fz(J3) = Fz(J3)-DT6
            Fx(K3) = Fx(K3)-DT7
            Fy(K3) = Fy(K3)-DT8
            Fz(K3) = Fz(K3)-DT9

          END DO

! Energy Calculations

          E = 0.0

          do i=1, nTA
             E = E + TK(i)*(ANTC(i)- ANT(i))**2
	!write(*,*) i,TK(i),(ANTC(i)- ANT(i))**2
          end do

          E = E/2.0

       RETURN
       END

!^^^^^^^^^^^^^^^^^^^^^^^^End of ANGL^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
