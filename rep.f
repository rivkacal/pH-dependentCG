!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* NonContacts computes the forces due to non native contacts         *
!**********************************************************************

      subroutine noncontacts(E)
      include 'MD.com'

      integer C1, C2,tt
      real  r2, rm2, rm14, f_over_r 

      E = 0.0

!      write(*,*) Npairnum

!       go over the NCset which is a subset of all the non-native contacts
!         the subset includs all non-native contacts that may be below the
!         repultion distance in this iteration.
	do i=1, NNCt
           tt = NCset(i)
           C1 = INC(tt)
           C2 = JNC(tt)

	dx = X(C1) - X(C2)

 	dy = Y(C1) - Y(C2)

	dz = Z(C1) - Z(C2)

	  r2 = dx**2 + dy**2 + dz**2

	  if(r2 .lt. 6.25*NCsigma(tt))then
             rm2 = 1/r2
             rm14 = rm2**7

! NNCsigma1 is actually 12*eps*sigma**12 (look at read.f and init.f)
              E = E + NNCSigma(tt)*rm14*r2

! f_over_r is the force over the magnitude of r so there is no need to resolve
! the dx, dy and dz into unit vectors
              f_over_r = (NNCSigma(tt))*rm14

! now add the acceleration 
	      FX(C1) = Fx(C1) + f_over_r * dx
	      FY(C1) = Fy(C1) + f_over_r * dy
	      FZ(C1) = Fz(C1) + f_over_r * dz

 	      Fx(C2) =  Fx(C2) - f_over_r * dx
	      Fy(C2) =  Fy(C2) - f_over_r * dy
	      Fz(C2) =  Fz(C2) - f_over_r * dz
	   endif

           end do

           E = E/12.

      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^End of NonContacts^^^^^^^^^^^^^^^^^^^^^^^^^

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* NonContactsNear computes the non-native contacts to be used for the*
!* next NCnear steps                                                  *
!*                                                                    *
!* every 50 steps this subroutine is called                           *
!* it checks each non-native contact,                                 *
!*if the contact length < 9*repultion_distance, it is inseted into the*
!*to check list (I guess the probabilty of an atom traveling          *
!*8*repultion_distance in 50 steps is very low )                      *
!**********************************************************************

      subroutine noncontactsNear
      include 'MD.com'

      integer C1, C2
      real  r2, rm2, rm14, f_over_r

	NNCt=0

         do i=1, NNC

           C1 = INC(i)
           C2 = JNC(i)

        dx = X(C1) - X(C2)

        dy = Y(C1) - Y(C2)

        dz = Z(C1) - Z(C2)

          r2 = dx**2 + dy**2 + dz**2

          if(r2 .lt. 9.0*NCsigma(i))then
	NNCt = NNCt+1
	NCset(NNCt) = i
           endif

           end do


	END

!^^^^^^^^^^^^^^^^^^^^^^^^End of nonContactsNear^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* NonContactscheck checks  the non-native contacts to be used for the*
!* next NCnear steps                                                  *
!*                                                                    *
!*                     NOT USED !!!!                                  *
!**********************************************************************

      subroutine noncontactscheck
      include 'MD.com'

      integer C1, C2,NNCtt,flag, flag1
      real  r2, rm2, rm14, f_over_r

	flag1 =0
!      NNC is the number of non native interation	
      do 4060  i=1, NNC

           C1 = INC(i)
           C2 = JNC(i)

        dx = X(C1) - X(C2)

        dy = Y(C1) - Y(C2)

        dz = Z(C1) - Z(C2)

          r2 = dx**2 + dy**2 + dz**2

          if(r2 .lt. 6.25*NCsigma(i))then

	flag  = 0

	do j=1, NNCt
        if(NCset(j) .eq. i)then
	flag =  1
	goto 4060
	endif

	enddo

	if(flag .eq. 0)then
	flag1 = flag1 + 1
	endif

       endif
4060   continue

        write(36,*) flag1

        END

!^^^^^^^^^^^^^^^^^^^^^^^^End of nonContactsCheck^^^^^^^^^^^^^^^^^^^^^^^^^^^^
