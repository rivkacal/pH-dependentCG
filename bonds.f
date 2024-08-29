!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Bonds  computes the hookean force between chosen atoms              *
!***********************************************************************

      subroutine Bonds(E)
      include 'MD.com'
      integer I2, J2,  outE
      real r2, f, RT,r1,g
      dimension RT(nBA)

      E = 0.0
	do 1 i=1, nBA
           I2 = Ib1(i)
           J2 = Ib2(i)

	dx = X(I2) - X(J2)
	dy = Y(I2) - Y(J2)
	dz = Z(I2) - Z(J2)
       
	  r2 = dx**2 + dy**2 + dz**2
          r1 = sqrt(r2)

! energy calculation
             E = E + bk(i)*(r1-Rb(i))**2

	if( bk(i)*(r1-Rb(i))**2 .gt. 50.0)then

!	write(46,*) 
!	write(46,*) bk(i)*(r1-Rb(i))**2
!	write(46,*) I2,J2
!	write(46,*) r1, Rb(i), bk(i)
	endif
! End energy calculation

! f_over_r is the force over the magnitude of r so there is no need to resolve
! the dx, dy and dz into unit vectors

! the index i indicates the interaction between particle i and i+1

            f = RBC(i)/r1 - bK(i)

            ! now add the force  
	      Fx(I2) = Fx(I2) + f * dx
	      Fy(I2) = Fy(I2) + f * dy
	      Fz(I2) = Fz(I2) + f * dz
! the negative sign is due to the computation of dx, dy and dz
	      Fx(J2) = Fx(J2) - f * dx
	      Fy(J2) = Fy(J2) - f * dy
	      Fz(J2) = Fz(J2) - f * dz

1         continue
             E = E/2.0

      END
      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF BONDS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* BONDSP does the same as ANGLP but for the bond lengths              *
!* computes RBC for each bond (this computation is done only once)     *
!* to compute the force: f = RBC(i)/r1 - bK(i)
!* where r1 is the delta from the optimal length 
!***********************************************************************

      subroutine bondsp
      include 'MD.com'

      do i=1, nBA

      RBC(i) = Rb(i)*bK(i)

      end do

      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF BONDSP^^^^^^^^^^^^^^^^^^^^^^^^^^^^
