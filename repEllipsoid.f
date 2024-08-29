!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* NonContacts computes the forces due to non native contacts         *
!**********************************************************************

      subroutine ellipsoidRepulsions(E)
      include 'MD.com'

      integer firstBeadIndex, secondBeadIndex,thirdBeadIndex, cerIter,
     Q        cerIndex
      real  dx1,dy1,dz1,dx2,dy2,dz2, r1,r2,one_over_r1_and_r2,f_over_r,
     Q      oor1a2m12, dr1_x1, dr1_x3, dr1_y1, dr1_y3, dr1_z1, dr1_z3,
     Q      dr2_x2, dr2_x3, dr2_y2, dr2_y3, dr2_z2, dr2_z3, r1ar2m12

      E = 0.0

	do cerIter=1, currentEllipsoidRepulsionsNum
           cerIndex = currentEllipsoidRepulsions(cerIter)
           firstBeadIndex = Iellipsoid(cerIndex)
           secondBeadIndex = Jellipsoid(cerIndex)
           thirdBeadIndex = Kellipsoid(cerIndex)

	dx1 = X(firstBeadIndex) - X(thirdBeadIndex)
	dy1 = Y(firstBeadIndex) - Y(thirdBeadIndex)
	dz1 = Z(firstBeadIndex) - Z(thirdBeadIndex)
	dx2 = X(secondBeadIndex) - X(thirdBeadIndex)
	dy2 = Y(secondBeadIndex) - Y(thirdBeadIndex)
	dz2 = Z(secondBeadIndex) - Z(thirdBeadIndex)

	  r1 = sqrt(dx1**2 + dy1**2 + dz1**2)
	  r2 = sqrt(dx2**2 + dy2**2 + dz2**2)

	  if((r1+r2) .lt. 3.0*ellipsoidSigma(cerIndex))then
             one_over_r1_and_r2 = 1.0/(r1+r2)
	     r1ar2m12 = (r1+r2)**12
             oor1a2m12 = one_over_r1_and_r2**12
              E = E + ellipsoidCoeff(cerIndex) 
     Q            * ellipsoidSigma(cerIndex)**12 * oor1a2m12

              f_over_r = 12.0*E*one_over_r1_and_r2
	      dr1_x1 = (r1**-1)*(dx1)
	      dr1_y1 = (r1**-1)*(dy1)
	      dr1_z1 = (r1**-1)*(dz1)
	      dr1_x3 = -1.0*(r1**-1)*(dx1)
	      dr1_y3 = -1.0*(r1**-1)*(dy1)
	      dr1_z3 = -1.0*(r1**-1)*(dz1)
	      dr2_x2 = (r2**-1)*(dx2)
	      dr2_y2 = (r2**-1)*(dy2)
	      dr2_z2 = (r2**-1)*(dz2)
	      dr2_x3 = -1.0*(r2**-1)*(dx2)
	      dr2_y3 = -1.0*(r2**-1)*(dy2)
	      dr2_z3 = -1.0*(r2**-1)*(dz2)

	      FX(firstBeadIndex) = Fx(firstBeadIndex) + f_over_r*dr1_x1
	      FY(firstBeadIndex) = Fy(firstBeadIndex) + f_over_r*dr1_y1
	      FZ(firstBeadIndex) = Fz(firstBeadIndex) + f_over_r*dr1_z1

	      FX(secondBeadIndex) = Fx(secondBeadIndex) + f_over_r*dr2_x2
	      FY(secondBeadIndex) = Fy(secondBeadIndex) + f_over_r*dr2_y2
	      FZ(secondBeadIndex) = Fz(secondBeadIndex) + f_over_r*dr2_z2

 	      Fx(thirdBeadIndex) =  Fx(thirdBeadIndex) 
     Q                              + f_over_r*(dr1_x3+dr2_x3)
	      Fy(thirdBeadIndex) =  Fy(thirdBeadIndex) 
     Q                              + f_over_r*(dr1_y3+dr2_y3)
	      Fz(thirdBeadIndex) =  Fz(thirdBeadIndex) 
     Q                              + f_over_r*(dr1_z3+dr2_z3)
	   endif

           end do

      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^End of ellipsoidRepulsions^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* ellipsoidRepulsionsCheck checks the ellipsoid repulsions to be     *
!* used in the next 50 steps                                          *
!*                                                                    *
!**********************************************************************

      subroutine ellipsoidRepulsionsCheck
      include 'MD.com'

       integer firstBeadIndex, secondBeadIndex,thirdBeadIndex,
     Q         cerCounter, erIter
       real  dx1,dy1,dz1,dx2,dy2,dz2,r1,r2

	cerCounter = 0

      do erIter=1, ellipsoidRepulsionsNum

        firstBeadIndex = Iellipsoid(erIter)
        secondBeadIndex = Jellipsoid(erIter)
        thirdBeadIndex = Kellipsoid(erIter)

        dx1 = X(thirdBeadIndex) - X(firstBeadIndex)
        dy1 = Y(thirdBeadIndex) - Y(firstBeadIndex)
        dz1 = Z(thirdBeadIndex) - Z(firstBeadIndex)

        dx2 = X(thirdBeadIndex) - X(secondBeadIndex)
        dy2 = Y(thirdBeadIndex) - Y(secondBeadIndex)
        dz2 = Z(thirdBeadIndex) - Z(secondBeadIndex)

	r1 = sqrt(dx1**2 + dy1**2 + dz1**2)
	r2 = sqrt(dx2**2 + dy2**2 + dz2**2)
   
        if((r1+r2) .lt. 3.0*ellipsoidSigma(erIter))then
	  cerCounter = cerCounter + 1
	  currentEllipsoidRepulsions(cerCounter) = erIter
        end if
      enddo

       currentEllipsoidRepulsionsNum = cerCounter
       END

!^^^^^^^^^^^^^^^^^^^^^^^^End of ellipsoidRepulsionsCheck^^^^^^^^^^^^^^^^^^^^^^^^^^^^
