!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! * subroutine RANDOM                                                  *
! * Random generator of variate on 0, 1                                *
! * Purpose: generates a pseudorandom number RAND                      *
! *    whose distribution is flat on the interval (0,1).               *
! *    Uses as seeds xrandom, yrandom, zrandom stored between runs     *
! *    in file RANDOM. DAT                                             *
! * Variables: rand,temprand,xrandom,yrandom,zrandom                   *
! * based on Wichman and Hill                                          *
! * begun GDJP 1/9/90                                                  *
! **********************************************************************

      subroutine random
      include 'MD.com'
      
      xrandom=171*(xrandom-177*NINT(xrandom/177))-2*NINT(xrandom/177)
      IF (xrandom .le. 0) xrandom = xrandom + 30269
      yrandom=172*(yrandom-176*NINT(yrandom/176))-35*NINT(yrandom/176)
      IF (yrandom .le. 0) yrandom = yrandom + 30307
      zrandom=170*(zrandom-178*NINT(zrandom/178))-63*NINT(zrandom/178)
      IF (zrandom .le. 0) zrandom = zrandom + 30323
      TEMPRAND = xrandom/30269.0+yrandom/30307.0+zrandom/30323.0
      rand = TEMPRAND - NINT(TEMPRAND)
	if (rand .lt. 0)then
	rand = rand + 1
	endif

      return

      end

! ^^^^^^^^^^^^^^^^^^^^^^^^^^END OF RANDOM^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! * subroutine gauss generates a random number from a gaussian         *
! * distribution with mean MM, and standard deviation SD               *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	subroutine gauss(MM,SD,RR)
        include 'MD.com'
	real tmp1, tmp2,MM,SD,RR,zeta
	dimension zeta(2)

10      do i=1,2
	call random
	zeta(i) = 2.0*(rand-0.5)
        end do

        tmp1 = zeta(1)*zeta(1) + zeta(2)*zeta(2)
        if( tmp1 .ge. 1.d0 .or. tmp1 .eq. 0.d0 ) goto 10
        tmp2 = sd*sqrt( -2.d0*log(tmp1)/tmp1 )
        RR =        zeta(1)*tmp2 + MM 

       end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
