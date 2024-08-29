!***********************************************************************
!*  Core Program by Paul Whitford <pcw@wpi.edu> April 2005             *
!*           Go.f                                                      *
!*                                                                     *
!*  This MD code runs simulations according to the values in an input  *
!*  that is written by Read.f.  settings.dat and random.dat are read in*
!*  for each simulation                                                *
!***********************************************************************
      PROGRAM DualGo
      include 'MD.com'
      real timee1, timee2, sum, sumc, DD1, DD2

      call start
      call init

      call integrate
      call stophere
     
      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF MAIN PROGRAM^^^^^^^^^^^^^^^^^^^^^^^^^




