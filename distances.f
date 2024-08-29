!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Distances calculate the distance between any two beads
!***********************************************************************
      subroutine distances()
      include 'MD.com'

      do i=1, AN
        do j=i+1,AN
  	  dx = X(i) - X(j)
	  dy = Y(i) - Y(j)
	  dz = Z(i) - Z(j)
	  write(93,*) dx**2+dy**2+dz**2
        enddo
      enddo
      end
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END OF DISTANCES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
