!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* THREEBODYINIT: determines all of the possible three-body            *
!* interactions in the system                                          *
!***********************************************************************
      subroutine ThreeBodyInit
      include 'MD.com'

! find the contacts that can make triples
	tripi = 0
	do i = 1, NC
	do j = i+1, NC
	do k = j+1, NC

	if(IC(i) .eq. IC(j))then

	if((JC(i) .eq. IC(k) .and. JC(j) .eq. JC(k))
     Q .or. (JC(i) .eq. JC(k) .and. JC(j) .eq. IC(k))   )then
	tripi = tripi + 1
	Trip(1,tripi) = i
        Trip(2,tripi) = j
        Trip(3,tripi) = k

!       write(11,*) IC(i), JC(i), IC(j), JC(j), IC(k), JC(k)

	endif
	endif


        if(IC(i) .eq. JC(j))then
                                                                                                                             
        if((JC(i) .eq. IC(k) .and. IC(j) .eq. JC(k))
     Q .or. (JC(i) .eq. JC(k) .and. IC(j) .eq. IC(k))   )then
        tripi = tripi + 1
        Trip(1,tripi) = i
        Trip(2,tripi) = j
        Trip(3,tripi) = k
                                                                                                     
!       write(11,*) IC(i), JC(i), IC(j), JC(j), IC(k), JC(k)
                                                                                                                             
        endif
        endif
                                                                                                                             

        if(JC(i) .eq. IC(j))then
                                                                                                                             
        if((IC(i) .eq. IC(k) .and. JC(j) .eq. JC(k))
     Q .or. (IC(i) .eq. JC(k) .and. JC(j) .eq. IC(k))   )then
        tripi = tripi + 1
        Trip(1,tripi) = i
        Trip(2,tripi) = j
        Trip(3,tripi) = k
                                                                                                                             
                                                                                                                            
        endif
        endif

        if(JC(i) .eq. JC(j))then
                                                                                                                             
        if((IC(i) .eq. IC(k) .and. IC(j) .eq. JC(k))
     Q .or. (IC(i) .eq. JC(k) .and. IC(j) .eq. IC(k)) )then
        tripi = tripi + 1
        Trip(1,tripi) = i
        Trip(2,tripi) = j
        Trip(3,tripi) = k
                                                                                                                             
        endif
        endif 

	enddo
	enddo
	enddo
        !if(ContactFile .ne. 'NO')then
        !write(12,*) NC, ' 2body '
	!endif

        !if(ThreeBodyFile .ne. 'NO')then
        !write(13,*) tripi , ' 3body '
	!endif
      End

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* THREEBODY: computes the number of three body interactions in the    *
!* system                                                              *
!***********************************************************************

      subroutine ThreeBody(Qn,Qi)
      include 'MD.com'

      integer Qn, Qi, Contact, DE
      dimension  Contact(NC), Qi(NC)
      DE = 0

	do i=1, Qn
	Contact(i) = 0
	enddo

	do i=1, Qn
	Contact(Qi(i)) = 1
        enddo

! see if the triplets are formed.  If they are, add an energy contribution to the perturbation

	do i=1,tripi
	if(contact(trip(1, i))+contact(trip(2, i))+
     Q contact(trip(3, i)) .eq. 3)then
	DE = DE + 1
	endif

	enddo

        if(ContactFile .ne. 'NO')then
        write(12,*) Qn
        endif
                                                                                                                             
        if(ThreeBodyFile .ne. 'NO')then
        write(13,*) DE
        endif

      end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of Contacts^^^^^^^^^^^^^^^^^^^^^^^^^^^
