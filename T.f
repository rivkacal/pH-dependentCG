!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Findtemp() finds the temperature of the system                       *
!  remark for DNA                                                                 *
! here there are lines which remove angular and center of mass moment  *
! to avoid that the system rotates of moves (for now I left it as it   *
! was)                                                                 *
!***********************************************************************

      subroutine Findtemp(KE)
      include 'MD.com'

	real  P2, pxsum, pysum, pzsum, Px, Py, Pz, Mass, LX, LY, LZ,
     Q IYZ, IXZ, IXY, WX, WY, WZ, CMX, CMY, CMZ
! NOA - changed form Px(lastDynamicAtom) => Px(DynamicAtomRange(DynLength)), same for Py and Pz
        dimension Px(DynamicAtomRange(DynLength)),  
     Q Py(DynamicAtomRange(DynLength)), Pz(DynamicAtomRange(DynLength))
	integer cnt

        real ICut

        data ICut /.1/

! NOA - the loop was: do i=1, lastDynamicAtom; i've changed it so it will run only on the Dynamic atom
	cnt = 1
	do while (cnt .le. DynLength)
		do i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
			Px(i) = Vx(i)*ms(i)
			Py(i) = Vy(i)*ms(i)
			Pz(i) = Vz(i)*ms(i)
		end do
		cnt = cnt + 2
	enddo
        
! don't reset center of mass momentums is ther are static atoms       
        if ((hasStaticAtoms .ne. 'YES') 
     Q      .and. (confineInBox .ne. 'YES')) then

!        goto 11111
! These next few steps reset center of mass momentum to zero
! This computes the total linear momentum
      
	pxsum=0.0
	pysum=0.0
	pzsum=0.0
        Mass = 0.0

! NOA - the loop was: do 1000 i=1, lastDynamicAtom; i've changed it so it will run only on the Dynamic atom
	cnt = 1
	do while (cnt .le. DynLength)
	   do 1000 i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
		pxsum = Px(i) + pxsum
		pysum = Py(i) + pysum
		pzsum = Pz(i) + pzsum
		Mass = Mass + ms(i)
1000	   continue
	   cnt = cnt + 2
	enddo

! This finds the average linear momementum
	pxsum = pxsum/Mass
	pysum = pysum/Mass
	pzsum = pzsum/Mass
! This brings the center of mass linear momentum to zero
! NOA - the loop was: do 2000 i=1, lastDynamicAtom; i've changed it so it will run only on the Dynamic atom
	cnt = 1
	do while (cnt .le. DynLength)
	   do 2000 i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
	 	Px(i) = Px(i) - pxsum*ms(i)
	 	Py(i) = Py(i) - pysum*ms(i)
	 	Pz(i) = Pz(i) - pzsum*ms(i)
                VX(i) = PX(i)/ms(i)
                VY(i) = PY(i)/ms(i)
                VZ(i) = PZ(i)/ms(i)
2000	   continue
	   cnt = cnt + 2
	enddo


! These next lines remove the angular momentum

! First center the molecule about the origin


        LX = 0.0
        LY = 0.0
        LZ = 0.0

        IXY = 0.0
        IYZ = 0.0
        IXZ = 0.0

        CMX = 0.0
        CMY = 0.0
        CMZ = 0.0

! NOA - the loop was: do i=1, lastDynamicAtom; i've changed it so it will run only on the Dynamic atom
	cnt = 1
	do while (cnt .le. DynLength)
           do i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
		CMX = CMX + X(i)*ms(i)
		CMY = CMY + Y(i)*ms(i)
		CMZ = CMZ + Z(i)*ms(i)
	   end do
	   cnt = cnt + 2
	enddo

        CMX = CMX/Mass
        CMY = CMY/Mass
        CMZ = CMZ/Mass

! NOA - the loop was: do i=1, lastDynamicAtom; i've changed it so it will run only on the Dynamic atom
	cnt = 1
	do while (cnt .le. DynLength)
           do i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
		X(i) = X(i) - CMX
		Y(i) = Y(i) - CMY
		Z(i) = Z(i) - CMZ
	   end do
	   cnt = cnt + 2
	enddo

! NOA - the loop was: do i=1, lastDynamicAtom; i've changed it so it will run only on the Dynamic atom
	cnt = 1
	do while (cnt .le. DynLength)
           do i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
		LX = LX + Y(i)*PZ(i)-Z(i)*PY(i)
		LY = LY + PX(i)*Z(i)-X(i)*PZ(i)
		LZ = LZ + X(i)*PY(i)-Y(i)*PX(i)

		IXY = IXY + ms(i)*(X(i)**2 + Y(i)**2)
		IXZ = IXZ + ms(i)*(X(i)**2 + Z(i)**2)
		IYZ = IYZ + ms(i)*(Y(i)**2 + Z(i)**2)

           end do
	   cnt = cnt + 2
	enddo

!        if(IYZ .lt. Icut)then
!           IYZ = Icut
!        endif

!        if(IXZ .lt. Icut)then
!           IXZ = Icut
!        endif

!        if(IXY .lt. Icut)then
!           IXY = Icut
!        endif

        WX = LX/IYZ
        WY = LY/IXZ
        WZ = LZ/IXY


! sutract out each L

! NOA - the loop was: do i=1, lastDynamicAtom; i've changed it so it will run only on the Dynamic atom
	cnt = 1
	do while (cnt .le. DynLength)
           do i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
		Vx(i) = VX(i) + WZ*Y(i) - WY*Z(i)
		VY(i) = VY(i) - X(i)*WZ + WX*Z(i)
		VZ(i) = VZ(i) + X(i)*WY - WX*Y(i)
	
		PX(i) = VX(i)*ms(i)
		PY(i) = VY(i)*ms(i)
		PZ(i) = VZ(i)*ms(i)
           end do
	   cnt = cnt + 2
	enddo

!End of L subtraction

	endif !has static atoms
	
	P2=0
! NOA - the loop was: do 33 i=1, lastDynamicAtom; i've changed it so it will run only on the Dynamic atom
	cnt = 1
	do while (cnt .le. DynLength)
	   do 33 i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
	       P2 =  P2 + (Px(i)**2 + Py(i)**2 + Pz(i)**2)/ms(i)
33	   continue
	   cnt = cnt + 2
	enddo

! temp is the temp of the system in reduced units
! there should be a 2 in the numerator and 2 in the demoninator.

! NOA - changed from lastDynamicAtom to numDyn (the number of Dynamic Atoms)
	temp = P2/((numDyn-2)*3)
        KE = P2/2.0

	return

	end

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of Findtemp^^^^^^^^^^^^^^^^^^^^^^^^^^
