!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* IntLD() : integrates over one time step for all particles           *
!***********************************************************************

      subroutine intLD(E,ET, Conts, outE, count, step )
! E is the total potential energy
! ET is the total potential and kinetic energy 
      include 'MD.com'
      integer outE, Qi, SC1, SC2, count, cnt , step
      real DumE, DumKE, Res,FR,SD
      dimension DumE(9), Res(AN), count(NC)
      real RsBy

      E = 0.0

      do i=1, 9
         DumE(i)=0.0
      end do
! these routines compute the forces
      call Fstart 

      call Bonds( DumE(1))

      call ANGL(DumE(2))
 
      call dihedral(DumE(3))


      if (useChirals .eq. 'YES') then
       call chiral(DumE(8))
      end if


      call contacts(DumE(4), Conts, count, outE,step)

      call NonContacts( DumE(5))

      if (useEllipsoidRepulsions .eq. 'YES') then
        call ellipsoidRepulsions(DumE(9))
      endif

! calculate the box force
      if (confineInBox .eq. 'YES') then
        call box(X,Y,Z,Fx,Fy,Fz ,DynamicAtomRange, DynLength 
     Q       ,DumE(6),boxMin,boxMax,boxCoeff)
      end if

! calculate electrostatic force
      if (useElectrostatics .eq. 'YES') then
      ! E1,2 - residue indexes
      ! Q102 - charge of bond
      ! X,Y,Z - residue position
      !  FX, FY, FZ - residue applied force
      ! NE - number of electrostatic contacts
        if (useDebyeHuckel .eq. 'YES') then
!       write(7,*) 'esCharge(3) ',esCharge(3)
!       write(7,*) 'esCharge(76) ',esCharge(76)
          call debyehuckel(useESCutoff,useDHTable,esFirstAtomIndex,
     Q  esSecondAtomIndex,esCharge,X,Y,Z,Fx,Fy,Fz,esPairsNum,
     Q  esDistanceCutoff,DumE(7),deConstant,screeningFactor, 
     Q  saltCoefficient,DebyeHuckelPotentials,DebyeHuckelForces)
	else
          call coulomb(esFirstAtomIndex,esSecondAtomIndex,
     Q  esCharge,X,Y,Z,Fx,Fy,Fz,esPairsNum,DumE(7), deConstant)
	endif
      end if	

      
      Call Findtemp(DumKE)
      ET = E + DumKE

      if(outE .eq. 0)then
      E = DumE(1) + DumE(2) + DumE(3) +  DumE(4) + 
     Q    DumE(5) + DumE(6) + DumE(7) + DumE(8) + DumE(9)
        if(EnergyTerm .ne. 'NO')then
      write(70,"(9F9.2)") DumE(1),DumE(2),DumE(3),DumE(4),DumE(5),
     Q DumE(6),DumE(7),DumE(8), DumE(9)
	endif
      endif

	cnt = 1
	do while (cnt .le. DynLength)
           do i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)
		SD = sqrt(2*ms(i)*gamma*T/tau)

		call gauss(0., SD, fr)
		Vx(i) = (Vx(i)*c_e + (Fx(i)+fr)*tau/ms(i))*c_i
		call gauss(0., SD, fr)
		Vy(i) = (Vy(i)*c_e + (Fy(i)+fr)*tau/ms(i))*c_i
		call gauss(0., SD, fr)
		Vz(i) = (Vz(i)*c_e + (Fz(i)+fr)*tau/ms(i))*c_i
	   enddo
	   cnt = cnt + 2
	enddo


!      do i=1, AN

!      Fx(i) = X(i)
!      Fy(i) = Y(i)
!      Fz(i) = Z(i)

!      enddo

      do i=1, AN
      X(i) = X(i) + Vx(i)*tau
      Y(i) = Y(i) + Vy(i)*tau
      Z(i) = Z(i) + Vz(i)*tau
      enddo

!      do i=1, AN
!      X(i) = X(i) + tau*Vx(i)
!      Y(i) = Y(i) + tau*Vy(i)
!      Z(i) = Z(i) + tau*Vz(i)
!      enddo

!      do i=1, AN

!      Vx(i) = (X(i)

!      enddo


      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of IntSymp^^^^^^^^^^^^^^^^^^^^^^^^^^^^
