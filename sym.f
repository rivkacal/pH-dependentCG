!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* IntSymp() : integrates over one time step for all particles         *
!***********************************************************************

      subroutine intsymp(E,ET, Conts, outE, count )
! E is the total potential energy
! ET is the total potential and kinetic energy 
      include 'MD.com'
      integer outE, Qi, SC1, SC2, count, cnt 
      real DumE, DumKE, Res
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


      call contacts(DumE(4), Conts, count, outE)

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
        if (useDebyeHuckel .eq. 'YES') then
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
     Q DumE(5) + DumE(6) + DumE(7) + DumE(8) + DumE(9)
         if(EnergyTerm .ne. 'NO')then
           write(70,"(9F9.2)") DumE(1),DumE(2),DumE(3),DumE(4),DumE(5),
     Q DumE(6),DumE(7), DumE(8), DumE(9)
         endif
       endif

! Thermostat
      RsBy =sqrt( 1+(T/temp-1)*tau/Rstau )   
!      RsBy = 1.0 ! this would turn off the thermostat
	cnt = 1
	do while (cnt .le. DynLength)
           do 2000 i=DynamicAtomRange(cnt), DynamicAtomRange(cnt+1)

! Momentum(V in reduced units) change of particle i
! Insert division by mass.
!	goto 2000
		Vx(i) = (Vx(i) + tau*Fx(i)/ms(i))*RsBy
		Vy(i) = (Vy(i) + tau*Fy(i)/ms(i))*RsBy
		Vz(i) = (Vz(i) + tau*Fz(i)/ms(i))*RsBy
! displacements of particle i
		X(i) = X(i) + tau*Vx(i)
		Y(i) = Y(i) + tau*Vy(i)
		Z(i) = Z(i) + tau*Vz(i)
2000       continue
	   cnt = cnt + 2
	enddo

      END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of IntSymp^^^^^^^^^^^^^^^^^^^^^^^^^^^^
