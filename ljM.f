!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* CONTACTS: computes the force on all atoms due to contacts via a     *
!* 10-12 potential                                                     *
!***********************************************************************

        subroutine contacts(E, Conts, count, outE, step)
        include 'MD.com'

        integer C1, C2, ConfID, ContN, Q, SC1, SC2, Cf1, cf2, count(NC),
     Q summc, outE, step
        real r2, rm2, rm10, f_over_r, dsig, deps, Energy, Res, ResCon,
     Q s1, s2, ep1, ep2, r1, rc, r, summm, currentEnergy
        dimension Energy(2), ContN(2), ResCon(AN), Res(AN)
        integer isContact, currentRangesIndex, currentContactInRange
        dimension currentRangesIndex(50)

        if (writeContactsRanges .eq. 'YES') then
        do i = 1, rangesNumber
            currentRangesIndex(i) = 1
            EtotalRanges(i) = 0.0
            TwoBodyRanges(i) = 0
        enddo
        end if

        Q = 0
        E = 0.0
        Conts = 0

        do i = 1, NC
            count(i) = 0
        end do

        do i = 1, NC
        C1 = IC(i)
        C2 = JC(i)
        dx = X(C1) - X(C2)
        dy = Y(C1) - Y(C2)
        dz = Z(C1) - Z(C2)
        r2 = dx**2 + dy**2 + dz**2
        rm2 = 1.0 / r2
        rm2 = rm2 * sigma(i)
        rm10 = rm2**5

        ! Determine the type of interaction based on residue IDs and charges
        if (changeHisCharge .eq. 'YES') then
             ! Histidine interactions will be updated soon :)
            ! update energy
        else
            ! Default case without histidine charge considerations
            currentEnergy = epsC(i) * rm10 * (5 * rm2 - 6)
            f_over_r = epsC(i) * rm10 * (rm2 - 1) * 60 / r2
        endif

        E = E + currentEnergy

        if (r2 <= sigma(i) * ConCut**2) then
        Conts = Conts + 1
        count(Conts) = i
        isContact = 1
        else
        isContact = 0
        endif


        if (writeAllContacts .eq. 'YES' .and. outE .eq. 0) then
            write(14,*) isContact
        endif

        if (writeContactsRanges .eq. 'YES') then
            do j = 1, rangesNumber
                currentContactInRange = rangeContacts(j, currentRangesIndex(j))
                if (i .eq. currentContactInRange) then
                    EtotalRanges(j) = EtotalRanges(j) + currentEnergy
                    TwoBodyRanges(j) = TwoBodyRanges(j) + isContact
                    currentRangesIndex(j) = currentRangesIndex(j) + 1
                endif
            enddo
        endif

        ! update forces
        FX(C1) = Fx(C1) + f_over_r * dx
        FY(C1) = Fy(C1) + f_over_r * dy
        FZ(C1) = Fz(C1) + f_over_r * dz

        Fx(C2) = Fx(C2) - f_over_r * dx
        Fy(C2) = Fy(C2) - f_over_r * dy
        Fz(C2) = Fz(C2) - f_over_r * dz
        enddo

        if (writeContactsRanges .eq. 'YES' .and. outE .eq. 0) then
        do j = 1, rangesNumber
            write(15, '(I5)') TwoBodyRanges(j)
            write(16, '(F9.2)') EtotalRanges(j)
        enddo
        write(15,*) ''
        write(16,*) ''
        end if

        end subroutine

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of Contacts^^^^^^^^^^^^^^^^^^^^^^^^^^^
