!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* Integrate performs the integration the desired repeated number of   *
!* times                                                               *
!***********************************************************************

      subroutine integrate
      include 'MD.com'
      
	integer writeout, writeoutT,Qdummy, Edum, writei, 
     Q QQ, SC1, SC2, Tcount, count, step, minimalExchFequency
	real  taucum, step1, step2, vxsum, vysum, vzsum, ttime, t1, t2,
     Q dummyy, Tsum, Res
        character(LEN=15) Form
        dimension Res(AN), count(NC)

        writei = 0
        
        form = "(I5,TR1,F12.5)"

	call CPU_Time(t1)

        if(Trajectory .ne. 'NO')then
        write(6,*) AN
	endif
        writecount = 0
        writeout = WO
        writeoutT = WOT
        Tsum = 0.0
        Tcount = 0
        Tcount = 0

!       call noncontactsnear to find the non native non bonded list 
	call noncontactsnear
         if (useEllipsoidRepulsions .eq. 'YES') then
	   call ellipsoidRepulsionsCheck
         endif
      do 999 step=1, stepstop

         writecount = writecount + 1

            Edum= writeout - step

	if(mod(step,50) .eq. 0)then
         call noncontactsnear
         if (useEllipsoidRepulsions .eq. 'YES') then
	   call ellipsoidRepulsionsCheck
         endif
	endif

! Call 1 timestep integration
	if(symtype .eq. 'MD')then
	  call intSymp(E ,ET, Conts, Edum, count)
	elseif(symtype .eq. 'LD')then
	  call intLD(E,ET,Conts,Edum,count,step)
	else
	write(*,*) 'Must select dynamics type!'
	call abort
	endif

	  call findtemp(dummyy)

          if(Edum .eq. 0)then
! Find and return the Three-body interactions
	  !call ThreeBody(conts, count)
        if(EnergyTot .ne. 'NO')then
	  write(99,*) E
	endif
! output temperature
          write(58,*) step, temp 

          Tsum = Tsum + temp
          Tcount = Tcount + 1
          writeout = writeout + WO
          
          endif

          if(writeoutT .eq. step)then
! writei indicates to write to restart file 1, or 2
             call writetofile(writei)

             if(writei .eq. 0)then
                writei = 1
             else
                writei = 0
             endif
         
        if(Trajectory .ne. 'NO')then
          
	  if (TrajDist .ne. 'NO') then
            if(mod(step,WOT) .eq. 0)then
	      write(93,*) step
	      call distances()
            endif
          end if

            if(stepstop-step .ge. WOT)then

            write(6,*) 'continue'

            else
            write(6,*) 'end'
            close(7, status = 'KEEP')  
            endif
        endif
            writeoutT = writeoutT + WOT

          endif

	  if (step - WO .eq. 0)then
	  call CPU_time(t2)
	  ttime = (t2-t1)*(stepstop-step)/step
      open(999, file = 'timeleft.dat', access = 'SEQUENTIAL',
     Q   status = 'unknown')
	  write(999,*) 'Expected time remaining'
	  write(999,*)  NINT(ttime/3600), 'hours', 
     Q ttime/60-60*nint(ttime/3600), 'mins'
	close(999, status= 'keep')
	  endif

999	  continue

          write(58,*) Tsum/Tcount

      	Call Cpu_time(t2)
        ttime = t2-t1
        open(91, file='time.dat', status='unknown')
	write(91,*) NINT(ttime/3600), 'hours', 
     Q ttime/60-60*nint(ttime/3600), 'mins'
        close(91)

	close(8, status = 'keep')
        ! call writetofile(writei) !, Res)

        close(10)
        close(11)


	END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end of Integrate^^^^^^^^^^^^^^^^^^^^^^^^
