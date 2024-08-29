	subroutine genDynRange(DynLength,DynamicStr,DynamicAtomRange)
	integer DynLength, DynamicAtomRange(DynLength), n, m, a2i
	character (1000) :: DynamicStr, modDynamicStr, tempStr

	do n=1,DynLength
		DynamicAtomRange(n)=0
	enddo

	modDynamicStr = DynamicStr
	do n=1,DynLength
		m = INDEX(modDynamicStr,',')
		if (m == 0) then
			DynamicAtomRange(n) = a2i(modDynamicStr)
		else
			tempStr = modDynamicStr(:m-1)
			modDynamicStr = modDynamicStr(m+1:)
			DynamicAtomRange(n) = a2i(tempStr)
		endif
	enddo
	end








	integer function isDyn(DynamicAtomRange,DynLength,AtomNum)
	integer DynLength, DynamicAtomRange(DynLength), AtomNum, cnt

	isDyn = 0
	cnt = 1
	do while (cnt .le. DynLength)
	   if (AtomNum .ge. DynamicAtomRange(cnt) .AND.
     Q         AtomNum .le. DynamicAtomRange(cnt+1)) then
		isDyn = 1
		return
	   endif
	   cnt = cnt + 2
	enddo

	return
	END

	integer function front_trim(buf)
		IMPLICIT NONE
		CHARACTER(1000) :: buf
		INTEGER n, lng
		lng = LEN_TRIM(buf)
		DO n = 1, lng
			IF (buf(n:n) .NE. ' ') THEN
				front_trim = n
				RETURN
			END IF
		END DO
	END

	integer function a2i(buf)
	IMPLICIT NONE
	CHARACTER(1000) :: buf
	INTEGER n, foffs, lng
	INTEGER front_trim
	LOGICAL :: neg
	neg = .FALSE.
	lng = LEN_TRIM(buf)
	foffs = front_trim(buf)
	a2i = 0
	DO n = foffs, lng
		IF(buf(n:n) .EQ. '-') THEN
			neg = .TRUE.
			CONTINUE
		END IF
		IF(buf(n:n) .GE. '0' .AND. buf(n:n) .LE. '9' .AND. 
     Q			buf(n:n) .NE. ' ') THEN
			a2i = a2i * 10
			a2i = a2i + (IACHAR(buf(n:n)) - 48)
		END IF
	END DO
	IF(neg) THEN
		a2i = a2i * (-1)
	END IF
	return
	END 

!      The function calculate the number of the Dynamic Atoms
	integer function numDynAtom(DynamicAtomRange,DynLength)
	integer cnt, DynLength, DynamicAtomRange(DynLength)

	cnt = 1
	numDynAtom = 0
	do while (cnt .le. DynLength)	
		numDynAtom = numDynAtom + 
     Q          (DynamicAtomRange(cnt+1)-DynamicAtomRange(cnt) +1)
		cnt = cnt + 2
	enddo
	return
	end
