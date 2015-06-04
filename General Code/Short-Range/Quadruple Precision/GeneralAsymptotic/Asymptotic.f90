
real*16 function AsymptoticPart(TValues, LambdaLower, N1, LambdaUpper)
	use LU
	implicit none
	integer N, N1
	real*16 TValues(N1)
	real*16 LambdaLower, Temp
	integer LambdaUpper, ZetaUpperBound, num
	
	external DGESV, ModifiedRiemannZeta
	real*16 ModifiedRiemannZeta

	real*16 AsymptoticSum, AsymptoticSum16
	integer q, i, j, Info
	real*8 A(LambdaUpper+1, LambdaUpper+1)
	real*8 BX(LambdaUpper+1)
	real*16 A16(LambdaUpper+1, LambdaUpper+1)
	real*16 BX16(LambdaUpper+1)
	real*16 AB16(LambdaUpper+1,LambdaUpper+2)
	integer IPIV(N1)
	integer ( kind = 4 ) nnew, rhs_num, infonew
	integer, allocatable, dimension(:) :: INDX
	integer rc, D
	
	N = N1

    ZetaUpperBound = 10000
    
	! Initialize the summation
	AsymptoticSum = 0.0q0
	AsymptoticSum16 = 0.0q0
	A = 0.0d0
	A16 = 0.0q0

	!// Fill in x vector  with the values of T(q) for q = N-LambdaUpper to N.
	do i = 1, LambdaUpper+1, 1
		BX(i) = TValues(i + (N - (LambdaUpper+1)))
		BX16(i) = TValues(i + (N - (LambdaUpper+1)))
	enddo

	!// This is our lower value for q, q = N-LambdaUpper.
	q = N - (LambdaUpper)

	! Fill in our matrix, A, with the coefficients.  These coefficients are the q^(i+LambdaLower) in equation (19).
	!  Also notice that i and j are swapped from how they normally are in nested for loops.
	do j = 1, LambdaUpper+1, 1
		do i = 1, LambdaUpper+1, 1
			A(j,i) = 1.0q0 / real(q,16) ** ((i-1)+real(LambdaLower,16))
			A16(j,i) = 1.0q0 / real(q,16) ** ((i-1)+real(LambdaLower,16))
		enddo

		q = q + 1  ! We want all integral values for q up to N
	enddo

	! LU decomposition
	num = LambdaUpper + 1
	A = A16
	allocate(INDX(num))
	!call LU decomposition routine
	call LUDCMP(A16,num,INDX,D,rc)
	!call appropriate solver if previous return code is ok
	if (rc.eq.0) then
		call LUBKSB(A16,num,INDX,BX16)
	endif
	deallocate(INDX)

	do i = 1, LambdaUpper+1, 1
		!AsymptoticSum = AsymptoticSum + BX(i) * ModifiedRiemannZeta(N-1, real((i-1)+LambdaLower,16), ZetaUpperBound)
		!Temp = ModifiedRiemannZeta(N, real((i-1)+LambdaLower,16), ZetaUpperBound)
		AsymptoticSum = AsymptoticSum + BX16(i) * ModifiedRiemannZeta(N, real((i-1)+LambdaLower,16), ZetaUpperBound)
		!AsymptoticSum16 = AsymptoticSum16 + AB16(i,LambdaUpper+2) * ModifiedRiemannZeta(N, real((i-1)+LambdaLower,16), ZetaUpperBound)
	enddo

	!write (*,*) 'Asymptotic part:', AsymptoticSum
	
	AsymptoticPart = AsymptoticSum
	!AsymptoticPart = AsymptoticSum16
	return
end
	

! This is the Riemann zeta function shown in (18) with its first N terms removed.
real*16 function ModifiedRiemannZeta(N, i, ZetaUpperBound)
	implicit none
	include 'ModifiedZeta.h'
	integer N
	real*16 i
	integer ZetaUpperBound
	real*16 Sum, Term, Zeta2
	integer j
	logical, save :: FirstCall = .true.

	select case (N)
		case (5)
			ModifiedRiemannZeta = Zeta5(int(2*i) - 2)
			return
		case (10)
			ModifiedRiemannZeta = Zeta10(int(2*i) - 2)
			return
		case (15)
			ModifiedRiemannZeta = Zeta15(int(2*i) - 2)
			return
		case (20)
			ModifiedRiemannZeta = Zeta20(int(2*i) - 2)
			return
		case (21)
			ModifiedRiemannZeta = Zeta21(int(2*i) - 2)
			return
		case (25)
			ModifiedRiemannZeta = Zeta25(int(2*i) - 2)
			return
		case (26)
			ModifiedRiemannZeta = Zeta26(int(2*i) - 2)
			return
		case (30)
			ModifiedRiemannZeta = Zeta30(int(2*i) - 2)
			return
		case (31)
			ModifiedRiemannZeta = Zeta31(int(2*i) - 2)
			return
		case (35)
			ModifiedRiemannZeta = Zeta35(int(2*i) - 2)
			return
		case (36)
			ModifiedRiemannZeta = Zeta36(int(2*i) - 2)
			return
		case (40)
			ModifiedRiemannZeta = Zeta40(int(2*i) - 2)
			return
		case (41)
			ModifiedRiemannZeta = Zeta41(int(2*i) - 2)
			return
		case (45)
			ModifiedRiemannZeta = Zeta45(int(2*i) - 2)
			return
		case (46)
			ModifiedRiemannZeta = Zeta46(int(2*i) - 2)
			return
		case (50)
			ModifiedRiemannZeta = Zeta50(int(2*i) - 2)
			return
		case (51)
			ModifiedRiemannZeta = Zeta51(int(2*i) - 2)
			return
		case (55)
			ModifiedRiemannZeta = Zeta55(int(2*i) - 2)
			return
		case (60)
			ModifiedRiemannZeta = Zeta60(int(2*i) - 2)
			return
		case (65)
			ModifiedRiemannZeta = Zeta65(int(2*i) - 2)
			return
		case (70)
			ModifiedRiemannZeta = Zeta70(int(2*i) - 2)
			return
		case (75)
			ModifiedRiemannZeta = Zeta75(int(2*i) - 2)
			return
	end select

	! This is to make it explicitly clear that the hardcoded values are more accurate.
	if (FirstCall .eqv. .true.)	then
		write (*,*) 'NOTE! Hardcoded values for the modified Riemann Zeta function are not being used.'
		write (*,*) ' The accuracy of the results will be less than if hardcoded values are used.'
		FirstCall = .false.
	end if
		
	Sum = 0.0q0

	! For the special case of i == 2, we have a p-series with a definite sum (Zeta2 below).
	!  This series converges *VERY* slowly, so we take it as a special case.
	if (int(i) == 2) then
		Zeta2 = 1.644934066848226436472415166646025189219q0
		do j = 1, N, 1
			Zeta2 = Zeta2 - 1.0q0 / real(j**i,16)
		enddo
		ModifiedRiemannZeta = Zeta2
		return
	endif

	j = N+1
	!do j = N+1, ZetaUpperBound, 1
	do
		Term = 1.0q0 / (real(j,16) ** real(i,16))
		Sum = Sum + Term
		!@TODO: How small do we really need to get this?
		if (Term < 1q-21) exit
		j = j + 1
	enddo

	ModifiedRiemannZeta = Sum
	return
end

