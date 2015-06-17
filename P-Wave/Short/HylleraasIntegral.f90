! Calculates the integral defined in the paper by Drake and Yan from 1995.

module FactPreComp
   implicit none
   real*16, dimension(0:401) :: FactArray  !@TODO: Any possible need to even go this high?
end module FactPreComp


real*16 function HylleraasIntegral(RunCalc, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, qmax, Method, WMatrix, CMatrix)
	use WCLimits
	implicit none
	common /globals/ PI
	logical RunCalc
	real*16 PI, T, AsymptoticPart
	integer j1, j2, j3, j12, j23, j31, qmax, Method
	real*16 alpha, beta, gamma
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*16, dimension(-2:Cjmax, 0:Cqmax, 0:Ckmax) :: CMatrix
	integer q, pmax, LambdaUpper, ZetaUpperBound
	real*16 LambdaLower
	real*16 Sum
	real*16, allocatable, dimension(:) :: TValues  ! We do not know the size of TValues yet (N+1).
	real*16 HylleraasIntegralRecursion

	! This is for the recursion relations from Pachucki's paper (an analytic solution).
	if (Method == 2) then
		HylleraasIntegral = HylleraasIntegralRecursion(j1, j2, j3, j12, j23, j31, &
				alpha, beta, gamma, qmax, Method)
		return
	endif
	
	! Initialize global variables
	!PI = 3.14159265358979323846q0
	!PI = 4.0q0 * datan(1.0q0)
	PI = 3.141592653589793238462643383279502884197q0

	! Initialize local variables
	pmax = 75
	LambdaUpper = 20
	ZetaUpperBound = 750
	Sum = 0.0q0

	! Allocate the TValues array
	allocate(TValues(qmax+1)) ! We have N+1 possibilities (0 to N).

	do q = 0, qmax, 1
		! We store these to calculate the asymptotic part.
		TValues(q+1) = T(RunCalc, q, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, pmax, WMatrix, CMatrix)
		! If the q summation terminates, then every term past this will be 0 as well.
!		if (TValues(q+1) == 0.0q0 .and. RunCalc == .TRUE.) then
!			exit
!		endif
		!write (*,*) 'T(', q, '): ', TValues(q+1)
	enddo

	! Do the direct summation
	Sum = 0.0q0
	do q = 0, qmax, 1
		! No need to sum more terms if this summation was finite.  Exiting early also
		!  allows us to do the check right after this do loop.
		if (TValues(q+1) == 0.0q0) then
			exit
		endif	
		Sum = Sum + TValues(q+1)
	enddo
	!write (*,*) 'N =', N, ':'
	! Output the direct sum.
	!write (*,*) 'Direct sum:', Sum

	! If the summation was finite, there is no need to do the asymptotic approximation.
	!  AsymptoticPart will return exactly 0 if it is called for a finite summation.
	if (q < qmax) then
		HylleraasIntegral = Sum
		return
	endif

	! Equation (13) for lowercase lambda
	LambdaLower = 0.5q0 * real(j12 + 1,16) + 0.5q0 * real(j23 + 1,16) + 0.5q0 * real(j31 + 1,16) + 4.0q0

	! Add on the asymptotic part
	if (Method == 1) then
		!Sum = Sum + AsymptoticPart(TValues, LambdaLower, qmax+1, LambdaUpper, ZetaUpperBound)
		Sum = Sum + AsymptoticPart(TValues, LambdaLower, qmax, LambdaUpper, ZetaUpperBound)
	endif

	! Output the final result, including the asymptotic part.
	!write (*,*) 'Final sum:', Sum
	!write (*,*)
	
	! Deallocate the TValues array before we end this routine.
	deallocate(TValues)

	HylleraasIntegral = Sum
	
	return
end


! This is a standard factorial function (not using recursion).
!  Only integer arguments are accepted.
!@TODO: Arguments past 400 will cause a crash. That is on the order of 10^868.
real*16 function Factorial(n)
	use FactPreComp
	implicit none
	integer n
	integer i

	! Right now, this is easier than erroring out.
	if (n < 0) then
		Factorial = 1.0q0  ! Just return a default value of 1
		write (*,*) 'n < 0 in factorial function'
		return
	endif
	
	Factorial = FactArray(n)
	return

	!Factorial = 1.0q0
	!do i = 1, n
	!	Factorial = Factorial * dble(i)
	!enddo
	!return
end


subroutine GenFactorial(MaxN)
	use FactPreComp
	implicit none
	real*16 Fact
	integer i, MaxN

	Fact = 1.0q0
	FactArray(0) = Fact
	do i = 1, MaxN
		Fact = Fact * real(i,16)
		FactArray(i) = Fact
	enddo
	
	return
end


! This is your standard "choose" function, or the binomial coefficient.
real*16 function Binomial(n, k)
	integer n, k
	real*16 Factorial

	Binomial = Factorial(n) / (Factorial(k) * Factorial(n-k))
	return
end


! S is given in equation (4) - just a simple minimum function.
!@TODO: Rewrite as returning an integer as in equation (13) of Perkins.
real*8 function S(q, j)
	integer q, j
	S = MIN(real(q-1,16), 0.5q0 * real(j+1,16))
	return
end

	
! C is given in equation (4).
real*16 function C(j, q, k)
	integer j, q, k
	real*16 Binomial
	real*8 S
	real*16 Multiplication, SVal, Bin
	integer t

	Bin = Binomial(j+2, 2*k + 1)
	Multiplication = 1.0q0
	SVal = S(q, j)
	
	do t = 0, int(SVal), 1
		Multiplication = Multiplication * (2*k + 2*t - j) / real(2*k + 2*q - 2*t + 1,16)
	enddo

	C = (2*q + 1) / real(j + 2,16) * Multiplication * Bin
	return
end
	

! This is the W function from equation (8).
real*16 function W(l, m, n, alpha, beta, gamma, pmax)
	implicit none
	integer l, m, n, b, c
	real*16 alpha, beta, gamma
	real*16 z, Hypergeometric2F1
	integer pmax
	real*16 Factorial
	real*16 WPartial
	real*16 LeadingFactor, Summation
	real*16 Hyper, HyperTest, Diff
	integer p

	LeadingFactor = Factorial(l) / ((alpha + beta + gamma) ** (l+m+n+3))
	Summation = 0.0q0

	! We do the summation backwards so that we can use the backwards recurrence relation in (9).
	!do p = 0, pmax, 1
	b = l+m+n+pmax+3
	c = l+m+pmax+3
	z = real((alpha + beta) / (alpha + beta + gamma),16)
	Hyper = Hypergeometric2F1(1.0q0, real(b,16), real(c,16), real(z,16))
	do p = pmax, 0, -1
		!write (*,*) p
		Summation = Summation + WPartial(l, m, n, alpha, beta, gamma, p, Hyper)
		b = b - 1;  c = c - 1;
		Hyper = 1.0q0 + real(b,16) / real(c,16) * z * Hyper
		!HyperTest = Hypergeometric2F1(1.0q0, real(l+m+n+p-1+3,8), real(l+m+p-1+3,8), real(z,8))
		!!HyperTest = Hypergeometric2F1(1.0q0, real(b+1,8), real(c+1,8), z)
		!Diff = abs((Hyper-HyperTest)/(Hyper+HyperTest)*2)
		!if (Diff > 1e-13) then
		!	write (*,*) Hyper, HyperTest, abs(Hyper-HyperTest), ":", b, c, z
		!	Hyper = HyperTest
		!end if
	enddo

	W = LeadingFactor * Summation
	
	return
end


! This is what's inside the summation of the W function from equation (8).
real*16 function WPartial(l, m, n, alpha, beta, gamma, p, hyper)
	integer l, m, n
	real*16 alpha, beta, gamma
	integer p
	real*16 Factorial
	real*16 Hypergeometric2F1, hyper
	
	WPartial = Factorial(l+m+n+p+2) / (Factorial(l+1+p) * (l+m+2+p)) * (alpha / (alpha + beta + gamma)) ** p  * &
		hyper
		!Hypergeometric2F1(1.0q0, real(l+m+n+p+3,8), real(l+m+p+3,8), real((alpha + beta) / (alpha + beta + gamma),8) )
		
	return
end


! The T function from equation (6)
real*16 function T(RunCalc, q, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, pmax, WMatrix, CMatrix)
	use WCLimits
	implicit none
	logical RunCalc
	integer q, j1, j2, j3, j12, j23, j31
	real*16 alpha, beta, gamma
	integer pmax
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*16, dimension(-2:Cjmax, 0:Cqmax, 0:Ckmax) :: CMatrix

	common /globals/ PI
	real*16 PI, W, C
	real*16 Summation
	real*8 L12, L23, L31
	integer k12, k23, k31

	real*16 Res1,Res2,Res3,Res4,Res5,Res6,Res
	real*16 Temp
	real*16 Rcoeff
	real*16 C12, C23, C31

	Summation = 0.0q0
	! L12, L23 and L31 are always even, because of the if statements following.
	L12 = (j12 + 1) / 2
	L23 = (j23 + 1) / 2
	L31 = (j31 + 1) / 2

	! If any of j12, j23 or j31 are even, then it changes the upper limit on its
	!  corresponding k12, k23 or k31 summation.  This also terminates the q summation
	!  as mentioned on page 1988 of the paper by J.F. Perkins (1967).
	if (mod(j12,2) == 0) then
		L12 = j12/2 - q
		if (q > j12 / 2) then
			T = 0.0q0
			return
		endif
	endif
	if (mod(j23,2) == 0) then
		L23 = j23/2 - q
		if (q > j23 / 2) then
			T = 0.0q0
			return
		endif			
	endif
	if (mod(j31,2) == 0) then
		L31 = j31/2 - q
		if (q > j31 / 2) then
			T = 0.0q0
			return
		endif			
	endif

	do k12 = 0, int(L12), 1  !@TODO: Determine if L12, L23 and L31 can be real.
		do k31 = 0, int(L31), 1
			do k23 = 0, int(L23), 1

				! This helps determine which elements of the WMatrix need to be calculated.
				if (RunCalc == .false.) then
					WMatrix(j1 + 2*q + 2*k12 + 2*k31 + 2, j2 + j12 - 2*k12 + 2*k23 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, 1) = 1000.0q0  ! The exact value doesn't matter, as
					WMatrix(j1 + 2*q + 2*k12 + 2*k31 + 2, j3 + j31 - 2*k31 + 2*k23 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, 2) = 1000.0q0  !  long as it is different from 0.
					WMatrix(j2 + 2*q + 2*k12 + 2*k23 + 2, j1 + j12 - 2*k12 + 2*k31 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, 3) = 1000.0q0
					WMatrix(j2 + 2*q + 2*k12 + 2*k23 + 2, j3 + j23 - 2*k23 + 2*k31 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, 4) = 1000.0q0
					WMatrix(j3 + 2*q + 2*k23 + 2*k31 + 2, j1 + j31 - 2*k31 + 2*k12 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, 5) = 1000.0q0
					WMatrix(j3 + 2*q + 2*k23 + 2*k31 + 2, j2 + j23 - 2*k23 + 2*k12 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, 6) = 1000.0q0
					CMatrix(j12, q, k12) = 1000.0q0
					CMatrix(j23, q, k23) = 1000.0q0
					CMatrix(j31, q, k31) = 1000.0q0
					Summation = 0.0q0
					cycle
				end if

				! The original version is preserved above.  Now this uses precalculated values for the W function via WMatrix.
				!  The ordering is the same as above, so that the 3 corresponds to beta, alpha, gamma and so on.
				Res1 = WMatrix(j1 + 2*q + 2*k12 + 2*k31 + 2, j2 + j12 - 2*k12 + 2*k23 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, 1)
				Res2 = WMatrix(j1 + 2*q + 2*k12 + 2*k31 + 2, j3 + j31 - 2*k31 + 2*k23 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, 2)
				Res3 = WMatrix(j2 + 2*q + 2*k12 + 2*k23 + 2, j1 + j12 - 2*k12 + 2*k31 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, 3)
				Res4 = WMatrix(j2 + 2*q + 2*k12 + 2*k23 + 2, j3 + j23 - 2*k23 + 2*k31 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, 4)
				Res5 = WMatrix(j3 + 2*q + 2*k23 + 2*k31 + 2, j1 + j31 - 2*k31 + 2*k12 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, 5)
				Res6 = WMatrix(j3 + 2*q + 2*k23 + 2*k31 + 2, j2 + j23 - 2*k23 + 2*k12 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, 6)

				C12 = CMatrix(j12, q, k12)
				C23 = CMatrix(j23, q, k23)
				C31 = CMatrix(j31, q, k31)
				Rcoeff = C12 * C31 * C23
				Summation = Summation + Rcoeff * (Res1+Res2+Res3+Res4+Res5+Res6)
			enddo
		enddo
	enddo

	T = ((4.0q0 * PI) ** 3) / real(2*q + 1,16) ** 2 * Summation

	return
end


real*16 function AsymptoticPart(TValues, LambdaLower, N, LambdaUpper, ZetaUpperBound)
	use LU
	implicit none
	integer N
	real*16 TValues(N)
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
	integer IPIV(N)
	integer ( kind = 4 ) nnew, rhs_num, infonew
	integer, allocatable, dimension(:) :: INDX
	integer rc, D

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
	q = N - (LambdaUpper+1)

	! Fill in our matrix, A, with the coefficients.  These coefficients are the q^(i+LambdaLower) in equation (19).
	!  Also notice that i and j are swapped from how they normally are in nested for loops.
	do j = 1, LambdaUpper+1, 1
		do i = 1, LambdaUpper+1, 1
			A(j,i) = 1.0q0 / real(q,16) ** ((i-1)+real(LambdaLower,16))
			A16(j,i) = 1.0q0 / real(q,16) ** ((i-1)+real(LambdaLower,16))
		enddo

		q = q + 1  ! We want all integral values for q up to N
	enddo

	!!!! Uses LAPACK to solve our system of LambdaUpper+1 linear equations for x.
	!!!!LaLinearSolve(A, x, b);
	!!!call dgesv(LambdaUpper+1, 1, A, LambdaUpper+1, IPIV, BX, LambdaUpper+1, Info)
	!!!do j = 1, LambdaUpper+1, 1
	!!!	do i = 1, LambdaUpper+1, 1
	!!!		AB16(i,j) = A16(i,j)
	!!!	end do
	!!!end do
	!!!do j = 1, LambdaUpper+1, 1
	!!!	AB16(j,LambdaUpper+2) = BX16(j)
	!!!end do
	!!!nnew = LambdaUpper+1
	!!!rhs_num = 1
	!!!call r16mat_solve(nnew, rhs_num, AB16, infonew)
	!!!!call NewDgesv(LambdaUpper+1, 1, A, LambdaUpper+1, IPIV, BX, LambdaUpper+1, Info)
	!!!!call NewDgesv(LambdaUpper+1, A, BX, Info)
	
	
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
		case (25)
			ModifiedRiemannZeta = Zeta25(int(2*i) - 2)
			return
		case (30)
			ModifiedRiemannZeta = Zeta30(int(2*i) - 2)
			return
		case (35)
			ModifiedRiemannZeta = Zeta35(int(2*i) - 2)
			return
		case (40)
			ModifiedRiemannZeta = Zeta40(int(2*i) - 2)
			return
		case (45)
			ModifiedRiemannZeta = Zeta45(int(2*i) - 2)
			return
		case (50)
			ModifiedRiemannZeta = Zeta50(int(2*i) - 2)
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

	! 
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
		if (Term < 1q-30) exit
		j = j + 1
	enddo

	ModifiedRiemannZeta = Sum
	return
end


real*16 function Hypergeometric2F1Sum(a, b, c, z)
	implicit none
	real*16 a, b, c, z
	real*16 Sum, Term, Prev, nd
	integer n, nUpper

	Sum = 1.0q0
	Prev = 0.0q0
	nUpper = 10000  ! Only for extremely singular 2F1 will it even reach close to this number.

	Term = 1.0q0
	do n = 1, nUpper, 1
		nd = n
		Term = Term * (a+nd-1.0q0) * (b+nd-1.0q0) / ((c+nd-1.0q0) * nd) * z
		Sum = Sum + Term
		!@TODO: Why did I have this here?
		!if (abs(Sum - Prev) <= 1d-17) then
		!!	write (*,*) 'Number of terms: ', n
		!	exit
		!endif
		Prev = Sum
	enddo

	Hypergeometric2F1Sum = Sum
end function Hypergeometric2F1Sum


!  Arguments b and c are assumed to be integers, and a is 1.
recursive function Hypergeometric2F1(a, b, c, z) result (hyper_result)
	implicit none
	real*16 a, b, c, z, s, t
	real*16 Hypergeometric2F1Sum, hyper_result
	real*16 Hypergeometric2F1First

	if (c > b) then
		hyper_result = Hypergeometric2F1Sum(a, b, c, z)
		!!!!hyper_result = Hypergeometric2F1Second(a, b, c, z)  !!! DO NOT USE !!!
		return
	end if
	
	!hyper_result = Hypergeometric2F1Sum(a, b, c, z)
	hyper_result = Hypergeometric2F1First(real(a,16), real(b,16), real(c,16), real(z,16))
	return
end function Hypergeometric2F1


! This is used when c <= b.
recursive function Hypergeometric2F1First(a, b, c, z) result (hyper_result)
	implicit none
	real*16 a, b, c, z, s, t, hyper_result
	
	if (c == 1) then
		hyper_result = (1-z)**(-b)
		return
	endif
	
	s = c-1
	t = b-1-s
	hyper_result = (Hypergeometric2F1First(a,b-1,c-1,z) - 1) * s/(z*(s+t))
end function Hypergeometric2F1First
