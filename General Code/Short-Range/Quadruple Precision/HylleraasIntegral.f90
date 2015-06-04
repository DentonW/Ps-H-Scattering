! Calculates the integral defined in the paper by Drake and Yan from 1995.
!  Started conversion to Fortran on 06/23/2009
!  Updated 08/02/2009 to use HYGFX subroutine for calculating 2F1 hypergeometric function.
!  Modified 08/02/2009 to be used in Ps-H integral calculations.
!  Modified 08/16/2009 to use precomputed values for the W function through WMatrix.
!  Modified 08/19/2009 to skip asymptotic expansion for any j_{mu nu} even.
!  Modified 09/05/2009 to move 1/(2q+1)**2 out of inner loop in T.
!  Somewhere in here, I got this to produce the same as PVR's code and then made some updates.
!  Updated 09/22/2009 to use ModifiedRiemannZeta from Drake and Yan '97 codebase.

! TODO: Rewrite array limits (TValues, etc.)
! TODO: Do we really need to reallocate TValues every time?
! TODO: Do we need pmax in HylleraasIntegral?

module FactPreComp
   implicit none
   real*16, dimension(0:401) :: FactArray  !@TODO: Any possible need to even go this high?
end module FactPreComp


real*16 function HylleraasIntegral(j1, j2, j3, j12, j23, j31, alpha, beta, gamma, qmax, WMatrix, Method)
	use WLimits
	implicit none
	common /globals/ PI, HyperB, HyperC, HyperVal
	real*16 PI, HyperB, HyperC, HyperVal
	external T, AsymptoticPart
	real*16 T, AsymptoticPart
	integer j1, j2, j3, j12, j23, j31, qmax, Method
	real*16 alpha, beta, gamma
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	integer q, pmax, LambdaUpper, ZetaUpperBound
	real*16 LambdaLower
	real*16 Sum
	real*16, allocatable, dimension(:) :: TValues  ! We do not know the size of TValues yet (N+1).
	real*16 HylleraasIntegralRecursion
	
	! This is for the recursion relations from Pachucki's paper (an analytic solution).
	if (Method == 2) then
		HylleraasIntegral = HylleraasIntegralRecursion(j1, j2, j3, j12, j23, j31, &
				real(alpha,16), real(beta,16), real(gamma,16), qmax, Method)
		return
	endif
	
	! Initialize global variables
	PI = 3.141592653589793238462643383279502884197q0
	HyperB = 0.0q0
	HyperC = 0.0q0
	HyperVal = 0.0q0

	! Initialize local variables
	pmax = 75
	!pmax = 65
	LambdaUpper = 6
	ZetaUpperBound = 500
	Sum = 0

	! Allocate the TValues array
	allocate(TValues(qmax+1)) ! We have N+1 possibilities (0 to N).

	do q = 0, qmax, 1
		! We store these to calculate the asymptotic part.
		TValues(q+1) = T(q, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, pmax, WMatrix)
		! If the q summation terminates, then every term past this will be 0 as well.
		if (TValues(q+1) == 0.0q0) then
			exit
		endif
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
		Sum = Sum + AsymptoticPart(TValues, LambdaLower, qmax+1, LambdaUpper)
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
! To keep the output of Calc3j correct, this needs to return 0 whenever we have a negative value for n.
! Note: Arguments past 400 will cause a crash.
real*16 function Factorial(n)
	use FactPreComp
	implicit none
	integer n
	integer i

	if (n < 0) then
		Factorial = 0.0q0
		!write (*,*) 'n < 0 in factorial function'
		return
	endif
	
	Factorial = FactArray(n)
	return

	!Factorial = 1.0q0
	!do i = 1, n
	!	Factorial = Factorial * real(i, 16)
	!enddo
	!return
end


subroutine GenFactorial(MaxN)
	use FactPreComp
	implicit none
	real*16 Fact
	integer i, MaxN

	Fact = 1.0_16
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
! TODO: Rewrite as returning an integer as in equation (13) of Perkins.
real*16 function S(q, j)
	integer q, j
	S = MIN(real(q-1, 16), 0.5q0 * real(j+1, 16))
	return
end

	
! C is given in equation (4).
real*16 function C(j, q, k)
	integer j, q, k
	real*16 S, Binomial
	real*16 Multiplication, SVal, Bin
	integer t

	Bin = Binomial(j+2, 2*k + 1)
	Multiplication = 1.0q0
	SVal = S(q, j)
	
	do t = 0, int(SVal), 1
		Multiplication = Multiplication * (2*k + 2*t - j) / real(2*k + 2*q - 2*t + 1, 16)
	enddo

	C = (2*q + 1) / real(j + 2, 16) * Multiplication * Bin
	return
end
	

! This is the W function from equation (8).
real*16 function W(l, m, n, alpha, beta, gamma, pmax)
	implicit none
	integer l, m, n, b, c
	real*16 alpha, beta, gamma
	real*16 z, Hypergeometric2F1
	integer pmax
	real*16 Factorial, WPartial
	real*16 LeadingFactor, Summation
	real*16 Hyper, HyperTest, Diff
	integer p
	integer, save :: LMin1 = 1000000, LMax1 = -1000000, MMin1 = 1000000, MMax1 = -1000000, NMin1 = 1000000, NMax1 = -1000000
	
	if (l < LMin1) LMin1 = l
	if (l > LMax1) LMax1 = l
	if (m < MMin1) MMin1 = m
	if (m > MMax1) MMax1 = m
	if (n < NMin1) NMin1 = n
	if (n > NMax1) NMax1 = n
	write (*,*) LMin1, LMax1, MMin1, MMax1, NMin1, NMax1
	
	LeadingFactor = Factorial(l) / ((alpha + beta + gamma) ** (l+m+n+3))
	Summation = 0.0q0

	! We do the summation backwards so that we can use the backwards recurrence relation in (9).
	!do p = 0, pmax, 1
	b = l+m+n+pmax+3
	c = l+m+pmax+3
	z = real((alpha + beta) / (alpha + beta + gamma),16)
	Hyper = Hypergeometric2F1(1.0_16, real(b,16), real(c,16), real(z,16))
	do p = pmax, 0, -1
		!write (*,*) p
		Summation = Summation + WPartial(l, m, n, alpha, beta, gamma, p, Hyper)
		b = b - 1;  c = c - 1;
		Hyper = 1.0_16 + real(b,16) / real(c,16) * z * Hyper
		!HyperTest = Hypergeometric2F1(1.0_16, real(l+m+n+p-1+3,16), real(l+m+p-1+3,16), real(z,16))
		!!HyperTest = Hypergeometric2F1(1.0_16, real(b+1,16), real(c+1,16), z)
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
		!Hypergeometric2F1(1.0_16, real(l+m+n+p+3,16), real(l+m+p+3,16), real((alpha + beta) / (alpha + beta + gamma),16) )
		
	return
end


! The T function from equation (6)
real*16 function T(q, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, pmax, WMatrix)
	use WLimits
	implicit none
	integer q, j1, j2, j3, j12, j23, j31
	real*16 alpha, beta, gamma
	integer pmax
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix

	common /globals/ PI, HyperB, HyperC, HyperVal
	real*16 PI, HyperB, HyperC, HyperVal
	real*16 W, C
	real*16 Summation
	real*16 L12, L23, L31
	integer k12, k23, k31

	real*16 RemoveMe1, RemoveMe2
	real*16 Res1,Res2,Res3,Res4,Res5,Res6,Res
	real*16 Temp
real*16 Rqt, Piterm
real*16 Rcoeff, C12, C23, C31

	Summation = 0.0q0
	! TODO: L12, L23 and L31 are always even, because of the if statements following.
!	L12 = .5 * (j12 + 1)  ! Need to convert each parentheses to real?
!	L23 = .5 * (j23 + 1)
!	L31 = .5 * (j31 + 1)
	L12 = (j12 + 1) / 2  ! Need to convert each parentheses to real?
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

	Rqt = (4.0q0*PI)**3 / ((2.0q0*q + 1.0q0)**2.0q0)

	do k12 = 0, int(L12), 1  ! TODO: Determine if L12, L23 and L31 can be real.
		do k31 = 0, int(L31), 1
			do k23 = 0, int(L23), 1

!		do k31 = 0, int(L31), 1
!			do k23 = 0, int(L23), 1

!					Summation = Summation + 1.0q0 / dble(2*q + 1) ** 2.0q0 * C(j12, q, k12) * C(j23, q, k23) * C(j31, q, k31) * &
!					( W(j1 + 2*q + 2*k12 + 2*k31 + 2, j2 + j12 - 2*k12 + 2*k23 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, &
!						alpha, beta, gamma, pmax) &
!					+ W(j1 + 2*q + 2*k12 + 2*k31 + 2, j3 + j31 - 2*k31 + 2*k23 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, &
!						alpha, gamma, beta, pmax) &
!					+ W(j2 + 2*q + 2*k12 + 2*k23 + 2, j1 + j12 - 2*k12 + 2*k31 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, &
!						beta, alpha, gamma, pmax) &
!					+ W(j2 + 2*q + 2*k12 + 2*k23 + 2, j3 + j23 - 2*k23 + 2*k31 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, &
!						beta, gamma, alpha, pmax) &
!					+ W(j3 + 2*q + 2*k23 + 2*k31 + 2, j1 + j31 - 2*k31 + 2*k12 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, &
!						gamma, alpha, beta, pmax) &
!					+ W(j3 + 2*q + 2*k23 + 2*k31 + 2, j2 + j23 - 2*k23 + 2*k12 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, &
!						gamma, beta, alpha, pmax))

				! The original version is preserved above.  Now this uses precalculated values for the W function via WMatrix.
				!  The ordering is the same as above, so that the 3 corresponds to beta, alpha, gamma and so on.
				Res1 = WMatrix(j1 + 2*q + 2*k12 + 2*k31 + 2, j2 + j12 - 2*k12 + 2*k23 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, 1)
				Res2 = WMatrix(j1 + 2*q + 2*k12 + 2*k31 + 2, j3 + j31 - 2*k31 + 2*k23 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, 2)
				Res3 = WMatrix(j2 + 2*q + 2*k12 + 2*k23 + 2, j1 + j12 - 2*k12 + 2*k31 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, 3)
				Res4 = WMatrix(j2 + 2*q + 2*k12 + 2*k23 + 2, j3 + j23 - 2*k23 + 2*k31 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, 4)
				Res5 = WMatrix(j3 + 2*q + 2*k23 + 2*k31 + 2, j1 + j31 - 2*k31 + 2*k12 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, 5)
				Res6 = WMatrix(j3 + 2*q + 2*k23 + 2*k31 + 2, j2 + j23 - 2*k23 + 2*k12 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, 6)

				C12 = C(j12, q, k12)
				C23 = C(j23, q, k23)
				C31 = C(j31, q, k31)
				Rcoeff = C12 * C31 * C23 !* Rqt
				Summation = Summation + Rcoeff * (Res1+Res2+Res3+Res4+Res5+Res6)

				! For some reason, this does not give the same results as above.  If I cast each of the WMatrix elements to real*16,
				!  I do get the correct answer.  I think I may be running out of registers when it is accessing the elements of the
				!  WMatrix, so it is not keeping the intermediate values in a larger-sized (80/96-bit) register.
!				Summation = Summation + &! C(j12, q, k12) * C(j23, q, k23) * C(j31, q, k31) * &
!				( WMatrix(j1 + 2*q + 2*k12 + 2*k31 + 2, j2 + j12 - 2*k12 + 2*k23 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, 1) &
!				+ WMatrix(j1 + 2*q + 2*k12 + 2*k31 + 2, j3 + j31 - 2*k31 + 2*k23 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, 2) &
!				+ WMatrix(j2 + 2*q + 2*k12 + 2*k23 + 2, j1 + j12 - 2*k12 + 2*k31 + 2, j3 + j23 - 2*q - 2*k23 + j31 - 2*k31 + 2, 3) &
!				+ WMatrix(j2 + 2*q + 2*k12 + 2*k23 + 2, j3 + j23 - 2*k23 + 2*k31 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, 4) &
!				+ WMatrix(j3 + 2*q + 2*k23 + 2*k31 + 2, j1 + j31 - 2*k31 + 2*k12 + 2, j2 + j12 - 2*q - 2*k12 + j23 - 2*k23 + 2, 5) &
!				+ WMatrix(j3 + 2*q + 2*k23 + 2*k31 + 2, j2 + j23 - 2*k23 + 2*k12 + 2, j1 + j12 - 2*q - 2*k12 + j31 - 2*k31 + 2, 6))
			enddo
		enddo
	enddo

	T = ((4.0q0 * PI) ** 3) / real(2*q + 1, 16) ** 2 * Summation
	!T = Summation
	!T = Summation / dble(2*q + 1)
	!T = Summation
	return
end


!! TODO: Description here
!real*16 function AsymptoticPart(TValues, LambdaLower, N, LambdaUpper, ZetaUpperBound)
!	integer N
!	real*16 TValues(N)
!	real*16 LambdaLower
!	integer LambdaUpper, ZetaUpperBound
!	
!	external DGESV, ModifiedRiemannZeta
!	real*16 ModifiedRiemannZeta
!
!	real*16 AsymptoticSum
!	integer q, i, j, Info
!	real*16 A(LambdaUpper+1, LambdaUpper+1)
!	real*16 BX(LambdaUpper+1);
!	real*16 IPIV(N)
!	
!
!	! Initialize the summation
!	AsymptoticSum = 0.0
!
!	!// Fill in x vector  with the values of T(q) for q = N-LambdaUpper to N.
!	do i = 1, LambdaUpper+1, 1
!		BX(i) = TValues(i + (N - (LambdaUpper+1)))
!	enddo
!
!	!// This is our lower value for q, q = N-LambdaUpper.
!	q = N - (LambdaUpper+1)
!
!	! Fill in our matrix, A, with the coefficients.  These coefficients are the q^(i+LambdaLower) in equation (19).
!	!  Also notice that i and j are swapped from how they normally are in nested for loops.
!	do j = 1, LambdaUpper+1, 1
!		do i = 1, LambdaUpper+1, 1
!			A(j, i) = 1.0q0 / dble(q) ** ((i-1)+dble(LambdaLower))
!		enddo
!
!		q = q + 1  ! We want all integral values for q up to N
!	enddo
!
!	! Uses LAPACK to solve our system of LambdaUpper+1 linear equations for x.
!	!LaLinearSolve(A, x, b);
!	CALL dgesv(LambdaUpper+1, 1, A, LambdaUpper+1, IPIV, BX, LambdaUpper+1, Info)
!
!	do i = 1, LambdaUpper+1, 1
!		AsymptoticSum = AsymptoticSum + BX(i) * ModifiedRiemannZeta(N-1, (i-1)+LambdaLower, ZetaUpperBound)
!	enddo
!
!	!write (*,*) 'Asymptotic part:', AsymptoticSum
!	
!	AsymptoticPart = AsymptoticSum
!	return
!end
!	
!
!! This is the Riemann zeta function shown in (18) with its first N terms removed.
!real*16 function ModifiedRiemannZeta(N, i, ZetaUpperBound)
!	integer N
!	real*16 i
!	integer ZetaUpperBound
!	
!	real*16 Sum, Term, Zeta2
!	integer j
!	
!	Sum = 0.0q0
!
!	! For the special case of i == 2, we have a p-series with a definite sum (Zeta2 below).
!	!  This series converges *VERY* slowly, so we take it as a special case.
!	if (int(i) == 2) then
!		Zeta2 = 1.6449340668482264365q0
!		do j = 1, N, 1
!			Zeta2 = Zeta2 - 1.0q0 / dble(j**i)
!		enddo
!		ModifiedRiemannZeta = Zeta2
!		return
!	endif
!
!	j = N+1
!	!do j = N+1, ZetaUpperBound, 1
!	do
!		Term = 1.0q0 / (dble(j) ** dble(i))
!		Sum = Sum + Term
!		! TODO: How small do we really need to get this?
!		if (Term < 1d-21) exit
!		j = j + 1
!	enddo
!
!	ModifiedRiemannZeta = Sum
!	return
!end


! 2F1 hypergeometric function (or Gauss's hypergeometric function) definition from
!  http://mathworld.wolfram.com/HypergeometricFunction.html
! @TODO: Do we need to worry about convergence for our purposes?
!real*16 function Hypergeometric2F1(a, b, c, z)
!	real*16 a, b, c, z
!
!	real*16 Sum
!	real*16 Term
!	real*16 s, t
!	integer nUpper
!	!integer n
!	!external Factorial, PochhammerSymbol
!	real*16 Hypergeometric
!	real*16 Factorial
!	real*16 PochhammerSymbol
!	
!	common /globals/ PI, HyperB, HyperC, HyperVal
!	real*16 PI, HyperB, HyperC, HyperVal
!
!	! @TODO: Error checking all in this section
!	! @TODO: Check how stable this is.
!	if (b == HyperB - 1.0q0) then
!		if (c == HyperC - 1.0q0) then
!			!write (*,*) 'Using recursion relation'
!			s = c
!			t = b - c
!			Hypergeometric2F1 = 1.0q0 + (s + t)/s * z * HyperVal
!			HyperVal = Hypergeometric2F1
!			HyperB = b
!			HyperC = c
!			return
!		endif
!	endif
!
!	Sum = 0.0_16
!	nUpper = 200 !100
!
!!	do n = 0, nUpper, 1
!!		Term = PochhammerSymbol(a, n) * PochhammerSymbol(b, n) / PochhammerSymbol(c, n) &
!!			* z ** n / Factorial(n)
!!		if (Term <= 1d-20) then
!!			Term = 0.0_16
!!			exit
!!		endif
!!		Sum = Sum + Term
!!	enddo
!	
!	!call HYGFX(a, b, c, z, Hypergeometric2F1)
!	
!!	Hypergeometric2F1 = Sum
!	Hypergeometric2F1 = Hypergeometric(a,b,c,z)
!	HyperVal = Hypergeometric2F1
!	HyperB = b
!	HyperC = c
!	return
!end


real*16 function Hypergeometric2F1Sum(a, b, c, z)
	implicit none
	real*16 a, b, c, z
	!integer a, b, c
	!real*16 z
	real*16 Sum, Term, Prev, nd
	integer n, nUpper

	Sum = 1.0_16
	Prev = 0.0_16
	nUpper = 10000  ! Only for extremely singular 2F1 will it even reach close to this number.
	!nUpper = 65

	Term = 1.0_16
	do n = 1, nUpper, 1
		nd = n
		Term = Term * (a+nd-1.0_16) * (b+nd-1.0_16) / ((c+nd-1.0_16) * nd) * z
		Sum = Sum + Term
		if (abs(Sum - Prev) <= 1d-30) then
		!	write (*,*) 'Number of terms: ', n
			exit
		endif
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
		!write (*,*) "This must have c <= b for the recursion relation to work properly.", b, c
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
    real*16 temp16
    real*8 temp8
	
	if (c == 1) then
		hyper_result = (1-z)**(-int(b))
		return
	endif
	
	s = c-1
	t = b-1-s
	hyper_result = (Hypergeometric2F1First(a,b-1,c-1,z) - 1) * s/(z*(s+t))
end function Hypergeometric2F1First



! Pochhammer symbol definition from http://mathworld.wolfram.com/PochhammerSymbol.html
!  Returns x(x+1)...(x+n-1)
real*16 function PochhammerSymbol(x, n)
	real*16 x
	integer n
	real*16 p

	p = 1.0_16;

	do i = 0, n - 1, 1
		p = p * (x + i)  ! Need to convert i to a double?
	enddo

	PochhammerSymbol = p  ! @TODO: Not any real need to use the intermediate p
	return
end

