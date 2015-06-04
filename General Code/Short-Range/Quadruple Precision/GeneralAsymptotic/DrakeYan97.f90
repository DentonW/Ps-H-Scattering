! Calculates Hylleraas-style integrals that have spherical harmonics.
!  This code is based on the algorithm in the paper by Drake and Yan:
!   J. Phys. B: At. Mol. Opt. Phys. 30 (1997) 4723-4750

! The constant 1984.4017075391884912 appears many times here.  This comes from the spherical harmonics.
!  In the Drake/Yan paper, only the l = 0 states are calculated, so this is the l=0, m=0 spherical
!  harmonic, which is 1/sqrt(4*pi).  Taking this to the sixth power and inverting gives the constant.


!module Limits
!   implicit none
!   integer, parameter :: linf = 20
!   integer, parameter :: minf = 20
!end module Limits


!! This is a standard factorial function using recursion.
!!  Only integer arguments are accepted.
!! To keep the output of Calc3j correct, this needs to return 0 whenever we have a negative value for n.
!real*16 function Factorial(n)
!	implicit none
!	integer n
!	integer i
!	
!	! Right now, this is easier than erroring out.
!	if (n < 0) then
!		Factorial = 0.0_16  ! Just return a default value of 1
!		!write (*,*) 'n < 0 in factorial function'
!		return
!	endif
!
!	Factorial = 1.0_16
!	do i = 1, n
!		Factorial = Factorial * real(i,16)
!	enddo
!
!	return
!end


real*16 function DoubleFactorial(n)
	implicit none
	integer n
	integer i
	real*16 Fact

	if (n < -1) then
		!write (*,*) "n < -1 in double factorial function"
		DoubleFactorial = 0.0q0
		return  ! TODO: Return some kind of error
	endif
	if (n == 0 .or. n == 1 .or. n == -1) then
		DoubleFactorial = 1.0q0
		return
	endif

	Fact = n;

	do i = n-2, 1, -2
		Fact = Fact * i
	enddo
	
	DoubleFactorial = Fact
	return
end


!! This is your standard "choose" function, or the binomial coefficient.
!real*16 function BinomialGeneral(n, k)
!	implicit none
!	integer n, k
!	real*16 Factorial
!	
!	! There really should not be cases of this. The binomial coefficient is only defined for k <= n.
!	if (k > n) then
!		BinomialGeneral = 0.0q0
!		return
!	end if
!	
!	BinomialGeneral = Factorial(n) / (Factorial(k) * Factorial(n-k))
!	return
!end


! S is given in equation () - just a simple minimum function.
real*16 function SGeneral(q, j)
	implicit none
	integer q, j
	SGeneral = min(real(q-1,16), 0.5_16 * real(j+1,16))
	return
end


! C is given in equation (4).
real*16 function CCoeff(j, q, k)
	implicit none
	integer j, q, k
	real*16 SGeneral, Binomial, DoubleFactorial
	real*16 LeadingFactor
	real*16 Multiplication
	real*16 SVal
	integer t

	! C for the case of j12 == -2 is given by equation (129).
	if (j == -2) then
		CCoeff = real(2*q + 1,16) * DoubleFactorial(2*q + 2*k) * DoubleFactorial(2*k - 1) &
			/ (DoubleFactorial(2*q + 2*k + 1) * DoubleFactorial(2*k))
		return
	endif

	LeadingFactor = real(2*q + 1,16) / real(j + 2,16) * Binomial(j+2, 2*k + 1)
	Multiplication = 1.0q0
	SVal = SGeneral(q, j)

	do t = 0, int(SVal), 1
		Multiplication = Multiplication * real(2*k + 2*t - j,16) / real(2*k + 2*q - 2*t + 1,16)
	enddo

	CCoeff = LeadingFactor * Multiplication
	return
end


! This is the W function from equation (8).
real*16 function WGeneral(l, m, n, alpha, beta, gamma, pmax)
	implicit none
	integer l, m, n
	real*16 alpha, beta, gamma
	integer pmax
	real*16 Factorial
	real*16 WPartialGeneral
	real*16 LeadingFactor
	real*16 Summation
	integer p

	LeadingFactor = Factorial(l) / ((alpha + beta + gamma) ** (l+m+n+3))
	Summation = 0.0

	! We do the summation backwards so that we can use the backwards recurrence relation in (9).
	!do p = 0, pmax, 1
	do p = pmax, 0, -1
		Summation = Summation + WPartialGeneral(l, m, n, alpha, beta, gamma, p)
	enddo

	WGeneral = LeadingFactor * Summation

	return
end


! This is what's inside the summation of the W function from equation (8).
real*16 function WPartialGeneral(l, m, n, alpha, beta, gamma, p)
	implicit none
	integer l, m, n
	real*16 alpha, beta, gamma
	integer p
	real*16 Factorial
	real*16 Hypergeometric2F1

	WPartialGeneral = Factorial(l+m+n+p+2) / (Factorial(l+1+p) * (l+m+2+p)) * (alpha / (alpha + beta + gamma)) ** p  * &
		Hypergeometric2F1(1.0_16, real(l+m+n+p+3, 16), real(l+m+p+3, 16), real((alpha + beta) / (alpha + beta + gamma), 16))
		
	return
end


!! 2F1 hypergeometric function (or Gauss's hypergeometric function) definition from
!!  http://mathworld.wolfram.com/HypergeometricFunction.html
!real*16 function Hypergeometric2F1(a, b, c, z)
!	implicit none
!	real*16 a, b, c, z
!	!integer a, b, c
!	!real*16 z
!	real*16 Sum, Term, Prev, nd
!	integer n, nUpper
!
!	Sum = 1.0_16
!	Prev = 0.0_16
!	nUpper = 10000  ! Only for extremely singular 2F1 will it even reach close to this number.
!	!nUpper = 65
!
!	Term = 1.0_16
!	do n = 1, nUpper, 1
!		nd = n
!		Term = Term * (a+nd-1.0_16) * (b+nd-1.0_16) / ((c+nd-1.0_16) * nd) * z
!		Sum = Sum + Term
!		if (abs(Sum - Prev) <= 1d-30) then
!		!	write (*,*) 'Number of terms: ', n
!			exit
!		endif
!		Prev = Sum
!	enddo
!
!	Hypergeometric2F1 = Sum
!end function Hypergeometric2F1


! TODO: Parameter order
real*16 function IAngular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, SphHarm, sl, sr)
	implicit none
	integer l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, SphHarm, sl, sr, bfunc, b
	integer Triangular
	real*16 ThreeJSymbol, SixJSymbol, CP12Angular, CP21Angular, CP23Angular, CP31Angular
	real*16 Ang, Term
	integer T1, T2, T3, T1min, T1max, T2min, T2max, T3min, T3max
	integer n1, n2, n3, n1min, n1max, n2min, n2max, n3min, n3max
	real*16 c, Removec1, Removec2

	Ang = 0.0_16

	select case (SphHarm)
		case (1,4)  ! This is the case where we are calculating equation 54 for C^P(r1.r2)
			T1min = abs(1 - l1)
			T1max = 1 + l1
			T2min = abs(1 - l2)
			T2max = 1 + l2
			do T1 = T1min, T1max, 1
				do T2 = T2min, T2max, 1
					if (sr == 2) then
						b = bfunc(l2, T2)
					else  ! Assuming sr == 1 here
						b = bfunc(l1, T1)
					end if
					c = 0
					if (b /= 0) then
						c = CP12Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T1, T2)
					end if
					!Removec1 = CP12Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T1, T2)
					!Removec2 = CP21Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T1, T2)
					!if (Removec1 /= Removec2) write (*,*) "Removec1 /= Removec2"
					Ang = Ang + b * c
				end do
			end do
		case (3,6)  ! This is the case where we are calculating equation 59 for C^P(r2.r3)
			T2min = abs(1 - l2)
			T2max = 1 + l2
			T3min = abs(1 - l3)
			T3max = 1 + l3
			do T2 = T2min, T2max, 1
				do T3 = T3min, T3max, 1
					if (sr == 3) then
						b = bfunc(l3, T3)
					else  ! Assuming sr == 2 here
						b = bfunc(l2, T2)
					end if
					c = 0
					if (b /= 0) then
						c = CP23Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T2, T3)
					end if
					Ang = Ang + b * c
				end do
			end do
		case (2,5)  ! This is the case where we are calculating equation 64 for C^P(r3.r1)
			T3min = abs(1 - l3)
			T3max = 1 + l3
			T1min = abs(1 - l1)
			T1max = 1 + l1
			do T3 = T3min, T3max, 1
				do T1 = T1min, T1max, 1
					if (sr == 1) then
						b = bfunc(l1, T1)
					else  ! Assuming sr == 3 here
						b = bfunc(l3, T3)
					end if
					c = 0
					if (b /= 0) then
						c = CP31Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T3, T1)
					end if
					Ang = Ang + b * CP31Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T3, T1)
				end do
			end do
		case default  ! This is when we are calculating only the basic integral.
			! Satisfies the triangular inequalities for the 3j symbol
			n1min = max(abs(q31-q12), abs(l1p-l1))
			n1max = min(q31+q12, l1p+l1)
			n2min = max(abs(q12-q23), abs(l2p-l2))
			n2max = min(q12+q23, l2p+l2)
			n3min = max(abs(q23-q31), abs(l3p-l3))
			n3max = min(q23+q31, l3p+l3)

			do n1 = n1min, n1max, 1
				do n2 = n2min, n2max, 1
					do n3 = n3min, n3max, 1
						if (Triangular(n1, n2, n3) /= 1) cycle
						Term = (2*n1+1)*(2*n2+1)*(2*n3+1) * SixJSymbol(n1,n2,n3,q23,q31,q12) &
							* ThreeJSymbol(n1,n2,n3,m1p-m1,m2p-m2,m3p-m3) * ThreeJSymbol(l1p,l1,n1,-m1p,m1,m1p-m1) &
							* ThreeJSymbol(l2p,l2,n2,-m2p,m2,m2p-m2) * ThreeJSymbol(l3p,l3,n3,-m3p,m3,m3p-m3) &
							* ThreeJSymbol(l1p,l1,n1,0,0,0) * ThreeJSymbol(l2p,l2,n2,0,0,0) &
							* ThreeJSymbol(l3p,l3,n3,0,0,0) * ThreeJSymbol(q31,q12,n1,0,0,0) &
							* ThreeJSymbol(q12,q23,n2,0,0,0) * ThreeJSymbol(q23,q31,n3,0,0,0)
						Ang = Ang + Term
						!if (Term == 0.0_16) exit  ! Every term past here should be 0.  @TODO: Verify this.
					enddo
				enddo
			enddo

			Ang = Ang * (-1)**(m1p+m2p+m3p+q12+q23+q31) * &
				sqrt(real((2*l1p+1)*(2*l2p+1)*(2*l3p+1)*(2*l1+1)*(2*l2+1)*(2*l3+1),16))
	end select

	IAngular = Ang
	return
end


integer function bfunc(l, T)
	implicit none
	integer l, T

	if (T == l-1) then
		bfunc = l + 1
	else if (T == l+1) then
		bfunc = -l
	else
		bfunc = 0
		!write (*,*) "b function = 0!"
	end if

	return
end


! Equation 54 from Drake and Yan's paper
real*16 function CP12Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T1, T2)
	implicit none
	integer Triangular, l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31
	integer T1, T2, T1min, T1max, T2min, T2max, bfunc
	integer n1, n2, n3, n1min, n1max, n2min, n2max, n3min, n3max
	real*16 ThreeJSymbol, SixJSymbol, CP12AngularTilde
	real*16 U, ThreeJs, CP, Mult

	! Satisfies the triangular inequalities for the 3j symbol
	n1min = max(abs(q31-q12), abs(l1p-T1))
	n1max = min(q31+q12, l1p+T1)
	n2min = max(abs(q12-q23), abs(l2p-T2))
	n2max = min(q12+q23, l2p+T2)
	n3min = max(abs(q23-q31), abs(l3p-l3))
	n3max = min(q23+q31, l3p+l3)

	! We do not include the l12p and l12 from their wavefunction, along with the two Omega functions.
	U = sqrt(real((2*l1p+1)*(2*l2p+1)*(2*l3p+1)*(2*l1+1)*(2*l2+1)*(2*l3+1), 16)) * (-1)**(q12+q23+q31)

	CP = 0.0q0
	do n1 = n1min, n1max, 1
		do n2 = n2min, n2max, 1
			do n3 = n3min, n3max, 1
				if (Triangular(n1, n2, n3) /= 1) cycle
				Mult = ((2*n1+1)*(2*n2+1)*(2*n3+1)*(2*T1+1)*(2*T2+1))
				ThreeJs = ThreeJSymbol(1, l1, T1, 0, 0, 0) * ThreeJSymbol(1, l2, T2, 0, 0, 0) &
					* ThreeJSymbol(l1p, T1, n1, 0, 0, 0) * ThreeJSymbol(l2p, T2, n2, 0, 0, 0) &
					* ThreeJSymbol(l3p, l3, n3, 0, 0, 0) * ThreeJSymbol(q31, q12, n1, 0, 0, 0) &
					* ThreeJSymbol(q12, q23, n2, 0, 0, 0) * ThreeJSymbol(q23, q31, n3, 0, 0, 0)
				if (ThreeJS == 0.0q0) then  ! No need to calculate the rest of the term, since it will be 0.
					continue
				end if
				CP = CP + U * Mult * ThreeJs * CP12AngularTilde(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, n1, n2, n3, T1, T2) &
					* SixJSymbol(n1, n2, n3, q23, q31, q12)
			end do
		end do
	end do

	CP12Angular = CP

	return
end


real*16 function CP12AngularTilde(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, n1, n2, n3, T1, T2)
	implicit none
	integer l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T1, T2, n1, n2, n3, mu
	real*16 ThreeJSymbol

	CP12AngularTilde = 0.0q0
	do mu = -1, 1, 1
		CP12AngularTilde = CP12AngularTilde + (-1)**(-m3+mu) * ThreeJSymbol(1, l1, T1, mu, m1, -mu-m1) * ThreeJSymbol(1, l2, T2, -mu, m2, mu-m2) &
							* ThreeJSymbol(l3p, l3, n3, -m3p, m3, m3p-m3) * ThreeJSymbol(l1p, T1, n1, -m1p, mu+m1, m1p-mu-m1) &
							* ThreeJSymbol(l2p, T2, n2, -m2p, -mu+m2, m2p+mu-m2) * ThreeJSymbol(n1, n2, n3, m1p-mu-m1, m2p+mu-m2, m3p-m3)
	end do

	return
end


! Equation 59 from Drake and Yan's paper
real*16 function CP23Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T2, T3)
	implicit none
	integer Triangular, l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31
	integer T2, T3, T1min, T1max, T2min, T2max, bfunc
	integer n1, n2, n3, n1min, n1max, n2min, n2max, n3min, n3max
	real*16 ThreeJSymbol, SixJSymbol, CP23AngularTilde
	real*16 U, ThreeJs, CP, Mult

	! Satisfies the triangular inequalities for the 3j symbol
	n1min = max(abs(q31-q12), abs(l1p-l1))
	n1max = min(q31+q12, l1p+l1)
	n2min = max(abs(q12-q23), abs(l2p-T2))
	n2max = min(q12+q23, l2p+T2)
	n3min = max(abs(q23-q31), abs(l3p-T3))
	n3max = min(q23+q31, l3p+T3)

	! We do not include the l12p and l12 from their wavefunction, along with the two Omega functions.
	U = sqrt(real((2*l1p+1)*(2*l2p+1)*(2*l3p+1)*(2*l1+1)*(2*l2+1)*(2*l3+1), 16)) * (-1)**(q12+q23+q31)

	CP = 0.0q0
	do n1 = n1min, n1max, 1
		do n2 = n2min, n2max, 1
			do n3 = n3min, n3max, 1
				if (Triangular(n1, n2, n3) /= 1) cycle
				Mult = ((2*n1+1)*(2*n2+1)*(2*n3+1)*(2*T2+1)*(2*T3+1))
				ThreeJs = ThreeJSymbol(1, l2, T2, 0, 0, 0) * ThreeJSymbol(1, l3, T3, 0, 0, 0) &
					* ThreeJSymbol(l1p, l1, n1, 0, 0, 0) * ThreeJSymbol(l2p, T2, n2, 0, 0, 0) &
					* ThreeJSymbol(l3p, T3, n3, 0, 0, 0) * ThreeJSymbol(q31, q12, n1, 0, 0, 0) &
					* ThreeJSymbol(q12, q23, n2, 0, 0, 0) * ThreeJSymbol(q23, q31, n3, 0, 0, 0)
				if (ThreeJS == 0.0q0) then  ! No need to calculate the rest of the term, since it will be 0.
					continue
				end if
				CP = CP + U * Mult * ThreeJs * CP23AngularTilde(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, n1, n2, n3, T2, T3) &
					* SixJSymbol(n1, n2, n3, q23, q31, q12)
			end do
		end do
	end do

	CP23Angular = CP

	return
end


real*16 function CP23AngularTilde(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, n1, n2, n3, T2, T3)
	implicit none
	integer l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T2, T3, n1, n2, n3, mu
	real*16 ThreeJSymbol

	CP23AngularTilde = 0.0q0
	do mu = -1, 1, 1
		CP23AngularTilde = CP23AngularTilde + (-1)**(-m1+mu) * ThreeJSymbol(1, l2, T2, mu, m2, -mu-m2) * ThreeJSymbol(1, l3, T3, -mu, m3, mu-m3) &
							* ThreeJSymbol(l1p, l1, n1, -m1p, m1, m1p-m1) * ThreeJSymbol(l2p, T2, n2, -m2p, mu+m2, m2p-mu-m2) &
							* ThreeJSymbol(l3p, T3, n3, -m3p, -mu+m3, m3p+mu-m3) * ThreeJSymbol(n1, n2, n3, m1p-m1, m2p-mu-m2, m3p+mu-m3)
	end do

	return
end


! Equation 64 from Drake and Yan's paper
real*16 function CP31Angular(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T3, T1)
	implicit none
	integer Triangular, l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31
	integer T3, T1, T1min, T1max, T2min, T2max, bfunc
	integer n1, n2, n3, n1min, n1max, n2min, n2max, n3min, n3max
	real*16 ThreeJSymbol, SixJSymbol, CP31AngularTilde
	real*16 U, ThreeJs, CP, Mult

	! Satisfies the triangular inequalities for the 3j symbol
	n1min = max(abs(q31-q12), abs(l1p-T1))
	n1max = min(q31+q12, l1p+T1)
	n2min = max(abs(q12-q23), abs(l2p-l2))
	n2max = min(q12+q23, l2p+l2)
	n3min = max(abs(q23-q31), abs(l3p-T3))
	n3max = min(q23+q31, l3p+T3)

	! We do not include the l12p and l12 from their wavefunction, along with the two Omega functions.
	U = sqrt(real((2*l1p+1)*(2*l2p+1)*(2*l3p+1)*(2*l1+1)*(2*l2+1)*(2*l3+1), 16)) * (-1)**(q12+q23+q31)

	CP = 0.0q0
	do n1 = n1min, n1max, 1
		do n2 = n2min, n2max, 1
			do n3 = n3min, n3max, 1
				if (Triangular(n1, n2, n3) /= 1) cycle
				Mult = ((2*n1+1)*(2*n2+1)*(2*n3+1)*(2*T3+1)*(2*T1+1))
				ThreeJs = ThreeJSymbol(1, l1, T1, 0, 0, 0) * ThreeJSymbol(1, l3, T3, 0, 0, 0) &
					* ThreeJSymbol(l2p, l2, n2, 0, 0, 0) * ThreeJSymbol(l1p, T1, n1, 0, 0, 0) &
					* ThreeJSymbol(l3p, T3, n3, 0, 0, 0) * ThreeJSymbol(q31, q12, n1, 0, 0, 0) &
					* ThreeJSymbol(q12, q23, n2, 0, 0, 0) * ThreeJSymbol(q23, q31, n3, 0, 0, 0)
				if (ThreeJS == 0.0q0) then  ! No need to calculate the rest of the term, since it will be 0.
					continue
				end if
				CP = CP + U * Mult * ThreeJs * CP31AngularTilde(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, n1, n2, n3, T3, T1) &
					* SixJSymbol(n1, n2, n3, q23, q31, q12)
			end do
		end do
	end do

	CP31Angular = CP

	return
end


real*16 function CP31AngularTilde(l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, n1, n2, n3, T3, T1)
	implicit none
	integer l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, q12, q23, q31, T3, T1, n1, n2, n3, mu
	real*16 ThreeJSymbol

	CP31AngularTilde = 0.0q0
	do mu = -1, 1, 1
		CP31AngularTilde = CP31AngularTilde + (-1)**(-m2+mu) * ThreeJSymbol(1, l1, T1, -mu, m1, mu-m1) * ThreeJSymbol(1, l3, T3, mu, m3, -mu-m3) &
							* ThreeJSymbol(l2p, l2, n2, -m2p, m2, m2p-m2) * ThreeJSymbol(l1p, T1, n1, -m1p, -mu+m1, m1p+mu-m1) &
							* ThreeJSymbol(l3p, T3, n3, -m3p, mu+m3, m3p-mu-m3) * ThreeJSymbol(n1, n2, n3, m1p+mu-m1, m2p-m2, m3p-mu-m3)
	end do

	return
end


integer function Lij(j, q, linf)
	implicit none
	integer j, q, linf
	
	if (j == -2) then
		Lij = linf  ! When j == -2, the summation is infinite.
		return
	endif
	if (mod(j, 2) == 0) then  ! j is even.
		Lij = j/2 - q
		return
	endif
	Lij = (j+1)/2  ! Odd
	
	return
end


integer function Mj(j, minf)
	implicit none
	integer j, minf

	if (j == -2) then  ! Infinite summation
		Mj = minf  ! When j == -2, the summation is infinite.
		return
	endif
	if (mod(j, 2) == 0) then  ! j is even.
		Mj = j/2
		return
	endif
	Mj = minf  ! Infinite summation when j is odd.
	
	return
end


! TODO: Ordering of parameters
real*16 function WR(RunCalc, UsePreCalc, q12, q23, q31, k12, k23, k31, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, pmax, WMatrix)
	use WLimits
	implicit none
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	integer q12, q23, q31, k12, k23, k31, j1, j2, j3, j12, j23, j31, pmax
	logical RunCalc, UsePreCalc
	real*16 alpha, beta, gamma
	real*16 WGeneral
	integer s12, s23, s31
	real*16 W1, W2, W3, W4, W5, W6
	
	s12 = q12 + 2*k12
	s23 = q23 + 2*k23
	s31 = q31 + 2*k31

	! This helps determine which elements of the WMatrix need to be calculated.
	if (RunCalc == .false.) then
		WMatrix(j1+2+s12+s31,j2+2+j12-s12+s23,j3+2+j23-s23+j31-s31,1) = 1000  ! The exact value doesn't matter, as
		WMatrix(j1+2+s12+s31,j3+2+s23+j31-s31,j2+2+j12-s12+j23-s23,2) = 1000  !  long as it is different from 0.
		WMatrix(j2+2+s12+s23,j1+2+j12-s12+s31,j3+2+j23-s23+j31-s31,3) = 1000
		WMatrix(j2+2+s12+s23,j3+2+j23-s23+s31,j1+2+j12-s12+j31-s31,4) = 1000
		WMatrix(j3+2+s23+s31,j1+2+s12+j31-s31,j2+2+j12-s12+j23-s23,5) = 1000
		WMatrix(j3+2+s23+s31,j2+2+s12+j23-s23,j1+2+j12-s12+j31-s31,6) = 1000
		WR = 0.0q0
		return
	end if

	! Use precalculated W functions in WMatrix.
	if (UsePreCalc == .true.) then
		W1 = WMatrix(j1+2+s12+s31,j2+2+j12-s12+s23,j3+2+j23-s23+j31-s31,1)
		W2 = WMatrix(j1+2+s12+s31,j3+2+s23+j31-s31,j2+2+j12-s12+j23-s23,2)
		W3 = WMatrix(j2+2+s12+s23,j1+2+j12-s12+s31,j3+2+j23-s23+j31-s31,3)
		W4 = WMatrix(j2+2+s12+s23,j3+2+j23-s23+s31,j1+2+j12-s12+j31-s31,4)
		W5 = WMatrix(j3+2+s23+s31,j1+2+s12+j31-s31,j2+2+j12-s12+j23-s23,5)
		W6 = WMatrix(j3+2+s23+s31,j2+2+s12+j23-s23,j1+2+j12-s12+j31-s31,6)
		WR = W1 + W2 + W3 + W4 + W5 + W6
		
								!W1 = WGeneral(j1+2+s12+s31,j2+2+j12-s12+s23,j3+2+j23-s23+j31-s31,alpha,beta,gamma,pmax)
								!W2 = WGeneral(j1+2+s12+s31,j3+2+s23+j31-s31,j2+2+j12-s12+j23-s23,alpha,gamma,beta,pmax)
								!W3 = WGeneral(j2+2+s12+s23,j1+2+j12-s12+s31,j3+2+j23-s23+j31-s31,beta,alpha,gamma,pmax)
								!W4 = WGeneral(j2+2+s12+s23,j3+2+j23-s23+s31,j1+2+j12-s12+j31-s31,beta,gamma,alpha,pmax)
								!W5 = WGeneral(j3+2+s23+s31,j1+2+s12+j31-s31,j2+2+j12-s12+j23-s23,gamma,alpha,beta,pmax)
								!W6 = WGeneral(j3+2+s23+s31,j2+2+s12+j23-s23,j1+2+j12-s12+j31-s31,gamma,beta,alpha,pmax)
								!!if (WR /= W1 + W2 + W3 + W4 + W5 + W6) then
								!if (abs((WR - (W1 + W2 + W3 + W4 + W5 + W6)) / (WR + (W1 + W2 + W3 + W4 + W5 + W6)) * 2) > 1D-12) then
								!	write (*,*) "WTF?!", abs((WR - (W1 + W2 + W3 + W4 + W5 + W6)) / (WR + (W1 + W2 + W3 + W4 + W5 + W6)) * 2)
								!	write (*,*) W1, W2, W3, W4, W5, W6
								!	write (*,*)
								!	write (*,*) WMatrix(j1+2+s12+s31,j2+2+j12-s12+s23,j3+2+j23-s23+j31-s31,1), WMatrix(j1+2+s12+s31,j3+2+s23+j31-s31,j2+2+j12-s12+j23-s23,2), &
								!		WMatrix(j2+2+s12+s23,j1+2+j12-s12+s31,j3+2+j23-s23+j31-s31,3), WMatrix(j2+2+s12+s23,j3+2+j23-s23+s31,j1+2+j12-s12+j31-s31,4), &
								!		WMatrix(j3+2+s23+s31,j1+2+s12+j31-s31,j2+2+j12-s12+j23-s23,5), WMatrix(j3+2+s23+s31,j2+2+s12+j23-s23,j1+2+j12-s12+j31-s31,6)
								!	stop
								!end if
		return
	end if	
	
	! Calculate W functions directly.
	W1 = WGeneral(j1+2+s12+s31,j2+2+j12-s12+s23,j3+2+j23-s23+j31-s31,alpha,beta,gamma,pmax)
	W2 = WGeneral(j1+2+s12+s31,j3+2+s23+j31-s31,j2+2+j12-s12+j23-s23,alpha,gamma,beta,pmax)
	W3 = WGeneral(j2+2+s12+s23,j1+2+j12-s12+s31,j3+2+j23-s23+j31-s31,beta,alpha,gamma,pmax)
	W4 = WGeneral(j2+2+s12+s23,j3+2+j23-s23+s31,j1+2+j12-s12+j31-s31,beta,gamma,alpha,pmax)
	W5 = WGeneral(j3+2+s23+s31,j1+2+s12+j31-s31,j2+2+j12-s12+j23-s23,gamma,alpha,beta,pmax)
	W6 = WGeneral(j3+2+s23+s31,j2+2+s12+j23-s23,j1+2+j12-s12+j31-s31,gamma,beta,alpha,pmax)
	WR = W1 + W2 + W3 + W4 + W5 + W6

	return
end


! TODO: Parameter order
real*16 function HylleraasIntegralGeneral(RunCalc, UsePreCalc, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, &
								l1p, l2p, l3p, m1p, m2p, m3p, l1, l2, l3, m1, m2, m3, M12max, M23max, M31max, pmax, SphHarm, sl, sr, WMatrix, Method)
								!l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, qmax)
	use WLimits
	implicit none
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	integer j1, j2, j3, j12, j23, j31, l1p, m1p, l2p, m2p, l3p, m3p, l1, m1, l2, m2, l3, m3, M12max, M23max, M31max, pmax, sl, sr
	logical RunCalc, UsePreCalc
	real*16, allocatable, dimension(:) :: Values
	real*16 Radial, Angular, Constants, Asymp, alpha, beta, gamma, AsymptoticPart
	real*16 IAngular, WR, CCoeff, InnerSum, OuterSum, RemoveMe, LambdaLower
	integer Mj, Lij, SphHarm, l, i, Method, LambdaUpper
	integer M12, M23, M31, L12, L23, L31, q12, q23, q31, k12, k23, k31

	OuterSum = 0.0_16
	M12 = Mj(j12, M12max)
	M23 = Mj(j23, M23max)
	M31 = Mj(j31, M31max)
	!L12 = linf  ! @TODO: Only valid for j == -2.
	allocate(Values(0:M12))

	! This is only for debugging purposes.
	!call omp_set_num_threads(1)
    
	LambdaUpper = 20
	!LambdaUpper = 10

	Radial = 0.0q0
	Angular = 0.0q0
	Constants = 0.0q0
    Values = 0.0q0

	!do q12 = 0, M12, 1
	!	InnerSum = 0.0_16
	!	L12 = Lij(j12,q12,linf)
	!	do q23 = 0, M23, 1
	!		L23 = Lij(j23,q23,linf)
	!		do q31 = 0, M31, 1
	!			L31 = Lij(j31,q31,linf)
	!			do k12 = 0, L12, 1
	!				do k23 = 0, L23, 1
	!					do k31 = 0, L31, 1
	!						if (RunCalc == .true.) then
	!							Angular = IAngular(l1p,m1p,l2p,m2p,l3p,m3p,l1,m1,l2,m2,l3,m3,q12,q23,q31,SphHarm,sl,sr)
	!							if (Angular == 0.0_16) cycle
	!							Constants = CCoeff(j12,q12,k12) * CCoeff(j23,q23,k23) * CCoeff(j31,q31,k31)
	!						end if
	!						Radial = WR(RunCalc, UsePreCalc, q12, q23, q31, k12, k23, k31, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, pmax, WMatrix)
 !
	!						InnerSum = InnerSum + Constants * Radial * Angular
	!					enddo
	!				enddo
	!			enddo
	!		enddo
	!	enddo
    
	do q12 = 0, M12, 1
		InnerSum = 0.0_16
		do q23 = 0, M23, 1
			do q31 = 0, M31, 1
				do k12 = 0, Lij(j12,q12,20), 1
					do k23 = 0, Lij(j23,q23,20), 1
						do k31 = 0, Lij(j31,q31,20), 1
							Angular = IAngular(l1p,m1p,l2p,m2p,l3p,m3p,l1,m1,l2,m2,l3,m3,q12,q23,q31,SphHarm,sl,sr)
							if (Angular == 0.0_16) cycle
							if (RunCalc .eqv. .true.) Constants = CCoeff(j12,q12,k12) * CCoeff(j23,q23,k23) * CCoeff(j31,q31,k31)
							Radial = WR(RunCalc, UsePreCalc, q12, q23, q31, k12, k23, k31, j1, j2, j3, j12, j23, j31, alpha, beta, gamma, pmax, WMatrix)

							!$omp critical(RemoveMe2)
							InnerSum = InnerSum + Constants * Radial * Angular
							!$omp end critical(RemoveMe2)
						enddo
					enddo
				enddo
			enddo
		enddo
	
		Values(q12) = InnerSum

		OuterSum = OuterSum + InnerSum
		Radial = 0.0q0  ! TODO: These three lines needed?
		Angular = 0.0q0
		Constants = 0.0q0
	enddo

	!write (*,'(d24.16)') 1984.4017075391884912_16*OuterSum

	!if (RunCalc == .false. .and. OuterSum /= 0) then
	!	write (*,*) "What?!"
	!	stop
	!end if
    
	! 1st case: Direct sum only (do not include asymptotic part)
	! 2nd case: Asymptotic expansion requested but summation is finite
    ! 3rd case: We are not performing the full computation at this step.
	if (Method == 0 .or. M12 /= M12max .or. RunCalc == .FALSE.) then
		HylleraasIntegralGeneral = 1984.401707539188491230484164294489292942_16*OuterSum
        deallocate(Values)
		return
	end if

	! @TODO: Skipping the asymptotic part right now
	!write (*,*) 'Doing asymptotic part'
!	Asymp = 1984.4017075391884912_16*OuterSum
!	!do l = 0, L12, 1
!	do l = 1, L12, 1
!		if (l <= 2) then
!			RemoveMe = AsymptoticPart(Values, 2.0q0, l, 0, 5000)
!		else
!			RemoveMe = AsymptoticPart(Values, 2.0q0, l, min(l-2,15), 5000)  ! TODO: Is 15 good enough? Maybe put as a parameter.
!		endif
!
!		OuterSum = 0.0q0
!		do i = 0, l-1, 1
!			OuterSum = OuterSum + Values(i)
!		enddo
!		Asymp = OuterSum + RemoveMe
!
!!		write (*,*) l, "Sum:", 1984.4017075391884912_16 * OuterSum, " Asymptotic:", 1984.4017075391884912_16 * Asymp
!	enddo

	l = M12+1
!	if (l <= 2) then
!		RemoveMe = AsymptoticPart(Values, 2.0q0, l, 0)
!	else
!		RemoveMe = AsymptoticPart(Values, 2.0q0, l, min(l-2, LambdaUpper))  ! TODO: Is 15 good enough? Maybe put as a parameter.
!	endif

    LambdaLower = 0.5q0 * real(j12 + 1,16) + 0.5q0 * real(j23 + 1,16) + 0.5q0 * real(j31 + 1,16) + 4.0q0
    !LambdaLower = 2.0q0  ! Only if we have r_{ij} == -2

    RemoveMe = AsymptoticPart(Values, LambdaLower, l, LambdaUpper)
	OuterSum = 0.0q0
	!do i = 0, M12-2, 1
	do i = 0, M12, 1
		OuterSum = OuterSum + Values(i)
	enddo
	Asymp = OuterSum + RemoveMe
	Asymp = Asymp * 1984.401707539188491230484164294489292942_16

	deallocate(Values)

	HylleraasIntegralGeneral = Asymp
	return
end

