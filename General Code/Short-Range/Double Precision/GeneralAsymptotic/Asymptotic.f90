! Asymptotic expansion based on the paper by Drake and Yan:
!  Phys. Rev. A 52, 3681–3685 (1995) 
!  DOI: 10.1103/PhysRevA.52.3681

real*8 function AsymptoticPartGeneral(TValues, LambdaLower, N, LambdaUpper, ZetaUpperBound)
	implicit none
	integer N
	real*8 TValues(N)
	real*8 LambdaLower
	integer qmax, LambdaUpper, ZetaUpperBound
	real*8 ModifiedRiemannZetaGeneral
	real*8 AsymptoticSum
	integer q, i, j, Info
	real*8 A(LambdaUpper+1, LambdaUpper+1)
	real*8 BX(LambdaUpper+1)
	real*8 IPIV(N)


	! Initialize the summation
	AsymptoticSum = 0.0

	!// Fill in x vector with the values of T(q) for q = N-LambdaUpper to N.
	do i = 1, LambdaUpper+1, 1
		BX(i) = TValues(i + (N - (LambdaUpper+1)))
	enddo

	!// This is our lower value for q, q = N-LambdaUpper.
	q = N - (LambdaUpper)

	! Fill in our matrix, A, with the coefficients.  These coefficients are the q^(i+LambdaLower) in equation (19).
	!  Also notice that i and j are swapped from how they normally are in nested for loops.
	do j = 0, LambdaUpper, 1
		do i = 0, LambdaUpper, 1
			A(j+1, i+1) = 1.0_8 / real(q,8) ** (real(i,8)+real(LambdaLower,8))
		enddo

		q = q + 1  ! We want all integral values for q up to N
	enddo

	! Uses LAPACK to solve our system of LambdaUpper+1 linear equations for x.
	!LaLinearSolve(A, x, b);
!	call dgesv(LambdaUpper+1, 1, A, LambdaUpper+1, IPIV, BX, LambdaUpper+1, Info)
	call GaussElim(A, BX, LambdaUpper+1)

	do i = 0, LambdaUpper, 1
		AsymptoticSum = AsymptoticSum + BX(i+1) * ModifiedRiemannZetaGeneral(N, i+LambdaLower, ZetaUpperBound)
	enddo

	write (*,*) 'Asymptotic part:', AsymptoticSum
	
	AsymptoticPartGeneral = AsymptoticSum
	return
end
	

! This is the Riemann zeta function shown in (18) with its first N terms removed.
real*8 function ModifiedRiemannZetaGeneral(N, i, ZetaUpperBound)
	implicit none
	integer N
	real*8 i
	integer ZetaUpperBound
	
	real*8 Sum, Term, Zeta2
	integer j
	
	Sum = 0.0d0

	! For the special case of i == 2, we have a p-series with a definite sum (Zeta2 below).
	!  This series converges *VERY* slowly, so we take it as a special case.
	if (int(i) == 2) then
		Zeta2 = 1.6449340668482264365_8
		do j = 1, N, 1
			Zeta2 = Zeta2 - 1.0_8 / real(j**i,8)
		enddo
		ModifiedRiemannZetaGeneral = Zeta2
		return
	endif

	j = N+1
	!do j = N+1, ZetaUpperBound, 1
	do
		Term = 1.0_8 / (real(j,8) ** real(i,8))
		Sum = Sum + Term
		! TODO: How small do we really need to get this?
		if (Term < 1d-21) exit
		j = j + 1
	enddo

	ModifiedRiemannZetaGeneral = Sum
	return
end


! The matrix A gets overwritten, along with b.
subroutine GaussElim(A, b, n)
	implicit none
	real*8 A(n,n), b(n), x(n)
	!real*8 A(n,n), b(n), x(n)
	integer n
	real*8 eps, mult
	integer i, j, k

	eps = 1d-25
	
	! Gaussian elimination
	do j = 1, n, 1
!		if (abs(A(j,j)) < eps) then
!			write (*,*) 'Zero pivot encountered in Gaussian elimination.'
!			return
!		endif
		do i = j+1, n, 1
			mult = a(i,j)/a(j,j)
			do k = j+1, n, 1
				a(i,k) = a(i,k) - mult*a(j,k)
			enddo
			b(i) = b(i) - mult*b(j)
		enddo
	enddo

	! Backward substitution
	do i = n, 1, -1
		do j = i+1, n
			b(i) = b(i) - a(i,j)*x(j)
		enddo
		x(i) = b(i)/a(i,i)
	enddo
	
	do i = 1, n, 1
		b(i) = x(i)
	enddo

	! x has the results at this point.    
	return
end




