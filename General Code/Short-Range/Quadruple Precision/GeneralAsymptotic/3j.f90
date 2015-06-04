
! Checks the triangular inequalities
integer function Triangular(a, b, c)
	implicit none
	integer a, b, c

	if (a + b < c) then
		Triangular = 0
	else if (abs(a - b) > c) then
		Triangular = 0
	else
		Triangular = 1
	end if

	return
end

! Calculates the 3j symbol using the Racah formula in Messiah, page 1058
! @TODO: Right now, this only takes integer arguments.  This just made calls from
!  IAngular much easier, since they are all integers.  Otherwise, I would have to
!  cast all functions to real*16.
!
! Calculates the Wigner 3J-symbol. The Racah formula is given in Albert Messiah's book on page 1058.
!  Had some help from David Terr's Matlab file (http://www.mathworks.com/matlabcentral/fileexchange/5420-physics-zip/content/Wigner3j.m)
real*16 function ThreeJSymbol(j1, j2, j3, m1, m2, m3)
	implicit none
	integer j1, j2, j3, m1, m2, m3
	integer Terms3j, t, tMin, tMax, Temp, NumTerms
	real*16 Factorial
	real*16 TriSym, r_3j, f

	! Just calculating the 3j right now.
	if (m1 + m2 + m3 /= 0) then
		!write(iwrite,*) 'm1 + m2 + m3 must be equal to 0'
		ThreeJSymbol = 0.0q0
		return
	end if
	if (abs(j1 - j2) > j3 .or. j1 + j2 < j3) then
		!write(iwrite,*) 'j1, j2 and j3 must follow the triangular inequality'
		ThreeJSymbol = 0.0q0
		return
	end if
	if (abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m3) > j3) then
		!write(iwrite,*) 'Not physical with |m_i| > j_i
		ThreeJSymbol = 0.0q0
		return
	end if
	if (j1-m1 /= int(j1-m1) .or. j2-m2 /= int(j2-m2) .or. j3-m3 /= int(j3-m3)) then
		!write (iwrite,*) 'Each j and m pair must have the same parity.'
		ThreeJSymbol = 0.0q0
		return
	end if	

	! Have to initialize, since we are summing into this.
	r_3j = 0.0q0

	! We do the sum first, so that we can just multiply by the first part afterward.
	tMin = max(0, j2-j3-m1, j1+m2-j3)
	tMax = min(j1-m1, j2+m2, j1+j2-j3)
	do t = tMin, tMax, 1
		f = Factorial(int(t)) * Factorial(int(j3-j2+t+m1)) * Factorial(int(j3-j1+t-m2)) &
			* Factorial(int(j1+j2-j3-t)) * Factorial(int(j1-t-m1)) * Factorial(int(j2-t+m2))
		if (f /= 0.0q0) then  ! We really do not want division by 0.
			r_3j = r_3j + (-1)**t / f
		endif
	end do

	r_3j = r_3j * (-1)**(j1 - j2 - m3) * sqrt(TriSym(j1, j2, j3)) &
			* sqrt(Factorial(int(j1+m1)) * Factorial(int(j1-m1)) * Factorial(int(j2+m2)) * Factorial(int(j2-m2)) &
			* Factorial(int(j3+m3)) * Factorial(int(j3-m3)))

	ThreeJSymbol = r_3j
return
end


! Triangle symbol
real*16 function TriSym(a, b, c)
	integer a, b, c
	real*16 Factorial
	TriSym = Factorial(a+b-c) * Factorial(b+c-a) * Factorial(c+a-b) / Factorial(a+b+c+1)
	return
end
