
! Calculates the 6j symbol using the Racah formula in Messiah, page 1065
! @TODO: Read the note on the ThreeJSymbol function.
!WHY WERE THE ARGUMENTS SWAPPED? IANGULAR HAS THEM SWAPPED. ANYWHERE ELSE?
real*16 function SixJSymbol(rj1, rj2, rj3, J1, J2, J3)
	implicit none
	integer J1, J2, J3, rj1, rj2, rj3
	integer iterm, NumTerms
	integer CheckTriad, Terms6j, t
	real*16 Factorial
	real*16 TriSym, r_6j, temp

	r_6j = 0.0_16
	t = 0; iterm = 0

	! Checks for the values based on information from http://mathworld.wolfram.com/Wigner6j-Symbol.html
	if (CheckTriad(rj1,rj2,rj3) /= 1 .OR. CheckTriad(rj1,J2,J3) /= 1 .OR. &
		CheckTriad(J1,rj2,J3) /= 1 .OR. CheckTriad(J1,J2,rj3) /= 1) then
		!write (iwrite,*) 'Each triad must sum together to an integer and satisfy the triangular inequalities.'
		!write (*,*) 'Each triad must sum together to an integer and satisfy the triangular inequalities.'
		SixJSymbol = 0.0q0
		return
	endif

	! TODO: Change this to use a do loop.  Also implement an upper limit to this summation like in ThreeJSymbol.
	! Again, we do the sum first
    NumTerms = Terms6j(rj1,rj2,rj3,J1,J2,J3)
10	if (iterm < NumTerms) then  ! Have to determine how many terms we have
		! The only valid terms are those where the arguments to the factorials are non-negative.
		temp = (Factorial(int(t-rj1-rj2-rj3)) * Factorial(int(t-rj1-J2-J3)) * Factorial(int(t-J1-rj2-J3)) &
				* Factorial(int(t-J1-J2-rj3))) * Factorial(int(rj1+rj2+J1+J2-t)) * Factorial(int(rj2+rj3+J2+J3-t)) &
				* Factorial(int(rj3+rj1+J3+J1-t))
		! If temp is 0, then 
		if (temp /= 0.0q0) then
			iterm = iterm + 1  ! This is one of our terms that Terms6j returned.
			r_6j = r_6j + (-1)**t * Factorial(t+1) / temp
		endif

		t = t + 1
		goto 10
	endif

	r_6j = r_6j * Sqrt(TriSym(rj1,rj2,rj3) * TriSym(rj1,J2,J3) * TriSym(J1,rj2,J3) * TriSym(J1,J2,rj3))
	SixJSymbol = r_6j
	return
end	


! Each triad has to add to an integer and must satisfy the the triangular inequalities.
!  Returns 1 if the conditions are satisfied and 0 if they are not.
integer function CheckTriad(a, b, c)
	integer a, b, c
	
	!if (a+b+c /= dble(int(a+b+c))) then
	!	CheckTriad = 0
	!	return
	!endif
	
	! Triangular inequalities
	if (c > a + b .OR. c < abs(a - b)) then
		CheckTriad = 0
		return
	endif
	
	CheckTriad = 1
	return
end


! Calculates the number of terms in the 6j sum.
!  The list comes from page 1065 of Messiah, volume 2.
integer function Terms6j(rj1, rj2, rj3, J1, J2, J3)
	integer rj1, rj2, rj3, J1, J2, J3
	integer s
	
	! Initialize the smallest value
	s = rj1 + rj2 - rj3

! TODO: Rewrite in terms of min function.
	if (rj1+J2-J3 < s) then
		s = rj1+J2-J3
	endif
	if (J1+rj2-J3 < s) then
		s = J1+rj2-J3
	endif
	if (J1+J2-rj3 < s) then
		s = J1+J2-rj3
	endif
	if (rj2+rj3-rj1 < s) then
		s = rj2+rj3-rj1
	endif
	if (J2+J3-rj1 < s) then
		s = J2+J3-rj1
	endif
	if (rj2+J3-J1 < s) then
		s = rj2+J3-J1
	endif
	if (J2+rj3-J1 < s) then
		s = J2+rj3-J1
	endif
	if (rj3+rj1-rj2 < s) then
		s = rj3+rj1-rj2
	endif
	if (J3+rj1-J2 < s) then
		s = J3+rj1-J2
	endif
	if (J3+J1-rj2 < s) then
		s = J3+J1-rj2
	endif
	if (rj3+J1-J2 < s) then
		s = rj3+J1-J2
	endif
	
	! Is this completely correct?
	Terms6j = s + 1
return
end
