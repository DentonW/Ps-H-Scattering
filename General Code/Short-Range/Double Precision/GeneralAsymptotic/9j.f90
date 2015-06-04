! This is a standard factorial function using recursion. Negative input values return 0.
real*8 function Fact(n)
	implicit none
	real*8 n
	integer i
	
	!write (*,*) n
	!if (n > 20) then
	!	write (*,*) n
	!end if
	
	! Right now, this is easier than erroring out.
	if (n < 0) then
		Fact = 0.0_8  ! Just return a default value of 0
		!write (*,*) 'n < 0 in Fact function'
		return
	endif

	Fact = 1.0_8
	do i = 1, int(n)
		Fact = Fact * real(i,8)
	enddo

	return
end


real*8 function Del(a, b, c)
	implicit none
	real*8 a, b, c, Fact

	Del = dsqrt(Fact(a - b + c) * Fact(a + b - c) * Fact(a + b + c + 1.0d0) / Fact(b + c - a))
	return
end


real*8 function NineJPart(j11, j21, j31, j12, j22, j32, j13, j23, j33, x, y, z)
	implicit none
	real*8 j11, j21, j31, j12, j22, j32, j13, j23, j33, Fact
	integer x, y, z

	
	if (-j11 + j12 + j33 - j23 + x + z < 0 .or. j21 - j12 + j32 - j23 + x + y < 0 &
		.or. j11 + j21 - j32 + j33 - y - z < 0) then
		!write (*,*) "Negatives in NineJPart", -j11 + j12 + j33 - j23 + x + z, j21 - j12 + j32 - j23 + x + y, j11 + j21 - j32 + j33 - y - z
		!write (*,*) "Cont.:", -j11 + j12 + j33 - j23
		NineJPart = 0
		return
	end if

	NineJPart = (-1)**(x + y + z) * Fact(2*j23 - x)*Fact(j21 + j22 - j23 + x) &
		/Fact(dble(x))/Fact(j22 - j21 + j23 - x)/Fact(j13 + j23 - j33 - x)/Fact(j21 - j12 + j32 - j23 + x + y) &
		*Fact(j13 - j23 + j33 + x)*Fact(j22 - j12 + j32 + y)/Fact(-j11 + j12 + j33 - j23 + x + z) &
		/Fact(dble(y))/Fact(j12 + j22 - j32 - y)*Fact(j31 + j32 - j33 + y)*Fact(j11 + j21 - j32 + j33 - y - z) &
		/Fact(j31 - j32 + j33 - y)/Fact(2*j32 + 1 + y)/Fact(dble(z))/Fact(j11 + j21 - j31 - z) &
		*Fact(2*j11 - z)*Fact(j12 - j11 + j13 + z)/Fact(j11 - j12 + j13 - z)/Fact(j11 + j21 + j31 + 1 - z)
	
	return
end


real*8 function NineJSymbol(j11, j21, j31, j12, j22, j32, j13, j23, j33)
	implicit none
	real*8 j11, j21, j31, j12, j22, j32, j13, j23, j33
	real*8 Del, NineJPart
	real*8 s0, Sum
	integer x, y, z, x1, y1, z1, CheckTriad

	!@TODO: Add in checks for integer or half-integers.
	!@TODO: Add in checks for triangle inequalities.
	if (j11 < 0 .or. j21 < 0 .or. j31 < 0 .or. j12 < 0 .or. j22 < 0 .or. j32 < 0 .or. j13 < 0 .or. j23 < 0 .or. j33 < 0) then
		NineJSymbol = 0
		!write (*,*) "Negative input value"
		return
	end if

	if (CheckTriad(j11, j21, j31) == 0 .or. CheckTriad(j12, j22, j32) == 0 .or. CheckTriad(j13, j23, j33) == 0 .or. &
		CheckTriad(j11, j12, j13) == 0 .or. CheckTriad(j21, j22, j23) == 0 .or. CheckTriad(j31, j32, j33) == 0) then
		!write (*,*) "Input values do not satisfy the triangle inequalities."
		NineJSymbol = 0
		return
	end if

	s0 = (-1)**(j13 + j23 - j33) * (Del(j21, j11, j31) * Del(j12, j22, j32) * Del(j33, j31, j32)) &
				/ (Del(j21, j22, j23) * Del(j12, j11, j13) * Del(j33, j13, j23))
	x1 = min(j22 - j21 + j23, j13 + j23 - j33)
	y1 = min(j31 - j32 + j33, j12 + j22 - j32)
	z1 = min(j11 - j12 + j13, j11 + j21 - j31)
	
	Sum = 0
	do x = 0, x1, 1
		do y = 0, y1, 1
			do z = 0, z1, 1
				Sum = Sum + NineJPart(j11, j21, j31, j12, j22, j32, j13, j23, j33, x, y, z)
			end do
		end do
	end do
	
	NineJSymbol = s0 * Sum

	return
end


!!real*8 function Binomial(n, k)
!!	implicit none
!	real*8 n, k, Fact
!
!	if (n < k .or. n < 0 .or. k < 0) then
!		Binomial = 0.0d0
!		return
!	end if
!
!	Binomial = Fact(n) / (Fact(k) * Fact(n-k))
!	return
!end

real*8 function Binomial(n, k)
	implicit none
	real*8 n, k, Fact, BinPrev
	integer knew

	if (n < k .or. n < 0 .or. k < 0) then
		Binomial = 0.0d0
		return
	end if

	BinPrev = 1.0d0
	do knew = 1, k, 1
		BinPrev = BinPrev * ((n+1.0d0) - knew) / knew
	end do
	Binomial = BinPrev

	return
end



real*8 function Del2(a, b, c)
	implicit none
	real*8 a, b, c, Fact

	Del2 = dsqrt(Fact(a + b - c) * Fact(a - b + c) * Fact(-a + b + c) / Fact(a + b + c + 1.0d0))

	return
end


real*8 function SixjBracket(m1, m2, m3, m4, m5, m6)
	implicit none
	real*8 m1, m2, m3, m4, m5, m6
	real*8 Sum, Binomial, RemoveMe
	integer p, q, n

	p = max(m1+m5+m6, m2+m4+m6, m3+m4+m5, m1+m2+m3)
	q = min(m1+m2+m4+m5, m1+m3+m4+m6, m2+m3+m5+m6)

	Sum = 0.0d0
	do n = p, q, 1
		RemoveMe = Binomial(-m1+m5+m6, dble(n)-m1-m2-m3)
		Sum = Sum + (-1.0d0)**n * Binomial(dble(n+1), dble(n)-m1-m5-m6) * Binomial(m1+m5-m6, dble(n)-m2-m4-m6) * Binomial(m1-m5+m6, dble(n)-m3-m4-m5) * &
								Binomial(-m1+m5+m6, dble(n)-m1-m2-m3)
	end do

	SixjBracket = Sum

	return
end


real*8 function NineJSymbol2(a, b, c, d, e, f, g, h, j)
	implicit none
	real*8 a, b, c, d, e, f, g, h, j
	real*8 s0, Sum, SixjBracket, Del2
	integer k, I1, I2, CheckTriad
	real*8 RemoveMe

	! Check whether triangle inequalities are satisfied.
	if (CheckTriad(a, b, c) == 0 .or. CheckTriad(d, e, f) == 0 .or. CheckTriad(g, h, j) == 0 .or. &
		CheckTriad(a, d, g) == 0 .or. CheckTriad(b, e, h) == 0 .or. CheckTriad(c, f, j) == 0) then
		NineJSymbol2 = 0.0d0
		return
	end if

	I1 = max(dabs(h-d), dabs(b-f), dabs(a-j))
	I2 = min(h+d, b+f, a+j)
	s0 = Del2(a, b, c) * Del2(d, e, f) * Del2(g, h, j) * Del2(a, d, g) * Del2(b, e, h) * Del2(c, f, j)
	Sum = 0.0d0
	
	do k = I1, I2, 1
		RemoveMe = SixjBracket(a, b, c, f, j, dble(k))
		Sum = Sum + (-1.0d0)**(2*k) * (2*k + 1) * SixjBracket(a, b, c, f, j, dble(k)) * SixjBracket(f, d, e, h, b, dble(k)) * SixjBracket(h, j, g, a, d, dble(k))
	end do

	NineJSymbol2 = s0 * Sum

	return
end


