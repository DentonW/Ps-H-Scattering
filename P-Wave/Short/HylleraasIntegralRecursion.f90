!       MODULE PRMTS1
!       INTEGER MAXF,MAXG,MAXFM,MAXFN,SHELL,SWITCH
!       PARAMETER(SHELL=24) 
!!      SHELL=SHELL+1 if grad. min. or hfs
!!      SHELL=SHELL+2*L, L - angular momentum
!!      max(SHELL)=2*6+6 = 18 in quad
!!      max(SHELL)=2*9+6 = 24 in sixr
!!      max(SHELL)=2*12+6= 30 in octr
!
!       PARAMETER(MAXG=(7+SHELL)*(8+SHELL)*(9+SHELL)/6-1)
!       PARAMETER(MAXFN=((1+SHELL)*(2+SHELL)*(3+SHELL)*(4+SHELL)*(5+SHELL)*(6+SHELL))/720)
!!      PARAMETER(MAXFM=((1+SHELL)*(2+SHELL)*(3+SHELL)*(4+SHELL)*(5+SHELL))/120)
!       PARAMETER(MAXF = MAXFN) ! if SWITCH=0
!!      PARAMETER(MAXF = MAXFN+6*MAXFM)  ! if SWITCH=1
!!      PARAMETER(MAXF = MAXFN+12*MAXFM) ! if SWITCH=2
!       END MODULE PRMTS1




module gam
	use PRMTS
	implicit none
	real*16 tg1(0:MAXG), tg2(0:MAXG), tg3(0:MAXG), tg4(0:MAXG), tg5(0:MAXG), tg6(0:MAXG)
	real*16 f0val, f1val, f2val, f3val, f12val, f13val, f23val, f123val
	real*16 fpre(0:MAXF)
end module gam


subroutine PreCalcGamma(alpha, beta, gamma, sm)
	use gam
	implicit none
	real*16 alpha, beta, gamma
	integer sm
	real*16 f0, f1, f2, f3, f12, f13, f23, f123
	integer n1,n2,n3,n4,n5,n6,i,j,k,l,m,n
	real*16 nd1,nd2,nd3,nd4,nd5,nd6
	real*16 expr,expr2
	real*16 w1, w2, w3
	integer idx,idx6,index

	call p_g(tg1,beta+gamma,alpha,0.0q0,sm) 
	call p_g(tg2,alpha+gamma,beta,0.0q0,sm)
	call p_g(tg3,alpha+beta,gamma,0.0q0,sm)
	call p_g(tg4,beta,gamma,0.0q0,sm)
	call p_g(tg5,gamma,alpha,0.0q0,sm)
	call p_g(tg6,alpha,beta,0.0q0,sm)

	! Temporarily trying precomputation of f0, f1, etc.
	f0val = f0(alpha, beta, gamma)
	f1val = f1(alpha, beta, gamma)
	f2val = f2(alpha, beta, gamma)
	f3val = f3(alpha, beta, gamma)
	f12val = f12(alpha, beta, gamma)
	f13val = f13(alpha, beta, gamma)
	f23val = f23(alpha, beta, gamma)
	f123val = f123(alpha, beta, gamma)

	! From Pachucki's code
	fpre(0) = 0.0q0
	w1 = alpha
	w2 = beta
	w3 = gamma

	if (sm>=0) then
		fpre(idx6(0,0,0,0,0,0)) = f0val
	end if

	if (sm>=1) then 
		fpre(idx6(1,0,0,0,0,0)) = f1val
		fpre(idx6(0,1,0,0,0,0)) = f2val
		fpre(idx6(0,0,1,0,0,0)) = f3val
	endif

	if (sm>=2) then
		fpre(idx6(1,1,0,0,0,0)) = f12val
		fpre(idx6(1,0,1,0,0,0)) = f13val 
		fpre(idx6(0,1,1,0,0,0)) = f23val
	endif

	if (sm>=3) then 
		fpre(idx6(1,1,1,0,0,0)) = f123val
	endif

	  do k=2,sm 
	   do j=0,k
		do i=0,j

		 n1 = i
		 n2 = j-i
		 n3 = k-j
		 nd1 = n1
		 nd2 = n2
		 nd3 = n3

		 IF(max(n1,n2,n3)<2) cycle

		 IF(n1.EQ.max(n1,n2,n3)) then

! expr=f(n1+2,n2,n3)

		  n1 = n1-2
		  nd1 = n1

		  expr2=nd2*(nd2-1)/(nd1+1)*fpre(idx6(n1+2,n2-2,n3,0,0,0))+ &
			  (nd1+2*nd2+nd3+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd1+1)*tg6(idx(n2-1,n1+1,n3))+ &
			  nd2*(nd2-1)/(nd3+1)*fpre(idx6(n1,n2-2,n3+2,0,0,0))+ &
			  1/(nd3+1)*tg4(idx(n3+1,n2-1,n1)) 

		  if (n2==0) then
			  expr2=expr2-1/(nd1+1)*tg2(idx(-1,0,n1+n3+1))- &
				  1/(nd3+1)*tg2(idx(-1,0,n1+n3+1))
		  endif

		  expr=1/(w3**2)*expr2

		  expr2=nd3*(nd3-1)/(nd1+1)*fpre(idx6(n1+2,n2,n3-2,0,0,0))+ &
			  (2*nd3+nd2+nd1+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd1+1)*tg5(idx(n1+1,n3-1,n2))+ &
			  nd3*(nd3-1)/(nd2+1)*fpre(idx6(n1,n2+2,n3-2,0,0,0))+ &
			  1/(nd2+1)*tg4(idx(n3-1,n2+1,n1)) 

		  if (n3==0) then
		   expr2=expr2-1/(nd1+1)*tg3(idx(-1,0,n1+n2+1))- &
			  1/(nd2+1)*tg3(idx(-1,0,n1+n2+1))
		  endif

		  expr=expr+1/(w2**2)*expr2  

		  expr2=nd1*(nd1-1)/(nd2+1)*fpre(idx6(n1-2,n2+2,n3,0,0,0))+ &
			  (nd3+nd2+2*nd1+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd2+1)*tg6(idx(n2+1,n1-1,n3))+ &
			  nd1*(nd1-1)/(n3+1)*fpre(idx6(n1-2,n2,n3+2,0,0,0))+ &
			  1/(nd3+1)*tg5(idx(n1-1,n3+1,n2))

		  if (n1==0) then
		   expr2=expr2-1/(nd2+1)*tg1(idx(-1,0,n2+n3+1))- &
			  1/(nd3+1)*tg1(idx(-1,0,n2+n3+1))
		  endif

		  expr=expr-w1**2/((w3*w2)**2)*expr2
		  expr=expr*(1+nd1)/2

		  fpre(idx6(n1+2,n2,n3,0,0,0))=expr 

		 CYCLE
		 ENDIF

		 IF(n2.EQ.max(n1,n2,n3)) then

! expr=f(n1,n2+2,n3)

		  n2=n2-2
		  nd2=n2
 
		  expr2=nd3*(nd3-1)/(nd2+1)*fpre(idx6(n1,n2+2,n3-2,0,0,0))+ &
			  (nd1+2*nd3+nd2+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd2+1)*tg4(idx(n3-1,n2+1,n1))+ &
			  nd3*(nd3-1)/(nd1+1)*fpre(idx6(n1+2,n2,n3-2,0,0,0))+ &
			  1/(nd1+1)*tg5(idx(n1+1,n3-1,n2)) 

		  if (n3==0) then
			  expr2=expr2-1/(nd2+1)*tg3(idx(-1,0,n1+n2+1))- &
				  1/(nd1+1)*tg3(idx(-1,0,n1+n2+1))  
		  endif

		  expr=1/(w1**2)*expr2

		  expr2=nd1*(nd1-1)/(nd2+1)*fpre(idx6(n1-2,n2+2,n3,0,0,0))+ &
			  (2*nd1+nd2+nd3+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd2+1)*tg6(idx(n2+1,n1-1,n3))+ &
			  nd1*(nd1-1)/(nd3+1)*fpre(idx6(n1-2,n2,n3+2,0,0,0))+ &
			  1/(nd3+1)*tg5(idx(n1-1,n3+1,n2))

		  if (n1==0) then
		   expr2=expr2-1/(nd2+1)*tg1(idx(-1,0,n3+n2+1))- &
			  1/(nd3+1)*tg1(idx(-1,0,n3+n2+1))
		  endif 

		  expr=expr+1/(w3**2)*expr2
 
		  expr2=nd2*(nd2-1)/(nd3+1)*fpre(idx6(n1,n2-2,n3+2,0,0,0))+ &
			  (nd1+nd3+2*nd2+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd3+1)*tg4(idx(n3+1,n2-1,n1))+ &
			  nd2*(nd2-1)/(nd1+1)*fpre(idx6(n1+2,n2-2,n3,0,0,0))+ &
			  1/(nd1+1)*tg6(idx(n2-1,n1+1,n3)) 

		  if (n2==0) then
		  expr2=expr2-1/(nd3+1)*tg2(idx(-1,0,n1+n3+1))- &
			  1/(nd1+1)*tg2(idx(-1,0,n1+n3+1))
		  endif
		 
		  expr=expr-w2**2/((w1*w3)**2)*expr2
		  expr=expr*(1+nd2)/2

		  fpre(idx6(n1,n2+2,n3,0,0,0))=expr 

		 CYCLE
		 ENDIF

		 IF(n3.EQ.max(n1,n2,n3)) then

! expr=f(n1,n2,n3+2)

		  n3=n3-2
		  nd3=n3

		  expr2= nd2*(nd2-1)/(nd3+1)*fpre(idx6(n1,n2-2,n3+2,0,0,0))+ & 
			  (nd1+2*nd2+nd3+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd3+1)*tg4(idx(n3+1,n2-1,n1))+ &
			  nd2*(nd2-1)/(nd1+1)*fpre(idx6(n1+2,n2-2,n3,0,0,0))+ &
			  1/(nd1+1)*tg6(idx(n2-1,n1+1,n3)) 

			if (n2==0) then
			  expr2=expr2-1/(nd3+1)*tg2(idx(-1,0,n1+n3+1))- &
				  1/(nd1+1)*tg2(idx(-1,0,n1+n3+1))
			endif

		  expr=1/(w1**2)*expr2

		  expr2=nd1*(nd1-1)/(nd3+1)*fpre(idx6(n1-2,n2,n3+2,0,0,0))+ &
			  (2*nd1+nd2+nd3+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd3+1)*tg5(idx(n1-1,n3+1,n2))+ &
			  nd1*(nd1-1)/(nd2+1)*fpre(idx6(n1-2,n2+2,n3,0,0,0))+ &
			  1/(nd2+1)*tg6(idx(n2+1,n1-1,n3))

		  if (n1==0) then
		   expr2=expr2-1/(nd3+1)*tg1(idx(-1,0,n3+n2+1))- &
			  1/(nd2+1)*tg1(idx(-1,0,n3+n2+1))
		  endif

		  expr=expr+1/(w2**2)*expr2

		  expr2= nd3*(nd3-1)/(nd2+1)*fpre(idx6(n1,n2+2,n3-2,0,0,0))+ &
			  (nd1+nd2+2*nd3+2)*fpre(idx6(n1,n2,n3,0,0,0))+ &
			  1/(nd2+1)*tg4(idx(n3-1,n2+1,n1))+ &
			  nd3*(nd3-1)/(nd1+1)*fpre(idx6(n1+2,n2,n3-2,0,0,0))+ &
			  1/(nd1+1)*tg5(idx(n1+1,n3-1,n2))

		  if (n3==0) then
		   expr2=expr2-1/(nd2+1)*tg3(idx(-1,0,n1+n2+1))- &
			  1/(nd1+1)*tg3(idx(-1,0,n1+n2+1))
		  endif

		  expr=expr-w3**2/((w1*w2)**2)*expr2
		  expr=expr*(1+nd3)/2
 
		  fpre(idx6(n1,n2,n3+2,0,0,0))=expr 

		 CYCLE
		 ENDIF

		enddo
	   enddo
	  enddo


	  index=1
	  do n=1,sm
	   do m=0,n
		do l=0,m
		 do k=0,l
		  do j=0,k
		   do i=0,j

			 index = index+1

			 n1 = i
			 n2 = j-i
			 n3 = k-j
			 n4 = l-k
			 n5 = m-l
			 n6 = n-m

			 if(n4+n5+n6.EQ.0) cycle

			 nd1 = n1
			 nd2 = n2
			 nd3 = n3
			 nd4 = n4
			 nd5 = n5
			 nd6 = n6

			 if(n4.EQ.max(n4,n5,n6)) then

! f(n1,n2,n3,n4+1,n5,n6)
			 n4 = n4-1
			 nd4 = n4

		  expr= &
			+(-1+nd3)*nd3*nd6*fpre(idx6(n1,n2,-2+n3,n4,1+n5,-1+n6)) &
			+(-1+nd2)*nd2*nd5*fpre(idx6(n1,-2+n2,n3,n4,-1+n5,1+n6)) &
			-(-1+nd1)*nd1*nd6*fpre(idx6(-2+n1,n2,n3,n4,1+n5,-1+n6)) &
			-((-1+nd1)*nd1*nd5*fpre(idx6(-2+n1,n2,n3,n4,-1+n5,1+n6))) &
			+(+1-nd1+nd2+nd3+nd4)*nd5*nd6*fpre(idx6(n1,n2,n3,n4,-1+n5,-1+n6)) &
			-(-1+nd3)*nd3*fpre(idx6(n1,n2,-2+n3,n4,1+n5,n6))*w3 &
			+(-1+nd1)*nd1*fpre(idx6(-2+n1,n2,n3,n4,1+n5,n6))*w3 &
			-(1-nd1+nd2+nd3+nd4)*nd5*fpre(idx6(n1,n2,n3,n4,-1+n5,n6))*w3 &
			-(-1+nd2)*nd2*fpre(idx6(n1,-2+n2,n3,n4,n5,1+n6))*w2 &
			+(-1+nd1)*nd1*fpre(idx6(-2+n1,n2,n3,n4,n5,1+n6))*w2 &
			-(1-nd1+nd2+nd3+nd4)*nd6*fpre(idx6(n1, n2, n3, n4,n5,-1+n6))*w2+ &
			-nd5*nd6*fpre(idx6(n1,n2,n3,1+n4,-1+n5,-1+n6))*w1 &
			+(1-nd1+nd2+nd3+nd4)*fpre(idx6(n1,n2,n3,n4,n5,n6))*w2*w3 &
			+nd6*fpre(idx6(n1,n2,n3,1+n4,n5,-1+ n6))*w1*w2 &
			+nd5*fpre(idx6(n1,n2,n3,1+n4,-1+n5,n6))*w1*w3



			 if (n2==0)then
			  expr=expr-nd5*tg2(idx(n4+n6,n5-1,n1+n3-1))+ &
			 tg2(idx(n4+n6,n5,n1+n3-1))*w2
			 endif 

			 if (n3==0)then
			  expr=expr-nd6*tg3(idx(n4+n5,n6-1,n1+n2-1))+ &
			  tg3(idx(n4+n5,n6,n1+n2-1))*w3
			 endif  

			 if (n1==0)then
			  expr=expr+nd5*tg1(idx(n5+n6-1,n4,n2+n3-1))+ &
			  nd6*tg1(idx(n5+n6-1,n4,n2+n3-1))- &
			  tg1(idx(n5+n6,n4,n2+n3-1))*w2- &
			  tg1(idx(n5+n6,n4,n2+n3-1))*w3
			 endif  

			 expr=expr/(w1*w2*w3)
			 fpre(index)=expr          

			 cycle
			 endif

			 if(n5.EQ.max(n4,n5,n6)) then         

! fpre(n1,n2,n3,n4,n5+1,n6)

			 n5 = n5-1
			 nd5 = n5

			 expr=(-1+nd1)*nd1*nd4*fpre(idx6(-2+n1,n2,n3,-1+n4,n5,1+n6))- &
			  (-1+nd2)*nd2*nd4*fpre(idx6(n1,-2+n2,n3,-1+n4,n5,1+n6))- &
			  (-1+nd2)*nd2*nd6*fpre(idx6(n1,-2+n2,n3,1+n4,n5,-1+n6))+ &
			  (-1+nd3)*nd3*nd6*fpre(idx6(n1,n2,-2+n3,1+n4,n5,-1+n6))+ &
			  nd4*(1+nd1-nd2+nd3+nd5)*nd6* &
			  fpre(idx6(n1,n2,n3,-1+n4,n5,-1+n6))- &
			  (-1+nd1)*nd1*fpre(idx6(-2+n1,n2,n3,n4,n5,1+n6))*w1+ &
			  (-1+nd2)*nd2*fpre(idx6(n1,-2+n2,n3,n4,n5,1+n6))*w1- &
			  (1+nd1-nd2+nd3+nd5)*nd6* &
			  fpre(idx6(n1,n2,n3,n4,n5,-1+n6))*w1- &
			  nd4*nd6*fpre(idx6(n1,n2,n3,-1+n4,1+n5,-1+n6))*w2+ &
			  nd6*fpre(idx6(n1,n2,n3,n4,1+n5,-1+n6))*w1*w2+ &
			  (-1+nd2)*nd2*fpre(idx6(n1,-2+n2,n3,1+n4,n5,n6))*w3- &
			  (-1+nd3)*nd3*fpre(idx6(n1,n2,-2+n3,1+n4,n5,n6))*w3- &
			  nd4*(1+nd1-nd2+nd3+nd5)* &
			  fpre(idx6(n1,n2,n3,-1+n4,n5,n6))*w3+ &
			  (1+nd1-nd2+nd3+nd5)*fpre(idx6(n1,n2,n3,n4,n5,n6))*w1*w3+ &
			  nd4*fpre(idx6(n1,n2,n3,-1+n4,1+n5,n6))*w2*w3

			if (n1==0)then
			 expr=expr-nd4*tg1(idx(n5+n6,n4-1,n3+n2-1))+ &
			  tg1(idx(n5+n6,n4,n3+n2-1))*w1
			endif 

			if (n3==0)then
			 expr=expr-nd6*tg3(idx(n4+n5,n6-1,n1+n2-1))+ &
			  tg3(idx(n4+n5,n6,n1+n2-1))*w3
			endif  

			if (n2==0)then
			 expr=expr+nd4*tg2(idx(n4+n6-1,n5,n1+n3-1))+ &
			   nd6*tg2(idx(n4+n6-1,n5,n1+n3-1))- &
			   tg2(idx(n4+n6,n5,n1+n3-1))*w1- &
			   tg2(idx(n4+n6,n5,n1+n3-1))*w3
			endif


			expr=expr/(w1*w2*w3)
			fpre(index)=expr

			cycle
			endif

			if(n6.EQ.max(n4,n5,n6)) then 

! f(n1,n2,n3,n4,n5,n6+1)

			n6 = n6-1
			nd6 = n6
			
			expr=(-1+nd1)*nd1*nd4* &
			fpre(idx6(-2+n1,n2,n3,-1+n4,1+n5,n6))+ &
			(-1+nd2)*nd2*nd5*fpre(idx6(n1,-2+n2,n3,1+n4,-1+n5,n6))- &
			(-1+nd3)*nd3*nd4*fpre(idx6(n1,n2,-2+n3,-1+n4,1+n5,n6))- &
			(-1+nd3)*nd3*nd5*fpre(idx6(n1,n2,-2+n3,1+n4,-1+n5,n6))+ &
			nd4*nd5*(1+nd1+nd2-nd3+nd6)* &
			fpre(idx6(n1,n2,n3,-1+n4,-1+n5,n6))- &
			(-1+nd1)*nd1*fpre(idx6(-2+n1,n2,n3,n4,1+n5,n6))*w1+ &
			(-1+nd3)*nd3*fpre(idx6(n1,n2,-2+n3,n4,1+n5,n6))*w1- &
			nd5*(1+nd1+nd2-nd3+nd6)* &
			fpre(idx6(n1,n2,n3,n4,-1+n5,n6))*w1- &
			(-1+nd2)*nd2*fpre(idx6(n1,-2+n2,n3,1+n4,n5,n6))*w2+ &
			(-1+nd3)*nd3*fpre(idx6(n1,n2,-2+n3,1+n4,n5,n6))*w2- &
			nd4*(1+nd1+nd2-nd3+nd6)* &
			fpre(idx6(n1,n2,n3,-1+n4,n5,n6))*w2+ &
			(1+nd1+nd2-nd3+nd6)* &
			fpre(idx6(n1,n2,n3,n4,n5,n6))*w1*w2- &
			nd4*nd5*fpre(idx6(n1, n2, n3,-1+n4,-1+n5,1+n6))*w3+ &
			nd5*fpre(idx6(n1,n2,n3,n4,-1+n5,1+n6))*w1*w3+ &
			nd4*fpre(idx6(n1,n2,n3,-1+n4,n5,1+n6))*w2*w3

			if (n1==0)then
			expr=expr-nd4*tg1(idx(n5+n6,n4-1,n3+n2-1))+ &
			tg1(idx(n5+n6,n4,n3+n2-1))*w1
			endif 


			if (n2==0)then
			expr=expr-nd5*tg2(idx(n4+n6,n5-1,n1+n3-1))+ &
			tg2(idx(n4+n6,n5,n1+n3-1))*w2
			endif  

			if (n3==0)then
			expr=expr+nd4*tg3(idx(n4+n5-1,n6,n1+n2-1))+ &
			nd5*tg3(idx(n4+n5-1,n6,n1+n2-1))- &
			tg3(idx(n4+n5,n6,n1+n2-1))*w1- &
			tg3(idx(n4+n5,n6,n1+n2-1))*w2
			endif  

			expr=expr/(w1*w2*w3)
			fpre(index)=expr 

			cycle
			endif         

		   enddo
		  enddo
		 enddo
		enddo
	   enddo
	  enddo


	return
end

real*16 function HylleraasIntegralRecursion(j1, j2, j3, j12, j23, j31, alpha, beta, gamma, qmax, Method)
	use gam
	implicit none
	integer j1, j2, j3, j12, j23, j31
	integer qmax, Method  ! These two are not used but are kept 
	real*16 alpha, beta, gamma
	real*16 f
	real*16 Sum
	integer sm, idx6

	if (j1 < -1 .or. j2 < -1 .or. j3 < -1 .or. j12 < -1 .or. j23 < -1 .or. j31 < -1) then
		write (*,*) "Argument to recursion relations is too singular. Exiting."
		stop
	end if

	!Sum = f(j23+1, j31+1, j12+1, j1+1, j2+1, j3+1, real(alpha,16), real(beta,16), real(gamma,16))
	Sum = fpre(idx6(j23+1, j31+1, j12+1, j1+1, j2+1, j3+1))
	Sum = 1984.4017075391884912304841642945q0 * Sum  ! (4*pi)^3

	HylleraasIntegralRecursion = Sum

	return
end


real*16 function f1star(n2, n3, n4, n5, n6, w1, w2, w3)
	use gam
	implicit none
	integer n1, n2, n3, n4, n5, n6
	real*16 w1, w2, w3
	real*16 Gamma1
	!f1star = Gamma1(n5 + n6 - 1, n4, n3 + n2 - 1, w2 + w3, w1, 0.0q0)
	f1star = Gamma1(n5 + n6 - 1, n4, n3 + n2 - 1, w2 + w3, w1, 0.0q0, tg1)
	return
end

real*16 function f2star(n1, n3, n4, n5, n6, w1, w2, w3)
	use gam
	implicit none
	integer n1, n2, n3, n4, n5, n6
	real*16 w1, w2, w3
	real*16 Gamma1
	f2star = Gamma1(n4 + n6 - 1, n5, n1 + n3 - 1, w1 + w3, w2, 0.0q0, tg2)
	return
end

real*16 function f3star(n1, n2, n4, n5, n6, w1, w2, w3)
	use gam
	implicit none
	integer n1, n2, n3, n4, n5, n6
	real*16 w1, w2, w3
	real*16 Gamma1
	f3star = Gamma1(n4 + n5 - 1, n6, n1 + n2 - 1, w1 + w2, w3, 0.0q0, tg3)
	return
end

real*16 function f4star(n1, n2, n3, n5, n6, w1, w2, w3)
	use gam
	implicit none
	integer n1, n2, n3, n4, n5, n6
	real*16 w1, w2, w3
	real*16 Gamma1
	if (n5 /= 0 .or. n6 /= 0) then
		write (*,*) "Error in input values for f4*"
		f4star = 0.0q0
		stop
	end if
	f4star = Gamma1(n3 - 1, n2 - 1, n1, w2, w3, 0.0q0, tg4)
	return
end

real*16 function f5star(n1, n2, n3, n4, n6, w1, w2, w3)
	use gam
	implicit none
	integer n1, n2, n3, n4, n5, n6
	real*16 w1, w2, w3
	real*16 Gamma1
	if (n4 /= 0 .or. n6 /= 0) then
		write (*,*) "Error in input values for f5"
		f5star = 0.0q0
		stop
	end if
	f5star = Gamma1(n1 - 1, n3 - 1, n2, w3, w1, 0.0q0, tg5)
	return
end

real*16 function f6star(n1, n2, n3, n4, n5, w1, w2, w3)
	use gam
	implicit none
	integer n1, n2, n3, n4, n5, n6
	real*16 w1, w2, w3
	real*16 Gamma1
	if (n4 /= 0 .or. n5 /= 0) then
		write (*,*) "Error in input values for f6*"
		f6star = 0.0q0
		stop
	end if
	f6star = Gamma1(n2 - 1, n1 - 1, n3, w1, w2, 0.0q0, tg6)
	return
end


real*16 function f0(w1, w2, w3)
	implicit none
	real*16 w1, w2, w3
	real*16 dilog

!	f0 = -(1.0q0/(2.0q0*w1*w2*w3)) * (qlog(w3/(w1 + w2)) * qlog(1.0q0 + w3/(w1 + w2)) + &
!		dilog(real(-(w3/(w1 + w2)),16)) + dilog(real(1.0q0 - w3/(w1 + w2),16)) + &
!		qlog(w2/(w3 + w1)) * qlog(1.0q0 + w2/(w1 + w3)) + &
!		dilog(real(-(w2/(w1 + w3)),16)) + dilog(real(1.0q0 - w2/(w1 + w3),16)) + &
!		qlog(w1/(w2 + w3)) * qlog(1.0q0 + w1/(w2 + w3)) + &
!		dilog(-real(w1/(w2 + w3),16)) + dilog(real(1.0q0 - w1/(w2 + w3),16)))

	f0 = -(1.0q0/(2.0q0*w1*w2*w3)) * (qlog(w3/(w1 + w2)) * qlog(1.0q0 + w3/(w1 + w2)) + &
		dilog(real(1.0q0+(w3/(w1 + w2)),16)) + dilog(real(w3/(w1 + w2),16)) + &
		qlog(w2/(w3 + w1)) * qlog(1.0q0 + w2/(w1 + w3)) + &
		dilog(real(1.0q0+(w2/(w1 + w3)),16)) + dilog(real(w2/(w1 + w3),16)) + &
		qlog(w1/(w2 + w3)) * qlog(1.0q0 + w1/(w2 + w3)) + &
		dilog(real(1.0q0+w1/(w2 + w3),16)) + dilog(real(w1/(w2 + w3),16)))

	return
end

real*16 function f1(w1, w2, w3) 
	implicit none
	real*16 w1, w2, w3
	f1 = -(1.0q0/(w2**2 * w3**2)) * qlog((w1*(w1 + w2 + w3))/((w1 + w2)*(w1 + w3)))
	return
end

real*16 function f2(w1, w2, w3)
	implicit none
	real*16 w1, w2, w3
	f2 = -(1.0q0/(w1**2 * w3**2)) * qlog((w2*(w1 + w2 + w3))/((w2 + w3)*(w2 + w1)))
	return
end

real*16 function f3(w1, w2, w3)
	implicit none
	real*16 w1, w2, w3
	f3 = -(1.0q0/(w1**2 * w2**2)) * qlog((w3*(w1 + w2 + w3))/((w3 + w1)*(w3 + w2)))
	return
end

real*16 function f12(w1, w2, w3)
	implicit none
	real*16 w1, w2, w3
	f12 = 1.0q0/(w1*w2*(w1 + w2) * w3**2)
	return
end

real*16 function f13(w1, w2, w3)
	implicit none
	real*16 w1, w2, w3
	f13 = 1.0q0/(w1*w3*(w1 + w3) * w2**2)
	return
end

real*16 function f23(w1, w2, w3)
	implicit none
	real*16 w1, w2, w3
	f23 = 1.0q0/(w2*w3*(w2 + w3) * w1**2)
	return
end

real*16 function f123(w1, w2, w3)
	implicit none
	real*16 w1, w2, w3
	f123 = 1.0q0/(w1**2 * w2**2 * w3**2)
	return
end


real*16 function Gamma1(n1, n2, n3, a1, a2, a3, tg)
	use PRMTS
	implicit none
	integer n1, n2, n3
	real*16 a1, a2, a3
	integer idx, sm
	real*16 tg(0:MAXG)

	sm = 10

	!call p_g(tg, real(a1,16), real(a2,16), real(a3,16), sm)

	Gamma1 = tg(idx(n1, n2, n3))
	
	!Gamma1 = 0.0q0

	return
end


integer function Delta(x)
	integer x
	if (x == 0) then
		Delta = 1
		return
	end if
	Delta = 0
	return
end


integer function Max(m1, m2, m3)
	integer m1, m2, m3

	Max = m1

	if (m2 > Max) then
		Max = m2
	end if
	if (m3 > Max) then
		Max = m3
	end if

	return
end


real*16 function fn1plus2(n11, n2, n3, n4, n5, n6, w1, w2, w3)
	implicit none
	integer n1, n11, n2, n3, n4, n5, n6
	real*16 w1, w2, w3
	real*16 f, f1star, f2star, f3star, f4star, f5star, f6star
	integer Delta
	real*16 r1, r2, r3

	n1 = n11 - 2

	r1 = ((n2*(n2 - 1))/real(n1 + 1,16) * f(n1 + 2, n2 - 2, n3, 0, 0, 0, w1, w2, w3) &
			+ (n3 + 2*n2 + n1 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
			+ 1.0q0/real(n1 + 1,16) * f6star(n1 + 2, n2, n3, 0, 0, w1, w2, w3) &
			+ (n2*(n2 - 1))/real(n3 + 1,16) * f(n1, n2 - 2, n3 + 2, 0, 0, 0, w1, w2, w3) &
			+ 1.0q0/real(n3 + 1,16) * f4star(n1, n2, n3 + 2, 0, 0, w1, w2, w3) &
			- Delta(n2)/real(n1 + 1,16) * f2star(n1 + 2, n3, 0, 0, 0, w1, w2, w3) &
			- Delta(n2)/real(n3 + 1,16) * f2star(n1, n3 + 2, 0, 0, 0, w1, w2, w3))

	r2 = ((n3*(n3 - 1))/real(n1 + 1,16) * f(n1 + 2, n2, n3 - 2, 0, 0, 0, w1, w2, w3) &
			+ (2*n3 + n2 + n1 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
			+ 1.0q0/real(n1 + 1,16) * f5star(n1 + 2, n2, n3, 0, 0, w1, w2, w3) &
			+ (n3*(n3 - 1))/real(n2 + 1,16) * f(n1, n2 + 2, n3 - 2, 0, 0, 0, w1, w2, w3) &
			+ 1.0q0/real(n2 + 1,16) * f4star(n1, n2 + 2, n3, 0, 0, w1, w2, w3) &
			- Delta(n3)/real(n1 + 1,16) * f3star(n1 + 2, n2, 0, 0, 0, w1, w2, w3) &
			- Delta(n3)/real(n2 + 1,16) * f3star(n1, n2 + 2, 0, 0, 0, w1, w2, w3))

	r3 = ((n1*(n1 - 1))/real(n2 + 1,16) * f(n1 - 2, n2 + 2, n3, 0, 0, 0, w1, w2, w3) &
			+ (n3 + n2 + 2*n1 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
			+ 1.0q0/real(n2 + 1,16) * f6star(n1, n2 + 2, n3, 0, 0, w1, w2, w3) &
			+ (n1 * (n1 - 1))/real(n3 + 1,16) * f(n1 - 2, n2, n3 + 2, 0, 0, 0, w1, w2, w3) &
			+ 1.0q0/real(n3 + 1,16) * f5star(n1, n2, n3 + 2, 0, 0, w1, w2, w3) &
			- Delta(n1)/real(n2 + 1,16) * f1star(n2 + 2, n3, 0, 0, 0, w1, w2, w3) &
			- Delta(n1)/real(n3 + 1,16) * f1star(n2, n3 + 2, 0, 0, 0, w1, w2, w3))

	fn1plus2 = (1.0q0 + n1)/2.0q0 * (1.0q0/w3**2 * r1 + 1.0q0/w2**2 * r2 - w1**2/(w3**2 * w2**2) * r3)
	return
end


real*16 function fn2plus2(n1, n21, n3, n4, n5, n6, w1, w2, w3)
	implicit none
	integer n1, n2, n21, n3, n4, n5, n6
	real*16 w1, w2, w3
	real*16 f, f1star, f2star, f3star, f4star, f5star, f6star
	integer Delta
	real*16 r1, r2, r3

	n2 = n21 - 2

	r1 = ((n3*(n3 - 1))/real(n2 + 1,16) * f(n1, n2 + 2, n3 - 2, 0, 0, 0, w1, w2, w3) &
		+ (n1 + 2*n3 + n2 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n2 + 1,16) * f4star(n1, n2 + 2, n3, 0, 0, w1, w2, w3) &
		+ (n3*(n3 - 1))/real(n1 + 1,16) * f(n1 + 2, n2, n3 - 2, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n1 + 1,16) * f5star(n1 + 2, n2, n3, 0, 0, w1, w2, w3) &
		- Delta(n3)/real(n2 + 1,16) * f3star(n1, n2 + 2, 0, 0, 0, w1, w2, w3) &
		- Delta(n3)/real(n1 + 1,16) * f3star(n1 + 2, n2, 0, 0, 0, w1, w2, w3))

	r2 = ((n1*(n1 - 1))/real(n2 + 1,16) * f(n1 - 2, n2 + 2, n3, 0, 0, 0, w1, w2, w3) &
		+ (2*n1 + n3 + n2 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n2 + 1,16) * f6star(n1, n2 + 2, n3, 0, 0, w1, w2, w3) &
		+ (n1*(n1 - 1))/real(n3 + 1,16) * f(n1 - 2, n2, n3 + 2, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n3 + 1,16) * f5star(n1, n2, n3 + 2, 0, 0, w1, w2, w3) &
		- Delta(n1)/real(n2 + 1,16) * f1star(n2 + 2, n3, 0, 0, 0, w1, w2, w3) &
		- Delta(n1)/real(n3 + 1,16) * f1star(n2, n3 + 2, 0, 0, 0, w1, w2, w3))

	r3 = ((n2*(n2 - 1))/real(n3 + 1,16) * f(n1, n2 - 2, n3 + 2, 0, 0, 0, w1, w2, w3) &
		+ (n1 + n3 + 2*n2 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n3 + 1,16) * f4star(n1, n2, n3 + 2, 0, 0, w1, w2, w3) &
		+ (n2*(n2 - 1))/real(n1 + 1,16) * f(n1 + 2, n2 - 2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n1 + 1,16) * f6star(n1 + 2, n2, n3, 0, 0, w1, w2, w3) &
		- Delta(n2)/real(n3 + 1,16) * f2star(n1, n3 + 2, 0, 0, 0, w1, w2, w3) &
		- Delta(n2)/real(n1 + 1,16) * f2star(n1 + 2, n3, 0, 0, 0, w1, w2, w3))

	fn2plus2 = (1.0q0 + n2)/2.0q0 * (1.0q0/w1**2 * r1 + 1.0q0/w3**2 * r2 - w2**2/(w1**2 * w3**2) * r3)

	return
end


real*16 function fn3plus2(n1, n2, n31, n4, n5, n6, w1, w2, w3)
	implicit none
	integer n1, n2, n3, n31, n4, n5, n6
	real*16 w1, w2, w3
	real*16 f, f1star, f2star, f3star, f4star, f5star, f6star
	integer Delta
	real*16 r1, r2, r3
	real*16 RemoveMe

	n3 = n31 - 2

	RemoveMe = f4star(n1, n2, n3 + 2, 0, 0, w1, w2, w3)

	r1 = ((n2*(n2 - 1))/real(n3 + 1,16) * f(n1, n2 - 2, n3 + 2, 0, 0, 0, w1, w2, w3) &
		+ (n1 + 2*n2 + n3 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n3 + 1,16) * f4star(n1, n2, n3 + 2, 0, 0, w1, w2, w3) &
		+ (n2*(n2 - 1))/real(n1 + 1,16) * f(n1 + 2, n2 - 2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n1 + 1,16) * f6star(n1 + 2, n2, n3, 0, 0, w1, w2, w3) &
		- Delta(n2)/real(n3 + 1,16) * f2star(n1, n3 + 2, 0, 0, 0, w1, w2, w3) &
		- Delta(n2)/real(n1 + 1,16) * f2star(n1 + 2, n3, 0, 0, 0, w1, w2, w3))

	r2 = ((n1*(n1 - 1))/real(n3 + 1,16) * f(n1 - 2, n2, n3 + 2, 0, 0, 0, w1, w2, w3) &
		+ (2*n1 + n2 + n3 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n3 + 1,16) * f5star(n1, n2, n3 + 2, 0, 0, w1, w2, w3) &
		+ (n1*(n1 - 1))/real(n2 + 1,16) * f(n1 - 2, n2 + 2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n2 + 1,16) * f6star(n1, n2 + 2, n3, 0, 0, w1, w2, w3) &
		- Delta(n1)/real(n3 + 1,16) * f1star(n2, n3 + 2, 0, 0, 0, w1, w2, w3) &
		- Delta(n1)/real(n2 + 1,16) * f1star(n2 + 2, n3, 0, 0, 0, w1, w2, w3))

	r3 = ((n3*(n3 - 1))/real(n2 + 1,16) * f(n1, n2 + 2, n3 - 2, 0, 0, 0, w1, w2, w3) &
		+ (n1 + n2 + 2*n3 + 2) * f(n1, n2, n3, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n2 + 1,16) * f4star(n1, n2 + 2, n3, 0, 0, w1, w2, w3) &
		+ (n3*(n3 - 1))/real(n1 + 1,16) * f(n1 + 2, n2, n3 - 2, 0, 0, 0, w1, w2, w3) &
		+ 1.0q0/real(n1 + 1,16) * f5star(n1 + 2, n2, n3, 0, 0, w1, w2, w3) &
		- Delta(n3)/real(n2 + 1,16) * f3star(n1, n2 + 2, 0, 0, 0, w1, w2, w3) &
		- Delta(n3)/real(n1 + 1,16) * f3star(n1 + 2, n2, 0, 0, 0, w1, w2, w3))

	fn3plus2 = (1 + n3)/2.0q0 * (1.0q0/w1**2 * r1 + 1.0q0/w2**2 * r2 - w3**2/(w1**2 * w2**2) * r3)

	return
end


real*16 function fn4plus1(n1, n2, n3, n41, n5, n6, w1, w2, w3)
	implicit none
	integer n1, n2, n3, n4, n41, n5, n6
	real*16 w1, w2, w3
	real*16 f, f1star, f2star, f3star
	integer Delta
	real*16 r1, r2, r3, r4, r5

	n4 = n41 - 1

	r1 = (n3 - 1)*n3*n6*f(n1, n2, n3 - 2, n4, n5 + 1, n6 - 1, w1, w2, w3) &
		+ (n2 - 1)*n2*n5*f(n1, n2 - 2, n3, n4, n5 - 1, n6 + 1, w1, w2, w3) &
		- (n1 - 1)*n1*n6*f(n1 - 2, n2, n3, n4, n5 + 1, n6 - 1, w1, w2, w3) &
		- (n1 - 1)*n1*n5*f(n1 - 2, n2, n3, n4, n5 - 1, n6 + 1, w1, w2, w3) &
		+ n6*n5*(n3 + n2 - n1 + n4 + 1) * f(n1, n2, n3, n4, n5 - 1, n6 - 1, w1, w2, w3)

	r2 = -(n3 - 1)*n3*f(n1, n2, n3 - 2, n4, n5 + 1, n6, w1, w2, w3) * w3 &
		+ (n1 - 1)*n1*f(n1 - 2, n2, n3, n4, n5 + 1, n6, w1, w2, w3) * w3 &
		- n5*(n3 + n2 - n1 + n4 + 1) * f(n1, n2, n3, n4, n5 - 1, n6, w1, w2, w3) * w3 &
		- (n2 - 1)*n2*f(n1, n2 - 2, n3, n4, n5, n6 + 1, w1, w2, w3) * w2 &
		+ (n1 - 1)*n1*f(n1 - 2, n2, n3, n4, n5, n6 + 1, w1, w2, w3) * w2

	r3 = -n6 * (n3 + n2 - n1 + n4 + 1) * f(n1, n2, n3, n4, n5, n6 - 1, w1, w2, w3) * w2 &
		- n6*n5*f(n1, n2, n3, n4 + 1, n5 - 1, n6 - 1, w1, w2, w3) * w1 &
		+ (n3 + n2 - n1 + n4 + 1) * f(n1, n2, n3, n4, n5, n6, w1, w2, w3) * w3 * w2 &
		+ n6 * f(n1, n2, n3, n4 + 1, n5, n6 - 1, w1, w2, w3) * w2 * w1 & 
		+ n5 * f(n1, n2, n3, n4 + 1, n5 - 1, n6, w1, w2, w3) * w1 * w3

	r4 = -Delta(n3) * n6 * f3star(n1, n2, n4, n5 + 1, n6 - 1, w1, w2, w3) &
		+ Delta(n3) * f3star(n1, n2, n4, n5 + 1, n6, w1, w2, w3) * w3 &
		- Delta(n2) * n5 * f2star(n1, n3, n4, n5 - 1, n6 + 1, w1, w2, w3) &
		+ Delta(n2) * f2star(n1, n3, n4, n5, n6 + 1, w1, w2, w3) * w2 &
		+ Delta(n1) * n6 * f1star(n2, n3, n4, n5 + 1, n6 - 1, w1, w2, w3)

	r5 = Delta(n1) * n5 * f1star(n2, n3, n4, n5 - 1, n6 + 1, w1, w2, w3) &
		- Delta(n1) * f1star(n2, n3, n4, n5 + 1, n6, w1, w2, w3) * w3 &
		- Delta(n1) * f1star(n2, n3, n4, n5, n6 + 1, w1, w2, w3) * w2

	fn4plus1 = 1.0q0/(w3*w2*w1) * (r1 + r2 + r3 + r4 + r5)

	return
end


real*16 function fn5plus1(n1, n2, n3, n4, n51, n6, w1, w2, w3)
	implicit none
	integer n1, n2, n3, n4, n5, n51, n6
	real*16 w1, w2, w3
	real*16 f, f1star, f2star, f3star
	integer Delta
	real*16 r1, r2, r3, r4, r5

	n5 = n51 - 1

	r1 = (n1 - 1)*n1*n4*f(n1 - 2, n2, n3, n4 - 1, n5, n6 + 1, w1, w2, w3) &
		+ (n3 - 1)*n3*n6*f(n1, n2, n3 - 2, n4 + 1, n5, n6 - 1,  w1, w2, w3) &
		- (n2 - 1)*n2*n4*f(n1, n2 - 2, n3, n4 - 1, n5, n6 + 1, w1, w2, w3) &
		- (n2 - 1)*n2*n6*f(n1, n2 - 2, n3, n4 + 1, n5, n6 - 1, w1, w2, w3) &
		+ n4*n6*(n1 + n3 - n2 + n5 + 1) * f(n1, n2, n3, n4 - 1, n5, n6 - 1, w1, w2, w3)

	r2 = -(n1 - 1)*n1*f(n1 - 2, n2, n3, n4, n5, n6 + 1, w1, w2, w3) * w1 &
		+ (n2 - 1)*n2*f(n1, n2 - 2, n3, n4, n5, n6 + 1, w1, w2, w3) * w1 &
		- n6*(n1 + n3 - n2 + n5 + 1)*f(n1, n2, n3, n4, n5, n6 - 1, w1, w2, w3) * w1 &
		- (n3 - 1)*n3*f(n1, n2, n3 - 2, n4 + 1, n5, n6, w1, w2, w3) * w3 &
		+ (n2 - 1)*n2*f(n1, n2 - 2, n3, n4 + 1, n5, n6, w1, w2, w3) * w3

	r3 = -n4 * (n1 + n3 - n2 + n5 + 1) * f(n1, n2, n3, n4 - 1, n5, n6, w1, w2, w3) * w3 &
		- n4*n6*f(n1, n2, n3, n4 - 1, n5 + 1, n6 - 1, w1, w2, w3) * w2 &
		+ (n1 + n3 - n2 + n5 + 1) * f(n1, n2, n3, n4, n5, n6, w1, w2, w3)*w1*w3 &
		+ n4*f(n1, n2, n3, n4 - 1, n5 + 1, n6, w1, w2, w3)*w3*w2 &
		+ n6*f(n1, n2, n3, n4, n5 + 1, n6 - 1, w1, w2, w3)*w2*w1

	r4 = -Delta(n1) * n4 * f1star(n2, n3, n4 - 1, n5, n6 + 1, w1, w2, w3) &
		+ Delta(n1) * f1star(n2, n3, n4, n5, n6 + 1, w1, w2, w3) * w1 &
		- Delta(n3) * n6 * f3star(n1, n2, n4 + 1, n5, n6 - 1, w1, w2, w3) &
		+ Delta(n3) * f3star(n1, n2, n4 + 1, n5, n6, w1, w2, w3) * w3 &
		+ Delta(n2) * n4 * f2star(n1, n3, n4 - 1, n5, n6 + 1, w1, w2, w3)

	r5 = Delta(n2) * n6 * f2star(n1, n3, n4 + 1, n5, n6 - 1, w1, w2, w3) &
		- Delta(n2) * f2star(n1, n3, n4, n5, n6 + 1, w1, w2, w3) * w1 &
		- Delta(n2) * f2star(n1, n3, n4 + 1, n5, n6, w1, w2, w3) * w3

	fn5plus1 = 1.0q0/(w1*w3*w2) * (r1 + r2 + r3 + r4 + r5)

	return
end


real*16 function fn6plus1(n1, n2, n3, n4, n5, n61, w1, w2, w3)
	implicit none
	integer n1, n2, n3, n4, n5, n6, n61
	real*16 w1, w2, w3
	real*16 f, f1star, f2star, f3star
	integer Delta
	real*16 r1, r2, r3, r4, r5

	n6 = n61 - 1

	r1 = (n1 - 1)*n1*n4*f(n1 - 2, n2, n3, n4 - 1, n5 + 1, n6, w1, w2, w3) &
		+ (n2 - 1)*n2*n5*f(n1, n2 - 2, n3, n4 + 1, n5 - 1, n6, w1, w2, w3) &
		- (n3 - 1)*n3*n4*f(n1, n2, n3 - 2, n4 - 1, n5 + 1, n6, w1, w2, w3) &
		- (n3 - 1)*n3*n5*f(n1, n2, n3 - 2, n4 + 1, n5 - 1, n6, w1, w2, w3) &
		+ n4*n5*(n1 + n2 - n3 + n6 + 1) * f(n1, n2, n3, n4 - 1, n5 - 1, n6, w1, w2, w3)

	r2 = -(n1 - 1)*n1*f(n1 - 2, n2, n3, n4, n5 + 1, n6, w1, w2, w3) * w1 &
		+ (n3 - 1)*n3*f(n1, n2, n3 - 2, n4, n5 + 1, n6, w1, w2, w3) * w1 &
		- n5 * (n1 + n2 - n3 + n6 + 1) * f(n1, n2, n3, n4, n5 - 1, n6, w1, w2, w3) * w1 &
		- (n2 - 1) * n2 * f(n1, n2 - 2, n3, n4 + 1, n5, n6, w1, w2, w3) * w2 &
		+ (n3 - 1) * n3 * f(n1, n2, n3 - 2, n4 + 1, n5, n6, w1, w2, w3) * w2

	r3 = -n4 * (n1 + n2 - n3 + n6 + 1) * f(n1, n2, n3, n4 - 1, n5, n6, w1, w2, w3) * w2 &
		- n4 * n5 * f(n1, n2, n3, n4 - 1, n5 - 1, n6 + 1, w1, w2, w3) * w3 &
		+ (n1 + n2 - n3 + n6 + 1) * f(n1, n2, n3, n4, n5, n6, w1, w2, w3) * w1 * w2 &
		+ n4 * f(n1, n2, n3, n4 - 1, n5, n6 + 1, w1, w2, w3) * w2 * w3 &
		+ n5 * f(n1, n2, n3, n4, n5 - 1, n6 + 1, w1, w2, w3) * w3 * w1

	r4 = -Delta(n1) * n4 * f1star(n2, n3, n4 - 1, n5 + 1, n6, w1, w2, w3) &
		+ Delta(n1) * f1star(n2, n3, n4, n5 + 1, n6, w1, w2, w3) * w1 &
		- Delta(n2) * n5 * f2star(n1, n3, n4 + 1, n5 - 1, n6, w1, w2, w3) &
		+ Delta(n2) * f2star(n1, n3, n4 + 1, n5, n6, w1, w2, w3) * w2 & 
		+ Delta(n3) * n4 * f3star(n1, n2, n4 - 1, n5 + 1, n6, w1, w2, w3)

	r5 = Delta(n3) * n5 * f3star(n1, n2, n4 + 1, n5 - 1, n6, w1, w2, w3) &
		- Delta(n3) * f3star(n1, n2, n4, n5 + 1, n6, w1, w2, w3) * w1 &
		- Delta(n3) * f3star(n1, n2, n4 + 1, n5, n6, w1, w2, w3) * w2

	fn6plus1 = 1.0q0/(w1*w2*w3) * (r1 + r2 + r3 + r4 + r5)

	return
end


real*16 function f(n1, n2, n3, n4, n5, n6, w1, w2, w3)
	use gam
	implicit none
	integer n1, n2, n3, n4, n5, n6, verbose
	real*16 w1, w2, w3
	integer Max
	real*16 f0, f1, f2, f3, f12, f13, f23, f123, fn4plus1, fn5plus1, fn6plus1, fn1plus2, fn2plus2, fn3plus2

	if (n1 < 0 .or. n2 < 0 .or. n3 < 0 .or. n4 < 0 .or. n5 < 0 .or. n6 < 0) then
		!if (verbose == 1) then
		!	write (*,*) "Input to f out of range: ", n1, n2, n3, n4, n5, n6
			! Doesn't really matter what the output here is, because a 
			!  coefficient will always be 0 in this case
			f = 1.0q0
			return
		!end if
	end if
 
	! At the end of the recursion - initial terms
	if (n1 == 0 .and. n2 == 0 .and. n3 == 0 .and. n4 == 0 .and. n5 == 0 .and. n6 == 0) then
		!f = f0(w1, w2, w3)
		f = f0val
		return
	end if
	if (n1 == 1 .and. n2 == 0 .and. n3 == 0 .and. n4 == 0 .and. n5 == 0 .and. n6 == 0) then
		!f = f1(w1, w2, w3)
		f = f1val
		return
	endif
	if (n1 == 0 .and. n2 == 1 .and. n3 == 0 .and. n4 == 0 .and. n5 == 0 .and. n6 == 0) then
		!f = f2(w1, w2, w3)
		f = f2val
		return
	endif
	if (n1 == 0 .and. n2 == 0 .and. n3 == 1 .and. n4 == 0 .and. n5 == 0 .and. n6 == 0) then
		!f = f3(w1, w2, w3)
		f = f3val
		return
	endif
	if (n1 == 1 .and. n2 == 1 .and. n3 == 0 .and. n4 == 0 .and. n5 == 0 .and. n6 == 0) then
		!f = f12(w1, w2, w3)
		f = f12val
		return
	endif
	if (n1 == 1 .and. n2 == 0 .and. n3 == 1 .and. n4 == 0 .and. n5 == 0 .and. n6 == 0) then
		!f = f13(w1, w2, w3)
		f = f13val
		return
	endif
	if (n1 == 0 .and. n2 == 1 .and. n3 == 1 .and. n4 == 0 .and. n5 == 0 .and. n6 == 0) then
		!f = f23(w1, w2, w3)
		f = f23val
		return
	endif
	if (n1 == 1 .and. n2 == 1 .and. n3 == 1 .and. n4 == 0 .and. n5 == 0 .and. n6 == 0) then
		!f = f123(w1, w2, w3)
		f = f123val
		return
	endif

	! Reduce n4, n5 and n6 down to 0 first
	if (Max(n4, n5, n6) > 0 .and. n4 == Max(n4, n5, n6)) then
		f = fn4plus1(n1, n2, n3, n4, n5, n6, w1, w2, w3)
		return
	endif
	if (Max(n4, n5, n6) > 0 .and. n5 == Max(n4, n5, n6)) then
		f = fn5plus1(n1, n2, n3, n4, n5, n6, w1, w2, w3)
		return
	endif
	if (Max(n4, n5, n6) > 0 .and. n6 == Max(n4, n5, n6)) then
		f = fn6plus1(n1, n2, n3, n4, n5, n6, w1, w2, w3)
		return
	endif
	! Then reduce n1, n2 and n3 down to the initial terms
	if (Max(n1, n2, n3) > 1 .and. n1 == Max(n1, n2, n3)) then
		f = fn1plus2(n1, n2, n3, n4, n5, n6, w1, w2, w3)
		return
	endif
	if (Max(n1, n2, n3) > 1 .and. n2 == Max(n1, n2, n3)) then
		f = fn2plus2(n1, n2, n3, n4, n5, n6, w1, w2, w3)
		return
	endif
	if (Max(n1, n2, n3) > 1 .and. n3 == Max(n1, n2, n3)) then
		f = fn3plus2(n1, n2, n3, n4, n5, n6, w1, w2, w3)
		return
	endif

	! Shouldn't get here - could print a warning message.
	f = 0.0q0

	return
end
