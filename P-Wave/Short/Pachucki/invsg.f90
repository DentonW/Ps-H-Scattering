!subroutine Newdsygv(1, 'V', 'L', NumTerms, Phi2HPhi8, NumTerms, PhiPhi8, NumTerms, Energies, Workspace, 3*NumTerms-1, Info)
subroutine Newdsygv(ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, INFO)
	implicit none
	CHARACTER          JOBZ, UPLO
	INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
	DOUBLE PRECISION   W( * ), WORK( * )
	real*16 A( LDA, * ), B( LDB, * )
	!real*16 A( LDA, * ), B( LDB, * )
	!DOUBLE PRECISION   W( * ), WORK( * )
!	real*16 AP(N*(N+1)/2), BP(N*(N+1)/2)
	real*16, allocatable, dimension(:) :: AP, BP, WA
	real*16 EE, EPS
	real*16 WF(N)
!	real*16 WA(3*N+1)
	integer, PARAMETER :: ITMAX = 40
	integer i, j, k
	
	allocate(AP(N*(N+1)/2))
	allocate(BP(N*(N+1)/2))
	allocate(WA(3*N+1))

	!write (*,*) "Entering Newdsygv"

	WA = 0.0q0

	EE = -0.7891967147q0  ! Guess for the first eigenvalue

	! Put into symmetric storage mode - from
	!  http://publib.boulder.ibm.com/infocenter/clresctr/vxrx/index.jsp?topic=%2Fcom.ibm.cluster.essl.v5r1.essl100.doc%2Fam501_upsm.html
	K = 0
	DO J=1,N
		DO I=1,J
			K = K+1
			AP(K)=A(I,J)
			BP(K)=B(I,J)
		end do
	end do

	DO I=1,N
		WF(I) = 1.0q0
	ENDDO
	EPS = 0.0q0

	call invsg(AP, BP, N, EE, WF, EPS, 0, ITMAX, WA)

	W(1) = EE

	Info = 0

	deallocate(AP, BP, WA)

	!write (*,*) "Exiting Newdsygv"

	return
end


subroutine invsg(a,b,n,eig,v,eps,ijob,maxit,wa)
!-------------------------------------------------------------
!
!   latest revision     - april 23, 1992
!
!   purpose             - inverse iteration for the
!                           generalized eigenvalue problem
!                           a*x = x*b*x - symmetric storage
!                           mode.
!
!   usage               - call invsg(a,b,n,eig,v,eps,ijob,maxit,wa)
!
!   arguments    a      - linear array, on input a contains
!                           the elements of matrix a in
!                           symmetric storage mode, the length
!                           of a is n*(n+1)/2, on output a is
!                           destroyed.
!                b      - linear array, on input b contains the
!                           elements of matrix b in symmetric
!                           storage mode.
!                n      - dimension of matrices a and b.
!                eig    - on input, eig contains the initial
!                           approximation to the eigenvalue,
!                         on output, eig contains renewed eigenvalue.
!                v      - on input,  v  contains the initial
!                           approximation to the eigenvector,
!                         on output,  v  contains renewed eigenvector
!                           with unit b-norm, (x,b*x)=1.
!                eps    - parameter of regularization,
!                           if eig is one of the lowest
!                             eigenvalues, then eps must be
!                             greater than or equal to zero,
!                           if eig is one of the uppest
!                             eigenvalues, then eps must be
!                             less than or equal to zero,
!                           if the problem is well defined,
!                             then it is better to set the
!                             value of eps to zero.
!                ijob   - job option parameter,
!                           if ijob = 0, then form matrix (a-xb)
!                             and prepare it for iterations,
!                           if ijob = 1, the matrix (a-xb)
!                             was previously factored by routine
!                             leq1s.
!                maxit  - maximal number of iterations (input).
!                wa     - work array of length 3*n+1
!                           on output, wa(1) contains
!                           the number of eigenvalues below
!                           the initial value of eig.
!
!   reqd. routines      - leq1s, avms
!
!-------------------------------------------------------------
      implicit real (16) (a-h,o-z)
      dimension           a(*),b(*),v(*),wa(*)
      real (16)      :: zero=0.0q0, one=1.0q0
      real (16)      :: sixtn=16.0q0
      real (16)      :: cscale=256.q0
      save                zero,one,sixtn,cscale
!
      ndet = 2*n+1
      aln256 = log(cscale)
      reps = epsilon(reps)
      errest = max(sixtn*n*reps,abs(eps))
!
!                         replace matrix a with matrix
!                           (a-eig*b)
!
      eold = eig
      nsc = n
      if (ijob /= 1) then
         ii = 0
         do i=1,n
            do j=1,i
               ii = ii+1
               a(ii) = a(ii)-eig*b(ii)
            end do
            temp = max(abs(a(ii)),abs(eig*b(ii)))
            isc = nint(log(temp)/aln256)
            wa(nsc+i) = sixtn**(-isc)
            do j=0,i-1
               a(ii-j) = a(ii-j)*wa(nsc+i)*wa(nsc+i-j)
               b(ii-j) = b(ii-j)*wa(nsc+i)*wa(nsc+i-j)
            end do
            a(ii) = a(ii)+eps*abs(a(ii))
         end do
!
         call leq1s(a,n,v,1,n,1,wa(ndet),ier)
!
      end if
!
!                         begin iterations
!
      do it=1,maxit
!
!                         multiply on matrix b and normalize
!                           vector v to unity
!
         call avms(b,n,v,wa)
         sm = dot_product(v(1:n),v(1:n))
         sm = one/sqrt(sm)
         do i=1,n
            t = wa(i)*sm
            wa(i) = v(i)*sm
            v(i) = t
         end do
!
!                                                     -1
!                         multiply on matrix (a-eig*b)
!
         call leq1s(a,n,v,1,n,2,wa(ndet),ier)
!
!                         calculate refined eigenvalue
!
         sm = dot_product(v(1:n),wa(1:n))
         eprev = eig
         eig = eold+one/sm
         if (abs(eprev-eig).lt.errest*abs(eig)) goto 100
!
      end do
!
!                         normalize eigenvector
!
  100 call avms(b,n,v,wa)
      sm = dot_product(wa(1:n),v(1:n))
      sm = one/sqrt(sm)
      v(1:n) = v(1:n)*wa(nsc+1:nsc+n)*sm
      wa(1) = wa(ndet+n)
!
!                        exit
!
      return
      end
      subroutine avms(a,n,v,vm)
!------------------------------------------------------------
!
!   purpose             - multiply a vector v on matrix a
!                           stored in symmetric mode
!                           - nucleus for invsg.
!
!   reqd. routines      - none required.
!
!-------------------------------------------------------------
      implicit real (16) (a-h,o-z)
      dimension           v(*),vm(*),a(*)
      real (16)      :: zero=0.0q0
      save                zero
!
      vm(1) = v(1)*a(1)
      ii = 1
      do i=2,n
         im1 = i-1
         vm(i) = zero
         do j=1,im1
            vm(i) = vm(i)+a(ii+j)*v(j)
            vm(j) = vm(j)+a(ii+j)*v(i)
         end do
         ii = ii+i
         vm(i) = vm(i)+a(ii)*v(i)
      end do
!
      return
      end
