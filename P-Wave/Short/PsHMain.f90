! This program calculates the matrix elements of the positronium-hydrogen P-wave
!  matrices (similar to equation 2.15 of the Armour and Humberston article).  Specifically,
!  we are calculating elements of the form (phi_i, L phi_j) as in equation (3.22).
! This does use OpenMP to speed up computation on multicore processors, so some familiarity
!  with parallel programming is required to understand the !$omp lines.

module WCLimits
   implicit none
   integer lmin, lmax
   integer mmin, mmax
   integer nmin, nmax
   integer Cjmax, Cqmax, Ckmax
end module WCLimits


program PsHMain
	implicit none
	
	interface
		integer(c_int) function optimizewavefn(MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, alpha, beta, gamma) bind(c, name='optimizewavefn')
		use iso_c_binding
		integer(c_int), value :: MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method
		real(c_double), value :: alpha, beta, gamma
		end function optimizewavefn
	end interface

	real*8 omp_get_wtime
	integer Omega  ! This sets the limits on the terms.
	integer NumTerms  ! Number of terms in our Hylleraas expansion
	integer iread, iwrite  ! Our input and output files
	real*16, allocatable, dimension(:,:) :: PhiPhi, Phi2HPhi
	real*16 Alpha, Beta, Gamma
	logical ReadParamFile, ReadResult
	integer, allocatable, dimension(:) :: UsedTerms
	integer AllEnergies, Optimize, EigenNum, Ordering, Iter, MaxIter, qmax, pmax
	integer iargc, IsTriplet, Method, LowerEigen, UpperEigen
	integer CalcPowerTableSize, NumUsed, i, j, k, n, fa
	real*8 StartTime, EndTime
	character *100 IOBuffer
	
	!call omp_set_num_threads(3)

	iread = 9
	iwrite = 10
	call GenFactorial(400)

	! This allows the possibility of using different input and output files than the defaults.
	!  This mainly lets me run this program on multiple computers at the same time without
	!  one overwriting the results from another.  If the program is called without any
	!  command-line arguments, it uses the default filenames.
	if (iargc() == 2) then
		call getarg(1, IOBuffer)
		open (iread, FILE=IOBuffer)
		call getarg(2, IOBuffer)
		open (iwrite, FILE=IOBuffer)
	else
		open (iread, FILE='input.txt')
		open (iwrite, FILE='output.xml')
	endif

	!@TODO: Should probably check the output
	ReadResult = ReadParamFile(iread, Omega, Alpha, Beta, Gamma, qmax, Method, pmax, IsTriplet, LowerEigen, UpperEigen, &
								AllEnergies, Optimize, EigenNum, MaxIter, Ordering)

	! Get the start time to determine the duration of the program.
	!call cpu_time(StartTime)
	StartTime = omp_get_wtime()

	NumTerms = CalcPowerTableSize(Omega)
	allocate(PhiPhi(NumTerms,NumTerms))
	allocate(Phi2HPhi(NumTerms,NumTerms))

	call WriteHeader(iwrite, Optimize, EigenNum, IsTriplet, Ordering, Method, Omega, Alpha, Beta, Gamma, NumTerms)
	
	if (Optimize == 1) then
		fa = optimizewavefn(MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, dble(Alpha), dble(Beta), dble(Gamma))
		call CalcMatrices(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, qmax, pmax, &
					IsTriplet, Ordering, Method)
	else
		! If we are not optimizing eigenvalues, we can still call the same function, which exits early with PhiPhi and Phi2HPhi filled.
		call CalcMatrices(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, qmax, pmax, &
					IsTriplet, Ordering, Method)
		write (*,*) fa
	endif
	
	call WriteData(iwrite, NumTerms, PhiPhi, Phi2HPhi)

	! Divide every entry in Phi2HPhi by 2, since we calculated <phi_i|2H|phi_j> in equation (3.22).
	Phi2HPhi = Phi2HPhi / 2

	call CalcEigenvalues(iwrite, NumTerms, AllEnergies, PhiPhi, Phi2HPhi, LowerEigen, UpperEigen)
	
	! Get the end time to find the duration of the program.
	!call cpu_time(EndTime)
	EndTime = omp_get_wtime()
	write (*,*) 'Time taken (s):', EndTime - StartTime
	write (iwrite,'(a)') "<runtime>"
	write (iwrite,*) '    (s):', EndTime - StartTime, '  (min):', (EndTime - StartTime) / 60.0
	write (iwrite,'(a)') "</runtime>"
	write (iwrite,'(a)') "</psh_data>"

	! Clean up memory before exiting
	deallocate(PhiPhi, Phi2HPhi)

	! Close file handles
	close(iwrite)
stop
end


! Read the input parameters.  It is assumed that these parameters are all valid (no real need to write parameter-checking code).
logical function ReadParamFile(iread, Omega, Alpha, Beta, Gamma, qmax, Method, pmax, IsTriplet, LowerEigen, UpperEigen, &
								AllEnergies, Optimize, EigenNum, MaxIter, Ordering)
	implicit none
	integer iread, Omega, qmax, Method, pmax, IsTriplet, LowerEigen, UpperEigen
	integer AllEnergies, Optimize, EigenNum, MaxIter, Ordering
	real*16 Alpha, Beta, Gamma
	
	read (iread,*) ! Description line
	read (iread,*) Omega
	read (iread,*) ! Description line
	read (iread,*) Alpha
	read (iread,*) Beta
	read (iread,*) Gamma
	read (iread,*) ! Description line
	read (iread,*) qmax
	read (iread,*) ! Description line
	read (iread,*) Method
	read (iread,*) ! Description line
	read (iread,*) pmax
	read (iread,*) ! Description line
	read (iread,*) IsTriplet
	read (iread,*) ! Description line
	read (iread,*) LowerEigen, UpperEigen
	read (iread,*)
	read (iread,*) AllEnergies
	read (iread,*) 
	read (iread,*) Optimize
	read (iread,*)
	read (iread,*) EigenNum
	read (iread,*)
	read (iread,*) MaxIter
	read (iread,*)
	read (iread,*) Ordering

	close(iread)

	if (IsTriplet /= 0 .and. IsTriplet /= 1) then
		write (*,*) 'Must choose between singlet and triplet calculation - exiting.'
		stop
	endif
	
	ReadParamFile = .true.
	return
end


! Outputs header information to the output file and to the screen
subroutine WriteHeader(iwrite, Optimize, EigenNum, IsTriplet, Ordering, Method, Omega, Alpha, Beta, Gamma, NumTerms)
	implicit none
	integer iwrite, Optimize, EigenNum, IsTriplet, Ordering, Method, Omega, NumTerms
	real*16 Alpha, Beta, Gamma
	character (len=8) cdate
	character (len=8) ctime

	if (Optimize == 1) write (iwrite,*) 'Optimizing eigenvalue - not a valid XML file!', EigenNum

	write (iwrite,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
	write (iwrite,'(a)') "<psh_data>"
	write (iwrite,'(a)') "<header>"
	
	if (IsTriplet == 0) then
		write (*,*) "P-Wave Singlet Ps-H: 1st formalism"
		write (iwrite,'(a)') "    <problem>P-Wave Singlet Ps-H: 1st formalism</problem>"
	else
		write (*,*) "P-Wave Triplet Ps-H: 1st formalism"
		write (iwrite,'(a)') "    <problem>P-Wave Triplet Ps-H: 1st formalism</problem>"
	endif

	write (iwrite,'(a)') "    <lvalue>1</lvalue>"  ! Indicates P-wave
	if (IsTriplet == 0) then
		write (iwrite,'(a,i0,a)') "    <spin>Singlet</spin>"
	else
		write (iwrite,'(a,i0,a)') "    <spin>Triplet</spin>"
	end if

	if (Ordering == 1) then
		write (*,*) "Using Peter Van Reeth's ordering"
		write (iwrite,'(a)') "    <ordering>Peter</ordering>"
	else
		write (*,*) "Using my ordering"
		write (iwrite,'(a)') "    <ordering>Denton</ordering>"
	endif

	write (iwrite,'(a)') "    <gradlapl>Gradient-Gradient</gradlapl>"  ! Only version of this code
	write (iwrite,'(a)') "    <formalism>1</formalism>"  ! Only does the first formalism
	
	if (Method == 0) then
		write (*,*) "Integration Technique: Direct Summation"
		write (iwrite,'(a)') "    <shortint>Direct summation</shortint>"
	else if (Method == 1) then
		write (*,*) "Integration Technique: Asymptotic Expansion"
		write (iwrite,'(a)') "    <shortint>Asymptotic Expansion</shortint>"
	else if (Method == 2) then
		write (*,*) "Integration Technique: Recursion Relations"
		write (iwrite,'(a)') "    <shortint>Recursion Relations</shortint>"
	else
		write (*,*) "Method parameter in input file must be 0, 1 or 2."
		stop
	end if
	
	write (*,*) 'Number of terms:', NumTerms
	write (*,*) 'Omega =', Omega
	write (iwrite,'(a,i0,a)') "    <omega>", Omega, "</omega>"
	write (iwrite,'(a,i0,a)') "    <numterms>", NumTerms, "</numterms>"
	write (iwrite,'(a)') "    <numsets>1</numsets>"  ! Can only have one set for the S-wave
	
	!write (iwrite,'(a)') "    <extraexp>false</extraexp>"  ! Not doing extra exponentials here
	write (iwrite,'(a)') "    <nonlinear>"
	write (iwrite,'(a,f12.10,a)') "        <alpha>", Alpha, "</alpha>"
	write (iwrite,'(a,f12.10,a)') "        <beta>", Beta, "</beta>"
	write (iwrite,'(a,f12.10,a)') "        <gamma>", Gamma, "</gamma>"
	write (iwrite,'(a)') "    </nonlinear>"
	
	write (iwrite,'(a)') "    <eigenroutine>LAPACK</eigenroutine>"

	! Note that this is not Y2K compliant.
	call date(cdate)
	call time(ctime)
	write (iwrite,'(5a)') "    <datetime>", cdate, " ", ctime, "</datetime>"

	write (iwrite,'(a)') "</header>"

	return
end


! Outputs only the matrix elements (before division of Phi2HPhi by 2)
subroutine WriteData(iwrite, NumTerms, PhiPhi, Phi2HPhi)
	implicit none
	integer iwrite, NumTerms, i, j
	real*16, dimension(NumTerms,NumTerms) :: PhiPhi, Phi2HPhi

	write (iwrite,'(a)') "<data>"
	do i = 1, NumTerms, 1
		!do j = 1, i, 1  ! Use this instead if only the lower triangle is required.
		do j = 1, NumTerms, 1
			write (iwrite,"(a,i7,i7,d38.30,d38.30)") '    ', i, j, PhiPhi(i,j), Phi2HPhi(i,j)
		enddo
	enddo
	write (iwrite,'(a)') "</data>"

	return
end


! Solves the generalized eigenvalue problem and outputs results
subroutine CalcEigenvalues(iwrite, NumTerms, AllEnergies, PhiPhi, Phi2HPhi, LowerEigen, UpperEigen)
	implicit none
	integer iwrite, NumTerms, AllEnergies, StartEnergy, LowerEigen, UpperEigen
	real*16, dimension(NumTerms,NumTerms) :: PhiPhi, Phi2HPhi
	real*8, allocatable, dimension(:,:) :: PhiPhi8, Phi2HPhi8
	real*8, allocatable, dimension(:) :: Energies, Workspace
	integer i, j, Info
	
	allocate(PhiPhi8(NumTerms,NumTerms))
	allocate(Phi2HPhi8(NumTerms,NumTerms))
	allocate(Energies(NumTerms))
	allocate(Workspace(3*NumTerms-1))
	
	write (iwrite,'(a)') '<energies>'
	
	if (AllEnergies == 0) then
		StartEnergy = NumTerms
	else
		StartEnergy = 1
	end if

	if (UpperEigen > NumTerms) UpperEigen = NumTerms

	do j = StartEnergy, NumTerms, 1
	!do j = 1, 1, 1
		! Copies our real*8 matrices to real*8 matrices so that LAPACK can operate on them.
		Phi2HPhi8 = Phi2HPhi
		PhiPhi8 = PhiPhi

		! This calculates the energy eigenvalues of the generalized eigenvalue problem generated from
		!  the Rayleigh-Ritz variational method:
		!  det(<phi_i|H|phi_j> - E <phi_i|phi_j>) = 0
		! We calculated the lower and upper triangles of the PhiPhi and Phi2HPhi matrices above, though
		!  they are symmetric.  LAPACK does not require us to fill in both halves; we can just
		!  specify 'L' in the third parameter to denote that the lower triangle (and diagonal) is filled.
		! Explanation of choices for the parameters are found at the Netlib site:
		!  http://www.netlib.org/lapack/double/dsygv.f
		call dsygv(1, 'N', 'L', j, Phi2HPhi8, NumTerms, PhiPhi8, NumTerms, Energies, Workspace, 3*NumTerms-1, Info)
		!call nag_sym_gen_eig_all('L', Phi2HPhi8, PhiPhi8, Energies)  ! NAG equivalent
		if (Info /= 0) then
			write (*,*)
			write (*,*) 'Energy eigenvalues could not be determined!'
			write (*,*) 'dsygv error code:', Info
			write (iwrite,*)
			write (iwrite,*) 'Energy eigenvalues could not be determined!'
			write (iwrite,*) 'dsygv error code:', Info
			exit  ! Exit the loop, because every successive j will give the same error, since we have a troublesome term.
		elseif (Info == 0) then  ! dsygv successfully completed.
			! Writes the eigenvalues.
			! This output is formatted properly for inclusion into Excel as a comma-separated value (.csv) file.
			write (iwrite,"(a,i4)",advance='no') '    ', j
			do i = LowerEigen, UpperEigen, 1
				write (iwrite,"(a,d21.14)",advance='no') ', ', Energies(i)
			enddo
			write (iwrite,'(a)') ' '  ! Finish the line

			write (*,"(i4)",advance='no') j
			do i = LowerEigen, UpperEigen, 1
				write (*,"(a,d21.14)",advance='no') ', ', Energies(i)
			enddo
			write (*,*) ' '  ! Finish the line
		endif
	enddo

	write (iwrite,'(a)') '</energies>'
	
	deallocate(PhiPhi8, Phi2HPhi8)
	deallocate(Energies, Workspace)
	
	return
end

	
! This subroutine calculates the number of terms for the function GenPowerTable (below).  This
!  is needed to determine the size of the dynamically allocated table before it is produced.
integer function CalcPowerTableSize(Omega)
	implicit none
	integer Omega  ! This sets the limits on the terms.
	integer NumTerms  ! The total number of terms
	integer om, ki, li, mi, ni, pi, qi  ! These are the exponents we are determining.

	NumTerms = 0
	do om = 0, Omega, 1
		do ki = 0, Omega, 1
			do li = 0, Omega, 1
				do mi = 0, Omega, 1
					do ni = 0, Omega, 1
						do pi = 0, Omega, 1
							do qi = 0, Omega, 1
								if (ki + li + mi + ni + pi + qi == om) then
									NumTerms = NumTerms + 1
								endif
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo
	enddo
	
	CalcPowerTableSize = NumTerms * 2  ! Returns elements of both symmetries.
	return
end


! This subroutine calculates values of k_i, l_i, m_i, n_i, p_i and q_i of
!  (3.14).  The summation of these is given by equation (3.15), but we do not have
!  the restriction on q of being even.
! This is written to always have the terms in order of increasing omega, i.e. the terms
!  for omega = 0 are first, etc.
!@TODO: Use logical for Ordering
subroutine GenOmegaPowerTable(Omega, PowerTable, NumTerms, Ordering)
	!implicit none
	integer Omega  ! This sets the limits on the terms.
	integer ArraySize
	integer, dimension(NumTerms,6) :: PowerTable
	integer NumTerm, i  ! The number of the current term
	integer om, ki, li, mi, ni, pi, qi  ! These are the exponents we are determining.
	integer Ordering  ! Whether to use Peter's ordering or mine

	if (Ordering == 0) then  ! Use my ordering
		NumTerm = 0
		do om = 0, Omega, 1
			do ki = 0, Omega, 1
				do li = 0, Omega, 1
					do mi = 0, Omega, 1
						do ni = 0, Omega, 1
							do pi = 0, Omega, 1
								do qi = 0, Omega, 1
									if (ki + li + mi + ni + pi + qi == om) then
										NumTerm = NumTerm + 1
										PowerTable(NumTerm,1) = ki
										PowerTable(NumTerm,2) = li
										PowerTable(NumTerm,3) = mi
										PowerTable(NumTerm,4) = ni
										PowerTable(NumTerm,5) = pi
										PowerTable(NumTerm,6) = qi
									endif
								enddo
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo

	else  ! Use Peter Van Reeth's ordering instead.
		! I pulled this (modified) snippet from Peter Van Reeth's code so that
		!  my indices match up with his to debug mine easier.
		MEGA=OMEGA
			IHDPP1 = MEGA +1
			INX=0
			DO 65 I=1,IHDPP1
			  DO 66 I23P1=1,I,1
				I23=I23P1 -1
				I12P1M = I-I23
				DO 67 I12P1=1,I12P1M
				  I12=I12P1-1
				  I2P1M =I12P1M -I12
				  DO 68 I2P1 =1,I2P1M
					 I2 = I2P1-1
					 I13P1M = I2P1M-I2
					 DO 69 I13P1=1,I13P1M
						I13 = I13P1 -1
						I3P1M = I13P1M -I13
						DO 70 I3P1 = 1, I3P1M
						   I3=I3P1-1
						   I1 = I3P1M -I3P1
						   
						   INX =INX +1
						   PowerTable(INX,1) = I1
						   PowerTable(INX,2) = I2
						   PowerTable(INX,3) = I12
						   PowerTable(INX,4) = I3
						   PowerTable(INX,5) = I13
						   PowerTable(INX,6) = I23
						   NumTerm = INX

	70                  CONTINUE
	69                CONTINUE
	68              CONTINUE
	67            CONTINUE
	66          CONTINUE
	65        CONTINUE
	endif

	do i = 1, NumTerm, 1
		! Duplicate these terms to have both r1 and r2 symmetries.
		PowerTable(NumTerm+i,1) = PowerTable(i,1)
		PowerTable(NumTerm+i,2) = PowerTable(i,2) + 1
		PowerTable(NumTerm+i,3) = PowerTable(i,3)
		PowerTable(NumTerm+i,4) = PowerTable(i,4)
		PowerTable(NumTerm+i,5) = PowerTable(i,5)
		PowerTable(NumTerm+i,6) = PowerTable(i,6)

		! Increase r1 by 1 for the first set only
		PowerTable(i,1) = PowerTable(i,1) + 1
	end do

	return
end


! This just cleans up the following subroutines.
subroutine AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj

	ki = PowerTablei(i,1)
	li = PowerTablei(i,2)
	mi = PowerTablei(i,3)
	ni = PowerTablei(i,4)
	pi = PowerTablei(i,5)
	qi = PowerTablei(i,6)
	kj = PowerTablej(j,1)
	lj = PowerTablej(j,2)
	mj = PowerTablej(j,3)
	nj = PowerTablej(j,4)
	pj = PowerTablej(j,5)
	qj = PowerTablej(j,6)

	return
end subroutine


! Sets all WMatrix elements equal to 0. 
subroutine ClearWCMatrices(WMatrix, CMatrix)
	use WCLimits
	implicit none
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*16, dimension(-2:Cjmax, 0:Cqmax, 0:Ckmax) :: CMatrix
	
	WMatrix = 0.0q0  ! Just a chosen value to show whether this has been used.
	CMatrix = 0.0q0
	return
end


! This precalculates the W matrices, since the W function is computationally expensive.
subroutine CalcWCMatrices(Omega, qmax, WMatrix, CMatrix, Alpha, Beta, Gamma, pmax)
	use WCLimits
	implicit none
	integer Omega, qmax
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*16, dimension(-2:Cjmax, 0:Cqmax, 0:Ckmax) :: CMatrix
	real*16 Alpha, Beta, Gamma
	integer pmax
	real*16 W, C
	integer j, k, l, m, n, q
	integer Used, Total
	
	Used = 0; Total = 0;

	!$omp parallel do shared(lmin,lmax,mmin,mmax,nmin,nmax,pmax,Alpha,Beta,Gamma,WMatrix) private(m,n) schedule(dynamic,5)
	do l = lmin, lmax, 1
		write (*,*) "wmatrix l:", l, "/", lmax
		!write (*,"(a)",advance='no') '.'
		do m = mmin, mmax, 1
			do n = nmin, nmax, 1
				if (l + m + n + 2 >= 0 .and. l + m + 1 >= 0) then  ! TODO: Needed?
					do k = 1, 6, 1
						! Note that the ordering of alpha, beta and gamma is the same as the ordering
						!  in the terms of equation (6) of Drake/Yan '95.
						Total = Total + 1
						if (WMatrix(l, m, n, k) /= 0.0d0) then
							Used = Used + 1
							select case(k)
								case (1)
									WMatrix(l, m, n, 1) = W(l, m, n, Alpha, Beta, Gamma, pmax)
								case (2)
									WMatrix(l, m, n, 2) = W(l, m, n, Alpha, Gamma, Beta, pmax)
								case (3)
									WMatrix(l, m, n, 3) = W(l, m, n, Beta, Alpha, Gamma, pmax)
								case (4)
									WMatrix(l, m, n, 4) = W(l, m, n, Beta, Gamma, Alpha, pmax)
								case (5)
									WMatrix(l, m, n, 5) = W(l, m, n, Gamma, Alpha, Beta, pmax)
								case (6)
									WMatrix(l, m, n, 6) = W(l, m, n, Gamma, Beta, Alpha, pmax)
							end select  ! There is no need for a default case.
						end if
					end do
				end if
			end do
		end do
	end do
	
	write (*,*) Used, "of", Total, ":", dble(Used)/Total*100
	
	do j = -2, Cjmax, 1
		do q = 0, Cqmax, 1
			do k = 0, Ckmax, 1
				if (CMatrix(j, q, k) /= 0.0q0) then
					CMatrix(j, q, k) = C(j, q, k)
				end if
			end do
		end do
	end do
	
	return
end


! Called by optimizewavefn in Minimize.c
subroutine calcmatriceswrapper(Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, Alpha, Beta, Gamma, Energy) bind(C)
	implicit none
	real*16, allocatable, dimension(:,:) :: PhiPhi, Phi2HPhi
	real*8, allocatable, dimension(:,:) :: PhiPhi8, Phi2HPhi8
	real*8, allocatable, dimension(:) :: Energies, Workspace
	integer Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, Info
	integer, save :: Iter = 1
	real*8 Alpha, Beta, Gamma, Energy

	allocate(PhiPhi(NumTerms,NumTerms))
	allocate(Phi2HPhi(NumTerms,NumTerms))
	allocate(PhiPhi8(NumTerms,NumTerms))
	allocate(Phi2HPhi8(NumTerms,NumTerms))
	allocate(Energies(NumTerms))
	allocate(Workspace(3*NumTerms-1))

	Energies = 0.0d0

	call CalcMatrices(Omega, NumTerms, real(Alpha,16), real(Beta,16), real(Gamma,16), 1, PhiPhi, Phi2HPhi, &
					qmax, pmax, IsTriplet, Ordering, Method)
	PhiPhi8 = PhiPhi
	Phi2HPhi8 = Phi2HPhi / 2
	call dsygv(1, 'N', 'L', NumTerms, Phi2HPhi8, NumTerms, PhiPhi8, NumTerms, Energies, Workspace, 3*NumTerms-1, Info)
	Energy = Energies(1)

	write (*,"(I4,A,3F8.5)") Iter, "; Alpha, Beta, Gamma: ", Alpha, Beta, Gamma
	write (10,"(I4,A,3F8.5)") Iter, "; Alpha, Beta, Gamma: ", Alpha, Beta, Gamma  ! I dislike hardcoding the 10, but it's easier at this point.
	write (10,"(A,F16.12)") "Energy: ", Energy
	write (10,*)
	flush(10)
	
	deallocate(PhiPhi, Phi2HPhi, PhiPhi8, Phi2HPhi8, Energies, Workspace)
	Iter = Iter + 1
end


! Calculates the derivative of the energy function.  If the Deriv parameter = 0, then it exits early with
!  PhiPhi and Phi2HPhi filled.
subroutine CalcMatrices(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, qmax, pmax, IsTriplet, Ordering, Method)
	use WCLimits
	!use nag_nsym_gen_eig
	implicit none
	real*16 Alpha, Beta, Gamma, Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer, allocatable, dimension(:,:) :: PowerTablei, PowerTablej  ! This contains the powers of the Hylleraas terms.
	integer CalcPowerTableSize  ! Function definition
	integer EigenNum  ! This is the energy eigenvalue number we are dealing with.
	integer Omega, NumTerms, qmax, pmax, IsTriplet, Method, Ordering, i, j, n, TempJ, Info
	real*8, allocatable, dimension(:) :: CoeffExpan
	real*8, allocatable, dimension(:) :: Energies, Workspace
	real*16, allocatable, dimension(:,:,:,:) :: WMatrix
	real*16, allocatable, dimension(:,:,:) :: CMatrix
	real*16, dimension(NumTerms,NumTerms) :: PhiPhi, Phi2HPhi
	real*8 PhiHPhiSum, PhiPhiSum, Norm
	
	allocate(PowerTablei(NumTerms, 6))
	call GenOmegaPowerTable(Omega, PowerTablei, NumTerms, Ordering)
	allocate(PowerTablej(NumTerms, 6))
	call GenOmegaPowerTable(Omega, PowerTablej, NumTerms, Ordering)
	allocate(Energies(NumTerms))
	allocate(Workspace(3*NumTerms-1))
	allocate(CoeffExpan(NumTerms))

	! Set this matrix equal to 0 because of the Phi2HPhi(i,j) = Phi2HPhi(i,j) + Sum line in CalcMatrices
	!  (and equivalent for PhiPhi).  This allows us to use the same subroutine for the direct and exchange elements.
	Phi2HPhi = 0.0q0
	PhiPhi = 0.0q0
	
	! For the direct integration, f_i and f_j have beta and gamma in the same place.
	Alphai = Alpha; Alphaj = Alpha; Betai = Beta; Betaj = Beta; Gammai = Gamma; Gammaj = Gamma;

	! Precalculate the W matrices.  lmin, lmax, etc. are the upper and lower limits for each of the three
	!  parameters to the W function (in HylleraasIntegral.f90).  These limits were determined by me, and
	!  there is a document detailing these limits (W Function Limits.pdf).
	lmin = 0;  lmax = 2*(Omega+1) + 2*qmax + 5;
	mmin = -1;  mmax = 2*(Omega+1) + 4;
	nmin = -2 - 2*qmax;  nmax = 2*(Omega+1) + 3;
	Cjmax = 2*Omega+4;  Cqmax = qmax;  Ckmax = 2*Omega+5;
	allocate(WMatrix(lmin:lmax, mmin:mmax, nmin:nmax, 6))
	allocate(CMatrix(-2:Cjmax, 0:Cqmax, 0:Ckmax))

	if (Method == 0 .or. Method == 1) then
		write (*,*) "Precomputing W matrix"
		call ClearWCMatrices(WMatrix, CMatrix)
		! First call it with RunCalc = false to determine what W functions to calculate.
		call CalcDirectMatrices(.false., PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
							Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, Omega, pmax, qmax, 0, Method)
		call CalcWCMatrices(Omega+1, qmax, WMatrix, CMatrix, Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, pmax)
		write (*,*) "Finished precomputing W matrix"
	else  ! Method == 2 (recursion relations)
		write (*,*) "Precomputing gamma functions..."
		call PreCalcGamma(Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, 25)
		write (*,*) "Finished precomputing gamma functions"
	end if

	! This subroutine call calculates the PhiPhi and Phi2HPhi matrices for the direct integrals.
	call CalcDirectMatrices(.true., PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, Omega, pmax, qmax, 0, Method)

	! Now we have to do the same thing for the exchange operator, P_23.
	call Permute(PowerTablej, NumTerms, Alphai, Alphaj, Alpha, Betai, Betaj, Beta, Gammai, Gammaj, Gamma)
	write (*,*)

	if (Method == 0 .or. Method == 1) then
		write (*,*) "Precomputing W matrix"
		call ClearWCMatrices(WMatrix, CMatrix)
		! First call it with RunCalc = false to determine what W functions to calculate.
		call CalcExchangeMatrices(.false., PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
							Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, Omega, pmax, qmax, IsTriplet, Method)
		call CalcWCMatrices(Omega+1, qmax, WMatrix, CMatrix, Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, pmax)
		write (*,*) "Finished precomputing W matrix"
	else  ! Method == 2 (recursion relations)
		write (*,*) "Precomputing gamma functions..."
		call PreCalcGamma(Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, 25)
		write (*,*) "Finished precomputing gamma functions"
	end if
	
	! Calculate the matrix elements for the exchanged terms, and add them to the direct terms.
	call CalcExchangeMatrices(.true., PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, Omega, pmax, qmax, IsTriplet, Method)

	! Clean up memory before exiting
	deallocate(WMatrix, CMatrix)
	deallocate(PowerTablei, PowerTablej)
	
	return
end

						
! We operate it on f_j by switching the
!  2<->3 powers in PowerTablej and recalculating, along with beta and gamma.
subroutine Permute(PowerTablej, NumTerms, Alphai, Alphaj, Alpha, Betai, Betaj, Beta, Gammai, Gammaj, Gamma)
	implicit none
	integer NumTerms, n, TempJ
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alpha, Alphai, Alphaj, Betai, Betaj, Beta, Gammai, Gammaj, Gamma
	
	do n = 1, NumTerms, 1
		! This exchanges r2 and r3.
		TempJ = PowerTablej(n,2)
		PowerTablej(n,2) = PowerTablej(n,4)  ! l with n
		PowerTablej(n,4) = TempJ
		! This exchanges r12 and r13 (r31).
		TempJ = PowerTablej(n,3)
		PowerTablej(n,3) = PowerTablej(n,5)  ! m with p
		PowerTablej(n,5) = TempJ
	enddo
	
	Alphai = Alpha; Alphaj = Alpha; Betaj = Gamma; Betai = Beta; Gammaj = Beta; Gammai = Gamma;
	return
end
