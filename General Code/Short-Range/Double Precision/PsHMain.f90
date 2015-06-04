! This program calculates the matrix elements of the positronium-hydrogen
!  matrices (similar to equation 2.15 of the Armour and Humberston article).  Specifically,
!  we are calculating elements of the form (phi_i, L phi_j) as in equation (3.22).
! This does use OpenMP to speed up computation on multicore processors, so some familiarity
!  with parallel programming is required to understand the !$omp lines.
! A future speed-up could be done with Open MPI (as is done for the energy code), but this
!  computation is usually done in a few hours.

! TODO: Make MaxIter and the 4 come from the configuration file.

module WLimits
   implicit none
   integer lmin, lmax
   integer mmin, mmax
   integer nmin, nmax
end module WLimits


program PsHMain
	use WLimits
	implicit none

	interface
		integer(c_int) function optimizewavefn(MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, alpha, beta, gamma) bind(c, name='optimizewavefn')
		use iso_c_binding
		integer(c_int), value :: MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine
		real(c_double), value :: alpha, beta, gamma
		end function optimizewavefn
	end interface
	
	real*8 omp_get_wtime
	integer Omega  ! This sets the limits on the terms.
	integer NumTerms  ! Number of terms in our Hylleraas expansion
	integer iread, iwrite  ! Our input and output files
	integer CalcPowerTableSize
	real*8, allocatable, dimension(:) :: Energies, Workspace
	real*16, allocatable, dimension(:,:) :: PhiPhi, Phi2HPhi
	real*8, allocatable, dimension(:,:) :: PhiPhi8, Phi2HPhi8
	real*8 Alpha, Beta, Gamma, Prev
	integer, allocatable, dimension(:) :: UsedTerms
	integer NumUsed, i, j, k, n
	real*8 StartTime, EndTime
	integer iargc, IsTriplet, Method, LowerEigen, UpperEigen, Info, AllEnergies, EigenRoutine
	character *100 IOBuffer
	integer Optimize, EigenNum, Ordering, Iter, MaxIter, MaxIterOuter, qmax, pmax, linf, minf, StartEnergy, LValue
	real*8 Tol, Err, h, EDeriv, fa, fh, y, Divisor, ederivwrapper
	real*8, dimension(3) :: Del, Delta, f1, f0, x1, x0, x, C31, C13
	real*8, dimension(3,3) :: B, BTemp, C33, Jacob
	real*8 ThreeJSymbol


	!! TODO: Remove this!
	!real*8, allocatable, dimension(:,:,:,:) :: WMatrix
	!real*8 HylleraasIntegralGeneral
	!lmin=0; lmax=0; mmin=0; mmax=0; nmin=0; nmax=0;
	!allocate(WMatrix(lmin:lmax, mmin:mmax, nmin:nmax, 6))
	!Tol = HylleraasIntegralGeneral(.true., .false., 1, 1, 1, -2, 1, 1, 2.7d0, 2.7d0, 2.7d0, 0, 0, 0, 0, 0, 0, &
	!			0, 0, 0, 0, 0, 0, 15, 50, 15, 15, 0, 0, 0, WMatrix, 1)


	iread = 9
	iwrite = 10
	
	! This allows the possibility of using different input and output files than the defaults.
	!  This mainly lets me run this program on multiple computers at the same time without
	!  one overwriting the results from another.  If the program is called without any
	!  command-line arguments, it uses the default filenames.
	if (iargc() == 2) then
		call getarg(1, IOBuffer)
		open(iread, FILE=IOBuffer)
		call getarg(2, IOBuffer)
		open(iwrite, FILE=IOBuffer)
	else
		open(iread, FILE='input.txt')
		open(iwrite, FILE='output.txt')
	endif

	! Read the input parameters.
	call ReadParamFile(iread, Omega, LValue, Alpha, Beta, Gamma, qmax, pmax, linf, minf, Method, IsTriplet, LowerEigen, UpperEigen, &
						 AllEnergies, Optimize, EigenNum, MaxIter, Ordering, EigenRoutine)
	
	if (IsTriplet /= 0 .and. IsTriplet /= 1) then
		write (*,*) 'Must choose between singlet and triplet calculation - exiting.'
		stop
	endif

	if (Optimize == 1) write (iwrite,*) 'Optimizing eigenvalue', EigenNum

	call WriteHeader(iwrite, LValue, IsTriplet, Ordering, Method, EigenRoutine, Omega, Alpha, Beta, Gamma)
	
	Tol = 1e-3
	Err = 1  ! Just has to be larger than tol.
	h = 1e-5
	
	! Get the start time to determine the duration of the program.
	!call cpu_time(StartTime)
	StartTime = omp_get_wtime()

	NumTerms = CalcPowerTableSize(Omega) * 2  ! First and second symmetries
	write (*,*) 'NumTerms:', NumTerms
	allocate(PhiPhi(NumTerms,NumTerms))
	allocate(Phi2HPhi(NumTerms,NumTerms))
	allocate(PhiPhi8(NumTerms,NumTerms))
	allocate(Phi2HPhi8(NumTerms,NumTerms))
	allocate(Energies(NumTerms))
	allocate(Workspace(3*NumTerms-1))

	write (*,*) 'Number of terms:', NumTerms
	write (*,*) 'Omega =', Omega
	write (iwrite,*) 'Number of terms:', NumTerms	
	write (iwrite,*)
	
	if (Omega == -1) then
		goto 200  ! No short-range terms, so no need to do these calculations (just a header)
	end if

	if (Optimize == 1) then
		! TODO: Add these back in.
		!fa = optimizewavefn(MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, Alpha, Beta, Gamma)
		!fa = EDeriv(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, PhiPhi8, Phi2HPhi8, qmax, pmax, 0, &
		!			IsTriplet, Ordering, Method)
	else
		! If we are not optimizing eigenvalues, we can still call the same function, which exits early with PhiPhi and Phi2HPhi filled.
		fa = EDeriv(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, PhiPhi8, Phi2HPhi8, qmax, pmax, linf, minf, 0, &
					IsTriplet, Ordering, Method, EigenRoutine, LValue)
		write (*,*) fa
	endif

	do i = 1, NumTerms, 1
		!do j = 1, i, 1  ! Use this instead if only the lower triangle is required.
		do j = 1, NumTerms, 1
			write (iwrite,"(i6,i6,d28.20,d28.20)") i, j, PhiPhi(i,j), Phi2HPhi(i,j)
		enddo
	enddo
	
	! Divide every entry in Phi2HPhi by 2, since we calculated <phi_i|2H|phi_j> in equation (3.22).
	Phi2HPhi = Phi2HPhi / 2.0_16

	write (iwrite,*)
	write (iwrite,*) 'Energies:'
	
	if (AllEnergies == 0) then
		StartEnergy = NumTerms
	else
		StartEnergy = 1
	end if

	do j = StartEnergy, NumTerms, 1
	!do j = 1, 1, 1
		! Copies our real*16 matrices to real*8 matrices so that LAPACK can operate on them.
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
		!call Newdsygv(1, 'N', 'L', j, Phi2HPhi, NumTerms, PhiPhi, NumTerms, Energies, Workspace, 3*NumTerms-1, Info)
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
			write (iwrite,"(i4)",advance='no') j
			do i = LowerEigen, UpperEigen, 1
				write (iwrite,"(a,d21.14)",advance='no') ', ', Energies(i)
			enddo
			write (iwrite,*) ' '  ! Finish the line

			write (*,"(i4)",advance='no') j
			do i = LowerEigen, UpperEigen, 1
				write (*,"(a,d21.14)",advance='no') ', ', Energies(i)
			enddo
			write (*,*) ' '  ! Finish the line
		endif
	enddo


	! Get the end time to find the duration of the program.
	!call cpu_time(EndTime)
200	EndTime = omp_get_wtime()
	write (*,*) 'Time taken (s):', EndTime - StartTime
	write (iwrite,*)
	write (iwrite,*) 'Time taken (s):', EndTime - StartTime, '  (min):', (EndTime - StartTime) / 60.0

	! Clean up memory before exiting
	deallocate(Workspace)
	deallocate(Energies)
	deallocate(PhiPhi8)
	deallocate(Phi2HPhi8)
	deallocate(PhiPhi)
	deallocate(Phi2HPhi)
	
	! Close file handles
	close(iwrite)

	stop
end


subroutine ReadParamFile(iread, Omega, LValue, Alpha, Beta, Gamma, qmax, pmax, linf, minf, Method, IsTriplet, LowerEigen, UpperEigen, &
						 AllEnergies, Optimize, EigenNum, MaxIter, Ordering, EigenRoutine)
	implicit none
	integer iread, LValue, Omega, Method, qmax, pmax, linf, minf, IsTriplet, LowerEigen, UpperEigen, AllEnergies
	integer Optimize, EigenNum, MaxIter, Ordering, EigenRoutine
	real*8 Alpha, Beta, Gamma

	read (iread,*) ! Description line
	read (iread,*) Omega
	read (iread,*) ! Description line
	read (iread,*) LValue
	read (iread,*) 
	read (iread,*) Alpha
	read (iread,*) Beta
	read (iread,*) Gamma
	read (iread,*) 
	read (iread,*) Method
	read (iread,*) 
	read (iread,*) qmax
	read (iread,*) 
	read (iread,*) pmax
	read (iread,*) 
	read (iread,*) linf
	read (iread,*) minf
	read (iread,*) 
	read (iread,*) IsTriplet
	read (iread,*) 
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
	read (iread,*)
	read (iread,*) EigenRoutine
	close (iread)
	
	return
end

subroutine WriteHeader(iwrite, LValue, IsTriplet, Ordering, Method, EigenRoutine, Omega, Alpha, Beta, Gamma)
	implicit none
	integer iwrite, LValue, IsTriplet, Ordering, Method, EigenRoutine, Omega
	real*8 Alpha, Beta, Gamma
	
	select case (LValue)
		case (0)  ! S-Wave
			if (IsTriplet == 0) then
				write (*,*) "S-Wave Singlet Ps-H"
				write (iwrite,*) "S-Wave Singlet Ps-H"
			else
				write (*,*) "S-Wave Triplet Ps-H"
				write (iwrite,*) "S-Wave Triplet Ps-H"
			endif
		case (1)
			if (IsTriplet == 0) then
				write (*,*) "P-Wave Singlet Ps-H: 1st formalism"
				write (iwrite,*) "P-Wave Singlet Ps-H: 1st formalism"
			else
				write (*,*) "P-Wave Triplet Ps-H: 1st formalism"
				write (iwrite,*) "P-Wave Triplet Ps-H: 1st formalism"
			endif
		case (2)
			if (IsTriplet == 0) then
				write (*,*) "D-Wave Singlet Ps-H: 1st formalism"
				write (iwrite,*) "D-Wave Singlet Ps-H: 1st formalism"
			else
				write (*,*) "D-Wave Triplet Ps-H: 1st formalism"
				write (iwrite,*) "D-Wave Triplet Ps-H: 1st formalism"
			endif
		case (3)  ! F-Wave
			if (IsTriplet == 0) then
				write (*,*) "F-Wave Singlet Ps-H"
				write (iwrite,*) "F-Wave Singlet Ps-H"
			else
				write (*,*) "F-Wave Triplet Ps-H"
				write (iwrite,*) "F-Wave Triplet Ps-H"
			endif
		case (4)  ! G-Wave
			if (IsTriplet == 0) then
				write (*,*) "G-Wave Singlet Ps-H"
				write (iwrite,*) "G-Wave Singlet Ps-H"
			else
				write (*,*) "G-Wave Triplet Ps-H"
				write (iwrite,*) "G-Wave Triplet Ps-H"
			endif
		case (5)  ! H-Wave
			if (IsTriplet == 0) then
				write (*,*) "H-Wave Singlet Ps-H"
				write (iwrite,*) "H-Wave Singlet Ps-H"
			else
				write (*,*) "H-Wave Triplet Ps-H"
				write (iwrite,*) "H-Wave Triplet Ps-H"
			endif
		case (6)  ! I-Wave
			if (IsTriplet == 0) then
				write (*,*) "I-Wave Singlet Ps-H"
				write (iwrite,*) "I-Wave Singlet Ps-H"
			else
				write (*,*) "I-Wave Triplet Ps-H"
				write (iwrite,*) "I-Wave Triplet Ps-H"
			endif
		case (7)  ! K-Wave
			if (IsTriplet == 0) then
				write (*,*) "K-Wave Singlet Ps-H"
				write (iwrite,*) "K-Wave Singlet Ps-H"
			else
				write (*,*) "K-Wave Triplet Ps-H"
				write (iwrite,*) "K-Wave Triplet Ps-H"
			endif
		case (8)  ! L-Wave
			if (IsTriplet == 0) then
				write (*,*) "L-Wave Singlet Ps-H"
				write (iwrite,*) "L-Wave Singlet Ps-H"
			else
				write (*,*) "L-Wave Triplet Ps-H"
				write (iwrite,*) "L-Wave Triplet Ps-H"
			endif
		case default
			write (*,*) "Higher partial waves are not supported yet...exiting."
			stop
	end select

	if (Ordering == 1) then
		write (*,*) "Using Peter Van Reeth's ordering"
		write (iwrite,*) "Using Peter Van Reeth's ordering"
	else
		write (*,*) "Using my ordering"
		write (iwrite,*) "Using my ordering"
	endif

	if (Method == 0) then
		write (*,*) "Integration Technique: Direct Summation"
		write (iwrite,*) "Integration Technique: Direct Summation"
	else if (Method == 1) then
		write (*,*) "Integration Technique: Asymptotic Expansion"
		write (iwrite,*) "Integration Technique: Asymptotic Expansion"
	else if (Method == 2) then
		write (*,*) "Integration Technique: Recursion Relations"
		write (iwrite,*) "Integration Technique: Recursion Relations"
	else
		write (*,*) "Method parameter in input file must be 0, 1 or 2."
		stop
	end if
	
	if (EigenRoutine == 1) then
		write (*,*) "Eigenvalue routine: dsygv from LAPACK"
		write (iwrite,*) "Eigenvalue routine: dsygv from LAPACK"
	else if (EigenRoutine == 2) then
		write (*,*) "Eigenvalue routine: Pachucki eigenproblem solver"
		write (iwrite,*) "Eigenvalue routine: Pachucki eigenproblem solver"
	else
		write (*,*) "Eigenvalue routine parameter in input file must be 1 or 2."
		stop
	end if

	write (iwrite,*) 'Omega =', Omega
	write (iwrite,*) 'Alpha =', Alpha
	write (iwrite,*) 'Beta =', Beta
	write (iwrite,*) 'Gamma =', Gamma

	
	return
end

	
! This subroutine calculates the number of terms for the function GenOmegaPowerTable (below).  This
!  is needed to determine the size of the dynamically allocated table before it is produced.
integer function CalcPowerTableSize(Omega)
	implicit none
	integer Omega  ! This sets the limits on the terms.
	integer NumTerms  ! The total number of terms
	integer om, ki, li, mi, ni, pi, qi  ! These are the exponents we are determining.

	if (Omega == -1) then  ! Special case where we don't want any short-range terms
		CalcPowerTableSize = 0
		return
	end if
	
	NumTerms = 0
	do om = 0, Omega, 1
		do ki = 0, Omega, 1
			do li = 0, Omega, 1
				do mi = 0, Omega, 1
					do ni = 0, Omega, 1
						do pi = 0, Omega, 1
							do qi = 0, Omega, 1
								if (ki + li + mi + ni + pi + qi == om) then
									!if (ki>li) cycle
									!if (ki==li .and. qi>pi) cycle

!									if (om > 3 .and. qi > 0) cycle  ! Restrict r23 to max of 1
!									if (om > 3 .and. pi > 0) cycle  ! Restrict r13 to max of 3
!									if (om > 2 .and. ni > 0) cycle  ! Restrict r3 to max of 0
!									if (ki >= li .and. qi > pi) cycle  ! Check if r1>r2 and r23 > r13 (powers)
									NumTerms = NumTerms + 1
								endif
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo
	enddo
	
	CalcPowerTableSize = NumTerms
	return
end


! This subroutine calculates values of k_i, l_i, m_i, n_i, p_i and q_i of
!  (3.14).  The summation of these is given by equation (3.15), but we do not have
!  the restriction on q of being even.
! This is written to always have the terms in order of increasing omega, i.e. the terms
!  for omega = 0 are first, etc.
! TODO: Use logical for Ordering
subroutine GenOmegaPowerTable(Omega, PowerTable, ArraySize, Ordering, l1, l2)
	!implicit none
	integer Omega  ! This sets the limits on the terms.
	integer ArraySize
	integer, dimension(ArraySize,6) :: PowerTable
	integer NumTerm  ! The number of the current term
	integer om, ki, li, mi, ni, pi, qi  ! These are the exponents we are determining.
	integer Ordering  ! Whether to use Peter's ordering or mine
	integer l1, l2, l3  ! l1 = 0 for S-wave, 1 for P-wave, etc.

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
										!if (ki + li + mi + ni + pi + qi <= Omega) then
										!if (ki>li) cycle
										!if (ki==li .and. qi>pi) cycle


	!									if (om > 3 .and. qi > 0) cycle  ! Restrict r23 to max of 1
	!									if (om > 3 .and. pi > 0) cycle  ! Restrict r13 to max of 3
	!									if (om > 2 .and. ni > 0) cycle  ! Restrict r3 to max of 0
	!									if (ki >= li .and. qi > pi) cycle  ! Check if r1>r2 and r23 > r13 (powers)
									
										NumTerm = NumTerm + 1
										PowerTable(NumTerm,1) = ki
										PowerTable(NumTerm,2) = li
										PowerTable(NumTerm,3) = mi
										PowerTable(NumTerm,4) = ni
										PowerTable(NumTerm,5) = pi
										PowerTable(NumTerm,6) = qi
										
										write (*,*) ki, li, mi, ni, pi, qi!, ki + li + mi + ni + pi + qi
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

						   !if (I1>I2) cycle
						   !if (I1==I2 .and. I23>I13) cycle
						   
						   INX =INX +1
						   PowerTable(INX,1) = I1
						   PowerTable(INX,2) = I2
						   PowerTable(INX,3) = I12
						   PowerTable(INX,4) = I3
						   PowerTable(INX,5) = I13
						   PowerTable(INX,6) = I23

						   !write (*,*) I1, I2, I12, I3, I13, I23

	70                  CONTINUE
	69                CONTINUE
	68              CONTINUE
	67            CONTINUE
	66          CONTINUE
	65        CONTINUE
		NumTerm = INX
	endif

! REMOVE THIS !
	if (Ordering == 2) then  ! Use my ordering
		NumTerm = 0
		do om = 0, Omega, 1
			do qi = 0, Omega, 1
				do pi = 0, Omega, 1
					do ni = 0, Omega, 1
						do mi = 0, Omega, 1
							do li = 0, Omega, 1
								do ki = 0, Omega, 1
									if (ki + li + mi + ni + pi + qi == om) then
										!if (ki + li + mi + ni + pi + qi <= Omega) then
										!if (ki>li) cycle
										!if (ki==li .and. qi>pi) cycle


	!									if (om > 3 .and. qi > 0) cycle  ! Restrict r23 to max of 1
	!									if (om > 3 .and. pi > 0) cycle  ! Restrict r13 to max of 3
	!									if (om > 2 .and. ni > 0) cycle  ! Restrict r3 to max of 0
	!									if (ki >= li .and. qi > pi) cycle  ! Check if r1>r2 and r23 > r13 (powers)
									
										NumTerm = NumTerm + 1
										PowerTable(NumTerm,1) = ki
										PowerTable(NumTerm,2) = li
										PowerTable(NumTerm,3) = ni
										PowerTable(NumTerm,4) = mi
										PowerTable(NumTerm,5) = qi
										PowerTable(NumTerm,6) = pi
										
										write (*,*) ki, li, ni, mi, qi, pi!, ki + li + mi + ni + pi + qi
									endif
								enddo
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo
	endif

	do i = 1, NumTerm, 1
		! Increase r1 by l1 and r2 by l2/l3
		PowerTable(i,1) = PowerTable(i,1) + l1
		PowerTable(i,2) = PowerTable(i,2) + l2
		write (*,'(i3,a,6i2)') i, ":", PowerTable(i,1), PowerTable(i,2), PowerTable(i,3), PowerTable(i,4), PowerTable(i,5), PowerTable(i,6)
	end do
	
	return
end


subroutine GenCoeffTable(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj, l1r, l2r, l3r)
	implicit none
	real*8, dimension(34) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*8 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer kj, lj, mj, nj, pj, qj, l1r, l2r, l3r

	! We do not really have to do this, but it makes reading the assignments later
	!  easier to follow instead of assigning index numbers to each of k, l, m, etc.
	! In the power tables, we have the following indices:
	!  k = 1, l = 2, m = 3, n = 4, p = 5, q = 6
	kj = PowerTablej(j,1)
	lj = PowerTablej(j,2)
	mj = PowerTablej(j,3)
	nj = PowerTablej(j,4)
	pj = PowerTablej(j,5)
	qj = PowerTablej(j,6)

	CoeffTable(1) = -Alphaj*Alphaj - Betaj*Betaj - Gammaj*Gammaj
	CoeffTable(2) = -kj - kj*kj - kj*mj - kj*pj + l1r + l1r*l1r
	CoeffTable(3) = -2.0d0*mj - kj*mj - 2.0d0*mj*mj - mj*lj - mj*qj - mj*pj
	CoeffTable(4) = -2.0d0
	CoeffTable(5) = -lj - mj*lj - lj*lj - lj*qj + l2r + l2r*l2r
	CoeffTable(6) = mj*lj
	CoeffTable(7) = kj*mj
	CoeffTable(8) = -2.0d0*qj - mj*qj - lj*qj - 2.0d0*qj*qj - qj*nj - qj*pj
	CoeffTable(9) = 2.0d0
	CoeffTable(10) = -nj - qj*nj - nj*nj - nj*pj + l3r + l3r*l3r
	CoeffTable(11) = qj*nj
	CoeffTable(12) = lj*qj
	CoeffTable(13) = -2.0d0*pj - kj*pj - mj*pj - qj*pj - nj*pj - 2.0d0*pj*pj
	CoeffTable(14) = qj*pj
	CoeffTable(15) = mj*pj
	CoeffTable(16) = nj*pj
	CoeffTable(17) = kj*pj
	CoeffTable(18) = -2.0d0
	CoeffTable(19) = mj*qj
	CoeffTable(20) = mj*Alphaj
	CoeffTable(21) = -mj*Alphaj
	CoeffTable(22) = pj*Alphaj
	CoeffTable(23) = -pj*Alphaj
	CoeffTable(24) = 2.0d0 + 2.0d0*Alphaj + 2.0d0*kj*Alphaj + mj*Alphaj + pj*Alphaj
	CoeffTable(25) = -mj*Betaj
	CoeffTable(26) = mj*Betaj
	CoeffTable(27) = qj*Betaj
	CoeffTable(28) = -qj*Betaj
	CoeffTable(29) = -2.0d0 + 2.0d0*Betaj + mj*Betaj + 2.0d0*lj*Betaj + qj*Betaj
	CoeffTable(30) = -qj*Gammaj
	CoeffTable(31) = qj*Gammaj
	CoeffTable(32) = -pj*Gammaj
	CoeffTable(33) = pj*Gammaj
	CoeffTable(34) = -2.0d0 + 2.0d0*Gammaj + qj*Gammaj + 2.0d0*nj*Gammaj + pj*Gammaj

	return
end


! Note that the naming is a little different here.  I used rPowerTable here as the parameter name but
!  use rPowers elsewhere.
subroutine GenRPowerTable(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(34,6) :: rPowerTable
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(34,6) :: rPowers
	integer n

	! This is a 2-dimensional array expressed in a contiguous 1-D array.  I see no way to assign a 2-D
	!  array in Fortran (at least in G95).  I found several examples that would not compile with G95.
	!  These are taken directly from the table in rpowers.pdf.  Each power is covered by 2 lines, since
	!  there are 34 terms to take care of.  These can be read as going down the columns of rpowers.pdf.
	data rPowers / 0,-2, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1,-1, 1,-1,-1, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, &  ! r1
				   0, 0, 0, 0,-2,-2, 2, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-1, 1, 1,-1,-1, 2, 0, 0, 0, 0, &  ! r2
				   0, 0,-2,-1, 0,-2,-2, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0,-2,-2,-2, 0, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, &  ! r12
				   0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 2, 0, 0, 0,-2, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0,-1, 1,-1, 1,-1, &  ! r3
				   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2,-2,-1, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, &  ! r31
				   0, 0, 0, 0, 0, 0, 0,-2,-1, 0,-2,-2, 0,-2, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0,-2,-2, 0, 0, 0 /   ! r23
				  !1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34

	! We do not really have to do this, but it makes reading the assignments later
	!  easier to follow instead of assigning index numbers to each of k, l, m, etc.
	! In the power tables, we have the following indices:
	!  k = 1, l = 2, m = 3, n = 4, p = 5, q = 6
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

	! Combine powers of like r's from f_i, f_j and the powers in rpowers.pdf.
	do n = 1, 34, 1
		rPowerTable(n,1) = ki + kj + rPowers(n,1)
		rPowerTable(n,2) = li + lj + rPowers(n,2)
		rPowerTable(n,3) = mi + mj + rPowers(n,3)
		rPowerTable(n,4) = ni + nj + rPowers(n,4)
		rPowerTable(n,5) = pi + pj + rPowers(n,5)
		rPowerTable(n,6) = qi + qj + rPowers(n,6)
	enddo
	
	return
end


subroutine GenCoeffTableSH(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, l1r, l2r, l3r)
	implicit none
	real*8, dimension(6) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*8 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer kj, lj, mj, nj, pj, qj, l1r, l2r, l3r

	kj = PowerTablej(j,1)
	lj = PowerTablej(j,2)
	mj = PowerTablej(j,3)
	nj = PowerTablej(j,4)
	pj = PowerTablej(j,5)
	qj = PowerTablej(j,6)

	CoeffTable(1) = mj
	CoeffTable(2) = pj
	CoeffTable(3) = qj
	CoeffTable(4) = mj
	CoeffTable(5) = pj
	CoeffTable(6) = qj

	return
end


! Note that the naming is a little different here.  I used rPowerTable here as the parameter name but
!  use rPowers elsewhere.
subroutine GenRPowerTableSH(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(6,6) :: rPowerTable
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(6,6) :: rPowers
	integer n

	data rPowers / -1,-1, 0, 1, 1, 0, &  ! r1
				    1, 0,-1,-1, 0, 1, &  ! r2
				   -2, 0, 0,-2, 0, 0, &  ! r12
				    0, 1, 1, 0,-1,-1, &  ! r3
				    0,-2, 0, 0,-2, 0, &  ! r31
				    0, 0,-2, 0, 0,-2 /   ! r23
				   !1, 2, 3, 4, 5, 6

	! In the power tables, we have the following indices:
	!  k = 1, l = 2, m = 3, n = 4, p = 5, q = 6
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

	! Combine powers of like r's from f_i, f_j and the powers in rpowers.pdf.
	do n = 1, 6, 1
		rPowerTable(n,1) = ki + kj + rPowers(n,1)
		rPowerTable(n,2) = li + lj + rPowers(n,2)
		rPowerTable(n,3) = mi + mj + rPowers(n,3)
		rPowerTable(n,4) = ni + nj + rPowers(n,4)
		rPowerTable(n,5) = pi + pj + rPowers(n,5)
		rPowerTable(n,6) = qi + qj + rPowers(n,6)
	enddo
	
	return
end


! Sets all WMatrix elements equal to 0. 
! TODO: Is there an equivalent to the C++ memset?
subroutine ClearWMatrix(WMatrix)
	use WLimits
	implicit none
	real*8, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	integer l, m, n
	
	do l = lmin, lmax, 1
		do m = mmin, mmax, 1
			do n = nmin, nmax, 1
				WMatrix(l, m, n, 1) = 0.0d0
				WMatrix(l, m, n, 2) = 0.0d0
				WMatrix(l, m, n, 3) = 0.0d0
				WMatrix(l, m, n, 4) = 0.0d0
				WMatrix(l, m, n, 5) = 0.0d0
				WMatrix(l, m, n, 6) = 0.0d0
			end do
		end do
	end do
	
	return
end


! This precalculates the W matrices, since the W function is computationally expensive.
subroutine CalcWMatrices(Omega, qmax, WMatrix, Alpha, Beta, Gamma, pmax)
	use WLimits
	implicit none
	integer Omega, qmax
	real*8, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*8 Alpha, Beta, Gamma
	integer pmax
	real*8 W
	integer l, m, n, k
	
	!$omp parallel do shared(lmin,lmax,mmin,mmax,nmin,nmax,pmax,Alpha,Beta,Gamma,WMatrix) private(m,n) schedule(dynamic,5)
	do l = lmin, lmax, 1
		write (*,*) "wmatrix l:", l, "/", lmax
		!write (*,"(a)",advance='no') '.'
		do m = mmin, mmax, 1
			do n = nmin, nmax, 1
				!if (l + m + n + 2 >= 0 .and. l + m + 1 >= 0) then  ! TODO: Needed?
					do k = 1, 6, 1
						! Note that the ordering of alpha, beta and gamma is the same as the ordering
						!  in the terms of equation (6) of Drake/Yan '95.
						if (WMatrix(l, m, n, k) /= 0.0d0) then
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
				!end if
			end do
		end do
	end do
	write (*,*)

	return
end


subroutine CalcMatricesSub(RunCalc, UsePreCalc, CalcSH, PowerTablei, PowerTablej, WMatrix, PhiPhi, Phi2HPhi, NumTerms, OffI, OffJ, Alphai, Alphaj, Betai, Betaj, &
							Gammai, Gammaj, l1l, l2l, l3l, m1l, m2l, m3l, l1r, l2r, l3r, m1r, m2r, m3r, qmax, pmax, linf, minf, IsTriplet, Method)
	use WLimits
	implicit none
	real*8, dimension(34) :: CoeffTable
	real*8, dimension(6) :: CoeffTableSH
	integer, dimension(34,6) :: rPowers
	integer, dimension(6,6) :: rPowersSH
	integer, dimension(6) :: sldata, srdata
	integer, dimension(NumTerms,6) :: PowerTablei, PowerTablej
	real*8, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*16, dimension(NumTerms*2,NumTerms*2) :: PhiPhi, Phi2HPhi
	integer qmax, pmax, linf, minf, NumTerms, OffI, OffJ, Method
	logical RunCalc, UsePreCalc, CalcSH
	integer l1l, l2l, l3l, m1l, m2l, m3l, l1r, l2r, l3r, m1r, m2r, m3r
	real*8 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	real*8 HylleraasIntegral, HylleraasIntegralGeneral, CP12Angular
	real*16 Sum, SumPiHalf
	real*8 PI
	integer i, j, n, IsTriplet
	real*16 RemoveMe, SHPart, SHPartPiHalf
	real*16 PiHalf
	data sldata / 2, 3, 3, 1, 1, 2 /
	data srdata / 1, 1, 2, 2, 3, 3 /

	! Calculated in Mathematica
	PiHalf = 0.1591549430918953357689_16
	
	!call omp_set_num_threads(1)
	
	!$omp parallel do shared(NumTerms) private(j,Sum,SHPart,RemoveMe,CoeffTable,CoeffTableSH,rPowers,rPowersSH) schedule(dynamic,10)
	do i = 1, NumTerms, 1
		write (*,*) i
		!write (*,"(i7)",advance='no') i
		!do j = 1, i, 1
		do j = 1, NumTerms, 1
			Sum = 0.0d0
			call GenCoeffTable(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, &
								Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, l1r, l2r, l3r)
			call GenRPowerTable(rPowers, PowerTablei, i, PowerTablej, j, NumTerms)
			call GenCoeffTableSH(CoeffTableSH, PowerTablei, i, PowerTablej, j, NumTerms, &
								Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, l1r, l2r, l3r)
			call GenRPowerTableSH(rPowersSH, PowerTablei, i, PowerTablej, j, NumTerms)

			do n = 1, 34, 1
				if (CoeffTable(n) /= 0.0d0) then
					! Note that the order of the parameters to the HylleraasIntegral function looks strange. This is
					!  because the Drake and Yan paper (1995) expressed the integrand in the following order:
					!  r1 r2 r3 r12 r23 r31, and the Armour and Humberston article has them in this order:
					!  r1 r2 r12 r3 r13 r23 (equation 3.14). So we have to swap the third and fourth parameters, along
					!  with the fifth and sixth.
   
					! Do the PhiPhi inner product
					RemoveMe = HylleraasIntegralGeneral(RunCalc, UsePreCalc, rPowers(n,1), rPowers(n,2), rPowers(n,4), rPowers(n,3), rPowers(n,6), &
												  rPowers(n,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, &
												  l1l, l2l, l3l, m1l, m2l, m3l, l1r, l2r, l3r, m1r, m2r, m3r, qmax, pmax, linf, minf, 0, 0, 0, WMatrix, Method)
					Sum = Sum + RemoveMe * CoeffTable(n)
				end if
			enddo
   
			! TODO: Check for whether these terms should be 0 before calculating.
			! Now we need to take care of the terms that act on the spherical harmonics. There are 6 terms total from this.
			if (CalcSH == .true.) then
				SHPart = 0.0d0
				do n = 1, 6, 1
					if (CoeffTableSH(n) /= 0.0d0) then
						RemoveMe = HylleraasIntegralGeneral(RunCalc, UsePreCalc, rPowersSH(n,1), rPowersSH(n,2), rPowersSH(n,4), rPowersSH(n,3), rPowersSH(n,6), &
															  rPowersSH(n,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, &
															  l1l, l2l, l3l, m1l, m2l, m3l, l1r, l2r, l3r, m1r, m2r, m3r, qmax, pmax, linf, minf, n, sldata(n), srdata(n), WMatrix, Method)
						SHPart = SHPart + RemoveMe * CoeffTableSH(n)
					end if
				end do
				!SumPiHalf = Sum * PiHalf  ! TODO: Remove this.
				!SHPartPiHalf = SHPart * PiHalf  ! TODO: Remove this.
				Sum = Sum + SHPart * 2  ! The other terms have a 1/2 in front, and we are calculating 2H
			end if
   
			if (IsTriplet == 0) then
				Phi2HPhi(i+OffI,j+OffJ) = Phi2HPhi(i+OffI,j+OffJ) + Sum * PiHalf
			else
				Phi2HPhi(i+OffI,j+OffJ) = Phi2HPhi(i+OffI,j+OffJ) - Sum * PiHalf
			end if

			! Do the PhiPhi inner product
			RemoveMe = HylleraasIntegralGeneral(RunCalc, UsePreCalc, PowerTablei(i,1)+PowerTablej(j,1), PowerTablei(i,2)+PowerTablej(j,2), &
						PowerTablei(i,4)+PowerTablej(j,4), PowerTablei(i,3)+PowerTablej(j,3), &
						PowerTablei(i,6)+PowerTablej(j,6), PowerTablei(i,5)+PowerTablej(j,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, &
						l1l, l2l, l3l, m1l, m2l, m3l, l1r, l2r, l3r, m1r, m2r, m3r, qmax, pmax, linf, minf, 0, 0, 0, WMatrix, Method)
			if (IsTriplet == 0) then
				PhiPhi(i+OffI,j+OffJ) = PhiPhi(i+OffI,j+OffJ) + RemoveMe * PiHalf
			else
				PhiPhi(i+OffI,j+OffJ) = PhiPhi(i+OffI,j+OffJ) - RemoveMe * PiHalf
			endif
		enddo
	enddo
	
	write (*,*)
	return
end


						
subroutine CalcMatrices(RunCalc, UsePreCalc, PowerTabler1i, PowerTabler1j, PowerTabler2i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, LValue, qmax, pmax, linf, minf, IsTriplet, Method, Exchanged)
	use WLimits
	implicit none
	integer, dimension(NumTerms/2,6) :: PowerTabler1i, PowerTabler1j, PowerTabler2i, PowerTabler2j
	real*8, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*16, dimension(NumTerms,NumTerms) :: PhiPhi, Phi2HPhi
	integer NumTerms, LValue, qmax, pmax, linf, minf, IsTriplet, Method
	logical RunCalc, UsePreCalc
	real*8 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	logical Exchanged
	
	! 1-1 matrix elements
	write (*,*) "Calculating 1-1 matrix elements"
	call CalcMatricesSub(RunCalc, UsePreCalc, .false., PowerTabler1i, PowerTabler1j, WMatrix, PhiPhi, Phi2HPhi, NumTerms/2, 0, 0, Alphai, Alphaj, Betai, Betaj, &
						Gammai, Gammaj, LValue, 0, 0, 0, 0, 0, LValue, 0, 0, 0, 0, 0, qmax, pmax, linf, minf, IsTriplet, Method)
	
	! 2-2 matrix elements
	write (*,*) "Calculating 2-2 matrix elements"
	if (Exchanged == .true.) then  ! Swap l2r and l3r
		call CalcMatricesSub(RunCalc, UsePreCalc, .true., PowerTabler2i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms/2, NumTerms/2, NumTerms/2, Alphai, Alphaj, Betai, Betaj, &
							Gammai, Gammaj, 0, LValue, 0, 0, 0, 0, 0, 0, LValue, 0, 0, 0, qmax, pmax, linf, minf, IsTriplet, Method)
	else
		call CalcMatricesSub(RunCalc, UsePreCalc, .false., PowerTabler2i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms/2, NumTerms/2, NumTerms/2, Alphai, Alphaj, Betai, Betaj, &
							Gammai, Gammaj, 0, LValue, 0, 0, 0, 0, 0, LValue, 0, 0, 0, 0, qmax, pmax, linf, minf, IsTriplet, Method)
	end if
		
	! 1-2 matrix elements
	write (*,*) "Calculating 1-2 matrix elements"
	if (Exchanged == .true.) then  ! Swap l2r and l3r
		call CalcMatricesSub(RunCalc, UsePreCalc, .true., PowerTabler1i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms/2, 0, NumTerms/2, Alphai, Alphaj, Betai, Betaj, &
							Gammai, Gammaj, LValue, 0, 0, 0, 0, 0, 0, 0, LValue, 0, 0, 0, qmax, pmax, linf, minf, IsTriplet, Method)
	else
		call CalcMatricesSub(RunCalc, UsePreCalc, .true.,  PowerTabler1i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms/2, 0, NumTerms/2, Alphai, Alphaj, Betai, Betaj, &
							Gammai, Gammaj, LValue, 0, 0, 0, 0, 0, 0, LValue, 0, 0, 0, 0, qmax, pmax, linf, minf, IsTriplet, Method)
	end if
	
	! 2-1 matrix elements
	write (*,*) "Calculating 2-1 matrix elements"
	call CalcMatricesSub(RunCalc, UsePreCalc, .true., PowerTabler2i, PowerTabler1j, WMatrix, PhiPhi, Phi2HPhi, NumTerms/2, NumTerms/2, 0, Alphai, Alphaj, Betai, Betaj, &
						Gammai, Gammaj, 0, LValue, 0, 0, 0, 0, LValue, 0, 0, 0, 0, 0, qmax, pmax, linf, minf, IsTriplet, Method)
		
	return
end


!real*8 function ederivwrapper(Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, Alpha, Beta, Gamma, Iter, Energy)
!	implicit none
!	real*16, allocatable, dimension(:,:) :: PhiPhi, Phi2HPhi
!	real*8, allocatable, dimension(:,:) :: PhiPhi8, Phi2HPhi8
!	real*8, allocatable, dimension(:) :: Energies, Workspace
!	integer Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, Info, Iter
!	real*8 Alpha, Beta, Gamma, Energy
!	real*8 EDeriv
!
!	allocate(PhiPhi(NumTerms,NumTerms))
!	allocate(Phi2HPhi(NumTerms,NumTerms))
!	allocate(Energies(NumTerms))
!	allocate(Workspace(3*NumTerms-1))  ! Not really needed for Newdsygv?
!
!	Energies = 0.0d0
!
!	Energy = EDeriv(Omega, NumTerms, Alpha, Beta, Gamma, 1, PhiPhi, Phi2HPhi, PhiPhi8, Phi2HPhi8, &
!					qmax, pmax, 0, IsTriplet, Ordering, Method, EigenRoutine)
!
!	Phi2HPhi = Phi2HPhi / 2.0q0  ! TODO: Remove?
!	if (EigenRoutine == 1) then
!		allocate(PhiPhi8(NumTerms,NumTerms))
!		allocate(Phi2HPhi8(NumTerms,NumTerms))
!		PhiPhi8 = PhiPhi
!		Phi2HPhi8 = Phi2HPhi
!		call dsygv(1, 'N', 'L', NumTerms, Phi2HPhi8, NumTerms, PhiPhi8, NumTerms, Energies, Workspace, 3*NumTerms-1, Info)
!		! TODO: Test Info!
!		write (*,*) "Info: ", Info
!		deallocate(PhiPhi8, Phi2HPhi8)
!	else  ! EigenRoutine == 2
!		call Newdsygv(1, 'N', 'L', NumTerms, Phi2HPhi, NumTerms, PhiPhi, NumTerms, Energies, Workspace, 3*NumTerms-1, Info)
!	end if
!	
!	Energy = Energies(1)
!
!	write (*,"(A,3F8.5)") "Alpha, Beta, Gamma: ", Alpha, Beta, Gamma
!	write (10,"(I4,A,3F8.5)") Iter, "; Alpha, Beta, Gamma: ", Alpha, Beta, Gamma  ! I dislike hardcoding the 10, but it's easier at this point.
!	write (10,"(A,F16.12)") "Energy: ", Energy
!	write (10,*)
!	flush(10)
!	
!	deallocate(PhiPhi, Phi2HPhi, Energies, Workspace)
!	ederivwrapper = 1.0d0
!end


! Calculates the derivative of the energy function.  If the Deriv parameter = 0, then it exits early with
!  PhiPhi and Phi2HPhi filled.
real*8 function EDeriv(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, PhiPhi8, Phi2HPhi8,&
						qmax, pmax, linf, minf, Deriv, IsTriplet, Ordering, Method, EigenRoutine, LValue)
	use WLimits
	!use nag_nsym_gen_eig
	implicit none
	real*8 Alpha, Beta, Gamma, Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, e
	integer, allocatable, dimension(:,:) :: PowerTabler1i, PowerTabler1j, PowerTabler2i, PowerTabler2j  ! This contains the powers of the Hylleraas terms.
	integer CalcPowerTableSize, PreCalcGammaSize  ! Function definitions
	integer EigenNum  ! This is the energy eigenvalue number we are dealing with.
	integer Omega, NumTerms, qmax, pmax, linf, minf, Deriv, IsTriplet, Method, Ordering, i, j, n, TempJ, Info, EigenRoutine, LValue
	real*16, dimension(NumTerms,NumTerms) :: PhiPhi, Phi2HPhi
	real*8, dimension(NumTerms,NumTerms) :: PhiPhi8, Phi2HPhi8
	real*8, allocatable, dimension(:,:,:,:) :: WMatrix
	real*16 PhiHPhiSum, PhiPhiSum, Norm
	integer l, m, k
	logical UsePreCalc
	
	UsePreCalc = .true.
	!UsePreCalc = .false.

	allocate(PowerTabler1i(NumTerms/2, 6))
	call GenOmegaPowerTable(Omega, PowerTabler1i, NumTerms/2, Ordering, LValue, 0)
	allocate(PowerTabler1j(NumTerms/2, 6))
	call GenOmegaPowerTable(Omega, PowerTabler1j, NumTerms/2, Ordering, LValue, 0)
	allocate(PowerTabler2i(NumTerms/2, 6))
	call GenOmegaPowerTable(Omega, PowerTabler2i, NumTerms/2, Ordering, 0, LValue)
	allocate(PowerTabler2j(NumTerms/2, 6))
	call GenOmegaPowerTable(Omega, PowerTabler2j, NumTerms/2, Ordering, 0, LValue)

	! Set this matrix equal to 0 because of the Phi2HPhi(i,j) = Phi2HPhi(i,j) + Sum line in CalcMatrices
	!  (and equivalent for PhiPhi).  This allows us to use the same subroutine for the direct and exchange elements.
	Phi2HPhi = 0.0_16
	PhiPhi = 0.0_16
	
	! For the direct integration, f_i and f_j have beta and gamma in the same place.
	Alphai = Alpha; Alphaj = Alpha; Betai = Beta; Betaj = Beta; Gammai = Gamma; Gammaj = Gamma;

	! Precalculate the W matrices.  lmin, lmax, etc. are the upper and lower limits for each of the three
	!  parameters to the W function (in HylleraasIntegral.f90).  These limits were determined by me, and
	!  there is a document detailing these limits (W Function Limits.pdf).
	lmin = 1;  lmax = 2*(Omega+LValue) + 4*max(linf,minf) + 4;
	mmin = -1 - 3*max(linf,minf);  mmax = 2*(Omega+LValue) + 3*max(linf,minf) + 3;
	nmin = -2 - 4*max(linf,minf);  nmax = 2*(Omega+LValue) + 2;
	allocate(WMatrix(lmin:lmax, mmin:mmax, nmin:nmax, 6))
	
	call GenFactorial(199)

	write (*,*)
	write (*,*) "Starting calculation of direct-direct terms"
	if (UsePreCalc == .true.) then
		write (*,*) "Precomputing direct-direct W matrix"
		call ClearWMatrix(WMatrix)
		! First call it with RunCalc = false to determine what W functions to calculate.
		call CalcMatrices(.false., .true., PowerTabler1i, PowerTabler1j, PowerTabler2i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms, &
							Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, LValue, qmax, pmax, linf, minf, IsTriplet, Method, .false.)
		! Then precalculate all W functions that are used.
		call CalcWMatrices(Omega, qmax, WMatrix, Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, pmax)
		write (*,*) "Finished precomputing direct-direct W matrix"
	end if
	! This subroutine calculates the PhiPhi and Phi2HPhi matrices for the direct integrals.
	call CalcMatrices(.true., UsePreCalc, PowerTabler1i, PowerTabler1j, PowerTabler2i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, LValue, qmax, pmax, linf, minf, IsTriplet, Method, .false.)
	if (IsTriplet == 1) then  ! Only do this for the triplet case.
		PhiPhi = -PhiPhi
		Phi2HPhi = -Phi2HPhi
	endif

	! Now we have to do the same thing for the exchange operator, P_23.  We operate it on f_j by switching the
	!  2<->3 powers in PowerTablej and recalculating, along with beta and gamma.
	do n = 1, NumTerms/2, 1
		! This exchanges r2 and r3.
		TempJ = PowerTabler1j(n,2)
		PowerTabler1j(n,2) = PowerTabler1j(n,4)  ! l with n
		PowerTabler1j(n,4) = TempJ
		TempJ = PowerTabler2j(n,2)
		PowerTabler2j(n,2) = PowerTabler2j(n,4)  ! l with n
		PowerTabler2j(n,4) = TempJ
		! This exchanges r12 and r13 (r31).
		TempJ = PowerTabler1j(n,3)
		PowerTabler1j(n,3) = PowerTabler1j(n,5)  ! m with p
		PowerTabler1j(n,5) = TempJ
		TempJ = PowerTabler2j(n,3)
		PowerTabler2j(n,3) = PowerTabler2j(n,5)  ! m with p
		PowerTabler2j(n,5) = TempJ
	enddo
	! Swap beta and gamma for f_j.
	Alphai = Alpha; Alphaj = Alpha; Betaj = Gamma; Betai = Beta; Gammaj = Beta; Gammai = Gamma;

	!i = 0
	!do l = lmin, lmax, 1
	!	do m = mmin, mmax, 1
	!		do n = nmin, nmax, 1
	!			do k = 1, 6, 1
	!				if (WMatrix(l, m, n, k) == 1000.0d0) i = i + 1
	!			end do
	!		end do
	!	end do
	!end do
	!
	!write (*,*) "Total WMatrix used: ", i
	!write (*,*) "Total WMatrix size: ", (lmax-lmin+1) * (mmax-mmin+1) * (nmax-nmin+1) * 6
	!stop
	
	! We have to recalculate the W matrices again, since we had the swaps.	
	write (*,*) "Starting calculation of direct-exchange terms"
	if (UsePreCalc == .true.) then
		write (*,*) "Precomputing direct-exchange W matrix"
		call ClearWMatrix(WMatrix)
		! The W functions to calculate are not, in general, the same as that for the direct-direct case.
		call CalcMatrices(.false., UsePreCalc, PowerTabler1i, PowerTabler1j, PowerTabler2i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, LValue, qmax, pmax, linf, minf, IsTriplet, Method, .true.)
		call CalcWMatrices(Omega, qmax, WMatrix, Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, pmax)
		write (*,*) "Finished precomputing direct-exchange W matrix"
	end if
	! Calculate the matrix elements for the exchanged terms, and add them to the direct terms.
	call CalcMatrices(.true., UsePreCalc, PowerTabler1i, PowerTabler1j, PowerTabler2i, PowerTabler2j, WMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, LValue, qmax, pmax, linf, minf, IsTriplet, Method, .true.)

	! If we are not optimizing the nonlinear parameters, return since matrices have been calculated.
	EDeriv = 1.0d0  ! Random value, return not used in this case.

	! Clean up memory before exiting
	deallocate(WMatrix)
	deallocate(PowerTabler1i)
	deallocate(PowerTabler1j)
	deallocate(PowerTabler2i)
	deallocate(PowerTabler2j)

	return
end
