! This program calculates the matrix elements of the positronium-hydrogen
!  matrices (similar to equation 2.15 of the Armour and Humberston article).  Specifically,
!  we are calculating elements of the form (phi_i, L phi_j) as in equation (3.22).
! This does use OpenMP to speed up computation on multicore processors, so some familiarity
!  with parallel programming is required to understand the !$omp lines.
! A future speed-up could be done with Open MPI (as is done for the energy code), but this
!  computation is usually done in a few hours.
! Created: 08/02/2009

!@TODO: Make MaxIter and the 4 come from the configuration file.

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
		integer(c_int) function optimizewavefn(MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, GradGrad, alpha, beta, gamma) bind(c, name='optimizewavefn')
		use iso_c_binding
		integer(c_int), value :: MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, GradGrad
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
	integer iargc, IsTriplet, Method, LowerEigen, UpperEigen, AllEnergies, EigenRoutine
	integer Optimize, EigenNum, Ordering, GradGrad, Iter, MaxIter, MaxIterOuter, qmax, pmax
	integer CalcPowerTableSize, NumUsed, i, j, k, n
	real*8 StartTime, EndTime, fa
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
								AllEnergies, Optimize, EigenNum, MaxIter, Ordering, EigenRoutine, GradGrad)
	
	! Get the start time to determine the duration of the program.
	!call cpu_time(StartTime)
	StartTime = omp_get_wtime()

	NumTerms = CalcPowerTableSize(Omega)
	allocate(PhiPhi(NumTerms,NumTerms))
	allocate(Phi2HPhi(NumTerms,NumTerms))

	call WriteHeader(iwrite, Optimize, EigenNum, IsTriplet, Ordering, Method, Omega, Alpha, Beta, Gamma, NumTerms, GradGrad, EigenRoutine)
	
	if (Optimize == 1) then
		fa = optimizewavefn(MaxIter, Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, GradGrad, dble(Alpha), dble(Beta), dble(Gamma))
		call CalcMatrices(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, qmax, pmax, &
					IsTriplet, Ordering, Method, GradGrad)
	else
		! If we are not optimizing eigenvalues, we can still call the same function, which exits early with PhiPhi and Phi2HPhi filled.
		call CalcMatrices(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, qmax, pmax, &
					IsTriplet, Ordering, Method, GradGrad)
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

	! Close file handle
	close(iwrite)
stop
end


! Read the input parameters.
logical function ReadParamFile(iread, Omega, Alpha, Beta, Gamma, qmax, Method, pmax, IsTriplet, LowerEigen, UpperEigen, &
								AllEnergies, Optimize, EigenNum, MaxIter, Ordering, EigenRoutine, GradGrad)
	implicit none
	integer iread, Omega, qmax, Method, pmax, IsTriplet, LowerEigen, UpperEigen
	integer AllEnergies, Optimize, EigenNum, MaxIter, Ordering, EigenRoutine, GradGrad
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
	read (iread,*) ! Description line
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
	read (iread,*)
	read (iread,*) GradGrad
	
	close(iread)  ! We are done with this parameter file.
	
	if (IsTriplet /= 0 .and. IsTriplet /= 1) then
		write (*,*) 'Must choose between singlet and triplet calculation - exiting.'
		stop
	endif
	
	ReadParamFile = .true.
	return
end


! Outputs header information to the output file and to the screen				
subroutine WriteHeader(iwrite, Optimize, EigenNum, IsTriplet, Ordering, Method, Omega, Alpha, Beta, Gamma, NumTerms, GradGrad, EigenRoutine)
	implicit none
	integer iwrite, Optimize, EigenNum, IsTriplet, Ordering, Method, Omega, NumTerms, GradGrad, EigenRoutine
	real*16 Alpha, Beta, Gamma
	character (len=8) cdate
	character (len=8) ctime

	if (Optimize == 1) write (iwrite,*) 'Optimizing eigenvalue - not a valid XML file!', EigenNum

	write (iwrite,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
	write (iwrite,'(a)') "<psh_data>"
	write (iwrite,'(a)') "<header>"
	if (IsTriplet == 0) then
		if (GradGrad == 1) then
			write (*,*) "S-Wave Singlet Ps-H - Grad-Grad Formalism"
			write (iwrite,'(a)') "    <problem>S-Wave Singlet Ps-H - Grad-Grad Formalism</problem>"
		else
			write (*,*) "S-Wave Singlet Ps-H - Laplacian Formalism"
			write (iwrite,'(a)') "    <problem>S-Wave Singlet Ps-H - Laplacian Formalism</problem>"
		end if
	else
		if (GradGrad == 1) then
			write (*,*) "S-Wave Triplet Ps-H - Grad-Grad Formalism"
			write (iwrite,'(a)') "    <problem>S-Wave Triplet Ps-H - Grad-Grad Formalism</problem>"
		else
			write (*,*) "S-Wave Triplet Ps-H - Laplacian Formalism"
			write (iwrite,'(a)') "    <problem>S-Wave Triplet Ps-H - Laplacian Formalism</problem>"
		end if
	endif

	write (iwrite,'(a)') "    <lvalue>0</lvalue>"  ! Indicates S-wave
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

	if (GradGrad == 1) then
		write (iwrite,'(a)') "    <gradlapl>Gradient-Gradient</gradlapl>"
	else  ! Assuming Laplacian
		write (iwrite,'(a)') "    <gradlapl>Laplacian</gradlapl>"
	end if
	
	write (iwrite,'(a)') "    <formalism>1</formalism>"  ! Only one type that we can do for the S-wave
	
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

	if (EigenRoutine == 1) then
		write (*,*) "Eigenvalue routine: dsygv from LAPACK"
		write (iwrite,'(a)') "    <eigenroutine>LAPACK</eigenroutine>"
	else if (EigenRoutine == 2) then
		write (*,*) "Eigenvalue routine: Pachucki eigenproblem solver"
		write (iwrite,'(a)') "    <eigenroutine>Pachucki</eigenroutine>"
	else
		write (*,*) "Eigenvalue routine parameter in input file must be 1 or 2."
		stop
	end if

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
	    if (LValue == 0 .and. j > NumTerms / 2) then  ! S-wave only has one symmetry
            exit
        end if
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

	
! This subroutine calculates the number of terms for the function GenOmegaPowerTable (below).  This
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
!@TODO: Use logical for Ordering
subroutine GenOmegaPowerTable(Omega, PowerTable, ArraySize, Ordering)
	!implicit none
	integer Omega  ! This sets the limits on the terms.
	integer ArraySize
	integer, dimension(ArraySize,6) :: PowerTable
	integer NumTerm  ! The number of the current term
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

						   INX =INX +1
						   PowerTable(INX,1) = I1
						   PowerTable(INX,2) = I2
						   PowerTable(INX,3) = I12
						   PowerTable(INX,4) = I3
						   PowerTable(INX,5) = I13
						   PowerTable(INX,6) = I23
							write (11,'(i4,a,6i4)') INX, ": ", I1, I2, I12, I3, I13, I23
	70                  CONTINUE
	69                CONTINUE
	68              CONTINUE
	67            CONTINUE
	66          CONTINUE
	65        CONTINUE
	endif
	
	return
end


! I am using the same numbering as in Dr. Quintanilla's table (Tablecoeff.pdf).
!  This was modified to differentiate between the alpha, beta and gamma of f_i and f_j,
!  and the table has been changed to reflect this.
subroutine GenCoeffTableGradGrad(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
	implicit none
	real*16, dimension(34) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	real*8 ScatteringVal
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj

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
	
!	if (IsScattering == 0) ScatteringVal = 0.0q0
!	if (IsScattering == 1) ScatteringVal = 1.5q0
	ScatteringVal = 0.0q0

	! These coefficients are taken from equation (3.25) of Armour and Humberston's article,
	!  with one minor exception: we have included gamma in the derivation.  Gamma = Beta in the case
	!  of Ps-He scattering but does not for Ps-H.
	CoeffTable(1) =   Alphai * Alphaj + Betai * Betaj + Gammai * Gammaj + ScatteringVal
	CoeffTable(2) =   2.0q0 - (Alphaj*ki + Alphai*kj) - 0.5q0 * (Alphaj*mi + Alphai*mj) - 0.5q0 * (Alphaj*pi + Alphai*pj)
	CoeffTable(3) =  -2.0q0 - (Betaj*li + Betai*lj) - 0.5q0 * (Betaj*mi + Betai*mj) - 0.5q0 * (Betaj*qi + Betai*qj)
	CoeffTable(4) =  -2.0q0
	CoeffTable(5) =  -2.0q0 - (Gammaj*ni + Gammai*nj) - 0.5q0 * (Gammaj*qi + Gammai*qj) - 0.5q0 * (Gammaj*pi + Gammai*pj)
	CoeffTable(6) =  -2.0q0
	CoeffTable(7) =   2.0q0
	CoeffTable(8) =   0.5q0 * (ki*mj + kj*mi + ki*pj + kj*pi) + ki*kj
	CoeffTable(9) =   0.5q0 * (li*mj + lj*mi + li*qj + lj*qi) + li*lj
	CoeffTable(10) =  0.5q0 * (ki*mj + kj*mi + li*mj + lj*mi + mi*pj + mj*pi + mi*qj + mj*qi) + 2.0q0*mi*mj
	CoeffTable(11) =  0.5q0 * (ni*qj + nj*qi + ni*pj + nj*pi) + ni*nj
	CoeffTable(12) =  0.5q0 * (ki*pj + kj*pi + ni*pj + nj*pi + mi*pj + mj*pi + pi*qj + pj*qi) + 2.0q0*pi*pj
	CoeffTable(13) =  0.5q0 * (li*qj + lj*qi + ni*qj + nj*qi + mi*qj + mj*qi + pi*qj + pj*qi) + 2.0q0*qi*qj
	CoeffTable(14) = -0.5q0 * (Alphaj*mi + Alphai*mj)
	CoeffTable(15) = -0.5q0 * (Betaj*mi + Betai*mj)
	CoeffTable(16) =  0.5q0 * (Alphaj*mi + Alphai*mj)
	CoeffTable(17) = -0.5q0 * (ki*mj + kj*mi)
	CoeffTable(18) =  0.5q0 * (Betaj * mi + Betai*mj)
	CoeffTable(19) = -0.5q0 * (li*mj + lj*mi)
	CoeffTable(20) = -0.5q0 * (Alphaj*pi + Alphai*pj)
	CoeffTable(21) = -0.5q0 * (Gammaj*pi + Gammai*pj)
	CoeffTable(22) =  0.5q0 * (Alphaj*pi + Alphai*pj)
	CoeffTable(23) = -0.5q0 * (ki*pj + kj*pi)
	CoeffTable(24) =  0.5q0 * (Gammaj*pi + Gammai*pj)
	CoeffTable(25) = -0.5q0 * (ni*pj + nj*pi)
	CoeffTable(26) = -0.5q0 * (Gammaj*qi + Gammai*qj)
	CoeffTable(27) = -0.5q0 * (Betaj*qi + Betai*qj)
	CoeffTable(28) =  0.5q0 * (Gammaj*qi + Gammai*qj)
	CoeffTable(29) = -0.5q0 * (ni*qj + nj*qi)
	CoeffTable(30) =  0.5q0 * (Betaj*qi + Betai*qj)
	CoeffTable(31) = -0.5q0 * (li*qj + lj*qi)
	CoeffTable(32) = -0.5q0 * (pi*qj + pj*qi)
	CoeffTable(33) = -0.5q0 * (mi*pj + mj*pi)
	CoeffTable(34) = -0.5q0 * (mi*qj + mj*qi)
	
	return
end


! Note that the naming is a little different here.  I used rPowerTable here as the parameter name but
!  use rPowers elsewhere.
subroutine GenRPowerTableGradGrad(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
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
	data rPowers / 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 0,-1,-2, 2, 2, 1, 0,-1,-2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, &  ! r1
				   0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 0, 0, 0, &  ! r2
				   0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0,-2,-2,-2,-2,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2,-2, &  ! r12
				   0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 1, 0,-1,-2, 2, 2, 0, 0, 0, &  ! r3
				   0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2,-2,-2, 0, 0, 0, 0, 0, 0,-2,-2, 2, &  ! r13
				   0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2,-2,-2,-2, 2,-2 /  ! r23

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


! This calculates the Hamiltonian directly with the Laplacians, as in Drake and Yan's 1997 paper. 
subroutine GenCoeffTableLaplacian(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(34) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer kj, lj, mj, nj, pj, qj, l1r, l2r, l3r

	l1r = 0; l2r = 0; l3r = 0;  ! Always 0 for the S-wave

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

	! These coefficients are taken from equation (3.25) of Armour and Humberston's article,
	!  with one minor exception: we have included gamma in the derivation.  Gamma = Beta in the case
	!  of Ps-He scattering but does not for Ps-H.
	CoeffTable(1) = -Alphaj*Alphaj - Betaj*Betaj - Gammaj*Gammaj
	CoeffTable(2) = -kj - kj*kj - kj*mj - kj*pj + l1r + l1r*l1r
	CoeffTable(3) = -2.0q0*mj - kj*mj - 2.0q0*mj*mj - mj*lj - mj*qj - mj*pj
	CoeffTable(4) = -2.0q0
	CoeffTable(5) = -lj - mj*lj - lj*lj - lj*qj + l2r + l2r*l2r
	CoeffTable(6) = mj*lj
	CoeffTable(7) = kj*mj
	CoeffTable(8) = -2.0q0*qj - mj*qj - lj*qj - 2.0q0*qj*qj - qj*nj - qj*pj
	CoeffTable(9) = 2.0q0
	CoeffTable(10) = -nj - qj*nj - nj*nj - nj*pj + l3r + l3r*l3r
	CoeffTable(11) = qj*nj
	CoeffTable(12) = lj*qj
	CoeffTable(13) = -2.0q0*pj - kj*pj - mj*pj - qj*pj - nj*pj - 2.0q0*pj*pj
	CoeffTable(14) = qj*pj
	CoeffTable(15) = mj*pj
	CoeffTable(16) = nj*pj
	CoeffTable(17) = kj*pj
	CoeffTable(18) = -2.0q0
	CoeffTable(19) = mj*qj
	CoeffTable(20) = mj*Alphaj
	CoeffTable(21) = -mj*Alphaj
	CoeffTable(22) = pj*Alphaj
	CoeffTable(23) = -pj*Alphaj
	CoeffTable(24) = 2.0q0 + 2.0q0*Alphaj + 2.0q0*kj*Alphaj + mj*Alphaj + pj*Alphaj
	CoeffTable(25) = -mj*Betaj
	CoeffTable(26) = mj*Betaj
	CoeffTable(27) = qj*Betaj
	CoeffTable(28) = -qj*Betaj
	CoeffTable(29) = -2.0q0 + 2.0q0*Betaj + mj*Betaj + 2.0q0*lj*Betaj + qj*Betaj
	CoeffTable(30) = -qj*Gammaj
	CoeffTable(31) = qj*Gammaj
	CoeffTable(32) = -pj*Gammaj
	CoeffTable(33) = pj*Gammaj
	CoeffTable(34) = -2.0q0 + 2.0q0*Gammaj + qj*Gammaj + 2.0q0*nj*Gammaj + pj*Gammaj

	return
end


subroutine GenRPowerTableLaplacian(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
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
						if (WMatrix(l, m, n, k) /= 0.0q0) then
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
	write (*,*)
	
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


!@TODO: Rename this to something more descriptive.
subroutine CalcMatricesSub(RunCalc, CoeffTable, rPowers, PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, qmax, pmax, IsTriplet, Method, GradGrad)
	use WCLimits
	implicit none
	real*16, dimension(34) :: CoeffTable
	integer, dimension(34,6) :: rPowers
	integer, dimension(NumTerms,6) :: PowerTablei, PowerTablej
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*16, dimension(-2:Cjmax, 0:Cqmax, 0:Ckmax) :: CMatrix
	real*16, dimension(NumTerms,NumTerms) :: PhiPhi, Phi2HPhi
	integer qmax, pmax, NumTerms, Method, GradGrad
	logical RunCalc
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	real*16 HylleraasIntegral
	real*16 Sum
	real*8 PI
	integer i, j, n, IsTriplet
	real*16 RemoveMe
	real*16 PiHalf
	integer WUsed, WUnused, l, m
	! Calculated in Mathematica
	PiHalf = 0.1591549430918953357688837633725143620345q0
	
	!$omp parallel do shared(NumTerms) private(j,Sum,RemoveMe,CoeffTable,rPowers) schedule(dynamic,10)
	do i = 1, NumTerms, 1
		write (*,*) i
		!write (*,"(i7)",advance='no') i
		!do j = 1, i, 1
		do j = 1, NumTerms, 1
			Sum = 0.0q0
			if (GradGrad == 0) then  ! Use grad-grad formalism
				call GenCoeffTableGradGrad(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, &
									Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
				call GenRPowerTableGradGrad(rPowers, PowerTablei, i, PowerTablej, j, NumTerms)
			else  ! Use Laplacian formalism
				call GenCoeffTableLaplacian(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
				call GenRPowerTableLaplacian(rPowers, PowerTablei, i, PowerTablej, j, NumTerms)
			end if

			do n = 1, 34, 1
				if (CoeffTable(n) /= 0.0q0) then
					! Note that the order of the parameters to the HylleraasIntegral function looks strange.  This is
					!  because the Drake and Yan paper (1995) expressed the integrand in the following order:
					!  r1 r2 r3 r12 r23 r31, and the Armour and Humberston article has them in this order:
					!  r1 r2 r12 r3 r13 r23 (equation 3.14).  So we have to swap the third and fourth parameters, along
					!  with the fifth and sixth.
					
					RemoveMe = HylleraasIntegral(RunCalc, rPowers(n,1), rPowers(n,2), rPowers(n,4), rPowers(n,3), rPowers(n,6), &
												  rPowers(n,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, pmax, Method, WMatrix, CMatrix)
					Sum = Sum + RemoveMe * CoeffTable(n)
				endif
			enddo

			if (IsTriplet == 0) then
				Phi2HPhi(i,j) = Phi2HPhi(i,j) + Sum * PiHalf
			else
				Phi2HPhi(i,j) = Phi2HPhi(i,j) - Sum * PiHalf
			endif

			! Do the PhiPhi inner product
			RemoveMe = HylleraasIntegral(RunCalc, PowerTablei(i,1)+PowerTablej(j,1), PowerTablei(i,2)+PowerTablej(j,2), &
						PowerTablei(i,4)+PowerTablej(j,4), PowerTablei(i,3)+PowerTablej(j,3), &
						PowerTablei(i,6)+PowerTablej(j,6), PowerTablei(i,5)+PowerTablej(j,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, pmax, Method, WMatrix, CMatrix)
			if (IsTriplet == 0) then
				PhiPhi(i,j) = PhiPhi(i,j) + RemoveMe * PiHalf
			else
				PhiPhi(i,j) = PhiPhi(i,j) - RemoveMe * PiHalf
			endif
		enddo
	enddo

	write (*,*)
end


! Called by optimizewavefn in Minimize.c
subroutine calcmatriceswrapper(Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, GradGrad, Alpha, Beta, Gamma, Energy) bind(C)
	implicit none
	real*16, allocatable, dimension(:,:) :: PhiPhi, Phi2HPhi
	real*8, allocatable, dimension(:,:) :: PhiPhi8, Phi2HPhi8
	real*8, allocatable, dimension(:) :: Energies, Workspace
	integer Omega, NumTerms, qmax, pmax, IsTriplet, Ordering, Method, EigenRoutine, GradGrad, Info
	integer, save :: Iter = 1
	real*8 Alpha, Beta, Gamma
	real*8 Energy
	
	allocate(PhiPhi(NumTerms,NumTerms))
	allocate(Phi2HPhi(NumTerms,NumTerms))
	allocate(Energies(NumTerms))
	allocate(Workspace(3*NumTerms-1))  ! Not really needed for Newdsygv?
	Energies = 0.0q0

	call CalcMatrices(Omega, NumTerms, real(Alpha,16), real(Beta,16), real(Gamma,16), 1, PhiPhi, Phi2HPhi, &
					qmax, pmax, IsTriplet, Ordering, Method, GradGrad)

	Phi2HPhi = Phi2HPhi / 2
	if (EigenRoutine == 1) then
		allocate(PhiPhi8(NumTerms,NumTerms))
		allocate(Phi2HPhi8(NumTerms,NumTerms))
		PhiPhi8 = PhiPhi
		Phi2HPhi8 = Phi2HPhi
		call dsygv(1, 'N', 'L', NumTerms, Phi2HPhi8, NumTerms, PhiPhi8, NumTerms, Energies, Workspace, 3*NumTerms-1, Info)
		deallocate(PhiPhi8, Phi2HPhi8)
	else  ! EigenRoutine == 2
		call Newdsygv(1, 'N', 'L', NumTerms, Phi2HPhi, NumTerms, PhiPhi, NumTerms, Energies, Workspace, 3*NumTerms-1, Info)
	end if
	
	Energy = Energies(1)

	write (*,"(A,3F8.5)") "Alpha, Beta, Gamma: ", Alpha, Beta, Gamma
	write (10,"(I4,A,3F8.5)") Iter, "; Alpha, Beta, Gamma: ", Alpha, Beta, Gamma  ! I dislike hardcoding the 10, but it's easier at this point.
	write (10,"(A,F16.12)") "Energy: ", Energy
	write (10,*)
	flush(10)
	
	deallocate(PhiPhi, Phi2HPhi, Energies, Workspace)
	Iter = Iter + 1
end


! Calculates the derivative of the energy function.  If the Deriv parameter = 0, then it exits early with
!  PhiPhi and Phi2HPhi filled.
subroutine CalcMatrices(Omega, NumTerms, Alpha, Beta, Gamma, EigenNum, PhiPhi, Phi2HPhi, qmax, pmax, IsTriplet, Ordering, Method, GradGrad)
	use WCLimits
	!use nag_nsym_gen_eig
	implicit none
	real*16 Alpha, Beta, Gamma, Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	real*8 e
	integer, allocatable, dimension(:,:) :: PowerTablei, PowerTablej, rPowers  ! This contains the powers of the Hylleraas terms.
	integer CalcPowerTableSize, PreCalcGammaSize  ! Function definitions
	integer EigenNum  ! This is the energy eigenvalue number we are dealing with.
	integer Omega, NumTerms, qmax, pmax, IsTriplet, Method, GradGrad, Ordering, i, j, n, Info
	real*16, allocatable, dimension(:) :: CoeffTable
	real*16, dimension(NumTerms,NumTerms) :: PhiPhi, Phi2HPhi
	real*16, allocatable, dimension(:,:,:,:) :: WMatrix
	real*16, allocatable, dimension(:,:,:) :: CMatrix
	real*16 PhiHPhiSum, PhiPhiSum, Norm
	
	allocate(PowerTablei(NumTerms, 6))
	call GenOmegaPowerTable(Omega, PowerTablei, NumTerms, Ordering)
	allocate(PowerTablej(NumTerms, 6))
	call GenOmegaPowerTable(Omega, PowerTablej, NumTerms, Ordering)
	allocate(CoeffTable(34))
	allocate(rPowers(34,6))

	! Set this matrix equal to 0 because of the Phi2HPhi(i,j) = Phi2HPhi(i,j) + Sum line in CalcMatricesSub
	!  (and equivalent for PhiPhi).  This allows us to use the same subroutine for the direct and exchange elements.
	Phi2HPhi = 0.0_16
	PhiPhi = 0.0_16
	
	! For the direct integration, f_i and f_j have beta and gamma in the same place.
	Alphai = Alpha; Alphaj = Alpha; Betai = Beta; Betaj = Beta; Gammai = Gamma; Gammaj = Gamma;

	! Precalculate the W matrices.  lmin, lmax, etc. are the upper and lower limits for each of the three
	!  parameters to the W function (in HylleraasIntegral.f90).  These limits were determined by me, and
	!  there is a document detailing these limits (W Function Limits.pdf).
	lmin = 1;  lmax = 2*(Omega+1) + 2*qmax + 4;
	mmin = 0;  mmax = 2*(Omega+1) + 3;
	nmin = -1 - 2*qmax;  nmax = 2*(Omega+1) + 2;
	allocate(WMatrix(lmin:lmax, mmin:mmax, nmin:nmax, 6))
	Cjmax = 2*Omega+4;  Cqmax = qmax;  Ckmax = 2*Omega+5;
	allocate(CMatrix(-2:Cjmax, 0:Cqmax, 0:Ckmax))

	if (Method == 0 .or. Method == 1) then
		write (*,*) "Precomputing W matrix"
		call ClearWCMatrices(WMatrix, CMatrix)
		! First call it with RunCalc = false to determine what W functions to calculate.
		call CalcMatricesSub(.false., CoeffTable, rPowers, PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
					Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, qmax, pmax, IsTriplet, Method, GradGrad)
		! Then precalculate all W functions that are used.
		call CalcWCMatrices(Omega+1, qmax, WMatrix, CMatrix, Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, pmax)
		write (*,*) "Finished precomputing W matrix"
	else  ! Method == 2 (recursion relations)
		write (*,*) "Precomputing gamma functions"
		call PreCalcGamma(real(Alphai+Alphaj,16), real(Betai+Betaj,16), real(Gammai+Gammaj,16), PreCalcGammaSize(Omega))
		write (*,*) "Finished precomputing gamma functions"
	end if
 
	! This subroutine call calculates the PhiPhi and Phi2HPhi matrices for the direct integrals.
	call CalcMatricesSub(.true., CoeffTable, rPowers, PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, qmax, pmax, IsTriplet, Method, GradGrad)
	!goto 100
	if (IsTriplet == 1) then  ! Only do this for the triplet case.
		PhiPhi = -PhiPhi
		Phi2HPhi = -Phi2HPhi
	endif

	! Now we have to do the same thing for the exchange operator, P_23.
	call Permute(PowerTablej, NumTerms, Alphai, Alphaj, Alpha, Betai, Betaj, Beta, Gammai, Gammaj, Gamma)

	! We have to recalculate the W matrices again, since we had the swaps.  The limits are still the same as above,
	!  though, since the nonlinear parameters do not affect the limits.
	if (Method == 0 .or. Method == 1) then
		write (*,*) "Precomputing W matrix"
		call ClearWCMatrices(WMatrix, CMatrix)
		call CalcMatricesSub(.false., CoeffTable, rPowers, PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
							Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, qmax, pmax, IsTriplet, Method, GradGrad)
		call CalcWCMatrices(Omega+1, qmax, WMatrix, CMatrix, Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, pmax)
		write (*,*) "Finished precomputing W matrix"
	else  ! Method == 2
		write (*,*) "Precomputing gamma functions"
		call PreCalcGamma(real(Alphai+Alphaj,16), real(Betai+Betaj,16), real(Gammai+Gammaj,16), PreCalcGammaSize(Omega))
		write (*,*) "Finished precomputing gamma functions"
	end if

	! Calculate the matrix elements for the exchanged terms, and add them to the direct terms.
	call CalcMatricesSub(.true., CoeffTable, rPowers, PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, qmax, pmax, IsTriplet, Method, GradGrad)

	! Clean up memory before exiting
100	deallocate(WMatrix, CMatrix)
	deallocate(rPowers, CoeffTable)
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
