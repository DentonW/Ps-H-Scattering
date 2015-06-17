!
! This file calculates the direct-exchange parts of the matrix elements.  For the direct-direct
!  part, refer to the file Direct.f90.
!


! Called by GenCoeffTablePhi2Phi2Exchange - The array here gets multipled by the terms of cos23 in GenCoeffTablePhi2Phi2Exchange.
subroutine SubGenCoeffTablePhi2Phi2Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(34) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)

	! This section only includes the grad-grad and potential parts
	CoeffTable(1) =   Alphai * Alphaj + Betai * Betaj + Gammai * Gammaj
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


subroutine GenCoeffTablePhi2Phi2ExchangeTrig(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(36) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)

	! This section includes the rest of the terms (mixture of cos and sin)
	!  Refer to the Mathematica notebook to see how to get to the expanded version.	
	CoeffTable(1) = -0.75q0 * mj - 0.375q0 * pi - 1.125q0 * qi - 0.75q0 * qj
	CoeffTable(2) = -0.375q0 * pi
	CoeffTable(3) = 0.75q0 * pi
	CoeffTable(4) = 0.75q0 * mj
	CoeffTable(5) = -0.375q0 * qi - 0.375q0 * qj
	CoeffTable(6) = 0.75q0 * mj + 1.125q0 * qj
	CoeffTable(7) = -0.75q0 * mj
	CoeffTable(8) = -0.75q0 * pi
	CoeffTable(9) = -0.375q0 * pi - 1.125q0 * qi
	CoeffTable(10) = 0.375q0 * pi
	CoeffTable(11) = 0.375q0 * qi
	CoeffTable(12) = 0.75q0 * pi + 1.125q0 * qi
	CoeffTable(13) = -0.75q0 * pi
	CoeffTable(14) = -0.375q0 * pi - 0.375q0 * qi
	CoeffTable(15) = 0.375q0 * pi
	CoeffTable(16) = -0.375q0 * mj - 0.75q0 * pi - 0.75q0 * qi - 1.125q0 * qj
	CoeffTable(17) = -0.375q0 * mj
	CoeffTable(18) = 0.75q0 * pi
	CoeffTable(19) = 0.75q0 * mj
	CoeffTable(20) = 0.375q0 * mj
	CoeffTable(21) = -0.375q0 * pi
	CoeffTable(22) = -0.375q0 * qi + 0.375q0 * qj
	CoeffTable(23) = -0.75q0 * mj
	CoeffTable(24) = 0.75q0 * mj + 0.75q0 * pi + 1.125q0 * qi + 1.125q0 * qj
	CoeffTable(25) = -0.75q0 * pi
	CoeffTable(26) = -0.75q0 * mj
	CoeffTable(27) = -0.375q0 * mj - 0.375q0 * qj
	CoeffTable(28) = 0.375q0 * mj
	CoeffTable(29) = 0.375q0 * mj
	CoeffTable(30) = 0.375q0 * pi
	CoeffTable(31) = -0.375q0 * mj - 1.125q0 * qj
	CoeffTable(32) = 0.375q0 * mj
	CoeffTable(33) = -0.375q0 * mj
	CoeffTable(34) = 0.375q0 * pi
	CoeffTable(35) = 0.375q0 * qi - 0.375q0 * qj
	CoeffTable(36) = 0.375q0 * qj
	!CoeffTable = 0.0q0
	!CoeffTable(1) =  -1.125q0 * qi - 0.75q0 * qj
	!CoeffTable(2) =  -0.375q0 * qi - 0.375q0 * qj
	!CoeffTable(3) =   1.125q0 * qj
	!CoeffTable(4) =  -1.125q0 * qi
	!CoeffTable(5) =   0.375q0 * qi
	!CoeffTable(6) =   1.125q0 * qi
	!CoeffTable(7) =  -0.375q0 * qi
	!CoeffTable(8) =  -0.75q0 * qi - 1.125q0 * qj
	!CoeffTable(9) =  -0.375q0 * qi + 0.375q0 * qj
	!CoeffTable(10) =  1.125q0 * qi + 1.125q0 * qj
	!CoeffTable(11) = -0.375q0 * qj
	!CoeffTable(12) = -1.125q0 * qj
	!CoeffTable(13) =  0.375q0 * qi - 0.375q0 * qj
	!CoeffTable(14) =  0.375q0 * qj

	return
end


! Called by GenCoeffTablePhi1Phi2Exchange - The array here gets multipled by the terms of cos13 in GenCoeffTablePhi1Phi2Exchange.
subroutine SubGenCoeffTablePhi1Phi2Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(34) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)
	
	! This section only includes the grad-grad and potential parts
	CoeffTable(1) =   Alphai * Alphaj + Betai * Betaj + Gammai * Gammaj
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


subroutine GenCoeffTablePhi1Phi2ExchangeTrig(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(36) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)
	
	! This section includes the rest of the terms (mixture of cos and sin)
	!  Refer to the Mathematica notebook to see how to get to the expanded version.
	CoeffTable(1) =  -0.75q0 * mj - 1.125q0 * pi - 0.75q0 * pj - 0.375q0 * qi
	CoeffTable(2) =  -0.375q0 * pi - 0.375q0 * pj
	CoeffTable(3) =   0.75q0 * mj + 1.125q0 * pj
	CoeffTable(4) =  -0.75q0 * mj
	CoeffTable(5) =   0.75q0 * qi
	CoeffTable(6) =  -0.75q0 * qi
	CoeffTable(7) =  -0.375q0 * qi
	CoeffTable(8) =   0.75q0 * mj
	CoeffTable(9) =  -1.125q0 * pi - 0.375q0 * qi
	CoeffTable(10) =  0.375q0 * pi
	CoeffTable(11) =  1.125q0 * pi + 0.75q0 * qi
	CoeffTable(12) = -0.375q0 * pi - 0.375q0 * qi
	CoeffTable(13) =  0.375q0 * qi
	CoeffTable(14) = -0.75q0 * qi
	CoeffTable(15) =  0.375q0 * qi
	CoeffTable(16) = -0.375q0 * mj - 0.75q0 * pi - 1.125q0 * pj - 0.75q0 * qi
	CoeffTable(17) =  0.375q0 * mj
	CoeffTable(18) = -0.375q0 * pi + 0.375q0 * pj
	CoeffTable(19) =  0.75q0 * mj + 1.125q0 * pi + 1.125q0 * pj + 0.75q0 * qi
	CoeffTable(20) = -0.75q0 * mj
	CoeffTable(21) = -0.375q0 * mj - 0.375q0 * pj
	CoeffTable(22) =  0.375q0 * mj
	CoeffTable(23) = -0.375q0 * mj
	CoeffTable(24) =  0.375q0 * mj
	CoeffTable(25) = -0.375q0 * qi
	CoeffTable(26) =  0.75q0 * qi
	CoeffTable(27) = -0.75q0 * qi
	CoeffTable(28) =  0.375q0 * qi
	CoeffTable(29) =  0.75q0 * mj
	CoeffTable(30) = -0.75q0 * mj
	CoeffTable(31) = -0.375q0 * mj - 1.125q0 * pj
	CoeffTable(32) = -0.375q0 * mj
	CoeffTable(33) =  0.375q0 * pi - 0.375q0 * pj
	CoeffTable(34) =  0.375q0 * mj
	CoeffTable(35) =  0.375q0 * qi
	CoeffTable(36) =  0.375q0 * pj

	return
end


! Called by GenCoeffTablePhi2Phi1Exchange - The array here gets multipled by the terms of cos12 in GenCoeffTablePhi2Phi1Exchange.
subroutine SubGenCoeffTablePhi2Phi1Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(34) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)
	
	! This section only includes the grad-grad and potential parts
	CoeffTable(1) =   Alphai * Alphaj + Betai * Betaj + Gammai * Gammaj
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


subroutine GenCoeffTablePhi2Phi1ExchangeTrig(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(36) :: CoeffTable
	integer i, j, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)
	
	! This section includes the rest of the terms (mixture of cos and sin)
	!  Refer to the Mathematica notebook to see how to get to the expanded version.
	CoeffTable(1) =  -0.75q0 * mi - 1.125q0 * mj - 0.75q0 * pi - 0.375q0 * qj
	CoeffTable(2) =  -0.375q0 * mi - 0.375q0 * mj
	CoeffTable(3) =   1.125q0 * mi + 0.75q0 * pi
	CoeffTable(4) =  -1.125q0 * mj - 0.375q0 * qj
	CoeffTable(5) =   0.375q0 * mj
	CoeffTable(6) =   1.125q0 * mj + 0.75q0 * qj
	CoeffTable(7) =  -0.375q0 * mj - 0.375q0 * qj
	CoeffTable(8) =  -1.125q0 * mi - 0.75q0 * mj - 0.375q0 * pi - 0.75q0 * qj
	CoeffTable(9) =   0.375q0 * mi - 0.375q0 * mj
	CoeffTable(10) =  1.125q0 * mi + 1.125q0 * mj + 0.75q0 * pi + 0.75q0 * qj
	CoeffTable(11) = -0.375q0 * mi - 0.375q0 * pi
	CoeffTable(12) =  0.375q0 * pi
	CoeffTable(13) = -0.75q0 * pi
	CoeffTable(14) =  0.375q0 * pi
	CoeffTable(15) = -1.125q0 * mi - 0.375q0 * pi
	CoeffTable(16) = -0.375q0 * mi + 0.375q0 * mj
	CoeffTable(17) = -0.375q0 * pi
	CoeffTable(18) =  0.375q0 * mi
	CoeffTable(19) = -0.75q0 * qj
	CoeffTable(20) =  0.75q0 * qj
	CoeffTable(21) = -0.375q0 * qj
	CoeffTable(22) =  0.375q0 * qj
	CoeffTable(23) =  0.75q0 * qj
	CoeffTable(24) = -0.75q0 * qj
	CoeffTable(25) =  0.375q0 * qj
	CoeffTable(26) =  0.75q0 * pi
	CoeffTable(27) =  0.75q0 * pi
	CoeffTable(28) = -0.75q0 * pi
	CoeffTable(29) = -0.75q0 * pi
	CoeffTable(30) = -0.375q0 * pi
	CoeffTable(31) =  0.375q0 * pi
	CoeffTable(32) =  0.375q0 * pi
	CoeffTable(33) = -0.375q0 * qj
	CoeffTable(34) =  0.375q0 * qj
	CoeffTable(35) = -0.75q0 * qj
	CoeffTable(36) =  0.375q0 * qj

	return
end


subroutine GenRPowerTablePhi2Phi2ExchangeTrig(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(36,6) :: rPowerTable
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms, n
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(36,6) :: rPowers

	data rPowers /	 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, &  ! r1
					-2,-2,-2,-2, 0,-4,-4,-2, 2, 2, 4, 0, 0,-2,-2, 0, 0, 0, 0, 2, 2, 2, 0,-2,-2,-2,-4,-4,-2,-2,-4,-4,-2,-2,-2,-4, &  ! r2
					 0, 0, 2,-2, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 2,-2,-2, 0, 0,-2, 0, 2,-2, 0,-2,-2, 0, 0,-2,-2, 0, 0, 0, &  ! r12
					 0, 0, 0, 0, 0, 0, 0, 0,-4,-4,-4,-4,-4,-4,-4,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 4, &  ! r3
					 0,-2,-2, 2, 0, 0, 0,-2, 0,-2, 0, 0,-2, 0,-2, 0, 0,-2, 2, 0,-2, 0, 0, 0,-2, 2, 0, 0, 0,-2, 0, 0, 0,-2, 0, 0, &  ! r13
					 0, 0, 0, 0,-2, 2, 2, 2, 0, 0,-2, 2, 2, 4, 4, 0, 0, 0, 0, 0, 0,-2, 2, 2, 2, 2, 4, 4, 4, 4, 0, 0, 0, 0,-2,-2 /   ! r23
	!! sin^2 terms only
	!data rPowers /	 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &  ! r1
	!				-2, 0,-4, 2, 4, 0,-2, 0, 2,-2,-4,-4,-2,-4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &  ! r2
	!				 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &  ! r12
	!				 0, 0, 0,-4,-4,-4,-4,-2,-2,-2,-2, 2, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &  ! r3
	!				 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &  ! r13
	!				 0,-2, 2, 0,-2, 2, 4, 0,-2, 2, 4, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /   ! r23	
	
	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)

	! Combine powers of like r's from phi_i, phi_j and the powers in dissertation.
	do n = 1, 36, 1
		rPowerTable(n,1) = ki + kj + rPowers(n,1)
		rPowerTable(n,2) = li + lj + rPowers(n,2)
		rPowerTable(n,3) = mi + mj + rPowers(n,3)
		rPowerTable(n,4) = ni + nj + rPowers(n,4)
		rPowerTable(n,5) = pi + pj + rPowers(n,5)
		rPowerTable(n,6) = qi + qj + rPowers(n,6)
	enddo
	
	return
end


! Called by GenRPowerTablePhi1Phi2Exchange - The array here gets multipled by the terms of cos13 in GenRPowerTablePhi1Phi2Exchange.
subroutine SubGenRPowerTablePhi1Phi2Exchange(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(34,6) :: rPowerTable
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms, n
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(34,6) :: rPowers

	data rPowers /	0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 0,-1,-2, 2, 2, 1, 0,-1,-2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, &  !r1
					0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 0, 0, 0, &  !r2
					0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0,-2,-2,-2,-2,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2,-2, &  !r12
					0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 1, 0,-1,-2, 2, 2, 0, 0, 0, &  !r3
					0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2,-2,-2, 0, 0, 0, 0, 0, 0,-2,-2, 2, &  !r13
					0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2,-2,-2,-2, 2, -2 /  !r23

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)

	! Combine powers of like r's from phi_i, phi_j and the powers in dissertation.
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


subroutine GenRPowerTablePhi1Phi2ExchangeTrig(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(36,6) :: rPowerTable
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms, n
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(36,6) :: rPowers

	data rPowers /	-2, 0,-4,-4,-2,-2,-2,-2, 2, 4, 0,-2, 2, 0,-2, 0, 2, 2,-2, 0,-4,-2, 0,-4, 2, 0,-2,-2, 0,-2,-4,-2,-2,-4,-2,-4, &  ! r1
					 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, &  ! r2
					 0, 0, 0,-2, 2, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0,-2, 0,-2,-2,-2, 0, 2, 2, 0,-2,-2, 0,-2, 0,-2, 0, 0, &  ! r12
					 0, 0, 0, 0, 0, 0, 0, 0,-4,-4,-4,-4,-4,-4,-4,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 4, &  ! r3
					 0,-2, 2, 2, 0, 2, 0, 0, 0,-2, 2, 4, 0, 2, 4, 0, 0,-2, 2, 2, 4, 4, 0, 4, 0, 0, 2, 4, 0, 2, 0, 0,-2, 0, 0,-2, &  ! r13
					 0, 0, 0, 0,-2,-2,-2, 2, 0, 0, 0, 0,-2,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2, 2, 2, 0, 0, 0, 0,-2, 0 /   ! r23

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)

	! Combine powers of like r's from phi_i, phi_j and the powers in dissertation.
	do n = 1, 36, 1
		rPowerTable(n,1) = ki + kj + rPowers(n,1)
		rPowerTable(n,2) = li + lj + rPowers(n,2)
		rPowerTable(n,3) = mi + mj + rPowers(n,3)
		rPowerTable(n,4) = ni + nj + rPowers(n,4)
		rPowerTable(n,5) = pi + pj + rPowers(n,5)
		rPowerTable(n,6) = qi + qj + rPowers(n,6)
	enddo
	
	return
end


! Called by GenRPowerTablePh21Phi1Exchange - The array here gets multipled by the terms of cos12 in GenRPowerTablePhi2Phi1Exchange.
subroutine SubGenRPowerTablePhi2Phi1Exchange(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(34,6) :: rPowerTable
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms, n
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(34,6) :: rPowers

	data rPowers /	0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 0,-1,-2, 2, 2, 1, 0,-1,-2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, &  !r1
					0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 0, 0, 0, &  !r2
					0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0,-2,-2,-2,-2,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2,-2, &  !r12
					0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2,-1,-2, 1, 0,-1,-2, 2, 2, 0, 0, 0, &  !r3
					0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2,-2,-2, 0, 0, 0, 0, 0, 0,-2,-2, 2, &  !r13
					0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2,-2,-2,-2, 2,-2 /   !r23

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)

	! Combine powers of like r's from phi_i, phi_j and the powers in dissertation.
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


subroutine GenRPowerTablePhi2Phi1ExchangeTrig(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(36,6) :: rPowerTable
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms, n
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(36,6) :: rPowers

	data rPowers /	-2, 0,-4, 2, 4, 0,-2, 0, 2,-2,-4, 2, 0,-2,-4,-2,-2,-4,-2,-2, 2,-2, 0,-2,-2,-2, 0,-2,-4, 0,-4,-4,-2, 2, 0,-2, &  ! r1
					 0, 0, 0,-4,-4,-4,-4,-2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 4, 0, 0,-2,-2,-2,-2, 2, 0,-2,-2, 0,-2,-2, 2, 0,-4,-4,-4, &  ! r2
					 0,-2, 2, 0,-2, 2, 4, 0,-2, 2, 4, 0, 2, 4, 0,-2, 0,-2, 2, 0, 0, 4, 0, 2, 0, 0, 0, 2, 2, 0, 4, 0, 0, 0, 2, 4, &  ! r12
					 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, &  ! r3
					 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2, 0, 0,-2, 0, 0, 2, 0, 0, 2, 2, 0,-2,-2,-2,-2,-2,-2,-2, 0, 0, 0, 0, &  ! r13
					 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 0, 0, 0, 0,-2,-2,-2,-2 /   ! r23

	call AssignIndices(i, j, ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj, PowerTablei, PowerTablej, NumTerms)

	! Combine powers of like r's from phi_i, phi_j and the powers in dissertation.
	do n = 1, 36, 1
		rPowerTable(n,1) = ki + kj + rPowers(n,1)
		rPowerTable(n,2) = li + lj + rPowers(n,2)
		rPowerTable(n,3) = mi + mj + rPowers(n,3)
		rPowerTable(n,4) = ni + nj + rPowers(n,4)
		rPowerTable(n,5) = pi + pj + rPowers(n,5)
		rPowerTable(n,6) = qi + qj + rPowers(n,6)
	enddo
	
	return
end


!subroutine GenCoeffTablePhi2Phi2Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
!						 Gammai, Gammaj)
!	implicit none
!	real*16, dimension(240) :: CoeffTable
!	real*16, dimension(36) :: CoeffTable36
!	real*16, dimension(34) :: CoeffTable34
!	integer i, j, k, n, NumTerms
!	integer, dimension(NumTerms,6) :: PowerTablei
!	integer, dimension(NumTerms,6) :: PowerTablej
!	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
!	real*16, dimension(0:5) :: Cos23Coeff
!
!	!data Cos23Coeff / 0.25q0, -0.75q0, 0.375q0, -0.75q0, 0.375q0, 0.375q0 /
!	data Cos23Coeff / 0.25q0, 0.375q0, 0.375q0, -0.75q0, -0.75q0, 0.375q0 /
!
!	do k = 0, 5, 1
!		call SubGenCoeffTablePhi2Phi2Exchange(CoeffTable34, PowerTablei, i, PowerTablej, j, NumTerms, &
!							Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
!		do n = 1, 34, 1
!			CoeffTable(k*34+n) = CoeffTable34(n) * Cos23Coeff(k)
!		end do
!	end do
!
!	call GenCoeffTablePhi2Phi2ExchangeTrig(CoeffTable36, PowerTablei, i, PowerTablej, j, NumTerms, &
!						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
!	do n = 1, 36, 1
!		CoeffTable(204+n) = CoeffTable36(n)
!	end do
!
!	return
!end subroutine


subroutine GenCoeffTablePhi2Phi2Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(240) :: CoeffTable
	real*16, dimension(36) :: CoeffTable36
	real*16, dimension(0:33) :: CoeffTable34
	integer i, j, k, n, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	real*16, dimension(6) :: Cos23Coeff

	!data Cos23Coeff / 0.25q0, -0.75q0, 0.375q0, -0.75q0, 0.375q0, 0.375q0 /
	data Cos23Coeff / 0.25q0, 0.375q0, 0.375q0, -0.75q0, -0.75q0, 0.375q0 /

	do n = 0, 33, 1
		call SubGenCoeffTablePhi2Phi2Exchange(CoeffTable34, PowerTablei, i, PowerTablej, j, NumTerms, &
							Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
		do k = 1, 6, 1
			CoeffTable(n*6+k) = CoeffTable34(n) * Cos23Coeff(k)
		end do
	end do

	call GenCoeffTablePhi2Phi2ExchangeTrig(CoeffTable36, PowerTablei, i, PowerTablej, j, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
	do n = 1, 36, 1
		CoeffTable(204+n) = CoeffTable36(n)
	end do

	return
end subroutine


subroutine GenCoeffTablePhi1Phi2Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(240) :: CoeffTable
	real*16, dimension(36) :: CoeffTable36
	real*16, dimension(34) :: CoeffTable34
	integer i, j, k, n, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	real*16, dimension(0:5) :: Cos13Coeff

	data Cos13Coeff / 0.25q0, -0.75q0, 0.375q0, -0.75q0, 0.375q0, 0.375q0 /

	do k = 0, 5, 1
		call SubGenCoeffTablePhi1Phi2Exchange(CoeffTable34, PowerTablei, i, PowerTablej, j, NumTerms, &
							Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
		do n = 1, 34, 1
			CoeffTable(k*34+n) = CoeffTable34(n) * Cos13Coeff(k)
		end do
	end do

	call GenCoeffTablePhi1Phi2ExchangeTrig(CoeffTable36, PowerTablei, i, PowerTablej, j, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
	do n = 1, 36, 1
		CoeffTable(204+n) = CoeffTable36(n)
	end do

	return
end subroutine


subroutine GenCoeffTablePhi2Phi1Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, Alphai, Alphaj, Betai, Betaj, &
						 Gammai, Gammaj)
	implicit none
	real*16, dimension(240) :: CoeffTable
	real*16, dimension(36) :: CoeffTable36
	real*16, dimension(34) :: CoeffTable34
	integer i, j, k, n, NumTerms
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	real*16, dimension(0:5) :: Cos12Coeff

	data Cos12Coeff / 0.25q0, -0.75q0, 0.375q0, -0.75q0, 0.375q0, 0.375q0 /

	do k = 0, 5, 1
		call SubGenCoeffTablePhi2Phi1Exchange(CoeffTable34, PowerTablei, i, PowerTablej, j, NumTerms, &
							Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
		do n = 1, 34, 1
			CoeffTable(k*34+n) = CoeffTable34(n) * Cos12Coeff(k)
		end do
	end do

	call GenCoeffTablePhi2Phi1ExchangeTrig(CoeffTable36, PowerTablei, i, PowerTablej, j, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
	do n = 1, 36, 1
		CoeffTable(204+n) = CoeffTable36(n)
	end do

	return
end subroutine


!subroutine GenRPowerTablePhi2Phi2Exchange(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
!	integer, dimension(240,6) :: rPowerTable
!	integer, dimension(34,6) :: rPowerTable34
!	integer, dimension(36,6) :: rPowerTable36
!	integer, dimension(NumTerms,6) :: PowerTablei
!	integer, dimension(NumTerms,6) :: PowerTablej
!	integer i, j, NumTerms, n
!	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
!	integer, dimension(0:5,6) :: Cos23Power
!
!	!data Cos23Power / 0, 0, 0, 0, 0, 0, &  ! r1
!	!				  0,-2, 2, 0,-2,-2, &  ! r2
!	!				  0, 0, 0, 0, 0, 0, &  ! r12
!	!				  0, 0,-2,-2,-2, 2, &  ! r3
!	!				  0, 0, 0, 0, 0, 0, &  ! r13
!	!				  0, 2, 0, 2, 4, 0 /   ! r23
!	
!	data Cos23Power / 0, 0, 0, 0, 0, 0, &  ! r1
!					  0, 2,-2, 0,-2,-2, &  ! r2
!					  0, 0, 0, 0, 0, 0, &  ! r12
!					  0,-2, 2,-2, 0,-2, &  ! r3
!					  0, 0, 0, 0, 0, 0, &  ! r13
!					  0, 0, 0, 2, 2, 4 /   ! r23
!
!	do k = 0, 5, 1
!		call GenRPowerTablePhi2Phi2(rPowerTable34, PowerTablei, i, PowerTablej, j, NumTerms)
!
!		do n = 1, 34, 1
!			do m = 1, 6, 1
!				rPowerTable(k*34+n,m) = rPowerTable34(n,m) + Cos23Power(k,m)
!			end do
!		end do
!	end do
!
!	call GenRPowerTablePhi2Phi2ExchangeTrig(rPowerTable36, PowerTablei, i, PowerTablej, j, NumTerms)
!	do n = 1, 36, 1
!		do m = 1, 6, 1
!			rPowerTable(204+n,m) = rPowerTable36(n,m)
!		enddo
!	end do
!
!	return
!end subroutine


subroutine GenRPowerTablePhi2Phi2Exchange(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(240,6) :: rPowerTable
	integer, dimension(0:33,6) :: rPowerTable34
	integer, dimension(36,6) :: rPowerTable36
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms, n
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(6,6) :: Cos23Power

	!data Cos23Power / 0, 0, 0, 0, 0, 0, &  ! r1
	!				  0,-2, 2, 0,-2,-2, &  ! r2
	!				  0, 0, 0, 0, 0, 0, &  ! r12
	!				  0, 0,-2,-2,-2, 2, &  ! r3
	!				  0, 0, 0, 0, 0, 0, &  ! r13
	!				  0, 2, 0, 2, 4, 0 /   ! r23
	
	data Cos23Power / 0, 0, 0, 0, 0, 0, &  ! r1
					  0, 2,-2, 0,-2,-2, &  ! r2
					  0, 0, 0, 0, 0, 0, &  ! r12
					  0,-2, 2,-2, 0,-2, &  ! r3
					  0, 0, 0, 0, 0, 0, &  ! r13
					  0, 0, 0, 2, 2, 4 /   ! r23

	do n = 0, 33, 1
		do k = 1, 6, 1
		call GenRPowerTablePhi2Phi2(rPowerTable34, PowerTablei, i, PowerTablej, j, NumTerms)
			do m = 1, 6, 1
				rPowerTable(n*6+k,m) = rPowerTable34(n,m) + Cos23Power(k,m)
			end do
		end do
	end do

	call GenRPowerTablePhi2Phi2ExchangeTrig(rPowerTable36, PowerTablei, i, PowerTablej, j, NumTerms)
	do n = 1, 36, 1
		do m = 1, 6, 1
			rPowerTable(204+n,m) = rPowerTable36(n,m)
		enddo
	end do

	return
end subroutine


subroutine GenRPowerTablePhi1Phi2Exchange(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(240,6) :: rPowerTable
	integer, dimension(34,6) :: rPowerTable34
	integer, dimension(36,6) :: rPowerTable36
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms, n
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(0:5,6) :: Cos13Power

	data Cos13Power / 0,-2, 2, 0,-2,-2, &  ! r1
					  0, 0, 0, 0, 0, 0, &  ! r2
					  0, 0, 0, 0, 0, 0, &  ! r12
					  0, 0,-2,-2,-2, 2, &  ! r3
					  0, 2, 0, 2, 4, 0, &  ! r13
					  0, 0, 0, 0, 0, 0 /   ! r23

	do k = 0, 5, 1
		call SubGenRPowerTablePhi1Phi2Exchange(rPowerTable34, PowerTablei, i, PowerTablej, j, NumTerms)

		do n = 1, 34, 1
			do m = 1, 6, 1
				rPowerTable(k*34+n,m) = rPowerTable34(n,m) + Cos13Power(k,m)
			end do
		end do
	end do

	call GenRPowerTablePhi1Phi2ExchangeTrig(rPowerTable36, PowerTablei, i, PowerTablej, j, NumTerms)
	do n = 1, 36, 1
		do m = 1, 6, 1
			rPowerTable(204+n,m) = rPowerTable36(n,m)
		enddo
	end do

	return
end subroutine


subroutine GenRPowerTablePhi2Phi1Exchange(rPowerTable, PowerTablei, i, PowerTablej, j, NumTerms)
	integer, dimension(240,6) :: rPowerTable
	integer, dimension(34,6) :: rPowerTable34
	integer, dimension(36,6) :: rPowerTable36
	integer, dimension(NumTerms,6) :: PowerTablei
	integer, dimension(NumTerms,6) :: PowerTablej
	integer i, j, NumTerms, n
	integer ki, li, mi, ni, pi, qi, kj, lj, mj, nj, pj, qj
	integer, dimension(0:5,6) :: Cos12Power

	data Cos12Power / 0,-2, 2, 0,-2,-2, &  ! r1
					  0, 0,-2,-2,-2, 2, &  ! r2
					  0, 2, 0, 2, 4, 0, &  ! r12
					  0, 0, 0, 0, 0, 0, &  ! r3
					  0, 0, 0, 0, 0, 0, &  ! r13
					  0, 0, 0, 0, 0, 0 /   ! r23

	do k = 0, 5, 1
		call SubGenRPowerTablePhi2Phi1Exchange(rPowerTable34, PowerTablei, i, PowerTablej, j, NumTerms)

		do n = 1, 34, 1
			do m = 1, 6, 1
				rPowerTable(k*34+n,m) = rPowerTable34(n,m) + Cos12Power(k,m)
			end do
		end do
	end do

	call GenRPowerTablePhi2Phi1ExchangeTrig(rPowerTable36, PowerTablei, i, PowerTablej, j, NumTerms)
	do n = 1, 36, 1
		do m = 1, 6, 1
			rPowerTable(204+n,m) = rPowerTable36(n,m)
		enddo
	end do

	return
end subroutine


! This subroutine calculates the direct-exchange part of the matrix elements.
subroutine CalcExchangeMatrices(RunCalc, PowerTablei, PowerTablej, WMatrix, CMatrix, PhiPhi, Phi2HPhi, NumTerms, &
						Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj, Omega, pmax, qmax, UseMinus, Method)
	use WCLimits
	implicit none
	logical RunCalc
	real*16, dimension(lmin:lmax, mmin:mmax, nmin:nmax, 6) :: WMatrix
	real*16, dimension(-2:Cjmax, 0:Cqmax, 0:Ckmax) :: CMatrix
	real*16, dimension(240) :: CoeffTable
	integer, dimension(34,6) :: rPowers
	!integer, dimension(171,6) :: rPowers171
	integer, dimension(240,6) :: rPowers240
	integer, dimension(NumTerms,6) :: PowerTablei, PowerTablej
	real*16, dimension(NumTerms,NumTerms) :: PhiPhi, Phi2HPhi
	integer Omega, pmax, qmax, NumTerms, Method
	real*16 Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj
	real*16 HylleraasIntegral
	real*16 PhiHPhiSum, PhiPhiSum
	real*16 PI
	integer i, j, n, UseMinus, Offset
	real*16 RemoveMe
	real*16 PiHalf  ! 0.5 / pi
	real*16, dimension(0:5) :: Cos12Coeff, Cos13Coeff, Cos23Coeff
	integer, dimension(0:5,6) :: Cos12Power
	integer, dimension(0:5,6) :: Cos13Power
	integer, dimension(0:5,6) :: Cos23Power

	data Cos12Coeff / 0.25q0, -0.75q0, 0.375q0, -0.75q0, 0.375q0, 0.375q0 /
	data Cos13Coeff / 0.25q0, -0.75q0, 0.375q0, -0.75q0, 0.375q0, 0.375q0 /
	!data Cos23Coeff / 0.25q0, -0.75q0, 0.375q0, -0.75q0, 0.375q0, 0.375q0 /
	data Cos23Coeff / 0.25q0, 0.375q0, 0.375q0, -0.75q0, -0.75q0, 0.375q0 /
	data Cos12Power / 0,-2, 2, 0,-2,-2, &  ! r1
					  0, 0,-2,-2,-2, 2, &  ! r2
					  0, 2, 0, 2, 4, 0, &  ! r12
					  0, 0, 0, 0, 0, 0, &  ! r3
					  0, 0, 0, 0, 0, 0, &  ! r13
					  0, 0, 0, 0, 0, 0 /   ! r23
	data Cos13Power / 0,-2, 2, 0,-2,-2, &  ! r1
					  0, 0, 0, 0, 0, 0, &  ! r2
					  0, 0, 0, 0, 0, 0, &  ! r12
					  0, 0,-2,-2,-2, 2, &  ! r3
					  0, 2, 0, 2, 4, 0, &  ! r13
					  0, 0, 0, 0, 0, 0 /   ! r23
	data Cos23Power / 0, 0, 0, 0, 0, 0, &  ! r1
					  0, 2,-2, 0,-2,-2, &  ! r2
					  0, 0, 0, 0, 0, 0, &  ! r12
					  0,-2, 2,-2, 0,-2, &  ! r3
					  0, 0, 0, 0, 0, 0, &  ! r13
					  0, 0, 0, 2, 2, 4 /   ! r23

	PiHalf = 0.1591549430918953357688837633725143620345q0  ! Calculated in Mathematica
	Offset = NumTerms / 2  ! NumTerms should be divisible by 2 always.

	!if (Method == 0 .or. Method == 1) then
	!	write (*,*) "Precomputing W matrix..."
	!	! We have to recalculate the W matrices, since we had the swaps.  The limits are still the same as before
	!	!  though, since the nonlinear parameters do not affect the limits.
	!	!allocate(WMatrix(lmin:lmax, mmin:mmax, nmin:nmax, 6))
	!	call CalcWMatrices(Omega+2, qmax, WMatrix, Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, pmax)
	!	write (*,*) "Finished precomputing W matrix"
	!else  ! Method == 2
	!	write (*,*) "Precomputing gamma functions..."
	!	call PreCalcGamma(real(Alphai+Alphaj,16), real(Betai+Betaj,16), real(Gammai+Gammaj,16), 25)
	!	write (*,*) "Finished precomputing gamma functions"
	!end if

	! Calculate the phi1-phi1 elements
	write (*,*) "Starting calculation of phi1-phi1 direct-exchange matrix elements"
	!$omp parallel do shared(NumTerms,PiHalf,Offset,PowerTablei,PowerTablej) private(j,n,PhiHPhiSum,PhiPhiSum,RemoveMe,CoeffTable,rPowers) schedule(dynamic,10)
	do i = 1, Offset, 1
		!write (*,"(i7)",advance='no') i
		write (*,*) i
		do j = 1, Offset, 1
			PhiHPhiSum = 0.0q0
			call GenCoeffTablePhi1Phi1(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, &
								Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
			call GenRPowerTablePhi1Phi1(rPowers, PowerTablei, i, PowerTablej, j, NumTerms)

			do n = 1, 34, 1
				if (CoeffTable(n) /= 0.0q0) then
					! Note that the order of the parameters to the HylleraasIntegral function looks strange.  This is
					!  because the Drake and Yan paper (1995) expressed the integrand in the following order:
					!  r1 r2 r3 r12 r23 r31, and the Armour and Humberston article has them in this order:
					!  r1 r2 r12 r3 r13 r23 (equation 3.14).  So we have to swap the third and fourth parameters, along
					!  with the fifth and sixth.
					
					RemoveMe = HylleraasIntegral(RunCalc, rPowers(n,1), rPowers(n,2), rPowers(n,4), rPowers(n,3), rPowers(n,6), &
												  rPowers(n,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, Method, WMatrix, CMatrix)
					PhiHPhiSum = PhiHPhiSum + RemoveMe * CoeffTable(n)
				endif
			enddo

			! Do the PhiPhi inner product
			PhiPhiSum = HylleraasIntegral(RunCalc, PowerTablei(i,1)+PowerTablej(j,1), PowerTablei(i,2)+PowerTablej(j,2), &
						PowerTablei(i,4)+PowerTablej(j,4), PowerTablei(i,3)+PowerTablej(j,3), &
						PowerTablei(i,6)+PowerTablej(j,6), PowerTablei(i,5)+PowerTablej(j,5), &
						Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, Method, WMatrix, CMatrix)

			if (UseMinus == 0) then
				Phi2HPhi(i,j) = Phi2HPhi(i,j) + PhiHPhiSum * PiHalf
				PhiPhi(i,j) = PhiPhi(i,j) + PhiPhiSum * PiHalf
			else
				Phi2HPhi(i,j) = Phi2HPhi(i,j) - PhiHPhiSum * PiHalf
				PhiPhi(i,j) = PhiPhi(i,j) - PhiPhiSum * PiHalf
			endif
		enddo
	enddo
	
	! Calculate the phi2-phi2 elements
	write (*,*) "Starting calculation of phi2-phi2 direct-exchange matrix elements"
	!$omp parallel do shared(NumTerms,PiHalf,Offset,PowerTablei,PowerTablej) private(j,n,PhiHPhiSum,PhiPhiSum,RemoveMe,CoeffTable,rPowers240) schedule(dynamic,10)
	do i = Offset+1, NumTerms, 1
		!write (*,"(i7)",advance='no') i
		write (*,*) i
		do j = Offset+1, NumTerms, 1
			PhiHPhiSum = 0.0q0
 
			call GenCoeffTablePhi2Phi2Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, &
								Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
			call GenRPowerTablePhi2Phi2Exchange(rPowers240, PowerTablei, i, PowerTablej, j, NumTerms)
 
			do n = 1, 240, 1
			!do n = 1, 204, 1
			!do n = 205, 240, 1
				if (CoeffTable(n) /= 0.0q0) then
					RemoveMe = HylleraasIntegral(RunCalc, rPowers240(n,1), rPowers240(n,2), rPowers240(n,4), rPowers240(n,3), rPowers240(n,6), &
												 rPowers240(n,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, Method, WMatrix, CMatrix)
					PhiHPhiSum = PhiHPhiSum + RemoveMe * CoeffTable(n)
				endif
			enddo
 
			! Do the PhiPhi inner product
			PhiPhiSum = 0.0q0
			do n = 0, 5, 1
				RemoveMe = Cos23Coeff(n) * &
							HylleraasIntegral(RunCalc, PowerTablei(i,1)+PowerTablej(j,1), PowerTablei(i,2)+PowerTablej(j,2)+Cos23Power(n,2), &
							PowerTablei(i,4)+PowerTablej(j,4)+Cos23Power(n,4), PowerTablei(i,3)+PowerTablej(j,3), &
							PowerTablei(i,6)+PowerTablej(j,6)+Cos23Power(n,6), PowerTablei(i,5)+PowerTablej(j,5), &
							Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, Method, WMatrix, CMatrix)
				PhiPhiSum = PhiPhiSum + RemoveMe
			end do
 
			if (UseMinus == 0) then
				Phi2HPhi(i,j) = Phi2HPhi(i,j) + PhiHPhiSum * PiHalf
				PhiPhi(i,j) = PhiPhi(i,j) + PhiPhiSum * PiHalf
			else
				Phi2HPhi(i,j) = Phi2HPhi(i,j) - PhiHPhiSum * PiHalf
				PhiPhi(i,j) = PhiPhi(i,j) - PhiPhiSum * PiHalf
			endif
		enddo
	enddo	


	! Calculate the phi1-phi2 elements	
	write (*,*) "Starting calculation of phi1-phi2 direct-exchange matrix elements"
	!$omp parallel do shared(NumTerms,PiHalf,Offset,PowerTablei,PowerTablej) private(j,n,PhiHPhiSum,PhiPhiSum,RemoveMe,CoeffTable,rPowers240) schedule(dynamic,10)
	do i = 1, Offset, 1
		!write (*,"(i7)",advance='no') i
		write (*,*) i
		!do j = 1, i, 1
		do j = Offset+1, NumTerms, 1
			PhiHPhiSum = 0.0q0

			call GenCoeffTablePhi1Phi2Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, &
								Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
			call GenRPowerTablePhi1Phi2Exchange(rPowers240, PowerTablei, i, PowerTablej, j, NumTerms)

			do n = 1, 240, 1
				if (CoeffTable(n) /= 0.0q0) then
					RemoveMe = HylleraasIntegral(RunCalc, rPowers240(n,1), rPowers240(n,2), rPowers240(n,4), rPowers240(n,3), rPowers240(n,6), &
												rPowers240(n,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, Method, WMatrix, CMatrix)
					PhiHPhiSum = PhiHPhiSum + RemoveMe * CoeffTable(n)
				endif
			enddo

			! Do the PhiPhi inner product
			PhiPhiSum = 0.0q0
			do n = 0, 5, 1
				RemoveMe = Cos13Coeff(n) * &
							HylleraasIntegral(RunCalc, PowerTablei(i,1)+PowerTablej(j,1)+Cos13Power(n,1), PowerTablei(i,2)+PowerTablej(j,2), &
							PowerTablei(i,4)+PowerTablej(j,4)+Cos13Power(n,4), PowerTablei(i,3)+PowerTablej(j,3), &
							PowerTablei(i,6)+PowerTablej(j,6), PowerTablei(i,5)+PowerTablej(j,5)+Cos13Power(n,5), &
							Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, Method, WMatrix, CMatrix)
				PhiPhiSum = PhiPhiSum + RemoveMe
			end do

			if (UseMinus == 0) then
				Phi2HPhi(i,j) = Phi2HPhi(i,j) + PhiHPhiSum * PiHalf
				PhiPhi(i,j) = PhiPhi(i,j) + PhiPhiSum * PiHalf
			else
				Phi2HPhi(i,j) = Phi2HPhi(i,j) - PhiHPhiSum * PiHalf
				PhiPhi(i,j) = PhiPhi(i,j) - PhiPhiSum * PiHalf
			endif
		enddo
	enddo


	! Calculate the phi2-phi1 elements	
	write (*,*) "Starting calculation of phi2-phi1 direct-exchange matrix elements"
	!$omp parallel do shared(NumTerms,PiHalf,Offset,PowerTablei,PowerTablej) private(j,n,PhiHPhiSum,PhiPhiSum,RemoveMe,CoeffTable,rPowers240) schedule(dynamic,10)
	do i = Offset+1, NumTerms, 1
		!write (*,"(i7)",advance='no') i
		write (*,*) i
		do j = 1, Offset, 1
			PhiHPhiSum = 0.0q0

			call GenCoeffTablePhi2Phi1Exchange(CoeffTable, PowerTablei, i, PowerTablej, j, NumTerms, &
								Alphai, Alphaj, Betai, Betaj, Gammai, Gammaj)
			call GenRPowerTablePhi2Phi1Exchange(rPowers240, PowerTablei, i, PowerTablej, j, NumTerms)

			do n = 1, 240, 1
				if (CoeffTable(n) /= 0.0q0) then
					RemoveMe = HylleraasIntegral(RunCalc, rPowers240(n,1), rPowers240(n,2), rPowers240(n,4), rPowers240(n,3), rPowers240(n,6), &
												rPowers240(n,5), Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, Method, WMatrix, CMatrix)
					PhiHPhiSum = PhiHPhiSum + RemoveMe * CoeffTable(n)
				endif
			enddo

			! Do the PhiPhi inner product
			PhiPhiSum = 0.0q0
			do n = 0, 5, 1
				RemoveMe = Cos12Coeff(n) * &
							HylleraasIntegral(RunCalc, PowerTablei(i,1)+PowerTablej(j,1)+Cos12Power(n,1), PowerTablei(i,2)+PowerTablej(j,2)+Cos12Power(n,2), &
							PowerTablei(i,4)+PowerTablej(j,4), PowerTablei(i,3)+PowerTablej(j,3)+Cos12Power(n,3), &
							PowerTablei(i,6)+PowerTablej(j,6), PowerTablei(i,5)+PowerTablej(j,5), &
							Alphai+Alphaj, Betai+Betaj, Gammai+Gammaj, qmax, Method, WMatrix, CMatrix)
				PhiPhiSum = PhiPhiSum + RemoveMe
			end do

			if (UseMinus == 0) then
				Phi2HPhi(i,j) = Phi2HPhi(i,j) + PhiHPhiSum * PiHalf
				PhiPhi(i,j) = PhiPhi(i,j) + PhiPhiSum * PiHalf
			else
				Phi2HPhi(i,j) = Phi2HPhi(i,j) - PhiHPhiSum * PiHalf
				PhiPhi(i,j) = PhiPhi(i,j) - PhiPhiSum * PiHalf
			endif
		enddo
	enddo
	
	return
end
