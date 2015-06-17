//
// Phase Shift.cpp: This file calculates the phase shifts using several different methods.
//

#include <math.h>
#include <iostream>
#include "Ps-H Scattering.h"
#ifndef NO_MPI
	#include <mpi.h>
#endif
#include <mkl_lapack.h>
#include <complex>

typedef complex<double> dcmplx;

// Uses the LAPACK dgesv function to solve our system of linear equations to find the Kohn phaseshift.
//  An example of the usage of dgesv is available at http://software.intel.com/sites/products/documentation/hpc/mkl/lapack/mkl_lapack_examples/dgesv.htm.
double Kohn(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS)
{
	MKL_INT n, nrhs, lda, ldb, info;
	vector <double> A((NumShortTerms+1)*(NumShortTerms+1)), X(NumShortTerms+1);

	// Fill in first row and column of A
	for (int i = 0; i < NumShortTerms+1; i++) {
		A[i] = ARow[i];
		A[i*(NumShortTerms+1)] = ARow[i];
	}

	// Copy short-range terms to bottom-right NumShortTerms x NumShortTerms submatrix of A.
	for (int i = 0; i < NumShortTerms; i++) {
		for (int j = 0; j < NumShortTerms; j++) {
			A[(i+1)*(NumShortTerms+1) + (j+1)] = ShortTerms[i*NumShortTerms + j];
		}
	}

	// Operate on X, as LAPACK overwrites original.
	X = B;

	vector <int> ipiv(NumShortTerms+1);
	n = lda = ldb = NumShortTerms+1;
	nrhs = 1;
	dgesv(&n, &nrhs, &A[0], &lda, &ipiv[0], &X[0], &ldb, &info);
	if (info != 0) {
		cout << "LAPACK Error: " << info << endl;
		return 0.0;
	}

	double PsiLS = -X[0] * B[0];
	for (int i = 1; i < (NumShortTerms+1); i++) {
		PsiLS += -X[i] * B[i];
	}
	PsiLS = -(PsiLS + SLS);

	double PhaseShift = atan(PsiLS);
	return PhaseShift;
}


// Uses the LAPACK dgesv function to solve our system of linear equations to find the inverse Kohn (Rubinow) phaseshift.
double InverseKohn(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS)
{
	MKL_INT n, nrhs, lda, ldb, info;
	vector <double> A((NumShortTerms+1)*(NumShortTerms+1)), X(NumShortTerms+1), XCopy(NumShortTerms+1);
	double CLC = ARow[0];

	// Copy short-range terms to bottom-right NumShortTerms x NumShortTerms submatrix of A.
	for (int i = 0; i < NumShortTerms; i++) {
		for (int j = 0; j < NumShortTerms; j++) {
			A[(i+1)*(NumShortTerms+1) + (j+1)] = ShortTerms[i*NumShortTerms + j];
		}
	}

	// Swap the (phi,LS) terms with (phi,LC) terms
	X = ARow;
	XCopy = ARow;
	X[0] = B[0] + 1.0;  // Use (S,LC) = (C,LS) + 1
	XCopy[0] = X[0];

	for (int i = 1; i < NumShortTerms+1; i++) {
		A[i] = B[i];  // First row of A now has (S,Lphi) terms.
		A[i*(NumShortTerms+1)] = B[i];  // Put (S,Lphi) = (phi,LS) terms in first column.
	}
	A[0] = SLS;

	// LAPACK requires calls by reference, so we have to define all these variables.
	vector <int> ipiv(NumShortTerms+1);
	n = lda = ldb = NumShortTerms+1;
	nrhs = 1;
	dgesv(&n, &nrhs, &A[0], &lda, &ipiv[0], &X[0], &ldb, &info);
	if (info != 0) {
		cout << "LAPACK Error: " << info << endl;
		return 0.0;
	}

	double PsiLS = -X[0] * XCopy[0];
	for (int i = 1; i < (NumShortTerms+1); i++) {
		PsiLS += -X[i] * XCopy[i];
	}
	PsiLS = (PsiLS + CLC);

	double PhaseShift = atan(1.0/PsiLS);
	return PhaseShift;
}


// Uses the LAPACK dgesv function to solve our system of linear equations to find the complex Kohn phaseshift (T-matrix).
double ComplexKohnT(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS)
{
	MKL_INT n, nrhs, lda, ldb, info;
	//double *A = new double[(NumShortTerms+1)*(NumShortTerms+1)];
	vector <dcmplx> A((NumShortTerms+1)*(NumShortTerms+1));
	vector <dcmplx> X(NumShortTerms+1);
	double CLC = ARow[0], CLS = B[0];
	double SLC = CLS + 1.0;  // Use (S,LC) = (C,LS) + 1

	// Copy short-range terms to bottom-right NumShortTerms x NumShortTerms submatrix of A.
	for (int i = 0; i < NumShortTerms; i++) {
		for (int j = 0; j < NumShortTerms; j++) {
			A[(i+1)*(NumShortTerms+1) + (j+1)] = ShortTerms[i*NumShortTerms + j];
		}
	}

	// Fill in the rest of A
	A[0] = dcmplx(CLC - SLS, CLS + SLC);
	for (int i = 1; i < NumShortTerms+1; i++) {
		A[i] = dcmplx(ARow[i], B[i]);
		A[i*(NumShortTerms+1)] = A[i];
	}

	// Fill in B (or X)
	X[0] = dcmplx(CLS, SLS);
	for (int i = 1; i < NumShortTerms+1; i++) {
		X[i] = B[i];
	}
	vector <dcmplx> BComp(NumShortTerms+1);
	BComp = X;

	// LAPACK requires calls by reference, so we have to define all these variables.
	vector <int> ipiv(NumShortTerms+1);
	n = lda = ldb = NumShortTerms+1;
	nrhs = 1;
	zgesv(&n, &nrhs, (MKL_Complex16*)&A[0], &lda, &ipiv[0], (MKL_Complex16*)&X[0], &ldb, &info);
	if (info != 0) {
		cout << "LAPACK Error: " << info << endl;
		return 0.0;
	}

	// Equation () of notes
	dcmplx PsiLS = -X[0] * BComp[0];
	for (int i = 1; i < (NumShortTerms+1); i++) {
		PsiLS += -X[i] * BComp[i];
	}
	dcmplx T = -(PsiLS + SLS);
	// Go from T matrix element to K.
	dcmplx K = T / (1.0 + dcmplx(0,1) * T);

	//@TODO: Check K for imaginary part.
	double PhaseShift = atan(K.real());
	return PhaseShift;
}

