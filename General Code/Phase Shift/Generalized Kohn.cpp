// This calculates the generalized Kohn method of J. N. Cooper et al 2010 J. Phys. A: Math. Theor. 43 175302.

#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string.h>
#include <mkl_lapack.h>
#include <complex>

using namespace std;

typedef complex<double> dcmplx;
void	FixPhase(double &PhaseShift, int &LValue, int &IsTriplet);


// ARow and BVec are sent from the main program as the A and B from the Kohn method.  This function rearranges everything into
//  the matrix equation (7) of their paper and solves.
double GeneralizedKohn(int NumShortTerms, vector <double> &ARow, vector <double> &BVec, vector <double> &ShortTerms, double SLS, int LValue, int IsTriplet, double Tau)
{
	double CLC = ARow[0], CLS = BVec[0];
	double PhaseShift;
	//@TODO: Worth changing these to vectors?
	double *A = new double[(NumShortTerms+1)*(NumShortTerms+1)];
	double *B = new double[NumShortTerms+1];
	double *X = new double[NumShortTerms+1];
	int LWork;
	vector <double> Work(1);
	double CosTau = cos(Tau), SinTau = sin(Tau);
	double SLC = CLS + 1.0;  // Use (S,LC) = (C,LS) + 1

	if (Tau == 4.0 * atan(1.0)) {  // Special case of PI
		CosTau = -1.0;  SinTau = 0.0;
	}
	else if (Tau == 2.0 * atan(1.0)) {  // Special case of PI/2
		CosTau = 0.0;  SinTau = 1.0;
	}

	// Copy short-range - short-range terms into A
	for (int i = 0; i < NumShortTerms; i++) {
		for (int j = 0; j < NumShortTerms; j++) {
			A[(i+1)*(NumShortTerms+1) + (j+1)] = ShortTerms[i*NumShortTerms + j];
		}
	}

	// Fill in the rest of A
	A[0] = SinTau*SinTau*SLS - SinTau*CosTau*SLC - CosTau*SinTau*CLS + CosTau*CosTau*CLC;
	for (int i = 1; i < NumShortTerms+1; i++) {
		A[i] = -SinTau*BVec[i] + CosTau*ARow[i];
		A[i*(NumShortTerms+1)] = A[i];
	}

	// Fill in B vector
	B[0] = -SinTau*CosTau*SLS - SinTau*SinTau*SLC + CosTau*CosTau*CLS + CosTau*SinTau*CLC;
	for (int i = 1; i < NumShortTerms+1; i++) {
		B[i] = CosTau*BVec[i] + SinTau*ARow[i];
	}

	// Copy to X, since LAPACK destroys the original B.
	for (int i = 0; i < NumShortTerms+1; i++) {
		X[i] = B[i];
	}

	// LAPACK requires calls by reference, so we have to define all these variables.
	vector <int> ipiv(NumShortTerms+1);
	int n, lda, ldb, nrhs = 1, info;
	n = lda = ldb = NumShortTerms+1;
	dgesv(&n, &nrhs, A, &lda, &ipiv[0], X, &ldb, &info);
	/*Work.resize(1);
	char Upper = 'U';
	dsysv(&Upper, &n, &nrhs, A, &lda, &ipiv[0], X, &ldb, &Work[0], &LWork, &info);
	LWork = Work[0];
	Work.resize(LWork);
	dsysv(&Upper, &n, &nrhs, A, &lda, &ipiv[0], X, &ldb, &Work[0], &LWork, &info);*/
	if (info != 0) {
		cout << "LAPACK Error: " << info << endl;
		return 0.0;
	}

	// Equation () of notes
	double PsiLS = -X[0] * B[0];
	for (int i = 1; i < (NumShortTerms+1); i++) {
		PsiLS += -X[i] * B[i];
	}
	PsiLS = -(PsiLS + (CosTau*CosTau*SLS + CosTau*SinTau*SLC + CosTau*SinTau*CLS + SinTau*SinTau*CLC));
	PhaseShift = atan(PsiLS) + Tau;
	// Phase shift needs to be in the (-pi/2,pi/2) range, since the phase shifts are always mod pi.
	FixPhase(PhaseShift, LValue, IsTriplet);
	//@TODO: Do we use Pi anywhere else in this program?
	/*double Pi = 4.0 * atan(1.0);
	while (PhaseShift > Pi/2.0) {
		PhaseShift = PhaseShift - Pi;
	}*/

	delete [] A;
	delete [] B;
	delete [] X;

	return PhaseShift;
}





// ARow and BVec are sent from the main program as the A and B from the Kohn method.  This function rearranges everything into
//  the matrix equation (7) of their paper and solves.
double CombinedKohn(dcmplx (&u)[2][2], int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> &ShortTerms, double SLS, int LValue, int IsTriplet)
{
	MKL_INT n, nrhs, lda, ldb, info;
	//double *A = new double[(NumShortTerms+1)*(NumShortTerms+1)];
	vector <dcmplx> A((NumShortTerms+1)*(NumShortTerms+1));
	vector <dcmplx> X(NumShortTerms+1);
	double CLC = ARow[0], CLS = B[0];
	double SLC = CLS + 1.0;  // Use (S,LC) = (C,LS) + 1
	dcmplx SLSt, SLCt, CLSt, CLCt;

	//dcmplx u[2][2] = {{dcmplx(1,0),dcmplx(0,0)}, {dcmplx(0,0),dcmplx(1,0)}};  // Kohn
	//dcmplx u[2][2] = {{dcmplx(0,0),dcmplx(1,0)}, {dcmplx(-1,0),dcmplx(0,0)}};  // Inverse Kohn
	//dcmplx u[2][2] = {{dcmplx(0,1), dcmplx(-1,0)}, {dcmplx(0,1), dcmplx(1,0)}};  // S-matrix
	//dcmplx u[2][2] = {{dcmplx(1,0), dcmplx(0,0)}, {dcmplx(0,1), dcmplx(1,0)}};  // T-matrix
	//dcmplx u[2][2] = {{dcmplx(cos(Tau),0),dcmplx(sin(Tau),0)}, {dcmplx(-sin(Tau),0),dcmplx(cos(Tau),0)}};  // Generalized Kohn
	dcmplx detu = u[0][0]*u[1][1] - u[0][1]*u[1][0];

	SLSt = u[0][0]*u[0][0]*SLS + u[0][0]*u[0][1]*SLC + u[0][1]*u[0][0]*CLS + u[0][1]*u[0][1]*CLC;
	SLCt = u[0][0]*u[1][0]*SLS + u[0][0]*u[1][1]*SLC + u[0][1]*u[1][0]*CLS + u[0][1]*u[1][1]*CLC;
	CLSt = u[1][0]*u[0][0]*SLS + u[1][0]*u[0][1]*SLC + u[1][1]*u[0][0]*CLS + u[1][1]*u[0][1]*CLC;
	CLCt = u[1][0]*u[1][0]*SLS + u[1][0]*u[1][1]*SLC + u[1][1]*u[1][0]*CLS + u[1][1]*u[1][1]*CLC;

	// Copy short-range terms to bottom-right NumShortTerms x NumShortTerms submatrix of A.
	for (int i = 0; i < NumShortTerms; i++) {
		for (int j = 0; j < NumShortTerms; j++) {
			A[(i+1)*(NumShortTerms+1) + (j+1)] = ShortTerms[i*NumShortTerms + j];
		}
	}

	// Fill in the rest of A
	A[0] = CLCt;
	for (int i = 1; i < NumShortTerms+1; i++) {
		//A[i] = dcmplx(ARow[i], B[i]);
		A[i] = u[1][0]*B[i] + u[1][1]*ARow[i];
		A[i*(NumShortTerms+1)] = A[i];
	}

	// Fill in B (or X)
	X[0] = -CLSt;
	for (int i = 1; i < NumShortTerms+1; i++) {
		X[i] = - u[0][0]*B[i] - u[0][1]*ARow[i];
	}

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
	dcmplx PsiLS = X[0] * CLSt;
	for (int i = 1; i < (NumShortTerms+1); i++) {
		PsiLS += X[i] * (u[0][0]*B[i] + u[0][1]*ARow[i]);
	}
	dcmplx L = -(PsiLS + SLSt) / detu;
	// Go from general L matrix element to K.
	dcmplx K = (u[0][1] + u[1][1]*L) / (u[0][0] + u[1][0]*L);

	//@TODO: Check K for imaginary part.
	double PhaseShift = atan(K.real());
	FixPhase(PhaseShift, LValue, IsTriplet);
	return PhaseShift;
}