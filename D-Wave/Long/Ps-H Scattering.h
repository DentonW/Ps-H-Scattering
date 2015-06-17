#ifndef PSH_SCATTERING_H
#define PSH_SCATTERING_H

#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

//@TODO: Should this be extern?
//double M_PI;

// From http://software.intel.com/en-us/articles/intel-math-kernel-library-intel-mkl-blas-cblas-and-lapack-compilinglinking-functions-fortran-and-cc-calls/
typedef struct
{
	double re;
	double im;
} complex16;

class rPowers
{
	public:
		rPowers() { ki = li = mi = ni = pi = qi = 0; alpha = beta = gamma = 1.0; Index = 0; return; }
		rPowers(double a, double b, double g) { ki = li = mi = ni = pi = qi = 0; alpha = a; beta = b; gamma = g; Index = 0; return; }
		int ki, li, mi, ni, pi, qi;
		double alpha, beta, gamma;
		int Index;
};

typedef struct
{
	int LongLong_r1, LongLong_r2Leg, LongLong_r2Lag, LongLong_r3Leg, LongLong_r3Lag, LongLong_r12, LongLong_r13, LongLong_phi23;
	int LongLongr23_r1, LongLongr23_r2Leg, LongLongr23_r2Lag, LongLongr23_r3Leg, LongLongr23_r3Lag, LongLongr23_phi12, LongLongr23_r13, LongLongr23_r23;

	int ShortLong_r1, ShortLong_r2Leg, ShortLong_r2Lag, ShortLong_r3Leg, ShortLong_r3Lag, ShortLong_r12, ShortLong_r13, ShortLong_phi23;
	int ShortLongr23_r1, ShortLongr23_r2Leg, ShortLongr23_r2Lag, ShortLongr23_r3Leg, ShortLongr23_r3Lag, ShortLongr23_r12, ShortLongr23_phi13, ShortLongr23_r23;
	int ShortLongQiGt0_r1, ShortLongQiGt0_r2Leg, ShortLongQiGt0_r2Lag, ShortLongQiGt0_r3Leg, ShortLongQiGt0_r3Lag, ShortLongQiGt0_r12, ShortLongQiGt0_r13, ShortLongQiGt0_phi23;
} QuadPoints;

template <class T> string to_string(const T& t);
void	ReadParamFile(ifstream &ParameterFile, QuadPoints &q, double &Mu, int &ShPower, double &Lambda1, double &Lambda2, double &Lambda3, double &r2Cusp, double &r3Cusp);
bool	ReadShortHeader(ifstream &FileShortRange, int &Omega, int &IsTriplet, int &Ordering, int &NumShortTerms, double &Alpha, double &Beta, double &Gamma);
void	ShowDateTime(ofstream &OutFile);
string	ShowTime(void);
long double	PsWaveFn(double r12);
long double	HWaveFn(double r3);
void	CalcARowAndBVector(int Node, int NumTermsQi0, int NumTermsQiGt0, int Omega, vector <rPowers> &PowerTableQi0, vector <rPowers> &PowerTableQiGt0, vector <double> &AResultsQi0,
			  vector <double> &AResultsQiGt0, vector <double> &ARow, vector <double> &BResultsQi0, vector <double> &BResultsQiGt0, vector <double> &B, double &SLS, double &SLC, QuadPoints &q, double r2Cusp,
			  double r3Cusp, double alpha, double beta, double gamma, double kappa, double mu, double lambda1, double lambda2, double lambda3, int shpower, int sf);
void	CombineResults(int Omega, int Ordering, vector <double> &ResultsQi0, vector <double> &ResultsQiGt0, vector <double> &Results, int Start, int End);

typedef	double (*FuncPtr)(rPowers &, double, double, double, double, double, double, double, double, double, double, int);

double	ipow(double a, int ex);
bool	IsFiniteNumber(double x);
void	WriteProgress(string Desc, int &Prog, int i, int N);

#ifdef USE_MPI
	#define FinishMPI() MpiError = MPI_Finalize()
#else
	#define FinishMPI()
#endif

// Gaussian Integration.cpp
void	ChangeOfInterval(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, long double a, long double b);
void	ChangeOfIntervalNoResize(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, long double a, long double b);
int		GaussLegendre(vector <long double> &Abscissas, vector <long double> &Weights, int n);
int		GaussLaguerre(vector <long double> &Abscissas, vector <long double> &Weights, int n);

// Vector Gaussian Integration.cpp
void	VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, int NumPowers, vector <rPowers> &Powers, int Omega, double lambda1, double lambda2, double lambda3);
void	VecGaussIntegrationPhi13_PhiLCBar_PhiLSBar_R23Term(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nPhi13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, int NumPowers, vector <rPowers> &Powers, int Omega, double lambda1, double lambda2, double lambda3);
void	VecGaussIntegrationPhi12_PhiLCBar_PhiLSBar_R23Term(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nPhi12, int nR13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, int NumPowers, vector <rPowers> &Powers, int Omega, double lambda1, double lambda2, double lambda3);
void	VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar_Full(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, int NumPowers, vector <rPowers> &Powers, int Omega, double lambda1, double lambda2, double lambda3);

// Short-Range.cpp
double	Phi(rPowers &rp, double r1, double r2, double r3, double r12, double r13, double r23);
double	PhiP23(rPowers &rp, double r1, double r2, double r3, double r12, double r13, double r23);
double	PhiBar(rPowers &rp, double r1, double r2, double r3, double r12, double r13, double r23, int sf);
double	CalcWithPhi(rPowers &Powers, double r1, double r2, double r3, double r12, double r13, double r23, double kappa, double rho, double rhop, double mu, int sf);
double	CalcWithPhiP23(rPowers &Powers, double r1, double r2, double r3, double r12, double r13, double r23, double kappa, double rho, double rhop, double mu, int sf);
double	CalcWithPhiBar(rPowers &Powers, double r1, double r2, double r3, double r12, double r13, double r23, double kappa, double rho, double rhop, double mu, int sf);
int		CalcPowerTableSize(int Omega);
int		CalcPowerTableSizeQi0(int Omega, int Ordering, int Start, int End);
int		CalcPowerTableSizeQiGt0(int Omega, int Ordering, int Start, int End);
void	GenOmegaPowerTable(int Omega, int Ordering, vector<rPowers> &PowerTable);
void	GenOmegaPowerTableQi0(int Omega, int Ordering, vector<rPowers> &PowerTable, int Start, int End);
void	GenOmegaPowerTableQiGt0(int Omega, int Ordering, vector<rPowers> &PowerTable, int Start, int End);

// Long-Range.cpp
long double	sf_bessel_j1(long double x);
long double	sf_bessel_n1(long double x);
long double	sf_bessel_j2(long double x);
long double	sf_bessel_n2(long double x);
long double	fshielding(long double rho, long double mu, int power);
long double	fshielding1(long double rho, long double mu, int power);
long double	fshielding2(long double rho, long double mu, int power);
void	GaussIntegrationPhi23_LongLong(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, double &CLC, double &SLC, double &CLS, double &SLS);
void	GaussIntegrationPhi12_LongLong_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nPhi12, int nR13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, double &CLC, double &SLC, double &CLS, double &SLS);
void	GaussIntegrationPhi13_LongLong_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nPhi13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, double &CLC, double &SLC, double &CLS, double &SLS);

// Phase Shift.cpp
double	Kohn(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS);
double	InverseKohn(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS);
double	LuccheseKohn(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS);
double	ComplexKohnT(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS);


#endif//PSH_SCATTERING_H
