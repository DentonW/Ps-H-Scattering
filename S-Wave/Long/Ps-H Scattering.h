#ifndef PSH_SCATTERING_H
#define PSH_SCATTERING_H

#include <vector>
#include <fstream>
#include <math.h>
using namespace std;

//#define DENTON_ORDERING

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

		friend ostream &operator<< (ostream &output, const rPowers &r)
		{ 
			output << r.ki << " " << r.li << " " << r.mi << " " << r.ni << " " << r.pi << " " << r.qi << " | " << r.alpha << " " << r.beta << " " << r.gamma;
			return output;            
		}
};

typedef struct
{
	int LongLong_x, LongLong_y, LongLong_z, LongLong_r3Leg, LongLong_r3Lag, LongLong_r13;
	int LongLongr23_r1, LongLongr23_r2Leg, LongLongr23_r2Lag, LongLongr23_r3Leg, LongLongr23_r3Lag, LongLongr23_phi12, LongLongr23_r13, LongLongr23_r23;

	int ShortLong_r1, ShortLong_r2Leg, ShortLong_r2Lag, ShortLong_r3Leg, ShortLong_r3Lag, ShortLong_r12, ShortLong_r13;
	int ShortLongr23_r1, ShortLongr23_r2Leg, ShortLongr23_r2Lag, ShortLongr23_r3Leg, ShortLongr23_r3Lag, ShortLongr23_r12, ShortLongr23_phi13, ShortLongr23_r23;
	int ShortLongQiGt0_r1, ShortLongQiGt0_r2Leg, ShortLongQiGt0_r2Lag, ShortLongQiGt0_r3Leg, ShortLongQiGt0_r3Lag, ShortLongQiGt0_r12, ShortLongQiGt0_r13, ShortLongQiGt0_phi23;
} QuadPoints;

template <class T> string to_string(const T& t);
bool	ReadShortHeader(string FileShort, int &Omega, int &IsTriplet, int &Ordering, int &NumShortTerms, double &Alpha, double &Beta, double &Gamma);
void	ReadParamFile(ifstream &ParameterFile, QuadPoints &q, double &Mu, double &r2Cusp, double &r3Cusp);
void	ShowDateTime(ofstream &OutFile);
string	ShowTime(void);
long double	PsWaveFn(double r12);
long double	HWaveFn(double r3);
void	CalcARow(int Node, int NumTermsQi0, int NumTermsQiGt0, int Omega, vector <rPowers> &PowerTableQi0, vector <rPowers> &PowerTableQiGt0, vector <double> &ResultsQi0,
			  vector <double> &ResultsQiGt0, vector <double> &ARow, QuadPoints &q, double r2Cusp, double r3Cusp, double alpha, double beta, double gamma, double kappa, double mu, int sf);
void	CalcBVector(int Node, int NumTermsQi0, int NumTermsQiGt0, int Omega, vector <rPowers> &PowerTableQi0, vector <rPowers> &PowerTableQiGt0, vector <double> &ResultsQi0,
			  vector <double> &ResultsQiGt0, vector <double> &B, QuadPoints &q, double r2Cusp, double r3Cusp, double alpha, double beta, double gamma, double kappa, double mu, int sf);
void	CalcARowAndBVector(int Node, int NumTermsQi0, int NumTermsQiGt0, int Omega, vector <rPowers> &PowerTableQi0, vector <rPowers> &PowerTableQiGt0, vector <double> &AResultsQi0,
			  vector <double> &AResultsQiGt0, vector <double> &ARow, vector <double> &BResultsQi0, vector <double> &BResultsQiGt0, vector <double> &B, double &SLS, double &SLC, QuadPoints &q, double r2Cusp,
			  double r3Cusp, double alpha, double beta, double gamma, double kappa, double mu, int sf);
void	CombineResults(int Omega, vector <double> &ResultsQi0, vector <double> &ResultsQiGt0, vector <double> &Results, int Start, int End);

typedef	double (*FuncPtr)(rPowers &, double, double, double, double, double, double, double, double, double, double, int);

double	ipow(double a, int ex);
bool	IsFiniteNumber(double x);
void	WriteProgress(string Desc, int &Prog, int N);

#ifdef USE_MPI
	#define FinishMPI() MpiError = MPI_Finalize()
#else
	#define FinishMPI()
#endif

// Gaussian Integration.cpp
void	ChangeOfInterval(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, double a, double b);
void	ChangeOfIntervalNoResize(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, double a, double b);
int		GaussLegendre(vector <long double> &Abscissas, vector <long double> &Weights, int n);
int		GaussLaguerre(vector <long double> &Abscissas, vector <long double> &Weights, int n);
void	GaussIntegrationPhi23Perimetric_LongLong(int nX, int nY, int nZ, int nR3Leg, int nR3Lag, int nR13, double CuspR2, double CuspR3, double kappa, double mu, int sf, double &CLC, double &CLS, double &SLS);
void	GaussIntegrationPhi12_LongLong_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nPhi12, int nR13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int sf, double &CLC, double &CLS, double &SLS);

// Vector Gaussian Integration.cpp
void	VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int sf, int NumPowers, vector <rPowers> &Powers, int Omega);
void	VecGaussIntegrationPhi13_PhiLCBar_PhiLSBar_R23Term(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nPhi13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int sf, int NumPowers, vector <rPowers> &Powers, int Omega);
void	VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar_Full(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int sf, int NumPowers, vector <rPowers> &Powers, int Omega);

// Long-Range.cpp
double	S(double r3, double r12, double kappa, double rho);
double	SP23(double r2, double r13, double kappa, double rhop);
double	C(double r3, double r12, double kappa, double rho, double mu);
double	CP23(double r2, double r13, double kappa, double rhop, double mu);
double	LOnCBar(double r1, double r2, double r3, double r12, double r13, double r23, double kappa, double rho, double rhop, double mu, int sf);

// Short-Range.cpp
int		CalcPowerTableSize(int Omega);
int		CalcPowerTableSizeQi0(int Omega, int Start, int End);
int		CalcPowerTableSizeQiGt0(int Omega, int Start, int End);
void	GenOmegaPowerTable(int Ordering, int Omega, vector<rPowers> &PowerTable);
void	GenOmegaPowerTableQi0(int Ordering, int Omega, vector<rPowers> &PowerTable, int Start, int End);
void	GenOmegaPowerTableQiGt0(int Ordering, int Omega, vector<rPowers> &PowerTable, int Start, int End);

// Phase Shift.cpp
double	Kohn(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS);
double	InverseKohn(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS);
double	LuccheseKohn(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS);
double	ComplexKohnT(int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> ShortTerms, double SLS);


#endif//PSH_SCATTERING_H
