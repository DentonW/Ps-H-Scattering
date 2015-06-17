//
// Long-Range.cpp: This file contains the asymptotic long-range (S and C) functions.
//

#include <iostream>
#include <omp.h>
#include "Ps-H Scattering.h"


//@TODO: In "Ps-H Scattering.h"
using namespace std;
extern long double PI;


// Spherical Bessel function of the first kind
long double sf_bessel_j1(long double x)
{
	return (sinl(x)/x - cosl(x))/x;
};


// Spherical Bessel function of the second kind (spherical Neumann)
long double sf_bessel_n1(long double x)
{
	return -(cosl(x)/x + sinl(x))/x;
};


// Spherical Bessel function of the first kind
long double sf_bessel_j2(long double x)
{
	return (3.0L/(x*x*x) - 1.0L/x)*sinl(x) - 3.0L/(x*x)*cosl(x);
};


// Spherical Bessel function of the second kind (spherical Neumann)
long double sf_bessel_n2(long double x)
{
	return -(3.0L/(x*x*x) - 1.0L/x)*cosl(x) - 3.0L/(x*x)*sinl(x);
};


// // f_sh in C22 and C23
// long double fshielding(long double rho, long double mu)
// {
// 	long double sh = 1.0L - exp(-mu*rho) * (1.0L + mu*rho / 2.0L);
// 	return sh*sh*sh*sh*sh;
// };


// // First derivative of f_sh
// //@TODO: Would this and fshielding2 be better to do 1/ExpMuRho^5 instead of exp(-5.0L*MuRho)?
// long double fshielding1(long double rho, long double mu)
// {
// 	long double MuRho = mu*rho;
// 	long double Part = 2.0L - 2.0L * exp(MuRho) + MuRho;
// 	long double sh = 5.0L/32.0L * exp(-5.0L*MuRho) * mu * (1.0L + MuRho) * Part*Part*Part*Part;
// 	return sh;
// };


// // Second derivative of f_sh
// long double fshielding2(long double rho, long double mu)
// {
// 	long double MuRho = mu*rho;
// 	long double ExpMuRho = exp(MuRho);
// 	long double Part1 = 2.0L - 2.0L * ExpMuRho + MuRho;
// 	long double Part2 = 4.0L - 2.0L * (-5.0L + ExpMuRho) * MuRho + 5.0L * MuRho*MuRho;
// 	long double sh = -5.0L/32.0L * exp(-5.0L*MuRho) * mu*mu * Part1*Part1*Part1 * Part2;
// 	return sh;
// };


// f_sh in C22 and C23
long double fshielding(long double rho, long double mu, int n)
{
	long double sh = 1.0L - exp(-mu*rho) * (1.0L + mu*rho / 2.0L);
	long double shret = sh;
	for (int i = 1; i < n; i++) {
		shret *= sh;
	}
	return shret;
};


// First derivative of f_sh
//@TODO: Would this and fshielding2 be better to do 1/ExpMuRho^5 instead of exp(-5.0L*MuRho)?
long double fshielding1(long double rho, long double mu, int n)
{
	long double MuRho = mu*rho;
	long double Part1 = 2.0L - 2.0L * exp(MuRho) + MuRho;
	long double Part2 = 1.0L - 0.5L * exp(-MuRho) * (2.0L + MuRho);
	long double Part2Mult = Part2;
	for (int i = 1; i < n; i++) {
		Part2Mult *= Part2;
	}
	long double sh = -n * mu * (1.0L + MuRho) * Part2Mult / Part1;
	return sh;
};


// Second derivative of f_sh
long double fshielding2(long double rho, long double mu, int n)
{
	long double MuRho = mu*rho;
	long double ExpMuRho = exp(MuRho);
	long double ExpMuRhoInv = exp(-MuRho);
	long double Part1 = 2.0L - 2.0L * exp(MuRho) + MuRho;
	Part1 *= Part1;
	long double Part2 = 1.0L - 0.5L * ExpMuRhoInv * (2.0L + MuRho);
	long double Part2Mult = Part2;
	for (int i = 1; i < n; i++) {
		Part2Mult *= Part2;
	}
	long double Part3 = -1.0L - 2.0L * ExpMuRho * MuRho + n * (1.0L + MuRho) * (1.0L + MuRho);
	long double sh = n * mu * mu * Part3 * Part2Mult / Part1;
	return sh;
};


void GaussIntegrationPhi23_LongLong(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, double &CLC, double &SLC, double &CLS, double &SLS)
{
	long double r1, r2, r3;
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR12, LegendreWeightsR12;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r3Abscissas, r2Weights, r3Weights;
	vector <long double> r12Array(nR12), r13Array(nR13);
	int NumR2Points, NumR3Points, Prog = 0;
	long double SqrtKappa = sqrt(kappa);

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(r1Abscissas, r1Weights, nR1);
	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR12, LegendreWeightsR12, nR12);
	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);

	long double r1SumCLC = 0.0L, r1SumCLS = 0.0L, r1SumSLC = 0.0L, r1SumSLS = 0.0L;
	#pragma omp parallel for shared(r1SumCLC,r1SumCLS,r1SumSLC,r1SumSLS,r1Abscissas,r1Weights) private(r1,r2,r3,r12Array,r13Array,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		WriteProgress(string("Long-range - long-range"), Prog, i, nR1);

		// These are private, so they need to be initialized.
		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);
		r12Array.resize(nR12);
		r13Array.resize(nR13);

		r1 = r1Abscissas[i];

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR2) {
			for (int m = 0; m < nR2Lag; m++) {
				r2Weights[m] = LaguerreWeightsR2[m];
				r2Abscissas[m] = LaguerreAbscissasR2[m];
			}
			NumR2Points = nR2Lag;
		}
		// Do Gauss-Legendre for r2 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0L, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * exp(-r2Abscissas[m]) * (r1-0.0L)/2.0L;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = LaguerreAbscissasR2[m-nR2Leg] + r1;
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-r1);
			}
		}

		long double r2SumCLC = 0.0L, r2SumCLS = 0.0L, r2SumSLC = 0.0L, r2SumSLS = 0.0L;
		for (int j = 0; j < NumR2Points; j++) {  // r2 integration
			r2 = r2Abscissas[j];
			long double a12 = fabs(r1-r2);
			long double b12 = fabs(r1+r2);
			ChangeOfIntervalNoResize(LegendreAbscissasR12, r12Array, a12, b12);

			// Gauss-Laguerre only if r1 is large.
			if (r1 > CuspR3) {
				for (int m = 0; m < nR3Lag; m++) {
					r3Weights[m] = LaguerreWeightsR3[m];
					r3Abscissas[m] = LaguerreAbscissasR3[m];
				}
				NumR3Points = nR3Lag;
			}
			// Do Gauss-Legendre for r3 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0L, r1);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * exp(-r3Abscissas[m]) * (r1-0.0L)/2.0L;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = LaguerreAbscissasR3[m-nR3Leg] + r1;
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-r1);
				}
			}

			long double r3SumCLC = 0.0L, r3SumCLS = 0.0, r3SumSLC = 0.0, r3SumSLS = 0.0;
			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				r3 = r3Abscissas[g];
				long double TempCoeff = r3Weights[g] * r2Weights[j] * r1Weights[i];
				long double ExpR1R2R3 = expl(r1+r2+r3);
				long double a13 = fabs(r1-r3);
				long double b13 = fabs(r1+r3);
				ChangeOfIntervalNoResize(LegendreAbscissasR13, r13Array, a13, b13);

				long double r12SumCLC = 0.0L, r12SumCLS = 0.0L, r12SumSLC = 0.0L, r12SumSLS = 0.0L;
				for (int k = 0; k < nR12; k++) {  // r12 integration
					long double r12 = r12Array[k];
					long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0L*r1*r2);
					long double Sin12 = sqrt(1.0L - Cos12*Cos12);
					long double rho = 0.5 * sqrt(2.0L*(r1*r1 + r2*r2) - r12*r12);
					long double n1rho = sf_bessel_n1(kappa*rho);
					long double j2rho = sf_bessel_j2(kappa*rho);
					long double n2rho = sf_bessel_n2(kappa*rho);
					long double ExpR12R3 = exp(-(r12/2.0L + r3));
					long double S22 = ExpR12R3 * SqrtKappa * j2rho;
					long double fshrho = fshielding(rho, mu, shpower);
					long double fshterm = fshielding1(rho, mu, shpower) / rho * (2.0L * n2rho - kappa*rho * n1rho) - 0.5L * fshielding2(rho, mu, shpower) * n2rho;
					//long double Part1 = 2.0L * ExpR12R3 * n2rho * fshrho * fshterm;

					long double r13SumCLC = 0.0L, r13SumCLS = 0.0L, r13SumSLC = 0.0L, r13SumSLS = 0.0L;
					for (int p = 0; p < nR13; p++) {  // r13 integration
						long double r13 = r13Array[p];
						long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0L*r1*r3);
						long double Sin13 = sqrt(1.0L - Cos13*Cos13);
						long double rhop = 0.5 * sqrt(2.0L*(r1*r1 + r3*r3) - r13*r13);
						long double j2rhop = sf_bessel_j2(kappa*rhop);
						long double n2rhop = sf_bessel_n2(kappa*rhop);
						long double ExpR13R2 = exp(-(r13/2.0L + r2));
						long double fshrhop = fshielding(rhop, mu, shpower);

						long double S23 =  ExpR13R2 * SqrtKappa * j2rhop;
						long double C23 = -ExpR13R2 * SqrtKappa * n2rhop * fshrhop;

						long double Pot = (2.0L/r1 - 2.0L/r2 - 2.0L/r13);
						long double Part2 = n1rho * fshrho * Pot + fshterm;
						Part2 = Part2 * ExpR13R2 * n2rhop * fshrhop;
						long double dTau = r2 * r3 * r12 * r13;
						
						long double Phi23SumCLC = 0.0L, Phi23SumCLS = 0.0L, Phi23SumSLC = 0.0L, Phi23SumSLS = 0.0L;
						for (int m = 1; m <= nPhi23; m++) {  // phi_23 integration
							long double Phi23 = (2.0L*m - 1.0L)*PI/(2.0L*nPhi23);
							long double r23 = sqrt(r2*r2 + r3*r3 - 2.0L*r2*r3*(Sin12*Sin13*cosl(Phi23) + Cos12*Cos13));

							long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0L*r2*r3);
							//long double Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0L * rho * rhop);
							long double Ang1 = 4.0L*rho*rho + 4.0L*rhop*rhop - r23*r23;
							long double Ang = 3.0L/8.0L * Ang1*Ang1 / (16*rho*rho*rhop*rhop) - 0.5L;

							//@TODO: Remove.
							long double PotP = (2.0L/r1 - 2.0L/r3 - 2.0L/r12);
							long double C22 = -ExpR12R3 * SqrtKappa * n2rho * fshrho;
	
							// SLS
							//long double RetSLS = S23 * S22 * PotP * Ang;

							long double j1rho = sf_bessel_j1(kappa*rho);
							long double j1rhop = sf_bessel_j1(kappa*rhop);
							//Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0 * rho * rhop);
							//S22 = ExpR12R3 * SqrtKappa * j2rho;
							//S23 = ExpR13R2 * SqrtKappa * j2rhop;
							//long double RetSLS = S22 * S22 * Pot;
							//long double RetSLS = S23 * S23 * PotP;
							//long double RetSLS = S22 * S23 * PotP * Ang;
							long double RetSLS = S23 * S22 * Pot * Ang;
							Phi23SumSLS += sf * ExpR1R2R3 * RetSLS * dTau;

							// CLS
							//long double RetCLS = C23 * S22 * Pot * Ang;

							//long double RetCLS = C22 * S22 * Pot;
							//long double RetCLS = C23 * S23 * PotP;
							//long double RetCLS = C22 * S23 * PotP * Ang;
							long double RetCLS = C23 * S22 * Pot * Ang;
							Phi23SumCLS += sf * ExpR1R2R3 * RetCLS * dTau;

							// SLC
							long double LCPart1 = C22 * Pot * Ang;
							long double LCPart2 = SqrtKappa * ExpR12R3 * fshterm;
							long double RetSLC = sf * S23 * LCPart1 - (S22 + sf * Ang * S23) * LCPart2;
							Phi23SumSLC += ExpR1R2R3 * RetSLC * dTau;

							// CLC
							long double RetCLC = sf * C23 * LCPart1 - (C22 + sf * Ang * C23) * LCPart2;
							Phi23SumCLC += ExpR1R2R3 * RetCLC * dTau;

							//// SLC
							//long double fshtermp = fshielding1(rhop, mu, shpower) / rhop * (n2rhop + cosl(kappa*rhop)) - 0.5 * fshielding2(rhop, mu, shpower) * n2rhop;
							//long double RetSLC1 = -SqrtKappa * Ang * ExpR12R3 * (Pot * n2rho * fshrho + fshterm);
							//long double RetSLC2 = -SqrtKappa * 2.0L * ExpR13R2 * (PotP * n2rhop * fshrhop + fshtermp);
							//Phi23SumSLC += ExpR1R2R3 * S23 * (RetSLC1 + sf * RetSLC2) * dTau;
							//Phi23SumSLC = 0.0L;  //@TODO: Fix this!
							//
							//long double RetCLC = C23 * C22 * Pot * Ang;
							//RetCLC += 
							//
							//// CLC
							//long double RetCLC = kappa * ExpR12R3 * (Part1 + sf * Part2 * Ang);
							//Phi23SumCLC += RetCLC * ExpR1R2R3 * dTau;
							//Phi23SumSLC = 0.0L;  //@TODO: Fix this!
						}
						r13SumCLC += LegendreWeightsR13[p] * Phi23SumCLC / nPhi23 * (b13-a13)/2.0L;
						r13SumSLC += LegendreWeightsR13[p] * Phi23SumSLC / nPhi23 * (b13-a13)/2.0L;
						r13SumCLS += LegendreWeightsR13[p] * Phi23SumCLS / nPhi23 * (b13-a13)/2.0L;
						r13SumSLS += LegendreWeightsR13[p] * Phi23SumSLS / nPhi23 * (b13-a13)/2.0L;
					}
					r12SumCLC += LegendreWeightsR12[k] * r13SumCLC * (b12-a12)/2.0L;
					r12SumSLC += LegendreWeightsR12[k] * r13SumSLC * (b12-a12)/2.0L;
					r12SumCLS += LegendreWeightsR12[k] * r13SumCLS * (b12-a12)/2.0L;
					r12SumSLS += LegendreWeightsR12[k] * r13SumSLS * (b12-a12)/2.0L;
				}
				r3SumCLC += r3Weights[g] * r12SumCLC;
				r3SumSLC += r3Weights[g] * r12SumSLC;
				r3SumCLS += r3Weights[g] * r12SumCLS;
				r3SumSLS += r3Weights[g] * r12SumSLS;
			}
			r2SumCLC += r2Weights[j] * r3SumCLC;
			r2SumSLC += r2Weights[j] * r3SumSLC;
			r2SumCLS += r2Weights[j] * r3SumCLS;
			r2SumSLS += r2Weights[j] * r3SumSLS;
		}
		//@TODO: Do I really need the atomic?
		//#pragma omp atomic
		r1SumCLC += r1Weights[i] * r2SumCLC;
		r1SumSLC += r1Weights[i] * r2SumSLC;
		r1SumCLS += r1Weights[i] * r2SumCLS;
		r1SumSLS += r1Weights[i] * r2SumSLS;
	}

	//@TODO: 4x P-wave?
	//r1SumCLC /= 2.0L;
	//r1SumSLC /= 2.0L;
	//r1SumCLS /= 2.0L;
	//r1SumSLS /= 2.0L;

	CLC += r1SumCLC;
	SLC += r1SumSLC;
	CLS += r1SumCLS;
	SLS += r1SumSLS;

	return;
}


void GaussIntegrationPhi12_LongLong_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nPhi12, int nR13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, double &CLC, double &SLC, double &CLS, double &SLS)
{
	long double r1, r2, r3;
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR23, LegendreWeightsR23;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r3Abscissas, r2Weights, r3Weights;
	vector <long double> r13Array(nR13), r23Array(nR23);
	int NumR2Points, NumR3Points, Prog = 0;
	long double SqrtKappa = sqrt(kappa);

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(r1Abscissas, r1Weights, nR1);
	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR23, LegendreWeightsR23, nR23);
	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);

	long double r1SumCLC = 0.0L, r1SumCLS = 0.0L, r1SumSLC = 0.0L, r1SumSLS = 0.0L;
	//#pragma omp parallel for shared(r1Sum,r1Abscissas,r1Weights,Powers) private(r1,r2,r3,r13,r23,r2Sum,r3Sum,r23sum,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points)
	#pragma omp parallel for shared(r1SumCLC,r1SumCLS,r1SumSLC,r1SumSLS,r1Abscissas,r1Weights) private(r1,r2,r3,r13Array,r23Array,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		WriteProgress(string("Long-range - Long-range r23"), Prog, i, nR1);

		// These are private, so they need to be initialized.
		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);
		r13Array.resize(nR13);
		r23Array.resize(nR23);

		r1 = r1Abscissas[i];

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR3) {
			for (int m = 0; m < nR3Lag; m++) {
				r3Weights[m] = LaguerreWeightsR3[m];
				r3Abscissas[m] = LaguerreAbscissasR3[m];
			}
			NumR3Points = nR3Lag;
		}
		// Do Gauss-Legendre for r3 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0L, r1);
			NumR3Points = nR3Leg + nR3Lag;
			for (int m = 0; m < nR3Leg; m++) {
				r3Weights[m] = LegendreWeightsR3[m] * exp(-r3Abscissas[m]) * (r1-0.0L)/2.0L;
			}
			for (int m = nR3Leg; m < NumR3Points; m++) {
				r3Abscissas[m] = LaguerreAbscissasR3[m-nR3Leg] + r1;
				r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-r1);
			}
		}

		long double r3SumCLC = 0.0L, r3SumCLS = 0.0L, r3SumSLC = 0.0L, r3SumSLS = 0.0L;
		for (int j = 0; j < NumR3Points; j++) {  // r3 integration
			r3 = r3Abscissas[j];
			long double a13 = fabs(r1-r3);
			long double b13 = fabs(r1+r3);
			ChangeOfIntervalNoResize(LegendreAbscissasR13, r13Array, a13, b13);

			// Gauss-Laguerre only if r3 is large.
			if (r3 > CuspR2) {
				for (int m = 0; m < nR2Lag; m++) {
					r2Weights[m] = LaguerreWeightsR2[m];
					r2Abscissas[m] = LaguerreAbscissasR2[m];
				}
				NumR2Points = nR2Lag;
			}
			// Do Gauss-Legendre for r2 integration over [0,r3] and then Gauss-Laguerre over [r3,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0L, r3);
				NumR2Points = nR2Leg + nR2Lag;
				for (int m = 0; m < nR2Leg; m++) {
					r2Weights[m] = LegendreWeightsR2[m] * exp(-r2Abscissas[m]) * (r3-0.0L)/2.0L;
				}
				for (int m = nR2Leg; m < NumR2Points; m++) {
					r2Abscissas[m] = LaguerreAbscissasR2[m-nR2Leg] + r3;
					r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-r3);
				}
			}

			long double r2SumCLC = 0.0L, r2SumCLS = 0.0L, r2SumSLC = 0.0L, r2SumSLS = 0.0L;
			for (int g = 0; g < NumR2Points; g++) {  // r2 integration
				r2 = r2Abscissas[g];
				long double TempCoeff = r3Weights[j] * r2Weights[g] * r1Weights[i];
				long double ExpR1R2R3 = expl(r1+r2+r3);
				long double a23 = fabs(r2-r3);
				long double b23 = fabs(r2+r3);
				ChangeOfIntervalNoResize(LegendreAbscissasR23, r23Array, a23, b23);

				long double r23SumCLC = 0.0L, r23SumCLS = 0.0L, r23SumSLC = 0.0L, r23SumSLS = 0.0L;
				for (int k = 0; k < nR23; k++) {  // r23 integration
					long double r23 = r23Array[k];
					long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0L*r2*r3);
					long double Sin23 = sqrt(1.0L - Cos23*Cos23);
					long double Pot = 2.0L/r23;
					long double PotP = Pot;

					long double r13SumCLC = 0.0L, r13SumCLS = 0.0L, r13SumSLC = 0.0L, r13SumSLS = 0.0L;
					for (int p = 0; p < nR13; p++) {  // r13 integration
						long double r13 = r13Array[p];
						long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0L*r1*r3);
						long double Sin13 = sqrt(1.0L - Cos13*Cos13);
						long double rhop = 0.5L * sqrt(2.0L*(r1*r1 + r3*r3) - r13*r13);
						long double j2rhop = sf_bessel_j2(kappa*rhop);
						long double n2rhop = sf_bessel_n2(kappa*rhop);
						long double ExpR13R2 = exp(-(r13/2.0L + r2));
						long double S23 =  ExpR13R2 * SqrtKappa * j2rhop;
						long double C23 = -ExpR13R2 * SqrtKappa * n2rhop * fshielding(rhop, mu, shpower);
						long double dTau = r1 * r2 * r23 * r13;

						long double Phi12SumCLC = 0.0L, Phi12SumCLS = 0.0L, Phi12SumSLC = 0.0L, Phi12SumSLS = 0.0L;
						for (int m = 1; m <= nPhi12; m++) {  // phi_12 integration
							long double Phi12 = (2.0L*m - 1.0L)*PI/(2.0L*nPhi12);
							long double r12 = sqrt(r1*r1 + r2*r2 - 2.0L*r1*r2*(Sin13*Sin23*cosl(Phi12) + Cos13*Cos23));
							long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0L*r1*r2);
							long double rho = 0.5 * sqrt(2.0L*(r1*r1 + r2*r2) - r12*r12);
							long double j2rho = sf_bessel_j2(kappa*rho);
							long double n2rho = sf_bessel_n2(kappa*rho);
							long double ExpR12R3 = exp(-(r12/2.0L + r3));

							long double S22 =  ExpR12R3 * SqrtKappa * j2rho;
							long double C22 = -ExpR12R3 * SqrtKappa * n2rho * fshielding(rho, mu, shpower);
							//long double Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0L * rho * rhop);
							long double Ang1 = 4.0L*rho*rho + 4.0L*rhop*rhop - r23*r23;
							long double Ang = 3.0L/8.0L * Ang1*Ang1 / (16*rho*rho*rhop*rhop) - 0.5L;

							//long double RetSLS = S22 * S23 * PotP * Ang;
							//Phi12SumSLS += sf * ExpR1R2R3 * RetSLS * dTau;

							//long double RetCLS = C23 * S22 * Pot * Ang;
							//Phi12SumCLS += sf * ExpR1R2R3 * RetCLS * dTau;

							//long double RetSLC1 = -Ang * SqrtKappa * ExpR12R3 * Pot * n2rho * fshielding(rho, mu, shpower);
							//long double RetSLC2 = -2.0L * SqrtKappa * ExpR13R2 * Pot * n2rhop * fshielding(rhop, mu, shpower);
							//Phi12SumSLC += ExpR1R2R3 * S23 * (RetSLC1 + sf * RetSLC2) * dTau;
							//Phi12SumSLC = 0.0L;  //@TODO: Fix this!

							//long double RetCLC = C22 * C23 * Pot * Ang;
							//Phi12SumCLC += sf * ExpR1R2R3 * RetCLC * dTau;


							// SLS
							//long double RetSLS = S23 * S22 * Pot * Ang;

							long double j1rho = sf_bessel_j1(kappa*rho);
							long double j1rhop = sf_bessel_j1(kappa*rhop);
							//Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0 * rho * rhop);
							//S22 = ExpR12R3 * SqrtKappa * j2rho;
							//S23 = ExpR13R2 * SqrtKappa * j2rhop;
							//long double RetSLS = S22 * S22 * Pot;
							//long double RetSLS = S23 * S23 * PotP;
							//long double RetSLS = S22 * S23 * PotP * Ang;
							long double RetSLS = S23 * S22 * Pot * Ang;
							Phi12SumSLS += sf * ExpR1R2R3 * RetSLS * dTau;

							// CLS
							//long double RetCLS = C23 * S22 * Pot * Ang;

							//long double RetCLS = C22 * S22 * Pot;
							//long double RetCLS = C23 * S23 * PotP;
							//long double RetCLS = C22 * S23 * PotP * Ang;
							long double RetCLS = C23 * S22 * Pot * Ang;
							Phi12SumCLS += sf * ExpR1R2R3 * RetCLS * dTau;

							// SLC
							long double LCPart1 = C22 * Pot * Ang;
							long double RetSLC = sf * S23 * LCPart1;
							Phi12SumSLC += ExpR1R2R3 * RetSLC * dTau;

							// CLC
							long double RetCLC = sf * C23 * LCPart1;
							Phi12SumCLC += ExpR1R2R3 * RetCLC * dTau;
						}
						r13SumCLC += LegendreWeightsR13[p] * Phi12SumCLC / nPhi12 * (b13-a13)/2.0L;
						r13SumSLC += LegendreWeightsR13[p] * Phi12SumSLC / nPhi12 * (b13-a13)/2.0L;
						r13SumCLS += LegendreWeightsR13[p] * Phi12SumCLS / nPhi12 * (b13-a13)/2.0L;
						r13SumSLS += LegendreWeightsR13[p] * Phi12SumSLS / nPhi12 * (b13-a13)/2.0L;
					}
					r23SumCLC += LegendreWeightsR23[k] * r13SumCLC * (b23-a23)/2.0L;
					r23SumSLC += LegendreWeightsR23[k] * r13SumSLC * (b23-a23)/2.0L;
					r23SumCLS += LegendreWeightsR23[k] * r13SumCLS * (b23-a23)/2.0L;
					r23SumSLS += LegendreWeightsR23[k] * r13SumSLS * (b23-a23)/2.0L;
				}
				r2SumCLC += r2Weights[g] * r23SumCLC;
				r2SumSLC += r2Weights[g] * r23SumSLC;
				r2SumCLS += r2Weights[g] * r23SumCLS;
				r2SumSLS += r2Weights[g] * r23SumSLS;
			}
			r3SumCLC += r3Weights[j] * r2SumCLC;
			r3SumSLC += r3Weights[j] * r2SumSLC;
			r3SumCLS += r3Weights[j] * r2SumCLS;
			r3SumSLS += r3Weights[j] * r2SumSLS;
		}
		r1SumCLC += r1Weights[i] * r3SumCLC;
		r1SumSLC += r1Weights[i] * r3SumSLC;
		r1SumCLS += r1Weights[i] * r3SumCLS;
		r1SumSLS += r1Weights[i] * r3SumSLS;
	}

	//r1SumCLC /= 2.0L;
	//r1SumSLC /= 2.0L;
	//r1SumCLS /= 2.0L;
	//r1SumSLS /= 2.0L;

	CLC = r1SumCLC;
	SLC = r1SumSLC;
	CLS = r1SumCLS;
	SLS = r1SumSLS;

	return;
}


void GaussIntegrationPhi13_LongLong_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nPhi13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, double &CLC, double &SLC, double &CLS, double &SLS)
{
	long double r1, r2, r3;
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR23, LegendreWeightsR23;
	vector <long double> LegendreAbscissasR12, LegendreWeightsR12;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r3Abscissas, r2Weights, r3Weights;
	vector <long double> r12Array(nR12), r23Array(nR23);
	int NumR2Points, NumR3Points, Prog = 0;
	long double SqrtKappa = sqrt(kappa);

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(r1Abscissas, r1Weights, nR1);
	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR23, LegendreWeightsR23, nR23);
	GaussLegendre(LegendreAbscissasR12, LegendreWeightsR12, nR12);

	long double r1SumCLC = 0.0L, r1SumCLS = 0.0L, r1SumSLC = 0.0L, r1SumSLS = 0.0L;
	#pragma omp parallel for shared(r1SumCLC,r1SumCLS,r1SumSLC,r1SumSLS,r1Abscissas,r1Weights) private(r1,r2,r3,r12Array,r23Array,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		WriteProgress(string("Long-range - Long-range r23"), Prog, i, nR1);

		// These are private, so they need to be initialized.
		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);
		r12Array.resize(nR12);
		r23Array.resize(nR23);

		r1 = r1Abscissas[i];

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR2) {
			for (int m = 0; m < nR2Lag; m++) {
				r2Weights[m] = LaguerreWeightsR2[m];
				r2Abscissas[m] = LaguerreAbscissasR2[m];
			}
			NumR2Points = nR2Lag;
		}
		// Do Gauss-Legendre for r2 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0L, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * exp(-r2Abscissas[m]) * (r1-0.0L)/2.0L;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = LaguerreAbscissasR2[m-nR2Leg] + r1;
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-r1);
			}
		}

		long double r2SumCLC = 0.0L, r2SumCLS = 0.0L, r2SumSLC = 0.0L, r2SumSLS = 0.0L;
		for (int j = 0; j < NumR2Points; j++) {  // r2 integration
			r2 = r2Abscissas[j];
			long double a12 = fabs(r1-r2);
			long double b12 = fabs(r1+r2);
			ChangeOfIntervalNoResize(LegendreAbscissasR12, r12Array, a12, b12);

			// Gauss-Laguerre only if r2 is large.
			if (r2 > CuspR3) {
				for (int m = 0; m < nR3Lag; m++) {
					r3Weights[m] = LaguerreWeightsR3[m];
					r3Abscissas[m] = LaguerreAbscissasR3[m];
				}
				NumR3Points = nR3Lag;
			}
			// Do Gauss-Legendre for r3 integration over [0,r2] and then Gauss-Laguerre over [r2,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0L, r2);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * exp(-r3Abscissas[m]) * (r2-0.0L)/2.0L;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = LaguerreAbscissasR3[m-nR3Leg] + r2;
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-r2);
				}
			}

			long double r3SumCLC = 0.0L, r3SumCLS = 0.0L, r3SumSLC = 0.0L, r3SumSLS = 0.0L;
			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				r3 = r3Abscissas[g];
				long double TempCoeff = r3Weights[g] * r2Weights[j] * r1Weights[i];
				long double ExpR1R2R3 = expl(r1+r2+r3);
				long double a23 = fabs(r2-r3);
				long double b23 = fabs(r2+r3);
				ChangeOfIntervalNoResize(LegendreAbscissasR23, r23Array, a23, b23);

				long double r23SumCLC = 0.0L, r23SumCLS = 0.0L, r23SumSLC = 0.0L, r23SumSLS = 0.0L;
				for (int k = 0; k < nR23; k++) {  // r23 integration
					long double r23 = r23Array[k];
					long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0L*r2*r3);
					long double Sin23 = sqrt(1.0L - Cos23*Cos23);
					long double Pot = 2.0L/r23;

					long double r12SumCLC = 0.0L, r12SumCLS = 0.0L, r12SumSLC = 0.0L, r12SumSLS = 0.0L;
					for (int p = 0; p < nR12; p++) {  // r12 integration
						long double r12 = r12Array[p];
						long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0L*r1*r2);
						long double Sin12 = sqrt(1.0L - Cos12*Cos12);
						long double rho = 0.5L * sqrt(2.0L*(r1*r1 + r2*r2) - r12*r12);
						long double j2rho = sf_bessel_j2(kappa*rho);
						long double n2rho = sf_bessel_n2(kappa*rho);
						long double ExpR12R3 = exp(-(r12/2.0L + r3));
						long double S22 =  ExpR12R3 * SqrtKappa * j2rho;
						long double C22 = -ExpR12R3 * SqrtKappa * n2rho * fshielding(rho, mu, shpower);
						long double dTau = r1 * r3 * r23 * r12;

						long double Phi13SumCLC = 0.0L, Phi13SumCLS = 0.0L, Phi13SumSLC = 0.0L, Phi13SumSLS = 0.0L;
						for (int m = 1; m <= nPhi13; m++) {  // phi_13 integration
							long double Phi13 = (2.0L*m - 1.0L)*PI/(2.0L*nPhi13);
							long double r13 = sqrt(r1*r1 + r3*r3 - 2.0L*r1*r3*(Sin12*Sin23*cosl(Phi13) + Cos12*Cos23));
							long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0L*r1*r3);
							long double rhop = 0.5L * sqrt(2.0L*(r1*r1 + r3*r3) - r13*r13);
							long double j2rhop = sf_bessel_j2(kappa*rhop);
							long double n2rhop = sf_bessel_n2(kappa*rhop);
							long double ExpR13R2 = exp(-(r13/2.0L + r2));

							long double S23 =  ExpR13R2 * SqrtKappa * j2rhop;
							long double C23 = -ExpR13R2 * SqrtKappa * n2rhop * fshielding(rhop, mu, shpower);
							//long double Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0L * rho * rhop);
							long double Ang1 = 4.0L*rho*rho + 4.0L*rhop*rhop - r23*r23;
							long double Ang = 3.0L/8.0L * Ang1*Ang1 / (16*rho*rho*rhop*rhop) - 0.5L;

							long double PotP = Pot;  //@TODO: Remove
							//long double RetSLS = S22 * S23 * PotP * Ang;
							//Phi13SumSLS += sf * ExpR1R2R3 * RetSLS * dTau;

							//long double RetCLS = C23 * S22 * Pot * Ang;
							//Phi13SumCLS += sf * ExpR1R2R3 * RetCLS * dTau;

							//// THIS ONE WORKS FOR PHI13.
							//long double RetSLC1 = -2.0L * SqrtKappa * ExpR12R3 * Pot * n2rho * fshielding(rho, mu);
							//long double RetSLC2 = -Ang * SqrtKappa * ExpR13R2 * Pot * n2rhop * fshielding(rhop, mu);
							//Phi13SumSLC += ExpR1R2R3 * S22 * (RetSLC1 + sf * RetSLC2) * dTau;
							//Phi13SumSLC = 0.0L;  //@TODO: Fix this!

							//long double RetCLC = C22 * C23 * Pot * Ang;
							//Phi13SumCLC += sf * ExpR1R2R3 * RetCLC * dTau;


							// SLS
							//long double RetSLS = S23 * S22 * Pot * Ang;

							long double j1rho = sf_bessel_j1(kappa*rho);
							long double j1rhop = sf_bessel_j1(kappa*rhop);
							//Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0 * rho * rhop);
							//S22 = ExpR12R3 * SqrtKappa * j2rho;
							//S23 = ExpR13R2 * SqrtKappa * j2rhop;
							//long double RetSLS = S22 * S22 * Pot;
							//long double RetSLS = S23 * S23 * PotP;
							//long double RetSLS = S22 * S23 * PotP * Ang;
							long double RetSLS = S23 * S22 * Pot * Ang;
							Phi13SumSLS += sf * ExpR1R2R3 * RetSLS * dTau;

							// CLS
							//long double RetCLS = C23 * S22 * Pot * Ang;

							//long double RetCLS = C22 * S22 * Pot;
							//long double RetCLS = C23 * S23 * PotP;
							//long double RetCLS = C22 * S23 * PotP * Ang;
							long double RetCLS = C23 * S22 * Pot * Ang;
							Phi13SumCLS += sf * ExpR1R2R3 * RetCLS * dTau;

							// SLC
							long double LCPart1 = C22 * Pot * Ang;
							long double RetSLC = sf * S23 * LCPart1;
							Phi13SumSLC += ExpR1R2R3 * RetSLC * dTau;

							// CLC
							long double RetCLC = sf * C23 * LCPart1;
							Phi13SumCLC += ExpR1R2R3 * RetCLC * dTau;
						}
						r12SumCLC += LegendreWeightsR12[p] * Phi13SumCLC / nPhi13 * (b12-a12)/2.0L;
						r12SumSLC += LegendreWeightsR12[p] * Phi13SumSLC / nPhi13 * (b12-a12)/2.0L;
						r12SumCLS += LegendreWeightsR12[p] * Phi13SumCLS / nPhi13 * (b12-a12)/2.0L;
						r12SumSLS += LegendreWeightsR12[p] * Phi13SumSLS / nPhi13 * (b12-a12)/2.0L;
					}
					r23SumCLC += LegendreWeightsR23[k] * r12SumCLC * (b23-a23)/2.0L;
					r23SumSLC += LegendreWeightsR23[k] * r12SumSLC * (b23-a23)/2.0L;
					r23SumCLS += LegendreWeightsR23[k] * r12SumCLS * (b23-a23)/2.0L;
					r23SumSLS += LegendreWeightsR23[k] * r12SumSLS * (b23-a23)/2.0L;
				}
				r3SumCLC += r3Weights[g] * r23SumCLC;
				r3SumSLC += r3Weights[g] * r23SumSLC;
				r3SumCLS += r3Weights[g] * r23SumCLS;
				r3SumSLS += r3Weights[g] * r23SumSLS;
			}
			r2SumCLC += r2Weights[j] * r3SumCLC;
			r2SumSLC += r2Weights[j] * r3SumSLC;
			r2SumCLS += r2Weights[j] * r3SumCLS;
			r2SumSLS += r2Weights[j] * r3SumSLS;
		}
		r1SumCLC += r1Weights[i] * r2SumCLC;
		r1SumSLC += r1Weights[i] * r2SumSLC;
		r1SumCLS += r1Weights[i] * r2SumCLS;
		r1SumSLS += r1Weights[i] * r2SumSLS;
	}

	//r1SumCLC /= 2.0L;
	//r1SumSLC /= 2.0L;
	//r1SumCLS /= 2.0L;
	//r1SumSLS /= 2.0L;

	CLC += r1SumCLC;
	SLC += r1SumSLC;
	CLS += r1SumCLS;
	SLS += r1SumSLS;

	return;
}
