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
	long double sinx, cosx;
	sinx = sin(x);
	cosx = cos(x);
	//sincos(x, &sinx, &cosx);
	return (sinx/x - cosx)/x;
};


// Spherical Bessel function of the second kind (spherical Neumann)
long double sf_bessel_n1(long double x)
{
	return -(cosl(x)/x + sinl(x))/x;
};


// f_sh in C22 and C23
long double fshielding(long double rho, long double mu)
{
	long double sh = 1.0 - exp(-mu*rho) * (1.0 + mu*rho / 2.0);
	return sh*sh*sh;
};


// First derivative of f_sh
long double fshielding1(long double rho, long double mu)
{
	long double MuRho = mu*rho;
	long double ExpMuRho = exp(-MuRho);
	long double Part = (1.0 - ExpMuRho * (1.0 + MuRho/2.0));
	long double sh = 1.5 * ExpMuRho * mu * Part*Part * (1.0 + MuRho);
	return sh;
};


// Second derivative of f_sh
long double fshielding2(long double rho, long double mu)
{
	long double MuRho = mu*rho;
	long double ExpMuRho = exp(-MuRho);
	long double sh = 1.5 * ExpMuRho * mu*mu * (1.0 - ExpMuRho * (1.0 + MuRho/2.0));
	sh = sh * (-MuRho + ExpMuRho * (1.0 + 3.0*MuRho + 1.5*MuRho*MuRho));
	return sh;
};


void GaussIntegrationPhi23_LongLong(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int sf, double &CLC, double &SLC, double &CLS, double &SLS)
{
	double r1, r2, r3;
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

	long double r1SumCLC = 0.0, r1SumCLS = 0.0, r1SumSLC = 0.0, r1SumSLS = 0.0;
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
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * exp(-r2Abscissas[m]) * (r1-0.0)/2.0;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = LaguerreAbscissasR2[m-nR2Leg] + r1;
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-r1);
			}
		}

		long double r2SumCLC = 0.0, r2SumCLS = 0.0, r2SumSLC = 0.0, r2SumSLS = 0.0;
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
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0, r1);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * exp(-r3Abscissas[m]) * (r1-0.0)/2.0;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = LaguerreAbscissasR3[m-nR3Leg] + r1;
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-r1);
				}
			}

			long double r3SumCLC = 0.0, r3SumCLS = 0.0, r3SumSLC = 0.0, r3SumSLS = 0.0;
			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				r3 = r3Abscissas[g];
				long double TempCoeff = r3Weights[g] * r2Weights[j] * r1Weights[i];
				long double ExpR1R2R3 = expl(r1+r2+r3);
				long double a13 = fabs(r1-r3);
				long double b13 = fabs(r1+r3);
				ChangeOfIntervalNoResize(LegendreAbscissasR13, r13Array, a13, b13);

				long double r12SumCLC = 0.0, r12SumCLS = 0.0, r12SumSLC = 0.0, r12SumSLS = 0.0;
				for (int k = 0; k < nR12; k++) {  // r12 integration
					long double r12 = r12Array[k];
					long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0*r1*r2);
					long double Sin12 = sqrt(1.0 - Cos12*Cos12);
					long double rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12*r12);
					long double j1rho = sf_bessel_j1(kappa*rho);
					long double n1rho = sf_bessel_n1(kappa*rho);
					long double ExpR12R3 = exp(-(r12/2.0 + r3));
					long double S22 = ExpR12R3 * SqrtKappa * j1rho;
					long double fshrho = fshielding(rho, mu);
					long double fshterm = fshielding1(rho, mu) / rho * (n1rho + cosl(kappa*rho)) - 0.5 * fshielding2(rho, mu) * n1rho;
					long double Part1 = 2.0 * ExpR12R3 * n1rho * fshrho * fshterm;

					long double r13SumCLC = 0.0, r13SumCLS = 0.0, r13SumSLC = 0.0, r13SumSLS = 0.0;
					for (int p = 0; p < nR13; p++) {  // r13 integration
						long double r13 = r13Array[p];
						long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0*r1*r3);
						long double Sin13 = sqrt(1.0 - Cos13*Cos13);
						long double rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13*r13);
						long double j1rhop = sf_bessel_j1(kappa*rhop);
						long double n1rhop = sf_bessel_n1(kappa*rhop);
						long double ExpR13R2 = exp(-(r13/2.0 + r2));
						long double fshrhop = fshielding(rhop, mu);

						long double S23 =  ExpR13R2 * SqrtKappa * j1rhop;
						long double C23 = -ExpR13R2 * SqrtKappa * n1rhop * fshrhop;

						long double Pot = (2.0/r1 - 2.0/r2 - 2.0/r13);
						long double Part2 = n1rho * fshrho * Pot + fshterm;
						Part2 = Part2 * ExpR13R2 * n1rhop * fshrhop;
						long double dTau = r2 * r3 * r12 * r13;
						
						long double Phi23SumCLC = 0.0, Phi23SumCLS = 0.0, Phi23SumSLC = 0.0, Phi23SumSLS = 0.0;
						for (int m = 1; m <= nPhi23; m++) {  // phi_23 integration
							long double Phi23 = (2.0*m - 1.0)*PI/(2.0*nPhi23);
							long double r23 = sqrt(r2*r2 + r3*r3 - 2.0*r2*r3*(Sin12*Sin13*cosl(Phi23) + Cos12*Cos13));

							long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0*r2*r3);
							long double Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0 * rho * rhop);

							//@TODO: Remove.
							long double PotP = (2.0/r1 - 2.0/r3 - 2.0/r12);
							long double C22 = -ExpR12R3 * SqrtKappa * n1rho * fshrho;
	
							// SLS
							long double RetSLS = S23 * S22 * Pot * Ang;
							//long double RetSLS = S22 * S23 * PotP * Ang;
							//long double RetSLS = S23 * S23 * PotP * 2.0;
							//long double RetSLS = S22 * S22 * Pot * 2.0;
							Phi23SumSLS += sf * ExpR1R2R3 * RetSLS * dTau;
							//Phi23SumSLS += sf * expl(r1+r2) * RetSLS * dTau;

							// CLS
							long double RetCLS = C23 * S22 * Pot * Ang;
							//long double RetCLS = C22 * S23 * PotP * Ang;
							//long double RetCLS = C23 * S23 * PotP * 2.0;
							//long double RetCLS = C22 * S22 * Pot * 2.0;
							Phi23SumCLS += sf * ExpR1R2R3 * RetCLS * dTau;



							////long double RetSLC = -ExpR12R3 * SqrtKappa * (n1rho * fshrho * Pot + fshterm);
							////Phi23SumSLC += sf * S23 * Ang * RetSLC * ExpR1R2R3 * dTau;
							////long double RetSLC = -ExpR12R3 * SqrtKappa * (n1rho * fshrho * Pot + fshterm);
							////Phi23SumSLC += (S22 + sf * S23) * Ang * RetSLC * ExpR1R2R3 * dTau;
							//long double RetSLC = -ExpR12R3 * SqrtKappa * (n1rho * fshrho * Pot + fshterm);
							//Phi23SumSLC += (S22 * 2.0 + sf * S23 * Ang) * RetSLC * ExpR1R2R3 * dTau;
							////long double RetSLC = -ExpR12R3 * SqrtKappa * (n1rho * fshrho * Pot + fshterm);
							////Phi23SumSLC += (S22 * 2.0) * RetSLC * ExpR1R2R3 * dTau;


							//long double PotP = (2.0/r1 - 2.0/r3 - 2.0/r12);
							//long double fshtermp = fshielding1(rhop, mu) / rhop * (n1rhop + cosl(kappa*rhop)) - 0.5 * fshielding2(rhop, mu) * n1rhop;
							//long double RetSLC1 = -SqrtKappa * 2.0 * ExpR12R3 * (Pot * n1rho * fshrho + fshterm);
							//long double RetSLC2 = -SqrtKappa * Ang * ExpR13R2 * (PotP * n1rhop * fshrhop + fshtermp);
							//Phi23SumSLC += ExpR1R2R3 * S22 * (RetSLC1 + sf * RetSLC2) * dTau;

							//long double RetSLCBad = ExpR12R3 * SqrtKappa * (j1rho * Pot);
							////Phi23SumSLC += (S22 * 2.0 + sf * S23 * Ang) * RetSLCBad * ExpR1R2R3 * dTau;
							////Phi23SumSLC += (sf * S23 * Ang) * RetSLCBad * ExpR1R2R3 * dTau;
							//Phi23SumSLC += (S22 * 2.0) * RetSLCBad * ExpR1R2R3 * dTau;
							
							//long double PotP = (2.0/r1 - 2.0/r3 - 2.0/r12);
							////long double RetSLSBad = 2.0 * S22 * S22 * Pot + Ang * S23 * S22 * Pot + 2.0 * S23 * S23 * PotP + Ang * S22 * S23 * PotP;
							////long double RetSLSBad = 2.0 * S22 * S22 * Pot  + Ang * S23 * S22 * Pot;
							////long double RetSLSBad = 2.0 * S23 * S23 * PotP + Ang * S22 * S23 * PotP;
							//long double RetSLSBad = 2.0 * S22 * S22 * Pot;
							//Phi23SumSLC += RetSLSBad * ExpR1R2R3 * dTau;


							// THIS ONE WORKS.
							// SLC
							long double fshtermp = fshielding1(rhop, mu) / rhop * (n1rhop + cosl(kappa*rhop)) - 0.5 * fshielding2(rhop, mu) * n1rhop;
							long double RetSLC1 = -SqrtKappa * Ang * ExpR12R3 * (Pot * n1rho * fshrho + fshterm);
							long double RetSLC2 = -SqrtKappa * 2.0 * ExpR13R2 * (PotP * n1rhop * fshrhop + fshtermp);
							Phi23SumSLC += ExpR1R2R3 * S23 * (RetSLC1 + sf * RetSLC2) * dTau;

							// THIS ONE WORKS FOR PHI13.
							// SLC
							//long double fshtermp = fshielding1(rhop, mu) / rhop * (n1rhop + cosl(kappa*rhop)) - 0.5 * fshielding2(rhop, mu) * n1rhop;
							//long double RetSLC1 = -SqrtKappa * 2.0 * ExpR12R3 * (Pot * n1rho * fshrho + fshterm);
							//long double RetSLC2 = -SqrtKappa * Ang * ExpR13R2 * (PotP * n1rhop * fshrhop + fshtermp);
							//Phi23SumSLC += ExpR1R2R3 * S22 * (RetSLC1 + sf * RetSLC2) * dTau;

							// CLC
							long double RetCLC = kappa * ExpR12R3 * (Part1 + sf * Part2 * Ang);
							Phi23SumCLC += RetCLC * ExpR1R2R3 * dTau;
						}
						r13SumCLC += LegendreWeightsR13[p] * Phi23SumCLC / nPhi23 * (b13-a13)/2.0;
						r13SumSLC += LegendreWeightsR13[p] * Phi23SumSLC / nPhi23 * (b13-a13)/2.0;
						r13SumCLS += LegendreWeightsR13[p] * Phi23SumCLS / nPhi23 * (b13-a13)/2.0;
						r13SumSLS += LegendreWeightsR13[p] * Phi23SumSLS / nPhi23 * (b13-a13)/2.0;
					}
					r12SumCLC += LegendreWeightsR12[k] * r13SumCLC * (b12-a12)/2.0;
					r12SumSLC += LegendreWeightsR12[k] * r13SumSLC * (b12-a12)/2.0;
					r12SumCLS += LegendreWeightsR12[k] * r13SumCLS * (b12-a12)/2.0;
					r12SumSLS += LegendreWeightsR12[k] * r13SumSLS * (b12-a12)/2.0;
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

	r1SumCLC /= 2.0;
	r1SumSLC /= 2.0;
	r1SumCLS /= 2.0;
	r1SumSLS /= 2.0;

	CLC += r1SumCLC;
	SLC += r1SumSLC;
	CLS += r1SumCLS;
	SLS += r1SumSLS;

	return;
}


void GaussIntegrationPhi12_LongLong_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nPhi12, int nR13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int sf, double &CLC, double &SLC, double &CLS, double &SLS)
{
	double r1, r2, r3;
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

	long double r1SumCLC = 0.0, r1SumCLS = 0.0, r1SumSLC = 0.0, r1SumSLS = 0.0;
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
			ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0, r1);
			NumR3Points = nR3Leg + nR3Lag;
			for (int m = 0; m < nR3Leg; m++) {
				r3Weights[m] = LegendreWeightsR3[m] * exp(-r3Abscissas[m]) * (r1-0.0)/2.0;
			}
			for (int m = nR3Leg; m < NumR3Points; m++) {
				r3Abscissas[m] = LaguerreAbscissasR3[m-nR3Leg] + r1;
				r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-r1);
			}
		}

		long double r3SumCLC = 0.0, r3SumCLS = 0.0, r3SumSLC = 0.0, r3SumSLS = 0.0;
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
				ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0, r3);
				NumR2Points = nR2Leg + nR2Lag;
				for (int m = 0; m < nR2Leg; m++) {
					r2Weights[m] = LegendreWeightsR2[m] * exp(-r2Abscissas[m]) * (r3-0.0)/2.0;
				}
				for (int m = nR2Leg; m < NumR2Points; m++) {
					r2Abscissas[m] = LaguerreAbscissasR2[m-nR2Leg] + r3;
					r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-r3);
				}
			}

			long double r2SumCLC = 0.0, r2SumCLS = 0.0, r2SumSLC = 0.0, r2SumSLS = 0.0;
			for (int g = 0; g < NumR2Points; g++) {  // r2 integration
				r2 = r2Abscissas[g];
				long double TempCoeff = r3Weights[j] * r2Weights[g] * r1Weights[i];
				long double ExpR1R2R3 = expl(r1+r2+r3);
				long double a23 = fabs(r2-r3);
				long double b23 = fabs(r2+r3);
				ChangeOfIntervalNoResize(LegendreAbscissasR23, r23Array, a23, b23);

				long double r23SumCLC = 0.0, r23SumCLS = 0.0, r23SumSLC = 0.0, r23SumSLS = 0.0;
				for (int k = 0; k < nR23; k++) {  // r23 integration
					long double r23 = r23Array[k];
					long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0*r2*r3);
					long double Sin23 = sqrt(1.0 - Cos23*Cos23);
					long double Pot = 2.0/r23;
					long double PotP = Pot;

					long double r13SumCLC = 0.0, r13SumCLS = 0.0, r13SumSLC = 0.0, r13SumSLS = 0.0;
					for (int p = 0; p < nR13; p++) {  // r13 integration
						long double r13 = r13Array[p];
						long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0*r1*r3);
						long double Sin13 = sqrt(1.0 - Cos13*Cos13);
						long double rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13*r13);
						long double j1rhop = sf_bessel_j1(kappa*rhop);
						long double n1rhop = sf_bessel_n1(kappa*rhop);
						long double ExpR13R2 = exp(-(r13/2.0 + r2));
						long double S23 =  ExpR13R2 * SqrtKappa * j1rhop;
						long double C23 = -ExpR13R2 * SqrtKappa * n1rhop * fshielding(rhop, mu);
						long double dTau = r1 * r2 * r23 * r13;

						long double Phi12SumCLC = 0.0, Phi12SumCLS = 0.0, Phi12SumSLC = 0.0, Phi12SumSLS = 0.0;
						for (int m = 1; m <= nPhi12; m++) {  // phi_12 integration
							long double Phi12 = (2.0*m - 1.0)*PI/(2.0*nPhi12);
							long double r12 = sqrt(r1*r1 + r2*r2 - 2.0*r1*r2*(Sin13*Sin23*cosl(Phi12) + Cos13*Cos23));
							long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0*r1*r2);
							long double rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12*r12);
							long double j1rho = sf_bessel_j1(kappa*rho);
							long double n1rho = sf_bessel_n1(kappa*rho);
							long double ExpR12R3 = exp(-(r12/2.0 + r3));

							long double S22 =  ExpR12R3 * SqrtKappa * j1rho;
							long double C22 = -ExpR12R3 * SqrtKappa * n1rho * fshielding(rho, mu);
							long double Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0 * rho * rhop);
							//@TODO: Is Cos23 ever needed other than above?  If so, we can use this definition.
							//long double Ang = (rho*rho + rhop*rhop - 0.25*r23*r23) / (2.0 * rho * rhop);

							//long double RetSLS = S22 * S23 * PotP * Ang;
							//long double RetSLS = S23 * S22 * Pot * Ang;
							long double RetSLS = S23 * S23 * PotP * 2.0;
							//long double RetSLS = S22 * S22 * Pot * 2.0;
							Phi12SumSLS += sf * ExpR1R2R3 * RetSLS * dTau;
							//Phi12SumSLS += sf * expl(r1+r2) * RetSLS * dTau;

							long double RetCLS = C23 * S22 * Pot * Ang;
							//long double RetCLS = C22 * S23 * Pot * Ang;
							//long double RetCLS = C23 * S23 * Pot * 2.0;
							//long double RetCLS = C22 * S22 * Pot * 2.0;
							Phi12SumCLS += sf * ExpR1R2R3 * RetCLS * dTau;




							////long double RetSLC = S23 * ExpR1R2R3 * C22 * Pot * Ang;
							////Phi12SumSLC += sf * RetSLC * dTau;
							////long double RetSLC = ExpR1R2R3 * C22 * Pot * Ang;
							////Phi12SumSLC += (S22 + sf * S23) * RetSLC * dTau;
							//long double RetSLC = C22 * Pot;
							//Phi12SumSLC += (S22 * 2.0 + sf * S23 * Ang) * RetSLC * ExpR1R2R3 * dTau;
							////if (!IsFiniteNumber(Phi12SumSLC)) {
							////	Phi12SumSLC = Phi12SumSLC;
							////}
							////long double RetSLC = C22 * Pot;
							////Phi12SumSLC += (S22 * 2.0) * RetSLC * ExpR1R2R3 * dTau;


							//long double RetSLC1 = -2.0 * SqrtKappa * ExpR12R3 * Pot * n1rho * fshielding(rho, mu);
							//long double RetSLC2 = -Ang * SqrtKappa * ExpR13R2 * Pot * n1rhop * fshielding(rhop, mu);
							//Phi12SumSLC += ExpR1R2R3 * S22 * (RetSLC1 + sf * RetSLC2) * dTau;


							//long double RetSLCBad = S22 * Pot;
							////Phi12SumSLC += (S22 * 2.0 + sf * S23 * Ang) * RetSLCBad * ExpR1R2R3 * dTau;
							////Phi12SumSLC += (sf * S23 * Ang) * RetSLCBad * ExpR1R2R3 * dTau;
							//Phi12SumSLC += (S22 * 2.0) * RetSLCBad * ExpR1R2R3 * dTau;
							
							////long double RetSLSBad = 2.0 * S22 * S22 + Ang * S23 * S22 + 2.0 * S23 * S23 + Ang * S22 * S23;
							////long double RetSLSBad = 2.0 * S22 * S22 + Ang * S23 * S22;
							////long double RetSLSBad = 2.0 * S23 * S23 + Ang * S22 * S23;
							//long double RetSLSBad = 2.0 * S22 * S22;
							//Phi12SumSLC += Pot * RetSLSBad * ExpR1R2R3 * dTau;

							long double RetSLC1 = -Ang * SqrtKappa * ExpR12R3 * Pot * n1rho * fshielding(rho, mu);
							long double RetSLC2 = -2.0 * SqrtKappa * ExpR13R2 * Pot * n1rhop * fshielding(rhop, mu);
							Phi12SumSLC += ExpR1R2R3 * S23 * (RetSLC1 + sf * RetSLC2) * dTau;



							long double RetCLC = C22 * C23 * Pot * Ang;
							Phi12SumCLC += sf * ExpR1R2R3 * RetCLC * dTau;
						}
						r13SumCLC += LegendreWeightsR13[p] * Phi12SumCLC / nPhi12 * (b13-a13)/2.0;
						r13SumSLC += LegendreWeightsR13[p] * Phi12SumSLC / nPhi12 * (b13-a13)/2.0;
						r13SumCLS += LegendreWeightsR13[p] * Phi12SumCLS / nPhi12 * (b13-a13)/2.0;
						r13SumSLS += LegendreWeightsR13[p] * Phi12SumSLS / nPhi12 * (b13-a13)/2.0;
					}
					r23SumCLC += LegendreWeightsR23[k] * r13SumCLC * (b23-a23)/2.0;
					r23SumSLC += LegendreWeightsR23[k] * r13SumSLC * (b23-a23)/2.0;
					r23SumCLS += LegendreWeightsR23[k] * r13SumCLS * (b23-a23)/2.0;
					r23SumSLS += LegendreWeightsR23[k] * r13SumSLS * (b23-a23)/2.0;
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

	r1SumCLC /= 2.0;
	r1SumSLC /= 2.0;
	r1SumCLS /= 2.0;
	r1SumSLS /= 2.0;

	CLC = r1SumCLC;
	SLC = r1SumSLC;
	CLS = r1SumCLS;
	SLS = r1SumSLS;

	return;
}


void GaussIntegrationPhi13_LongLong_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nPhi13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int sf, double &CLC, double &SLC, double &CLS, double &SLS)
{
	double r1, r2, r3;
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

	long double r1SumCLC = 0.0, r1SumCLS = 0.0, r1SumSLC = 0.0, r1SumSLS = 0.0;
	//#pragma omp parallel for shared(r1Sum,r1Abscissas,r1Weights,Powers) private(r1,r2,r3,r13,r23,r2Sum,r3Sum,r23sum,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points)
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
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * exp(-r2Abscissas[m]) * (r1-0.0)/2.0;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = LaguerreAbscissasR2[m-nR2Leg] + r1;
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-r1);
			}
		}

		long double r2SumCLC = 0.0, r2SumCLS = 0.0, r2SumSLC = 0.0, r2SumSLS = 0.0;
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
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0, r2);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * exp(-r3Abscissas[m]) * (r2-0.0)/2.0;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = LaguerreAbscissasR3[m-nR3Leg] + r2;
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-r2);
				}
			}

			long double r3SumCLC = 0.0, r3SumCLS = 0.0, r3SumSLC = 0.0, r3SumSLS = 0.0;
			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				r3 = r3Abscissas[g];
				long double TempCoeff = r3Weights[g] * r2Weights[j] * r1Weights[i];
				long double ExpR1R2R3 = expl(r1+r2+r3);
				long double a23 = fabs(r2-r3);
				long double b23 = fabs(r2+r3);
				ChangeOfIntervalNoResize(LegendreAbscissasR23, r23Array, a23, b23);

				long double r23SumCLC = 0.0, r23SumCLS = 0.0, r23SumSLC = 0.0, r23SumSLS = 0.0;
				for (int k = 0; k < nR23; k++) {  // r23 integration
					long double r23 = r23Array[k];
					long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0*r2*r3);
					long double Sin23 = sqrt(1.0 - Cos23*Cos23);
					long double Pot = 2.0/r23;

					long double r12SumCLC = 0.0, r12SumCLS = 0.0, r12SumSLC = 0.0, r12SumSLS = 0.0;
					for (int p = 0; p < nR12; p++) {  // r12 integration
						long double r12 = r12Array[p];
						long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0*r1*r2);
						long double Sin12 = sqrt(1.0 - Cos12*Cos12);
						long double rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12*r12);
						long double j1rho = sf_bessel_j1(kappa*rho);
						long double n1rho = sf_bessel_n1(kappa*rho);
						long double ExpR12R3 = exp(-(r12/2.0 + r3));
						long double S22 =  ExpR12R3 * SqrtKappa * j1rho;
						long double C22 = -ExpR12R3 * SqrtKappa * n1rho * fshielding(rho, mu);
						long double dTau = r1 * r3 * r23 * r12;

						long double Phi13SumCLC = 0.0, Phi13SumCLS = 0.0, Phi13SumSLC = 0.0, Phi13SumSLS = 0.0;
						for (int m = 1; m <= nPhi13; m++) {  // phi_13 integration
							long double Phi13 = (2.0*m - 1.0)*PI/(2.0*nPhi13);
							long double r13 = sqrt(r1*r1 + r3*r3 - 2.0*r1*r3*(Sin12*Sin23*cosl(Phi13) + Cos12*Cos23));
							long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0*r1*r3);
							long double rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13*r13);
							long double j1rhop = sf_bessel_j1(kappa*rhop);
							long double n1rhop = sf_bessel_n1(kappa*rhop);
							long double ExpR13R2 = exp(-(r13/2.0 + r2));

							long double S23 =  ExpR13R2 * SqrtKappa * j1rhop;
							long double C23 = -ExpR13R2 * SqrtKappa * n1rhop * fshielding(rhop, mu);
							long double Ang = (r1*r1 + r1*r2*Cos12 + r1*r3*Cos13 + r2*r3*Cos23) / (2.0 * rho * rhop);
							//@TODO: Is Cos23 ever needed other than above?  If so, we can use this definition.
							//long double Ang = (rho*rho + rhop*rhop - 0.25*r23*r23) / (2.0 * rho * rhop);

							long double PotP = Pot;  //@TODO: Remove
							//long double RetSLS = S22 * S23 * PotP * Ang;
							long double RetSLS = S23 * S22 * Pot * Ang;
							//long double RetSLS = S23 * S23 * PotP * 2.0;
							//long double RetSLS = S22 * S22 * Pot * 2.0;
							Phi13SumSLS += sf * ExpR1R2R3 * RetSLS * dTau;
							//Phi12SumSLS += sf * expl(r1+r2) * RetSLS * dTau;

							long double RetCLS = C23 * S22 * Pot * Ang;
							//long double RetCLS = C22 * S23 * Pot * Ang;
							//long double RetCLS = C23 * S23 * Pot * 2.0;
							//long double RetCLS = C22 * S22 * Pot * 2.0;
							Phi13SumCLS += sf * ExpR1R2R3 * RetCLS * dTau;




							////long double RetSLC = S23 * ExpR1R2R3 * C22 * Pot * Ang;
							////Phi12SumSLC += sf * RetSLC * dTau;
							////long double RetSLC = ExpR1R2R3 * C22 * Pot * Ang;
							////Phi12SumSLC += (S22 + sf * S23) * RetSLC * dTau;
							//long double RetSLC = C22 * Pot;
							//Phi12SumSLC += (S22 * 2.0 + sf * S23 * Ang) * RetSLC * ExpR1R2R3 * dTau;
							////if (!IsFiniteNumber(Phi12SumSLC)) {
							////	Phi12SumSLC = Phi12SumSLC;
							////}
							////long double RetSLC = C22 * Pot;
							////Phi12SumSLC += (S22 * 2.0) * RetSLC * ExpR1R2R3 * dTau;


							//long double RetSLC1 = -2.0 * SqrtKappa * ExpR12R3 * Pot * n1rho * fshielding(rho, mu);
							//long double RetSLC2 = -Ang * SqrtKappa * ExpR13R2 * Pot * n1rhop * fshielding(rhop, mu);
							//Phi12SumSLC += ExpR1R2R3 * S22 * (RetSLC1 + sf * RetSLC2) * dTau;


							//long double RetSLCBad = S22 * Pot;
							////Phi12SumSLC += (S22 * 2.0 + sf * S23 * Ang) * RetSLCBad * ExpR1R2R3 * dTau;
							////Phi12SumSLC += (sf * S23 * Ang) * RetSLCBad * ExpR1R2R3 * dTau;
							//Phi12SumSLC += (S22 * 2.0) * RetSLCBad * ExpR1R2R3 * dTau;
							
							////long double RetSLSBad = 2.0 * S22 * S22 + Ang * S23 * S22 + 2.0 * S23 * S23 + Ang * S22 * S23;
							////long double RetSLSBad = 2.0 * S22 * S22 + Ang * S23 * S22;
							////long double RetSLSBad = 2.0 * S23 * S23 + Ang * S22 * S23;
							//long double RetSLSBad = 2.0 * S22 * S22;
							//Phi12SumSLC += Pot * RetSLSBad * ExpR1R2R3 * dTau;

							// THIS ONE WORKS FOR PHI12.
							//long double RetSLC1 = -Ang * SqrtKappa * ExpR12R3 * Pot * n1rho * fshielding(rho, mu);
							//long double RetSLC2 = -2.0 * SqrtKappa * ExpR13R2 * Pot * n1rhop * fshielding(rhop, mu);
							//Phi13SumSLC += ExpR1R2R3 * S23 * (RetSLC1 + sf * RetSLC2) * dTau;

							// THIS ONE WORKS FOR PHI13.
							long double RetSLC1 = -2.0 * SqrtKappa * ExpR12R3 * Pot * n1rho * fshielding(rho, mu);
							long double RetSLC2 = -Ang * SqrtKappa * ExpR13R2 * Pot * n1rhop * fshielding(rhop, mu);
							Phi13SumSLC += ExpR1R2R3 * S22 * (RetSLC1 + sf * RetSLC2) * dTau;




							long double RetCLC = C22 * C23 * Pot * Ang;
							Phi13SumCLC += sf * ExpR1R2R3 * RetCLC * dTau;
						}
						r12SumCLC += LegendreWeightsR12[p] * Phi13SumCLC / nPhi13 * (b12-a12)/2.0;
						r12SumSLC += LegendreWeightsR12[p] * Phi13SumSLC / nPhi13 * (b12-a12)/2.0;
						r12SumCLS += LegendreWeightsR12[p] * Phi13SumCLS / nPhi13 * (b12-a12)/2.0;
						r12SumSLS += LegendreWeightsR12[p] * Phi13SumSLS / nPhi13 * (b12-a12)/2.0;
					}
					r23SumCLC += LegendreWeightsR23[k] * r12SumCLC * (b23-a23)/2.0;
					r23SumSLC += LegendreWeightsR23[k] * r12SumSLC * (b23-a23)/2.0;
					r23SumCLS += LegendreWeightsR23[k] * r12SumCLS * (b23-a23)/2.0;
					r23SumSLS += LegendreWeightsR23[k] * r12SumSLS * (b23-a23)/2.0;
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

	r1SumCLC /= 2.0;
	r1SumSLC /= 2.0;
	r1SumCLS /= 2.0;
	r1SumSLS /= 2.0;

	CLC += r1SumCLC;
	SLC += r1SumSLC;
	CLS += r1SumCLS;
	SLS += r1SumSLS;

	return;
}
