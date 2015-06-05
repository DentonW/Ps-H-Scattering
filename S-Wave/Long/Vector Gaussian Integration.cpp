#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <omp.h>
//#include <gmpxx.h>
#include "Ps-H Scattering.h"
using namespace std;

extern double PI;


double ipow1(double a, int ex)
// Return a**ex
{
	if ( 0==ex )  return 1.0;
	else
	{
		double z = a;
		double y = 1.0;
		while ( 1 )
		{
			if ( ex & 1 )  y *= z;
			ex /= 2;
			if ( 0==ex )  break;
			z *= z;
		}
		return y;
	}
}


// Assumes that the LUT has already been allocated properly.
void CreateRPowerLUT(vector <long double> &LUT, long double r, int Omega)
{
	LUT[0] = 1.0;
	for (int i = 1; i <= Omega; i++) {
		LUT[i] = LUT[i-1] * r;
	}
}


struct KahanAccumulation
{
	long double sum;
	long double correction;
};

KahanAccumulation KahanSum(KahanAccumulation accumulation, long double value)
{
	KahanAccumulation result;
	long double y = value - accumulation.correction;
	long double t = accumulation.sum + y;
	result.correction = (t - accumulation.sum) - y;
	result.sum = t;
	return result;
}

void VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int sf, int NumPowers, vector <rPowers> &Powers, int Omega)
{
	double r1, r2, r3;
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR12, LegendreWeightsR12;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	vector <long double> r12, r13;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r2Weights, r3Abscissas, r3Weights;
	int NumR2Points, NumR3Points, Prog = 0;
	vector <double> r1Pow(Omega+1), r2Pow(Omega+1), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);
	double SqrtKappa = sqrt(kappa), MuCubed = mu*mu*mu;

	WriteProgress(string("PhiLS and PhiLC starting"), Prog, nR1);
	Prog = 0;

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(r1Abscissas, r1Weights, nR1);
	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR12, LegendreWeightsR12, nR12);
	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);

	r1Abscissas.resize(nR1);
	r1Weights.resize(nR1);
	for (int m = 0; m < nR1; m++) {  //@TODO: Can we just do this as a single line a la Fortran?
		r1Abscissas[m] = r1Abscissas[m] / Powers[0].alpha;
		r1Weights[m] = r1Weights[m] / Powers[0].alpha;
	}

	//KahanAccumulation KATest = { 0.0, 0.0 };

	//mpf_set_default_prec(1000);
	//vector <mpf_class> AResultsMP(NumPowers), BResultsMP(NumPowers);
	//for (int mp = 0; mp < NumPowers; mp++) {
	//	AResultsMP[mp] = 0;
	//	BResultsMP[mp] = 0;
	//}

	#pragma omp parallel for shared(/*KATest,*/r1Abscissas,r1Weights,Powers) private(r1,r2,r3,r12,r13,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		vector <long double> TempAResults(NumPowers, 0.0), TempBResults(NumPowers, 0.0);
		//vector <mpf_class> TempAResultsMP(NumPowers), TempBResultsMP(NumPowers);
		vector <long double> r1Pow(Omega+1), r2Pow(Omega+1), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);
		
		KahanAccumulation KA = { 0.0, 0.0 };
		//for (int mp = 0; mp < NumPowers; mp++) {
		//	TempAResultsMP[mp] = 0;
		//	TempBResultsMP[mp] = 0;
		//}

		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);
		r12.resize(nR12);
		r13.resize(nR13);

		r1 = r1Abscissas[i];
		CreateRPowerLUT(r1Pow, r1, Omega);

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR2) {
			for (int m = 0; m < nR2Lag; m++) {
				r2Weights[m] = LaguerreWeightsR2[m] / Powers[0].beta;
				r2Abscissas[m] = LaguerreAbscissasR2[m] / Powers[0].beta;
			}
			NumR2Points = nR2Lag;
		}
		// Do Gauss-Legendre for r2 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * exp(-Powers[0].beta*r2Abscissas[m]) * (r1-0.0)/2.0 /** 4.0*/;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r1*Powers[0].beta) / Powers[0].beta;
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-Powers[0].beta*r1) / Powers[0].beta;
			}
		}

		for (int j = 0; j < NumR2Points; j++) {  // r2 integration
			r2 = r2Abscissas[j];
			CreateRPowerLUT(r2Pow, r2, Omega);

			// Cusp unimportant, so only do Gauss-Laguerre integration
			if (r1 > CuspR3) {
				for (int m = 0; m < nR3Lag; m++) {
					r3Weights[m] = LaguerreWeightsR3[m] / Powers[0].gamma;
					r3Abscissas[m] = LaguerreAbscissasR3[m] / Powers[0].gamma;
				}
				NumR3Points = nR3Lag;
			}
			// Do Gauss-Legendre for r3 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0, r1);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * exp(-Powers[0].gamma*r3Abscissas[m]) * (r1-0.0)/2.0 /** 4.0*/;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r1*Powers[0].gamma) / Powers[0].gamma;
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-Powers[0].gamma*r1) / Powers[0].gamma;
				}
			}

			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				r3 = r3Abscissas[g];
				CreateRPowerLUT(r3Pow, r3, Omega);
				double Coeff = r3Weights[g] * r2Weights[j] * r1Weights[i] * SqrtKappa * 0.70710678118654752440 * PI/nPhi23;
				// Integrations over the 3 interparticle coordinates
				//Vecr12r13phi23Integration(Results, Coeff, NumPowers, Powers, r1, r2, r3, kappa, mu, sf, LegendreAbscissasR12, LegendreWeightsR12, LegendreAbscissasR13, LegendreWeightsR13, nR12, nR13, nPhi23, f1, f2);


				long double a12, b12, a13, b13;
				long double r13Sum;
				long double Phi23, Cos12, Sin12, Cos13, Sin13, rho, rhop;
				long double fOuterA, fOuterB;
				long double CoeffNew;

				a12 = fabs(r1-r2);
				b12 = fabs(r1+r2);
				a13 = fabs(r1-r3);
				b13 = fabs(r1+r3);
				//ChangeOfInterval(LegendreAbscissasR12, r12, a12, b12);
				//ChangeOfInterval(LegendreAbscissasR13, r13, a13, b13);
				ChangeOfIntervalNoResize(LegendreAbscissasR12, r12, a12, b12);
				ChangeOfIntervalNoResize(LegendreAbscissasR13, r13, a13, b13);

				for (int k = 0; k < nR12; k++) {  // r12 integration
					CreateRPowerLUT(r12Pow, r12[k], Omega);
					r13Sum = 0.0;
					Cos12 = (r1*r1 + r2*r2 - r12[k]*r12[k]) / (2.0*r1*r2);
					Sin12 = sqrt(1.0 - Cos12*Cos12);
					rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12[k]*r12[k]);
					long double ExpMuRho = exp(-(mu*rho));
					long double ExpR12R3 = exp(-(r12[k]/2.0 + r3));
					long double SinKRho, CosKRho;
					long double Pot2 = (2.0/r1 - 2.0/r3 - 2.0/r12[k]);
					SinKRho = sin(kappa*rho)/(kappa*rho);  CosKRho = cos(kappa*rho)/(kappa*rho);

					for (int p = 0; p < nR13; p++) {  // r13 integration
						CreateRPowerLUT(r13Pow, r13[p], Omega);
						Cos13 = (r1*r1 + r3*r3 - r13[p]*r13[p]) / (2.0*r1*r3);
						Sin13 = sqrt(1.0 - Cos13*Cos13);
						rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13[p]*r13[p]);
						CoeffNew = Coeff * LegendreWeightsR13[p] * (b13-a13) * LegendreWeightsR12[k] * (b12-a12) * r2 * r3 * r12[k] * r13[p];

						Phi23 = (2.0*1.0 - 1.0)*PI/(2.0*nPhi23);
						//r23 = sqrt(r2*r2 + r3*r3 - 2.0*r2*r3*(Sin12*Sin13*cos(Phi23) + Cos12*Cos13));
						//CreateRPowerLUT(r23Pow, r23, Omega);


						// PhiLC part (LOnCBar)
						long double CRet1, CRet2;
						long double ExpMuRhop = exp(-(mu*rhop));
						long double ExpR13R2 = exp(-(r13[p]/2.0 + r2));
						long double Pot1 = (2.0/r1 - 2.0/r2 - 2.0/r13[p]);
						long double SinKRhop, CosKRhop;
						SinKRhop = sin(kappa*rhop)/(kappa*rhop);  CosKRhop = cos(kappa*rhop)/(kappa*rhop);
	
						CRet1 = kappa*mu/2.0*ExpMuRho*(1.0+mu*rho)*SinKRho;
						CRet1 += MuCubed*rho/4.0*ExpMuRho*CosKRho;
						CRet1 += Pot1*CosKRho*(1.0-ExpMuRho*(1.0+mu*rho/2.0));
						CRet1 *= ExpR12R3;

						CRet2 = kappa*mu/2.0*ExpMuRhop*(1.0+mu*rhop)*SinKRhop;
						CRet2 += MuCubed*rhop/4.0*ExpMuRhop*CosKRhop;
						CRet2 += Pot2*CosKRhop*(1.0-ExpMuRhop*(1.0+mu*rhop/2.0));
						CRet2 *= ExpR13R2;

						fOuterA = (CRet1 + sf*CRet2) * CoeffNew;


						// PhiLS part (LOnSBar)
						long double SRet1 = Pot1 * ExpR12R3 * SinKRho;
						long double SRet2 = Pot2 * ExpR13R2 * SinKRhop;
						fOuterB = (SRet1 + sf*SRet2) * CoeffNew;


						for (int n = 0; n < NumPowers; n++) {
							rPowers *rp = &Powers[n];
							long double Phi1 = r1Pow[rp->ki] * r2Pow[rp->li] * r12Pow[rp->mi] * r3Pow[rp->ni] * r13Pow[rp->pi]; // * r23Pow[rp->qi];
							TempAResults[n] += Phi1 * fOuterA;
							TempBResults[n] += Phi1 * fOuterB;
							//if (n == 0)
							//	KA = KahanSum(KA, Phi1 * fOuterA);
							//TempAResultsMP[n] = TempAResultsMP[n] + double(Phi1 * fOuterA);
							//TempBResultsMP[n] = TempBResultsMP[n] + double(Phi1 * fOuterB);
						}
					}
				}

			}
		}

		//#pragma omp critical(build)
		for (int n = 0; n < NumPowers; n++) {
			AResults[n] += TempAResults[n];
			BResults[n] += TempBResults[n];

			//AResultsMP[n] += TempAResultsMP[n];
			//BResultsMP[n] += TempBResultsMP[n];

			//if (n == 0) {
			//	cout << "Direct vs. Kahan sums: " << TempAResults[n] << " " << KA.sum << endl;
			//	KATest = KahanSum(KATest, KA.sum);
			//	cout << AResults[0] << " " << KATest.sum << endl;
			//}

			//cout << "Double vs. MP sums (A): " << TempAResults[n] << " " << TempAResultsMP[n] << endl;
			//cout << "Double vs. MP sums (B): " << TempBResults[n] << " " << TempBResultsMP[n] << endl;
		}
		WriteProgress(string("PhiLS and PhiLC"), Prog, nR1);
	}

	//cout << "Results:" << endl;
	//for (int n = 0; n < NumPowers; n++) {
	//	cout << n << " AResults: " << AResults[n] << " " << AResultsMP[n] << endl;
	//	cout << n << " BResults: " << BResults[n] << " " << BResultsMP[n] << endl;
	//}

	return;
}


void VecGaussIntegrationPhi13_PhiLCBar_PhiLSBar_R23Term(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nPhi13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int sf, int NumPowers, vector <rPowers> &Powers, int Omega)
{
	double r1, r2, r3;
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR12, LegendreWeightsR12;
	vector <long double> LegendreAbscissasR23, LegendreWeightsR23;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r3Abscissas, r2Weights, r3Weights;
	int NumR2Points, NumR3Points, Prog = 0;
	double SqrtKappa = sqrt(kappa);

	WriteProgress(string("PhiLS and PhiLC R23 starting"), Prog, nR1);
	Prog = 0;

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(r1Abscissas, r1Weights, nR1);
	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR12, LegendreWeightsR12, nR12);
	GaussLegendre(LegendreAbscissasR23, LegendreWeightsR23, nR23);

	r1Abscissas.resize(nR1);
	r1Weights.resize(nR1);
	for (int m = 0; m < nR1; m++) {  //@TODO: Can we just do this as a single line a la Fortran?
		r1Abscissas[m] = r1Abscissas[m] / Powers[0].alpha;
		r1Weights[m] = r1Weights[m] / Powers[0].alpha;
	}

	#pragma omp parallel for shared(r1Abscissas,r1Weights,Powers) private(r1,r2,r3,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		vector <long double> TempAResults(NumPowers, 0.0), TempBResults(NumPowers, 0.0);
		vector <long double> r1Pow(Omega+1), r2Pow(Omega+1), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);

		//@TODO: Could move these further in
		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);

		r1 = r1Abscissas[i];
		CreateRPowerLUT(r1Pow, r1, Omega);

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR2) {
			for (int m = 0; m < nR2Lag; m++) {
				r2Weights[m] = LaguerreWeightsR2[m] / Powers[0].beta;
				r2Abscissas[m] = LaguerreAbscissasR2[m] / Powers[0].beta;
			}
			NumR2Points = nR2Lag;
		}
		// Do Gauss-Legendre for r2 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * exp(-Powers[0].beta*r2Abscissas[m]) * (r1-0.0)/2.0 /** 4.0*/;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r1*Powers[0].beta) / Powers[0].beta;
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-Powers[0].beta*r1) / Powers[0].beta;
			}
		}

		for (int j = 0; j < NumR2Points; j++) {  // r2 integration
			r2 = r2Abscissas[j];
			CreateRPowerLUT(r2Pow, r2, Omega);

			if (r2 > CuspR3) {
				for (int m = 0; m < nR3Lag; m++) {
					r3Weights[m] = LaguerreWeightsR3[m] / Powers[0].gamma;
					r3Abscissas[m] = LaguerreAbscissasR3[m] / Powers[0].gamma;
				}
				NumR3Points = nR3Lag;
			}
			// Do Gauss-Legendre for r3 integration over [0,r2] and then Gauss-Laguerre over [r2,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0, r2);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * exp(-Powers[0].gamma*r3Abscissas[m]) * (r2-0.0)/2.0 /** 4.0*/;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r2*Powers[0].gamma) / Powers[0].gamma;
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-Powers[0].gamma*r2) / Powers[0].gamma;
				}
			}

			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				r3 = r3Abscissas[g];
				CreateRPowerLUT(r3Pow, r3, Omega);
				double Coeff = r3Weights[g] * r2Weights[j] * r1Weights[i];
				// Integrations over the 3 interparticle coordinates
				//Vecr23r12phi13Integration(Results, Coeff, NumPowers, Powers, r1, r2, r3, kappa, mu, sf, LegendreAbscissasR12, LegendreWeightsR12, LegendreAbscissasR23, LegendreWeightsR23, nR12, nPhi13, nR23, f1, f2);

				long double a23, b23, a12, b12;
				vector <long double> r23(nR23), r12(nR12);
				long double r13, Phi13, Cos23, Sin23, Cos12, Sin12, rho, rhop;
				long double fOuterA, fOuterB;

				a23 = fabs(r2-r3);
				b23 = fabs(r2+r3);
				a12 = fabs(r1-r2);
				b12 = fabs(r1+r2);
				ChangeOfIntervalNoResize(LegendreAbscissasR23, r23, a23, b23);
				ChangeOfIntervalNoResize(LegendreAbscissasR12, r12, a12, b12);

				Coeff *= r1 * r3 * PI/nPhi13 * (b12-a12) * (b23-a23) * 0.70710678118654752440 * SqrtKappa;

				for (int k = 0; k < nR23; k++) {  // r23 integration
					//CreateRPowerLUT(r23Pow, r23[k], Omega);
					Cos23 = (r2*r2 + r3*r3 - r23[k]*r23[k]) / (2.0*r2*r3);
					Sin23 = sqrt(1.0 - Cos23*Cos23);
					for (int p = 0; p < nR12; p++) {  // r12 integration
						CreateRPowerLUT(r12Pow, r12[p], Omega);
						Cos12 = (r1*r1 + r2*r2 - r12[p]*r12[p]) / (2.0*r1*r2);
						Sin12 = sqrt(1.0 - Cos12*Cos12);
						rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12[p]*r12[p]);
						long double CoeffNew = Coeff * /*r23[k] **/ r12[p] * LegendreWeightsR12[p] * LegendreWeightsR23[k] * 2.0;  // All terms are multipled by 2/r23, so cancels with r23 here.
						long double ExpR12R3 = exp(-(r12[p]/2.0 + r3));
						long double Ret1S = ExpR12R3 * sin(kappa*rho) / (kappa*rho);
						long double Ret1C = ExpR12R3 * cos(kappa*rho)/(kappa*rho)*(1.0-exp(-(mu*rho))*(1.0+mu*rho/2.0));

						for (int m = 1; m <= nPhi13; m++) {  // phi_13 integration
							Phi13 = (2.0*m - 1.0)*PI/(2.0*nPhi13);
							r13 = sqrt(r1*r1 + r3*r3 - 2.0*r1*r3*(Sin12*Sin23*cos(Phi13) + Cos12*Cos23));
							CreateRPowerLUT(r13Pow, r13, Omega);
							rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13*r13);
							long double ExpR13R2 = exp(-(r13/2.0 + r2));

							//
							// PhiLC part (LOnCBar_R23Term)
							//
							long double Ret2C = ExpR13R2 * cos(kappa*rhop)/(kappa*rhop)*(1.0-exp(-mu*rhop)*(1.0+mu*rhop/2.0));
							fOuterA = (Ret1C + sf*Ret2C) * CoeffNew;


							//
							// PhiLS part (LOnSBar_R23Term)
							//
							long double Ret2S = ExpR13R2 * sin(kappa*rhop) / (kappa*rhop);
							fOuterB = (Ret1S + sf*Ret2S) * CoeffNew;


							for (int n = 0; n < NumPowers; n++) {
								rPowers *rp = &Powers[n];
								long double Phi1 = r1Pow[rp->ki] * r2Pow[rp->li] * r12Pow[rp->mi] * r3Pow[rp->ni] * r13Pow[rp->pi]; // * r23Pow[rp->qi];
								TempAResults[n] += Phi1 * fOuterA;
								TempBResults[n] += Phi1 * fOuterB;
								//#pragma omp atomic
								//Results[n] += Phi1 * fOuter * CoeffNew;
							}
						}
					}
				}
			}
		}

		//#pragma omp critical(build)
		for (int n = 0; n < NumPowers; n++) {
			AResults[n] += TempAResults[n];
			BResults[n] += TempBResults[n];
		}
		WriteProgress(string("PhiLS and PhiLC R23"), Prog, nR1);
	}

	return;
}


void VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar_Full(vector <double> &AResults, vector <double> &BResults, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int sf, int NumPowers, vector <rPowers> &Powers, int Omega)
{
	double r1, r2, r3;

	//vector <double> r12, r13;
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR12, LegendreWeightsR12;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r2Weights, r3Abscissas, r3Weights;
	vector <long double> r12, r13;
	double SqrtKappa = sqrt(kappa), MuCubed = mu*mu*mu;
	int NumR2Points, NumR3Points, Prog = 0;

	WriteProgress(string("PhiLS and PhiLC Full starting"), Prog, nR1);
	Prog = 0;

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(r1Abscissas, r1Weights, nR1);
	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR12, LegendreWeightsR12, nR12);
	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);

	r1Abscissas.resize(nR1);
	r1Weights.resize(nR1);
	for (int m = 0; m < nR1; m++) {  //@TODO: Can we just do this as a single line a la Fortran?
		r1Abscissas[m] = r1Abscissas[m] / Powers[0].alpha;
		r1Weights[m] = r1Weights[m] / Powers[0].alpha;
	}

	#pragma omp parallel for shared(r1Abscissas,r1Weights,Powers) private(r1,r2,r3,r12,r13,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		vector <long double> TempAResults(NumPowers, 0.0), TempBResults(NumPowers, 0.0);
		vector <long double> r1Pow(Omega+1), r2Pow(Omega+1), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);

		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);
		r12.resize(nR12);
		r13.resize(nR13);

		r1 = r1Abscissas[i];
		CreateRPowerLUT(r1Pow, r1, Omega);

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR2) {
			for (int m = 0; m < nR2Lag; m++) {
				r2Weights[m] = LaguerreWeightsR2[m] / Powers[0].beta;
				r2Abscissas[m] = LaguerreAbscissasR2[m] / Powers[0].beta;
			}
			NumR2Points = nR2Lag;
		}
		// Do Gauss-Legendre for r2 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * exp(-Powers[0].beta*r2Abscissas[m]) * (r1-0.0)/2.0 /** 4.0*/;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r1*Powers[0].beta) / Powers[0].beta;
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-Powers[0].beta*r1) / Powers[0].beta;
			}
		}

		for (int j = 0; j < NumR2Points; j++) {  // r2 integration
			r2 = r2Abscissas[j];
			CreateRPowerLUT(r2Pow, r2, Omega);

			// Cusp unimportant, so only do Gauss-Laguerre integration
			if (r1 > CuspR3) {
				for (int m = 0; m < nR3Lag; m++) {
					r3Weights[m] = LaguerreWeightsR3[m] / Powers[0].gamma;
					r3Abscissas[m] = LaguerreAbscissasR3[m] / Powers[0].gamma;
				}
				NumR3Points = nR3Lag;
			}
			// Do Gauss-Legendre for r3 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0, r1);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * exp(-Powers[0].gamma*r3Abscissas[m]) * (r1-0.0)/2.0 /** 4.0*/;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r1*Powers[0].gamma) / Powers[0].gamma;
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-Powers[0].gamma*r1) / Powers[0].gamma;
				}
			}

			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				r3 = r3Abscissas[g];
				CreateRPowerLUT(r3Pow, r3, Omega);
				long double Coeff = r3Weights[g] * r2Weights[j] * r1Weights[i] * SqrtKappa;
				// Integrations over the 3 interparticle coordinates
				//Vecr12r13phi23Integration(Results, Coeff, NumPowers, Powers, r1, r2, r3, kappa, mu, sf, LegendreAbscissasR12, LegendreWeightsR12, LegendreAbscissasR13, LegendreWeightsR13, nR12, nR13, nPhi23, f1, f2);


				long double a12, b12, a13, b13;
				long double r13Sum;
				long double r23, Phi23, Cos12, Sin12, Cos13, Sin13, rho, rhop;
				long double fOuterA, fOuterB;
				long double CoeffNew;
				long double SinKRho, CosKRho, SinKRhop, CosKRhop;
				long double ExpMuRho, ExpMuRhop;
				long double kr, krp;

				a12 = fabs(r1-r2);
				b12 = fabs(r1+r2);
				a13 = fabs(r1-r3);
				b13 = fabs(r1+r3);
				ChangeOfIntervalNoResize(LegendreAbscissasR12, r12, a12, b12);
				ChangeOfIntervalNoResize(LegendreAbscissasR13, r13, a13, b13);

				//Coeff = Coeff * (b12-a12)/2.0 / 4.0 * r2 * r3 * (b13-a13)/2.0 * PI/nPhi * sqrt(2.0) * 16.0 / 2.0;

				for (int k = 0; k < nR12; k++) {  // r12 integration
					CreateRPowerLUT(r12Pow, r12[k], Omega);
					r13Sum = 0.0;
					Cos12 = (r1*r1 + r2*r2 - r12[k]*r12[k]) / (2.0*r1*r2);
					Sin12 = sqrt(1.0 - Cos12*Cos12);
					rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12[k]*r12[k]);
					kr = kappa*rho;
					SinKRho = sin(kr)/kr;  CosKRho = cos(kr)/kr;
					ExpMuRho = exp(-(mu*rho));
					long double ExpR12R3 = exp(-(r12[k]/2.0 + r3));
					long double Ret1S = ExpR12R3 * SinKRho;
					//Coeff12 = Coeff * LegendreWeights[k] * r12[k];

					for (int p = 0; p < nR13; p++) {  // r13 integration
						CreateRPowerLUT(r13Pow, r13[p], Omega);
						Cos13 = (r1*r1 + r3*r3 - r13[p]*r13[p]) / (2.0*r1*r3);
						Sin13 = sqrt(1.0 - Cos13*Cos13);
						rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13[p]*r13[p]);
						krp = kappa*rhop;
						SinKRhop = sin(krp)/krp;  CosKRhop = cos(krp)/krp;  
						ExpMuRhop = exp(-(mu*rhop));
						CoeffNew = Coeff * 0.70710678118654752440 * PI/nPhi23 * LegendreWeightsR13[p] * (b13-a13) * LegendreWeightsR12[k] * (b12-a12) * r2 * r3 * r12[k] * r13[p];
						long double ExpR13R2 = exp(-(r13[p]/2.0 + r2));
						long double Ret2S = ExpR13R2 * SinKRhop;

						for (int m = 1; m <= nPhi23; m++) {  // phi_23 integration
							Phi23 = (2.0*m - 1.0)*PI/(2.0*nPhi23);
							r23 = sqrt(r2*r2 + r3*r3 - 2.0*r2*r3*(Sin12*Sin13*cos(Phi23) + Cos12*Cos13));
							CreateRPowerLUT(r23Pow, r23, Omega);
							//fOuter = f1(Powers[0], r1, r2, r3, r12[k], r13[p], r23, kappa, rho, rhop, mu, sf);  // Just use any Powers (not supposed to be used in f1)
							//fOuter = LOnCBarFull(r1, r2, r3, r12[k], r13[p], r23, kappa, rho, rhop, mu, sf);

							long double Pot1 = (2.0/r1 - 2.0/r2 - 2.0/r13[p] + 2.0/r23);
							long double Pot2 = (2.0/r1 - 2.0/r3 - 2.0/r12[k] + 2.0/r23);

							long double Ret1C, Ret2C;

							//
							// PhiLC part
							//
							Ret1C = kappa*mu/2.0*ExpMuRho*(1.0+mu*rho)*SinKRho;
							Ret1C += MuCubed*rho/4.0*ExpMuRho*CosKRho;
							Ret1C += Pot1*CosKRho*(1.0-ExpMuRho*(1.0+mu*rho/2.0));
							Ret1C *= ExpR12R3;

							Ret2C = kappa*mu/2.0*ExpMuRhop*(1.0+mu*rhop)*SinKRhop;
							Ret2C += MuCubed*rhop/4.0*ExpMuRhop*CosKRhop;
							Ret2C += Pot2*CosKRhop*(1.0-ExpMuRhop*(1.0+mu*rhop/2.0));
							Ret2C *= ExpR13R2;

							fOuterA = (Ret1C + sf*Ret2C) * CoeffNew;


							//
							// PhiLS part
							//
							long double Ret1T = Pot1 * Ret1S;
							long double Ret2T = Pot2 * Ret2S;
							fOuterB = (Ret1T + sf*Ret2T) * CoeffNew;


							for (int n = 0; n < NumPowers; n++) {
								//Results[n] += f2(Powers[n], r1, r2, r3, r12[k], r13[p], r23, kappa, rho, rhop, mu, sf) * fOuter * /*Coeff23*/ CoeffNew;

								rPowers *rp = &Powers[n];
								long double Phi1 = r1Pow[rp->ki] * r2Pow[rp->li] * r12Pow[rp->mi] * r3Pow[rp->ni] * r13Pow[rp->pi] * r23Pow[rp->qi];
								TempAResults[n] += Phi1 * fOuterA;
								TempBResults[n] += Phi1 * fOuterB;
							}
						}
					}
				}

			}
		}

		//#pragma omp critical(build)
		for (int n = 0; n < NumPowers; n++) {
			AResults[n] += TempAResults[n];
			BResults[n] += TempBResults[n];
		}
		WriteProgress(string("PhiLS and PhiLC full"), Prog, nR1);
	}

	return;
}
