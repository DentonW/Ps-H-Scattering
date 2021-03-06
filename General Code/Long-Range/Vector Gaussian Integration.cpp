#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <omp.h>
#include <gsl/gsl_sf_legendre.h>
#include "Ps-H Scattering.h"
using namespace std;

extern long double PI;


long double ipow1(long double a, int ex)
// Return a**ex
{
	if ( 0==ex ) return 1.0L;
	else
	{
		long double z = a;
		long double y = 1.0L;
		while(1)
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
	LUT[0] = 1.0L;
	for (int i = 1; i <= Omega; i++) {
		LUT[i] = LUT[i-1] * r;
	}
}


long double AngR1Rho(int l, long double r1, long double r2, long double Cos12, long double Sin12, long double rho)
{
	long double Ang;

	//switch (l)
	//{
	//	case 0:
	//		Ang = 1.0L;
	//		break;
	//	case 1:
	//		Ang = (r1 + r2 * Cos12) / (2.0L * rho);
	//		break;
	//	case 2:
	//		Ang = 1.0L - 3.0L/8.0L * r2*r2 * Sin12*Sin12 / (rho*rho);
	//		break;
	//	case 3:
	//		Ang = (r1*r1 - r2*r2 + 4.0L * rho*rho) * (5.0L * (r1*r1 - r2*r2)*(r1*r1 - r2*r2) - 8.0L * (r1*r1 + 5.0L * r2*r2) * rho*rho + 80.0L * rho*rho*rho*rho);
	//		Ang = Ang /(128.0L * r1*r1*r1 * rho*rho*rho);
	//		break;
	//	default:
	//		cout << "Internal error" << endl;
	//		exit(6);
	//}

	long double AngCos = (r1 + r2 * Cos12) / (2.0L * rho);
	Ang = gsl_sf_legendre_Plm(l, 0, AngCos);

	return Ang;
}


long double AngR1Rhop(int l, long double r1, long double r3, long double Cos13, long double Sin13, long double rhop)
{
	long double Ang;
	long double A1 = (r3*r3 - 4 * rhop*rhop);

	//switch (l)
	//{
	//	case 0:
	//		Ang = 1.0L;
	//		break;
	//	case 1:
	//		Ang = (r1 + r3 * Cos13) / (2.0L * rhop);
	//		break;
	//	case 2:
	//		Ang = 1.0L - 3.0L/8.0L * r3*r3 * Sin13*Sin13 / (rhop*rhop);
	//		break;
	//	case 3:
	//		Ang = 5.0L * r1*r1*r1*r1*r1*r1 - 5.0L * A1*A1*A1 - 3.0L * r1*r1*r1*r1 * (5.0L * r3*r3 - 4.0L * rhop*rhop) + 3.0L * r1*r1 * (5.0L * r3*r3*r3*r3 - 24.0L * r3*r3 * rhop*rhop + 16.0L * rhop*rhop*rhop*rhop);
	//		Ang = Ang / (128.0L * r1*r1*r1 * rhop*rhop*rhop);
	//		break;
	//	default:
	//		cout << "Internal error" << endl;
	//		exit(6);
	//}

	long double AngCos = (r1 + r3 * Cos13) / (2.0L * rhop);
	Ang = gsl_sf_legendre_Plm(l, 0, AngCos);

	return Ang;
}


long double AngR2Rho(int l, long double r1, long double r2, long double Cos12, long double Sin12, long double rho)
{
	long double Ang;

	//switch (l)
	//{
	//	case 0:
	//		Ang = 1.0L;
	//		break;
	//	case 1:
	//		Ang = (r2 + r1 * Cos12) / (2.0L * rho);
	//		break;
	//	case 2:
	//		Ang = 1.0L - 3.0L/8.0L * r1*r1 * Sin12*Sin12 / (rho*rho);
	//		break;
	//	case 3:
	//		Ang = -5.0L * r1*r1*r1*r1*r1*r1 + 5.0L * r2*r2*r2*r2*r2*r2 + 12.0L * r2*r2*r2*r2 * rho*rho + 48.0L * r2*r2 * rho*rho*rho*rho + 320.0L * rho*rho*rho*rho*rho*rho + 15.0L * r1*r1*r1*r1 * (r2*r2 + 4.0L * rho*rho) - 3.0L * r1*r1 * (5.0L * r2*r2*r2*r2 + 24.0L * r2*r2 * rho*rho + 80.0L * rho*rho*rho*rho);
	//		Ang = Ang / (128.0L * r2*r2*r2 * rho*rho*rho);
	//		break;
	//	default:
	//		cout << "Internal error" << endl;
	//		exit(6);
	//}

	long double AngCos = (r2 + r1 * Cos12) / (2.0L * rho);
	Ang = gsl_sf_legendre_Plm(l, 0, AngCos);

	return Ang;
}


long double AngR2Rhop(int l, long double r1, long double r3, long double Cos12, long double Cos13, long double Cos23, long double rhop)
{
	long double Ang;

	//switch (l)
	//{
	//	case 0:
	//		Ang = 1.0L;
	//		break;
	//	case 1:
	//		Ang = (r1 * Cos12 + r3 * Cos23) / (2.0L * rhop);
	//		break;
	//	case 2:
	//		Ang = r1 * Cos12 + r3 * Cos23;
	//		Ang = 3.0L/8.0L * Ang * Ang / (rhop*rhop) - 0.5L;
	//		break;
	//	case 3:
	//		Ang = 5.0L * Cos12*Cos12 * r1*r1 + 10.0L * Cos12 * Cos23 * r1 * r3 + 5.0L * Cos23*Cos23 * r3*r3 - 12.0L * rhop*rhop;
	//		Ang = Ang * (Cos12 * r1 + Cos23 * r3) / (16.0L * rhop*rhop*rhop);
	//		break;
	//	default:
	//		cout << "Internal error" << endl;
	//		exit(6);
	//}

	long double AngCos = (r1 * Cos12 + r3 * Cos23) / (2.0L * rhop);
	Ang = gsl_sf_legendre_Plm(l, 0, AngCos);

	return Ang;
}


// This is part of the Laplacian (plus kappa^2 term) acting on the parts of C22 depending on rho.
long double LaplacianC(int l, long double kappa, long double rho, long double mu, int shpower)
{
	long double LC, kapparho = kappa*rho;
	long double kapparho2 = kapparho*kapparho;
	long double kapparho3 = kapparho2*kapparho;
	long double kapparho4 = kapparho2*kapparho2;
	long double kapparho5 = kapparho3*kapparho2;
	long double kapparho6 = kapparho3*kapparho3;
	long double kapparho7 = kapparho4*kapparho3;
	long double kapparho8 = kapparho4*kapparho4;
	long double CosKR = cosl(kapparho), SinKR = sinl(kapparho);

	switch (l)
	{
		case 0:  // S-Wave
			LC = -2.0L * kappa * SinKR * fshielding1(rho, mu, shpower) + CosKR * fshielding2(rho, mu, shpower);
			LC /= 2.0L * kapparho;
			break;
		case 1:  // P-Wave
			LC = 2.0L * ((kapparho*kapparho - 1.0L) * CosKR - kapparho * SinKR) * fshielding1(rho, mu, shpower);
			LC += rho * (CosKR + kapparho * SinKR) * fshielding2(rho, mu, shpower);
			LC /= 2.0L * kapparho*kapparho * rho;
			break;
		case 2:  // D-Wave
			//cout << "What about the minus sign?" << endl;
			//cout << "CosKR / SinKR" << endl;
			//LC = fshielding1(rho, mu, shpower) / rho * (2.0L * n2rho - kappa*rho * n1rho) - 0.5L * fshielding2(rho, mu, shpower) * n2rho;
			LC = 2.0L * (3.0L * (kapparho*kapparho - 2.0L) * CosKR + kapparho * (kapparho*kapparho - 6.0L) * SinKR) * fshielding1(rho, mu, shpower);
			LC -= rho * ((kapparho*kapparho - 3.0L) * CosKR - 3.0L * kapparho * SinKR) * fshielding2(rho, mu, shpower);
			LC /= 2.0L * kapparho*kapparho*kapparho * rho;
			break;
		case 3:  // F-Wave
			LC = 2.0L * ((45.0L - 21.0L * kapparho*kapparho + kapparho*kapparho*kapparho*kapparho) * CosKR + 3.0L * kapparho * (15.0L - 2.0L * kapparho*kapparho) * SinKR) * fshielding1(rho, mu, shpower);
			LC += rho * (3.0L * (2.0L * kapparho*kapparho - 5.0L) * CosKR + kapparho * (kapparho*kapparho - 15.0L) * SinKR) * fshielding2(rho, mu, shpower);
			LC /= -2.0L * kapparho*kapparho*kapparho*kapparho * rho;
			break;
		case 4:  // G-Wave
			LC = 2.0L * (5.0L * (84.0L - 39.0L * kapparho2 + 2.0L * kapparho4) * CosKR + kapparho * (420.0L-55.0L * kapparho2+kapparho4) * SinKR) * fshielding1(rho, mu, shpower);
			LC -= rho * ((105.0L - 45.0L * kapparho2 + kapparho4) * CosKR+5.0L * kapparho * (21.0L-2.0L * kapparho2) * SinKR) * fshielding2(rho, mu, shpower);
			LC /= -2.0L * kapparho*kapparho*kapparho*kapparho*kapparho * rho;
			break;
		case 5:  // H-Wave
			LC = 2.0L * ((-4725.0L+2205.0L * kapparho2 - 120.0L * kapparho4+kapparho2*kapparho4) * CosKR-15.0L * kapparho * (315.0L-42.0L * kapparho2+kapparho4) * SinKR) * fshielding1(rho, mu, shpower);
			LC += rho * (15.0L * (63.0L-28.0L * kapparho2+kapparho4) * CosKR+kapparho * (945.0L-105.0L * kapparho2+kapparho4) * SinKR) * fshielding2(rho, mu, shpower);
			LC /= 2.0L * kapparho2*kapparho4 * rho;
			break;
		case 6:  // I-Wave
			/*LC = 2.0L * (21.0L * (-2970.0L+1395.0L * kapparho2-80.0L * kapparho4+kapparho2*kapparho4) * CosKR+kapparho * (-62370.0L+8505.0L * kapparho2-231.0L * kapparho4+kapparho2*kapparho4) * SinKR) * fshielding1(rho, mu, shpower);
			LC += rho * ((10395.0L-4725.0L * kapparho2+210.0L * kapparho4-kapparho2*kapparho4) * CosKR+21.0L * kapparho * (495.0L-60.0L * kapparho2+kapparho4) * SinKR) * fshielding2(rho, mu, shpower);
			LC /= -2.0L * kapparho2*kapparho4*kapparho * rho;*/
			/*LC = 2.0L * (21.0L * (-2970.0L + 1395.0L * kappa*kappa * rho*rho - 80.0L * kappa*kappa*kappa*kappa * rho*rho*rho*rho + kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho) * CosKR + kappa*rho * (-62370.0L + 8505.0L * kappa*kappa * rho*rho - 231.0L * kappa*kappa*kappa*kappa * rho*rho*rho*rho + kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho) * SinKR) * fshielding1(rho, mu, shpower);
			LC += rho * ((10395.0L - 4725.0L * kappa*kappa * rho*rho + 210.0L * kappa*kappa*kappa*kappa * rho*rho*rho*rho - kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho) * CosKR + 21.0L * kappa*rho * (495.0L - 60.0L * kappa*kappa * rho*rho + kappa*kappa*kappa*kappa * rho*rho*rho*rho) * SinKR) * fshielding2(rho, mu, shpower);
			LC /= (-2.0L * kappa*kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho*rho*rho);*/
			long double LC1 = -2970.0L + 1395.0L * kapparho2 - 80.0L * kapparho4 + kapparho6;
			long double LC2 = -62370.0L + 8505.0L * kapparho2 - 231.0L * kapparho4 + kapparho6;
			long double LC3 = 10395.0L - 4725.0L * kapparho2 + 210.0L * kapparho4 - kapparho6;
			long double LC4 = 495.0L - 60.0L * kapparho2 + kapparho4;
			LC = 2.0L * (21.0L * LC1 * CosKR + kapparho * LC2 * SinKR) * fshielding1(rho, mu, shpower);
			LC += rho * (LC3 * CosKR + 21.0L * kapparho * LC4 * SinKR) * fshielding2(rho, mu, shpower);
			LC = LC / (-2.0L * kapparho7 * rho);
			break;
		case 7:  // K-Wave
			LC = 2.0L * ((945945.0L-446985.0L * kapparho2+26775.0L * kapparho4-406.0L * kapparho2*kapparho4 + kapparho4*kapparho4) * CosKR+7.0L * kapparho * (135135.0L-18810.0L * kapparho2+558.0L * kapparho4-4.0L * kapparho2*kapparho4) * SinKR) * fshielding1(rho, mu, shpower);
			LC += rho * (7.0L * (-19305.0L+8910.0L * kapparho2-450.0L * kapparho4+4.0L * kapparho2*kapparho4) * CosKR+kapparho * (-135135.0L+17325.0L * kapparho2-378.0L * kapparho4+kapparho2*kapparho4) * SinKR) * fshielding2(rho, mu, shpower);
			LC /= 2.0L * kapparho4*kapparho4 * rho;
			break;
		case 8:  // L-Wave
			LC = 2.0L * (9.0L * (1801800.0L - 855855.0L * kapparho*kapparho + 53130.0L * kapparho*kapparho*kapparho*kapparho - 910.0L * kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho + 4.0L * kappa*kappa*kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho*rho*rho) * CosKR + kapparho * (16216200.0L - 2297295.0L * kapparho*kapparho + 72765.0L * kapparho*kapparho*kapparho*kapparho - 666.0L * kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho + kappa*kappa*kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho*rho*rho) * SinKR) * fshielding1(rho, mu, shpower);
			LC -= rho * ((2027025.0L - 945945.0L * kapparho * kapparho + 51975.0L * kapparho*kapparho*kapparho*kapparho - 630.0L * kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho + kappa*kappa*kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho*rho*rho) * CosKR - 9.0L * kapparho * (-225225.0L + 30030.0L * kapparho*kapparho - 770.0L * kapparho*kapparho*kapparho*kapparho + 4.0L * kappa*kappa*kappa*kappa*kappa*kappa * rho*rho*rho*rho*rho*rho) * SinKR) * fshielding2(rho, mu, shpower);
			LC /= 2.0L * kapparho*kapparho*kapparho*kapparho*kapparho*kapparho*kapparho*kapparho*kapparho * rho;
			break;
		default:
			cout << "Internal error" << endl;
			exit(6);
	}

	return LC;
}


void VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar(vector <double> &AResults, vector <double> &BResults, int l, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, int NumPowers, vector <rPowers> &Powers, int Omega, double Lambda1, double Lambda2, double Lambda3)
{
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR12, LegendreWeightsR12;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r2Weights, r3Abscissas, r3Weights;
	vector <long double> r12Array, r13Array;
	int NumR2Points, NumR3Points, Prog = 0;
	vector <long double> r1Pow(Omega+l+1), r2Pow(Omega+l+1), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);
	long double SqrtKappa = sqrtl(kappa);

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
	for (int m = 0; m < nR1; m++) {
		r1Abscissas[m] = r1Abscissas[m] / (Powers[0].alpha + Lambda1);
		r1Weights[m] = r1Weights[m] / (Powers[0].alpha + Lambda1);
	}

	#pragma omp parallel for shared(r1Abscissas,r1Weights,Powers) private(r12Array,r13Array,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		vector <double> TempAResults(2*NumPowers, 0.0L), TempBResults(2*NumPowers, 0.0L);
		vector <long double> r1Pow(Omega+l+1), r2Pow(Omega+l+1), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);

		WriteProgress(string("PhiLS and PhiLC"), Prog, i, nR1);

		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);
		r12Array.resize(nR12);
		r13Array.resize(nR13);

		long double r1 = r1Abscissas[i];
		CreateRPowerLUT(r1Pow, r1, Omega+l);
		/*vector <double> r1Pow2(r1Pow.size());
		for (int p = 0; p < r1Pow.size(); p++)
			r1Pow2[p] = r1Pow[p];*/
		long double ExpR1L2 = expl(-(Powers[0].beta + Lambda2)*r1);
		long double ExpR1L3 = expl(-(Powers[0].gamma + Lambda3)*r1);

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR2) {
			for (int m = 0; m < nR2Lag; m++) {
				r2Weights[m] = LaguerreWeightsR2[m] / (Powers[0].beta + Lambda2);
				r2Abscissas[m] = LaguerreAbscissasR2[m] / (Powers[0].beta + Lambda2);
			}
			NumR2Points = nR2Lag;
		}
		// Do Gauss-Legendre for r2 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0L, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * expl(-(Powers[0].beta + Lambda2)*r2Abscissas[m]) * (r1-0.0L)/2.0L;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r1*(Powers[0].beta + Lambda2)) / (Powers[0].beta + Lambda2);
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * ExpR1L2 / (Powers[0].beta + Lambda2);
			}
		}

		for (int j = 0; j < NumR2Points; j++) {  // r2 integration
			long double r2 = r2Abscissas[j];
			long double a12 = fabs(r1-r2);
			long double b12 = fabs(r1+r2);
			CreateRPowerLUT(r2Pow, r2, Omega+l);
			ChangeOfIntervalNoResize(LegendreAbscissasR12, r12Array, a12, b12);

			// Cusp unimportant, so only do Gauss-Laguerre integration
			if (r1 > CuspR3) {
				for (int m = 0; m < nR3Lag; m++) {
					r3Weights[m] = LaguerreWeightsR3[m] / (Powers[0].gamma + Lambda3);
					r3Abscissas[m] = LaguerreAbscissasR3[m] / (Powers[0].gamma + Lambda3);
				}
				NumR3Points = nR3Lag;
			}
			// Do Gauss-Legendre for r3 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0L, r1);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * expl(-(Powers[0].gamma + Lambda3)*r3Abscissas[m]) * (r1-0.0L)/2.0L;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r1*(Powers[0].gamma + Lambda3)) / (Powers[0].gamma + Lambda3);
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * ExpR1L3 / (Powers[0].gamma + Lambda3);
				}
			}

			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				long double r3 = r3Abscissas[g];
				long double a13 = fabs(r1-r3);
				long double b13 = fabs(r1+r3);
				long double Coeff = r3Weights[g] * r2Weights[j] * r1Weights[i] * SqrtKappa * 0.70710678118654752440L * PI/nPhi23 * expl(Lambda1*r1 + Lambda2*r2 + Lambda3*r3);
				CreateRPowerLUT(r3Pow, r3, Omega);
				ChangeOfIntervalNoResize(LegendreAbscissasR13, r13Array, a13, b13);

				for (int k = 0; k < nR12; k++) {  // r12 integration
					long double r12 = r12Array[k];
					long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0L*r1*r2);
					long double Sin12 = sqrtl(1.0L - Cos12*Cos12);
					long double rho = 0.5L * sqrtl(2.0L*(r1*r1 + r2*r2) - r12*r12);
					long double jlrho = sf_bessel_jl(l, kappa*rho);
					long double nlrho = sf_bessel_nl(l, kappa*rho);
					long double ExpR12R3 = expl(-(r12/2.0L + r3));
					long double PotP = 2.0L/r1 - 2.0L/r3 - 2.0L/r12;
					CreateRPowerLUT(r12Pow, r12, Omega);

					for (int p = 0; p < nR13; p++) {  // r13 integration
						long double r13 = r13Array[p];
						long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0L*r1*r3);
						long double Sin13 = sqrtl(1.0L - Cos13*Cos13);
						long double rhop = 0.5L * sqrtl(2.0L*(r1*r1 + r3*r3) - r13*r13);
						long double jlrhop = sf_bessel_jl(l, kappa*rhop);
						long double nlrhop = sf_bessel_nl(l, kappa*rhop);
						long double ExpR13R2 = expl(-(r13/2.0L + r2));
						long double CoeffFinal = Coeff * LegendreWeightsR13[p] * (b13-a13) * LegendreWeightsR12[k] * (b12-a12) * r2 * r3 * r12 * r13;
						long double Pot = 2.0L/r1 - 2.0L/r2 - 2.0L/r13;
						CreateRPowerLUT(r13Pow, r13, Omega);

						// Phi1LS part
						long double S22 = ExpR12R3 * jlrho;
						long double S23 = ExpR13R2 * jlrhop;
						//long double AngPhi1S22 = 1.0L - 3.0L/8.0L * r2*r2 * Sin12*Sin12 / (rho*rho);
						long double AngPhi1S22 = AngR1Rho(l, r1, r2, Cos12, Sin12, rho);
						//long double AngPhi1S23 = 1.0L - 3.0L/8.0L * r3*r3 * Sin13*Sin13 / (rhop*rhop);
						long double AngPhi1S23 = AngR1Rhop(l, r1, r3, Cos13, Sin13, rhop);
						long double fOuterS1 = AngPhi1S22 * Pot * S22 + sf * AngPhi1S23 * PotP * S23;
						// Phi1LC part
						//@TODO: Would this be better to do as the other form with phi and P_23 phi?
						long double fshrho = fshielding(rho, mu, shpower);
						long double fshrhop = fshielding(rhop, mu, shpower);
						////long double fsh1rho = fshielding1(rho, mu, shpower) / rho * (2.0L * n2rho - kappa*rho * n1rho) - 0.5L * fshielding2(rho, mu, shpower) * n2rho;
						//long double kapparho = kappa*rho;
						//long double fsh1rho = 2.0L * (3.0L * (kapparho*kapparho - 2.0L) * cosl(kapparho) + kapparho * (kapparho*kapparho - 6.0L) * sinl(kapparho)) * fshielding1(rho, mu, shpower);
						//fsh1rho -= rho * ((kapparho*kapparho - 3.0L) * cosl(kapparho) - 3.0L * kapparho * sinl(kapparho)) * fshielding2(rho, mu, shpower);
						//fsh1rho = fsh1rho / (2.0L * kappa*kappa*kappa * rho*rho*rho*rho);
						long double fsh1rho = LaplacianC(l, kappa, rho, mu, shpower);
						//long double fsh1rhop = fshielding1(rhop, mu, shpower) / rhop * (2.0L * n2rhop - kappa*rhop * n1rhop) - 0.5L * fshielding2(rhop, mu, shpower) * n2rhop;
						long double fsh1rhop = LaplacianC(l, kappa, rhop, mu, shpower);
						long double AngPhi1C22 = AngPhi1S22;
						long double AngPhi1C23 = AngPhi1S23;
						long double fOuterC1 = -AngPhi1C22 * ExpR12R3 * (Pot * nlrho * fshrho + fsh1rho);
						fOuterC1 -= sf * AngPhi1C23 * ExpR13R2 * (PotP * nlrhop * fshrhop + fsh1rhop);

						// Phi2LS part
						//long double AngPhi2S22 = 1.0L - 3.0L/8.0L * r1*r1 * Sin12*Sin12 / (rho*rho);
						long double AngPhi2S22 = AngR2Rho(l, r1, r2, Cos12, Sin12, rho);
						// Phi2LC part
						long double AngPhi2C22 = AngPhi2S22;

						for (int m = 1; m <= nPhi23; m++) {  // phi_23 integration
							long double Phi23 = (2.0L*m - 1.0L)*PI/(2.0L*nPhi23);
							long double r23 = sqrtl(r2*r2 + r3*r3 - 2.0L*r2*r3*(Sin12*Sin13*cosl(Phi23) + Cos12*Cos13));
							long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0L*r2*r3);
							CreateRPowerLUT(r23Pow, r23, Omega);

							// Phi2LS part
							//long double AngPhi2S23 = r1 * Cos12 + r3 * Cos23;
							//AngPhi2S23 = 3.0L/8.0L * AngPhi2S23 * AngPhi2S23 / (rhop*rhop) - 0.5L;
							long double AngPhi2S23 = AngR2Rhop(l, r1, r3, Cos12, Cos13, Cos23, rhop);
							long double fOuterS2 = AngPhi2S22 * Pot * S22 + sf * AngPhi2S23 * PotP * S23;

							// Phi2LC part
							long double AngPhi2C23 = AngPhi2S23;
							long double fOuterC2 = -AngPhi2C22 * ExpR12R3 * (Pot * nlrho * fshrho + fsh1rho);
							fOuterC2 -= sf * AngPhi2C23 * ExpR13R2 * (PotP * nlrhop * fshrhop + fsh1rhop);

							// Combine with phi for final values
							for (int n = 0; n < NumPowers; n++) {
								rPowers *rp = &Powers[n];
								long double Common = r12Pow[rp->mi] * r3Pow[rp->ni] * r13Pow[rp->pi] * r23Pow[rp->qi] * CoeffFinal;
								long double Phi1andCoeff = r1Pow[rp->ki] * r2Pow[rp->li] * Common;
								TempAResults[n] += Phi1andCoeff * fOuterC1;
								TempBResults[n] += Phi1andCoeff * fOuterS1;
								rp = &Powers[NumPowers+n];
								long double Phi2andCoeff = r1Pow[rp->ki] * r2Pow[rp->li] * Common;
								TempAResults[NumPowers+n] += Phi2andCoeff * fOuterC2;
								TempBResults[NumPowers+n] += Phi2andCoeff * fOuterS2;
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
			AResults[NumPowers+n] += TempAResults[NumPowers+n];
			BResults[NumPowers+n] += TempBResults[NumPowers+n];
		}
	}

	return;
}


void VecGaussIntegrationPhi13_PhiLCBar_PhiLSBar_R23Term(vector <double> &AResults, vector <double> &BResults, int l, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nPhi13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, int NumPowers, vector <rPowers> &Powers, int Omega, double Lambda1, double Lambda2, double Lambda3)
{
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR12, LegendreWeightsR12;
	vector <long double> LegendreAbscissasR23, LegendreWeightsR23;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r3Abscissas, r2Weights, r3Weights;
	vector <long double> r23Array, r12Array;
	int NumR2Points, NumR3Points, Prog = 0;
	long double SqrtKappa = sqrtl(kappa);

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
		r1Abscissas[m] = r1Abscissas[m] / (Powers[0].alpha + Lambda1);
		r1Weights[m] = r1Weights[m] / (Powers[0].alpha + Lambda1);
	}
	
	#pragma omp parallel for shared(r1Abscissas,r1Weights,Powers) private(r12Array,r23Array,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		vector <long double> TempAResults(2*NumPowers, 0.0L), TempBResults(2*NumPowers, 0.0L);
		vector <long double> r1Pow(Omega+l+1), r2Pow(Omega+l+1), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);

		WriteProgress(string("PhiLS and PhiLC R23"), Prog, i, nR1);

		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);
		r12Array.resize(nR12);
		r23Array.resize(nR23);

		long double r1 = r1Abscissas[i];
		CreateRPowerLUT(r1Pow, r1, Omega+l);
		long double ExpR1L2 = expl(-(Powers[0].beta + Lambda2)*r1);

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR2) {
			for (int m = 0; m < nR2Lag; m++) {
				r2Weights[m] = LaguerreWeightsR2[m] / (Powers[0].beta + Lambda2);
				r2Abscissas[m] = LaguerreAbscissasR2[m] / (Powers[0].beta + Lambda2);
			}
			NumR2Points = nR2Lag;
		}
		// Do Gauss-Legendre for r2 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0L, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * expl(-(Powers[0].beta + Lambda2)*r2Abscissas[m]) * (r1-0.0L)/2.0L;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r1*(Powers[0].beta + Lambda2)) / (Powers[0].beta + Lambda2);
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * ExpR1L2 / (Powers[0].beta + Lambda2);
			}
		}

		for (int j = 0; j < NumR2Points; j++) {  // r2 integration
			long double r2 = r2Abscissas[j];
			long double a12 = fabs(r1-r2);
			long double b12 = fabs(r1+r2);
			CreateRPowerLUT(r2Pow, r2, Omega+l);
			ChangeOfIntervalNoResize(LegendreAbscissasR12, r12Array, a12, b12);

			if (r2 > CuspR3) {
				for (int m = 0; m < nR3Lag; m++) {
					r3Weights[m] = LaguerreWeightsR3[m] / (Powers[0].gamma + Lambda3);
					r3Abscissas[m] = LaguerreAbscissasR3[m] / (Powers[0].gamma + Lambda3);
				}
				NumR3Points = nR3Lag;
			}
			// Do Gauss-Legendre for r3 integration over [0,r2] and then Gauss-Laguerre over [r2,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0L, r2);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * expl(-(Powers[0].gamma + Lambda3)*r3Abscissas[m]) * (r2-0.0L)/2.0L;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r2*(Powers[0].gamma + Lambda3)) / (Powers[0].gamma + Lambda3);
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * expl(-(Powers[0].gamma + Lambda3)*r2) / (Powers[0].gamma + Lambda3);
				}
			}

			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				long double r3 = r3Abscissas[g];
				long double a23 = fabs(r2-r3);
				long double b23 = fabs(r2+r3);
				long double Coeff = r3Weights[g] * r2Weights[j] * r1Weights[i] * SqrtKappa * 0.70710678118654752440L * PI/nPhi13 * expl(Lambda1*r1 + Lambda2*r2 + Lambda3*r3);
				CreateRPowerLUT(r3Pow, r3, Omega);
				ChangeOfIntervalNoResize(LegendreAbscissasR23, r23Array, a23, b23);

				for (int k = 0; k < nR23; k++) {  // r23 integration
					long double r23 = r23Array[k];
					CreateRPowerLUT(r23Pow, r23, Omega);
					long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0L*r2*r3);
					long double Sin23 = sqrtl(1.0L - Cos23*Cos23);
					long double Pot = 2.0L/r23;
					long double PotP = Pot;  // Easier comparison with VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar

					for (int p = 0; p < nR12; p++) {  // r12 integration
						long double r12 = r12Array[p];
						CreateRPowerLUT(r12Pow, r12, Omega);
						long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0L*r1*r2);
						long double Sin12 = sqrtl(1.0L - Cos12*Cos12);
						long double rho = 0.5L * sqrtl(2.0L*(r1*r1 + r2*r2) - r12*r12);
						long double jlrho = sf_bessel_jl(l, kappa*rho);
						long double nlrho = sf_bessel_nl(l, kappa*rho);
						long double CoeffFinal = Coeff * LegendreWeightsR12[p] * (b12-a12) * LegendreWeightsR23[k] * (b23-a23) * r1 * r3 * r12 * r23;
						long double ExpR12R3 = expl(-(r12/2.0L + r3));
						// Phi1LS
						long double S22 = ExpR12R3 * jlrho;
						//long double AngPhi1S22 = 1.0L - 3.0L/8.0L * r2*r2 * Sin12*Sin12 / (rho*rho);
						long double AngPhi1S22 = AngR1Rho(l, r1, r2, Cos12, Sin12, rho);
						// Phi2LS
						//long double AngPhi2S22 = r1 * Cos12 + r2;
						//AngPhi2S22 = 3.0L/8.0L * AngPhi2S22 * AngPhi2S22 / (rho*rho) - 0.5L;
						//long double AngPhi2S22 = 1.0L - 3.0L/8.0L * r1*r1 * Sin12*Sin12 / (rho*rho);
						long double AngPhi2S22 = AngR2Rho(l, r1, r2, Cos12, Sin12, rho);
						// Phi1LC
						long double fshrho = fshielding(rho, mu, shpower);
						long double AngPhi1C22 = AngPhi1S22;
						long double fOuterC1Part = -AngPhi1C22 * ExpR12R3 * (Pot * nlrho * fshrho);
						// Phi2LC
						long double AngPhi2C22 = AngPhi2S22;
						long double fOuterC2Part = -AngPhi2C22 * ExpR12R3 * (Pot * nlrho * fshrho);

						for (int m = 1; m <= nPhi13; m++) {  // phi_13 integration
							//long double Phi13 = (2.0L*m - 1.0L)*PI/(2.0L*nPhi13);
							long double Phi13 = (2*m - 1)*PI/(2*nPhi13);
							long double r13 = sqrtl(r1*r1 + r3*r3 - 2.0L*r1*r3*(Sin12*Sin23*cosl(Phi13) + Cos12*Cos23));
							long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0L*r1*r3);
							long double Sin13 = sqrtl(1.0L - Cos13*Cos13);
							CreateRPowerLUT(r13Pow, r13, Omega);
							//long double rhop = 0.5L * sqrtl(2.0L*(r1*r1 + r3*r3) - r13*r13);
							long double rhop = 0.5L * sqrtl(2*(r1*r1 + r3*r3) - r13*r13);
							long double jlrhop = sf_bessel_jl(l, kappa*rhop);
							long double nlrhop = sf_bessel_nl(l, kappa*rhop);
							//long double ExpR13R2 = expl(-(r13/2.0L + r2));
							long double ExpR13R2 = expl(-(r13/2.0L + r2));

							// Phi1LS part
							long double S23 = ExpR13R2 * jlrhop;
							//long double AngPhi1S23 = 1.0L - 3.0L/8.0L * r3*r3 * Sin13*Sin13 / (rhop*rhop);
							long double AngPhi1S23 = AngR1Rhop(l, r1, r3, Cos13, Sin13, rhop);
							long double fOuterS1 = AngPhi1S22 * Pot * S22 + sf * AngPhi1S23 * PotP * S23;

							// Phi2LS part
							//long double AngPhi2S23 = r1 * Cos12 + r3 * Cos23;
							//AngPhi2S23 = 3.0L/8.0L * AngPhi2S23 * AngPhi2S23 / (rhop*rhop) - 0.5L;
							long double AngPhi2S23 = AngR2Rhop(l, r1, r3, Cos12, Cos13, Cos23, rhop);
							long double fOuterS2 = AngPhi2S22 * Pot * S22 + sf * AngPhi2S23 * PotP * S23;

							// Phi1LC part
							//@TODO: Would this be better to do as the other form with phi and P_23 phi?
							long double fshrhop = fshielding(rhop, mu, shpower);
							long double AngPhi1C23 = AngPhi1S23;
							long double fOuterC1 = fOuterC1Part - sf * AngPhi1C23 * ExpR13R2 * (PotP * nlrhop * fshrhop);

							// Phi2LC part
							//@TODO: Same question as above
							long double AngPhi2C23 = AngPhi2S23;
							long double fOuterC2 = fOuterC2Part - sf * AngPhi2C23 * ExpR13R2 * (PotP * nlrhop * fshrhop);

							for (int n = 0; n < NumPowers; n++) {
								rPowers *rp = &Powers[n];
								long double Common = r12Pow[rp->mi] * r3Pow[rp->ni] * r13Pow[rp->pi] * r23Pow[rp->qi] * CoeffFinal;
								long double Phi1andCoeff = r1Pow[rp->ki] * r2Pow[rp->li] * Common;
								TempAResults[n] += Phi1andCoeff * fOuterC1;
								TempBResults[n] += Phi1andCoeff * fOuterS1;
								rp = &Powers[NumPowers+n];
								long double Phi2andCoeff = r1Pow[rp->ki] * r2Pow[rp->li] * Common;
								TempAResults[NumPowers+n] += Phi2andCoeff * fOuterC2;
								TempBResults[NumPowers+n] += Phi2andCoeff * fOuterS2;
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
			AResults[NumPowers+n] += TempAResults[NumPowers+n];
			BResults[NumPowers+n] += TempBResults[NumPowers+n];
		}
	}

	return;
}


//void VecGaussIntegrationPhi12_PhiLCBar_PhiLSBar_R23Term(vector <double> &AResults, vector <double> &BResults, int l, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nPhi12, int nR13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, int NumPowers, vector <rPowers> &Powers, int Omega, double Lambda1, double Lambda2, double Lambda3)
//{  THIS FUNCTION HAS NOT BEEN UPDATED WITH THE REST!
//	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
//	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
//	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
//	vector <long double> LegendreAbscissasR23, LegendreWeightsR23;
//	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r3Abscissas, r2Weights, r3Weights;
//	vector <long double> r23Array, r13Array;
//	int NumR2Points, NumR3Points, Prog = 0;
//	long double SqrtKappa = sqrtl(kappa);
//
//	// Create the abscissas and weights for the needed number of points.
//	GaussLaguerre(r1Abscissas, r1Weights, nR1);
//	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
//	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
//	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
//	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
//	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);
//	GaussLegendre(LegendreAbscissasR23, LegendreWeightsR23, nR23);
//
//	r1Abscissas.resize(nR1);
//	r1Weights.resize(nR1);
//	for (int m = 0; m < nR1; m++) {  //@TODO: Can we just do this as a single line a la Fortran?
//		r1Abscissas[m] = r1Abscissas[m] / Powers[0].alpha;
//		r1Weights[m] = r1Weights[m] / Powers[0].alpha;
//	}
//
//	#pragma omp parallel for shared(r1Abscissas,r1Weights,Powers) private(r13Array,r23Array,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
//	for (int i = 0; i < nR1; i++) {  // r1 integration
//		vector <long double> TempAResults(2*NumPowers, 0.0L), TempBResults(2*NumPowers, 0.0L);
//		vector <long double> r1Pow(Omega+2), r2Pow(Omega+2), r3Pow(Omega+2), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);
//
//		WriteProgress(string("PhiLS and PhiLC R23"), Prog, i, nR1);
//
//		r2Abscissas.resize(nR2Leg + nR2Lag);
//		r2Weights.resize(nR2Leg + nR2Lag);
//		r3Abscissas.resize(nR3Leg + nR3Lag);
//		r3Weights.resize(nR3Leg + nR3Lag);
//		r13Array.resize(nR13);
//		r23Array.resize(nR23);
//
//		long double r1 = r1Abscissas[i];
//		CreateRPowerLUT(r1Pow, r1, Omega+2);
//
//		// If r1 is large enough, we just do Gauss-Laguerre.
//		if (r1 > CuspR3) {
//			for (int m = 0; m < nR3Lag; m++) {
//				r3Weights[m] = LaguerreWeightsR3[m] / Powers[0].gamma;
//				r3Abscissas[m] = LaguerreAbscissasR3[m] / Powers[0].gamma;
//			}
//			NumR3Points = nR3Lag;
//		}
//		// Do Gauss-Legendre for r3 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
//		else {
//			ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0L, r1);
//			NumR3Points = nR3Leg + nR3Lag;
//			for (int m = 0; m < nR3Leg; m++) {
//				r3Weights[m] = LegendreWeightsR3[m] * expl(-Powers[0].gamma*r3Abscissas[m]) * (r1-0.0L)/2.0L;
//			}
//			for (int m = nR3Leg; m < NumR3Points; m++) {
//				r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r1*Powers[0].gamma) / Powers[0].gamma;
//				r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * expl(-Powers[0].gamma*r1) / Powers[0].gamma;
//			}
//		}
//
//		for (int j = 0; j < NumR3Points; j++) {  // r3 integration
//			long double r3 = r3Abscissas[j];
//			long double a13 = fabs(r1-r3);
//			long double b13 = fabs(r1+r3);
//			CreateRPowerLUT(r3Pow, r3, Omega+2);
//			ChangeOfIntervalNoResize(LegendreAbscissasR13, r13Array, a13, b13);
//
//			if (r3 > CuspR2) {
//				for (int m = 0; m < nR2Lag; m++) {
//					r2Weights[m] = LaguerreWeightsR2[m] / Powers[0].beta;
//					r2Abscissas[m] = LaguerreAbscissasR2[m] / Powers[0].beta;
//				}
//				NumR2Points = nR2Lag;
//			}
//			// Do Gauss-Legendre for r2 integration over [0,r3] and then Gauss-Laguerre over [r3,inf).
//			else {
//				ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0L, r3);
//				NumR2Points = nR2Leg + nR2Lag;
//				for (int m = 0; m < nR2Leg; m++) {
//					r2Weights[m] = LegendreWeightsR2[m] * expl(-Powers[0].beta*r2Abscissas[m]) * (r3-0.0L)/2.0L;
//				}
//				for (int m = nR2Leg; m < NumR2Points; m++) {
//					r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r3*Powers[0].beta) / Powers[0].beta;
//					r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * expl(-Powers[0].beta*r3) / Powers[0].beta;
//				}
//			}
//
//			for (int g = 0; g < NumR2Points; g++) {  // r2 integration
//				long double r2 = r2Abscissas[g];
//				long double a23 = fabs(r2-r3);
//				long double b23 = fabs(r2+r3);
//				long double Coeff = r3Weights[j] * r2Weights[g] * r1Weights[i] * SqrtKappa * 0.70710678118654752440L * PI/nPhi12;
//				CreateRPowerLUT(r2Pow, r2, Omega+2);
//				ChangeOfIntervalNoResize(LegendreAbscissasR23, r23Array, a23, b23);
//
//				for (int k = 0; k < nR23; k++) {  // r23 integration
//					long double r23 = r23Array[k];
//					CreateRPowerLUT(r23Pow, r23, Omega);
//					long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0L*r2*r3);
//					long double Sin23 = sqrtl(1.0L - Cos23*Cos23);
//					long double Pot = 2.0L/r23;
//					long double PotP = Pot;  // Easier comparison with VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar
//
//					for (int p = 0; p < nR13; p++) {  // r13 integration
//						long double r13 = r13Array[p];
//						CreateRPowerLUT(r13Pow, r13, Omega);
//						long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0L*r1*r3);
//						long double Sin13 = sqrtl(1.0L - Cos13*Cos13);
//						long double rhop = 0.5L * sqrt(2.0L*(r1*r1 + r3*r3) - r13*r13);
//						long double j2rhop = sf_bessel_jl(2, kappa*rhop);
//						long double n2rhop = sf_bessel_nl(2, kappa*rhop);
//						long double CoeffFinal = Coeff * LegendreWeightsR13[p] * (b13-a13) * LegendreWeightsR23[k] * (b23-a23) * r1 * r2 * r13 * r23;
//						long double ExpR13R2 = expl(-(r13/2.0L + r2));
//
//						for (int m = 1; m <= nPhi12; m++) {  // phi_12 integration
//							long double Phi12 = (2.0L*m - 1.0L)*PI/(2.0L*nPhi12);
//							long double r12 = sqrt(r1*r1 + r2*r2 - 2.0L*r1*r2*(Sin13*Sin23*cosl(Phi12) + Cos13*Cos23));
//							long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0L*r1*r2);
//							long double Sin12 = sqrtl(1.0L - Cos12*Cos12);
//							CreateRPowerLUT(r12Pow, r12, Omega);
//							long double rho = 0.5L * sqrt(2.0L*(r1*r1 + r2*r2) - r12*r12);
//							long double j2rho = sf_bessel_jl(2, kappa*rho);
//							long double n2rho = sf_bessel_nl(2, kappa*rho);
//							long double ExpR12R3 = expl(-(r12/2.0L + r3));
//
//							// Phi1LS
//							long double S22 = ExpR12R3 * jlrho;
//							long double AngPhi1S22 = 1.0L - 3.0L/8.0L * r2*r2 * Sin12*Sin12 / (rho*rho);
//							// Phi2LS
//							long double AngPhi2S22 = (r2 + r1 * Cos12) / rho;
//							// Phi1LC
//							long double fshrho = fshielding(rho, mu, shpower);
//							long double AngPhi1C22 = AngPhi1S22;
//							long double fOuterC1Part = -AngPhi1C22 * ExpR12R3 * (Pot * nlrho * fshrho);
//							// Phi2LC
//							long double AngPhi2C22 = AngPhi2S22;
//							long double fOuterC2Part = -AngPhi2C22 * ExpR12R3 * (Pot * nlrho * fshrho);
//
//							// Phi1LS part
//							long double S23 = ExpR13R2 * jlrhop;
//							long double AngPhi1S23 = 1.0L - 3.0L/8.0L * r3*r3 * Sin13*Sin13 / (rhop*rhop);
//							long double fOuterS1 = AngPhi1S22 * Pot * S22 + sf * AngPhi1S23 * PotP * S23;
//							// Phi2LS part
//							long double AngPhi2S23 = r1 * Cos13 + r2 * Cos23;
//							AngPhi2S23 = 3.0L/8.0L * AngPhi2S23 * AngPhi2S23 / (rhop*rhop) - 0.5L;
//							long double fOuterS2 = AngPhi2S22 * Pot * S22 + sf * AngPhi2S23 * PotP * S23;
//							// Phi1LC part
//							//@TODO: Would this be better to do as the other form with phi and P_23 phi?
//							long double fshrhop = fshielding(rhop, mu, shpower);
//							long double AngPhi1C23 = (r1 + r3 * Cos13) / rhop;
//							long double fOuterC1 = fOuterC1Part - sf * AngPhi1C23 * ExpR13R2 * (PotP * nlrhop * fshrhop);
//							// Phi2LC part
//							//@TODO: Same question as above
//							long double AngPhi2C23 = AngPhi2S23;
//							long double fOuterC2 = fOuterC2Part - sf * AngPhi2C23 * ExpR13R2 * (PotP * nlrhop * fshrhop);
//
//							for (int n = 0; n < NumPowers; n++) {
//								rPowers *rp = &Powers[n];
//								long double Common = r12Pow[rp->mi] * r3Pow[rp->ni] * r13Pow[rp->pi] * r23Pow[rp->qi] * CoeffFinal;
//								long double Phi1andCoeff = r1Pow[rp->ki] * r2Pow[rp->li] * Common;
//								TempAResults[n] += Phi1andCoeff * fOuterC1;
//								TempBResults[n] += Phi1andCoeff * fOuterS1;
//								rp = &Powers[NumPowers+n];
//								long double Phi2andCoeff = r1Pow[rp->ki] * r2Pow[rp->li] * Common;
//								TempAResults[NumPowers+n] += Phi2andCoeff * fOuterC2;
//								TempBResults[NumPowers+n] += Phi2andCoeff * fOuterS2;
//							}
//						}
//					}
//				}
//			}
//		}
//
//		//#pragma omp critical(build)
//		for (int n = 0; n < NumPowers; n++) {
//			AResults[n] += TempAResults[n];
//			BResults[n] += TempBResults[n];
//			AResults[NumPowers+n] += TempAResults[NumPowers+n];
//			BResults[NumPowers+n] += TempBResults[NumPowers+n];
//		}
//	}
//
//	return;
//}


void VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar_Full(vector <double> &AResults, vector <double> &BResults, int l, int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nR12, int nR13, int nPhi23, double CuspR2, double CuspR3, double kappa, double mu, int shpower, int sf, int NumPowers, vector <rPowers> &Powers, int Omega, double Lambda1, double Lambda2, double Lambda3)
{
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR12, LegendreWeightsR12;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r2Weights, r3Abscissas, r3Weights;
	vector <long double> r12Array, r13Array;
	int NumR2Points, NumR3Points, Prog = 0;
	vector <long double> r1Pow(Omega+2), r2Pow(Omega+2), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);
	long double SqrtKappa = sqrtl(kappa);

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
	for (int m = 0; m < nR1; m++) {
		r1Abscissas[m] = r1Abscissas[m] / (Powers[0].alpha + Lambda1);
		r1Weights[m] = r1Weights[m] / (Powers[0].alpha + Lambda1);
	}

	#pragma omp parallel for shared(r1Abscissas,r1Weights,Powers) private(r12Array,r13Array,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		vector <long double> TempAResults(2*NumPowers, 0.0L), TempBResults(2*NumPowers, 0.0L);
		vector <long double> r1Pow(Omega+l+1), r2Pow(Omega+l+1), r3Pow(Omega+1), r12Pow(Omega+1), r13Pow(Omega+1), r23Pow(Omega+1);

		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);
		r12Array.resize(nR12);
		r13Array.resize(nR13);

		long double r1 = r1Abscissas[i];
		CreateRPowerLUT(r1Pow, r1, Omega+l);
		long double ExpR1L2 = expl(-(Powers[0].beta + Lambda2)*r1);
		long double ExpR1L3 = expl(-(Powers[0].gamma + Lambda3)*r1);

		// If r1 is large enough, we just do Gauss-Laguerre.
		if (r1 > CuspR2) {
			for (int m = 0; m < nR2Lag; m++) {
				r2Weights[m] = LaguerreWeightsR2[m] / (Powers[0].beta + Lambda2);
				r2Abscissas[m] = LaguerreAbscissasR2[m] / (Powers[0].beta + Lambda2);
			}
			NumR2Points = nR2Lag;
		}
		// Do Gauss-Legendre for r2 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
		else {
			ChangeOfIntervalNoResize(LegendreAbscissasR2, r2Abscissas, 0.0L, r1);
			NumR2Points = nR2Leg + nR2Lag;
			for (int m = 0; m < nR2Leg; m++) {
				r2Weights[m] = LegendreWeightsR2[m] * expl(-(Powers[0].beta + Lambda2)*r2Abscissas[m]) * (r1-0.0L)/2.0L;
			}
			for (int m = nR2Leg; m < NumR2Points; m++) {
				r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r1*(Powers[0].beta + Lambda2)) / (Powers[0].beta + Lambda2);
				r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * ExpR1L2 / (Powers[0].beta + Lambda2);
			}
		}

		for (int j = 0; j < NumR2Points; j++) {  // r2 integration
			long double r2 = r2Abscissas[j];
			long double a12 = fabs(r1-r2);
			long double b12 = fabs(r1+r2);
			CreateRPowerLUT(r2Pow, r2, Omega+l);
			ChangeOfIntervalNoResize(LegendreAbscissasR12, r12Array, a12, b12);

			// Cusp unimportant, so only do Gauss-Laguerre integration
			if (r1 > CuspR3) {
				for (int m = 0; m < nR3Lag; m++) {
					r3Weights[m] = LaguerreWeightsR3[m] / (Powers[0].gamma + Lambda3);
					r3Abscissas[m] = LaguerreAbscissasR3[m] / (Powers[0].gamma + Lambda3);
				}
				NumR3Points = nR3Lag;
			}
			// Do Gauss-Legendre for r3 integration over [0,r1] and then Gauss-Laguerre over [r1,inf).
			else {
				ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0L, r1);
				NumR3Points = nR3Leg + nR3Lag;
				for (int m = 0; m < nR3Leg; m++) {
					r3Weights[m] = LegendreWeightsR3[m] * expl(-(Powers[0].gamma + Lambda3)*r3Abscissas[m]) * (r1-0.0L)/2.0L;
				}
				for (int m = nR3Leg; m < NumR3Points; m++) {
					r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r1*(Powers[0].gamma + Lambda3)) / (Powers[0].gamma + Lambda3);
					r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * ExpR1L3 / (Powers[0].gamma + Lambda3);
				}
			}

			for (int g = 0; g < NumR3Points; g++) {  // r3 integration
				long double r3 = r3Abscissas[g];
				long double a13 = fabs(r1-r3);
				long double b13 = fabs(r1+r3);
				long double Coeff = r3Weights[g] * r2Weights[j] * r1Weights[i] * SqrtKappa * 0.70710678118654752440L * PI/nPhi23 * expl(Lambda1*r1 + Lambda2*r2 + Lambda3*r3);
				CreateRPowerLUT(r3Pow, r3, Omega);
				ChangeOfIntervalNoResize(LegendreAbscissasR13, r13Array, a13, b13);

				for (int k = 0; k < nR12; k++) {  // r12 integration
					long double r12 = r12Array[k];
					long double Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0L*r1*r2);
					long double Sin12 = sqrtl(1.0L - Cos12*Cos12);
					long double rho = 0.5L * sqrtl(2.0L*(r1*r1 + r2*r2) - r12*r12);
					long double jlrho = sf_bessel_jl(l, kappa*rho);
					long double nlrho = sf_bessel_nl(l, kappa*rho);
					long double ExpR12R3 = expl(-(r12/2.0L + r3));
					CreateRPowerLUT(r12Pow, r12, Omega);

					for (int p = 0; p < nR13; p++) {  // r13 integration
						long double r13 = r13Array[p];
						long double Cos13 = (r1*r1 + r3*r3 - r13*r13) / (2.0L*r1*r3);
						long double Sin13 = sqrtl(1.0L - Cos13*Cos13);
						long double rhop = 0.5L * sqrtl(2.0L*(r1*r1 + r3*r3) - r13*r13);
						long double jlrhop = sf_bessel_jl(l, kappa*rhop);
						long double nlrhop = sf_bessel_nl(l, kappa*rhop);
						long double ExpR13R2 = expl(-(r13/2.0L + r2));
						long double CoeffFinal = Coeff * LegendreWeightsR13[p] * (b13-a13) * LegendreWeightsR12[k] * (b12-a12) * r2 * r3 * r12 * r13;
						CreateRPowerLUT(r13Pow, r13, Omega);

						// Phi1LS part
						long double S22 = ExpR12R3 * jlrho;  // S22
						long double S23 = ExpR13R2 * jlrhop;  // S23
						//long double AngPhi1S22 = 1.0L - 3.0L/8.0L * r2*r2 * Sin12*Sin12 / (rho*rho);
						long double AngPhi1S22 = AngR1Rho(l, r1, r2, Cos12, Sin12, rho);
						//long double AngPhi1S23 = 1.0L - 3.0L/8.0L * r3*r3 * Sin13*Sin13 / (rhop*rhop);
						long double AngPhi1S23 = AngR1Rhop(l, r1, r3, Cos13, Sin13, rhop);
						// Phi2LS part
						//long double AngPhi2S22 = 1.0L - 3.0L/8.0L * r1*r1 * Sin12*Sin12 / (rho*rho);
						long double AngPhi2S22 = AngR2Rho(l, r1, r2, Cos12, Sin12, rho);
						// Phi1LC part
						//@TODO: Would this be better to do as the other form with phi and P_23 phi?
						long double fshrho = fshielding(rho, mu, shpower);
						long double fshrhop = fshielding(rhop, mu, shpower);
						//long double fsh1rho = fshielding1(rho, mu, shpower) / rho * (2.0L * n2rho - kappa*rho * n1rho) - 0.5L * fshielding2(rho, mu, shpower) * n2rho;
						long double fsh1rho = LaplacianC(l, kappa, rho, mu, shpower);
						//long double fsh1rhop = fshielding1(rhop, mu, shpower) / rhop * (2.0L * n2rhop - kappa*rhop * n1rhop) - 0.5L * fshielding2(rhop, mu, shpower) * n2rhop;
						long double fsh1rhop = LaplacianC(l, kappa, rhop, mu, shpower);
						long double AngPhi1C22 = AngPhi1S22;
						long double AngPhi1C23 = AngPhi1S23;
						// Phi2LC part
						//@TODO: Same question as above
						long double AngPhi2C22 = AngPhi2S22;

						for (int m = 1; m <= nPhi23; m++) {  // phi_23 integration
							long double Phi23 = (2.0L*m - 1.0L)*PI/(2.0L*nPhi23);
							long double r23 = sqrtl(r2*r2 + r3*r3 - 2.0L*r2*r3*(Sin12*Sin13*cosl(Phi23) + Cos12*Cos13));
							long double Cos23 = (r2*r2 + r3*r3 - r23*r23) / (2.0L*r2*r3);
							CreateRPowerLUT(r23Pow, r23, Omega);

							long double PotP = 2.0L/r1 - 2.0L/r3 - 2.0L/r12 + 2.0L/r23;
							long double Pot = 2.0L/r1 - 2.0L/r2 - 2.0L/r13 + 2.0L/r23;

							// Phi1LS part
							long double fOuterS1 = AngPhi1S22 * Pot * S22 + sf * AngPhi1S23 * PotP * S23;

							// Phi1LC part
							long double fOuterC1 = -AngPhi1C22 * ExpR12R3 * (Pot * nlrho * fshrho + fsh1rho);
							fOuterC1 -= sf * AngPhi1C23 * ExpR13R2 * (PotP * nlrhop * fshrhop + fsh1rhop);

							// Phi2LS part
							//long double AngPhi2S23 = r1 * Cos12 + r3 * Cos23;
							//AngPhi2S23 = 3.0L/8.0L * AngPhi2S23 * AngPhi2S23 / (rhop*rhop) - 0.5L;
							long double AngPhi2S23 = AngR2Rhop(l, r1, r3, Cos12, Cos13, Cos23, rhop);
							long double fOuterS2 = AngPhi2S22 * Pot * S22 + sf * AngPhi2S23 * PotP * S23;

							// Phi2LC part
							long double AngPhi2C23 = AngPhi2S23;
							long double fOuterC2 = -AngPhi2C22 * ExpR12R3 * (Pot * nlrho * fshrho + fsh1rho);
							fOuterC2 -= sf * AngPhi2C23 * ExpR13R2 * (PotP * nlrhop * fshrhop + fsh1rhop);

							// Combine with phi for final values
							for (int n = 0; n < NumPowers; n++) {
								rPowers *rp = &Powers[n];
								long double Common = r12Pow[rp->mi] * r3Pow[rp->ni] * r13Pow[rp->pi] * r23Pow[rp->qi] * CoeffFinal;
								long double Phi1andCoeff = r1Pow[rp->ki] * r2Pow[rp->li] * Common;
								TempAResults[n] += Phi1andCoeff * fOuterC1;
								TempBResults[n] += Phi1andCoeff * fOuterS1;
								rp = &Powers[NumPowers+n];
								long double Phi2andCoeff = r1Pow[rp->ki] * r2Pow[rp->li] * Common;
								TempAResults[NumPowers+n] += Phi2andCoeff * fOuterC2;
								TempBResults[NumPowers+n] += Phi2andCoeff * fOuterS2;
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
			AResults[NumPowers+n] += TempAResults[NumPowers+n];
			BResults[NumPowers+n] += TempBResults[NumPowers+n];
		}
		WriteProgress(string("PhiLS and PhiLC Full"), Prog, i, nR1);
	}

	return;
}
