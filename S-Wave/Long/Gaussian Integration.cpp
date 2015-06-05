//
// Gaussian Integration.cpp: Has Gauss-Legendre and Gauss-Laguerre integration routines
//

#include <vector>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <float.h>
#include <cstdio>
#include "Ps-H Scattering.h"
#include "Gaussian Integration.h"

extern	double PI;

void	ChangeOfInterval(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, double a, double b);
void	ChangeOfIntervalNoResize(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, double a, double b);
int		GaussLegendre(vector <long double> &Abscissas, vector <long double> &Weights, int n);
int		GaussLaguerre(vector <long double> &Abscissas, vector <long double> &Weights, int n);


// Returns whether a double is a finite number.
//  Used only for debugging purposes to see if a number is NaN or otherwise
bool IsFiniteNumber(double x) 
{
	return (x <= DBL_MAX && x >= -DBL_MAX); 
}


void ChangeOfInterval(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, double a, double b)
{
	unsigned int Size = Abscissas.size();
	ChangedAbscissas.resize(Size);
	for (unsigned int i = 0; i < Size; i++) {
		ChangedAbscissas[i] = Abscissas[i] * (b-a)/2.0 + (a+b)/2.0;
	}
}

void ChangeOfIntervalNoResize(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, double a, double b)
{
	unsigned int Size = Abscissas.size();
	for (unsigned int i = 0; i < Size; i++) {
		ChangedAbscissas[i] = Abscissas[i] * (b-a)/2.0 + (a+b)/2.0;
	}
}


// Generates Gauss-Legendre quadrature abscissas and weights
//  This function is based on the gauleg function on pages 145-146 of
//  Numerical Recipes in Fortran, Second Edition.
//@TODO: We could also make use of the fact that for odd-degree polynomials, one root is always 0.
int GaussLegendre(vector <long double> & Abscissas, vector <long double> & Weights, int n)
{
	double Tolerance = 1e-16;  //@TODO: Adjust tolerance accordingly.

	//@TODO: Remove
	Abscissas.resize(n);
	Weights.resize(n);

	switch (n)
	{
		case 5:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs5[i];
				Weights[i] = LegWeight5[i];
			}
			return 0;
		case 10:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs10[i];
				Weights[i] = LegWeight10[i];
			}
			return 0;
		case 15:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs15[i];
				Weights[i] = LegWeight15[i];
			}
			return 0;
		case 20:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs20[i];
				Weights[i] = LegWeight20[i];
			}
			return 0;
		case 25:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs25[i];
				Weights[i] = LegWeight25[i];
			}
			return 0;
		case 30:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs30[i];
				Weights[i] = LegWeight30[i];
			}
			return 0;
		case 35:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs35[i];
				Weights[i] = LegWeight35[i];
			}
			return 0;
		case 40:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs40[i];
				Weights[i] = LegWeight40[i];
			}
			return 0;
		case 45:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs45[i];
				Weights[i] = LegWeight45[i];
			}
			return 0;
		case 50:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs50[i];
				Weights[i] = LegWeight50[i];
			}
			return 0;
		case 55:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs55[i];
				Weights[i] = LegWeight55[i];
			}
			return 0;
		case 60:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs60[i];
				Weights[i] = LegWeight60[i];
			}
			return 0;
		case 65:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs65[i];
				Weights[i] = LegWeight65[i];
			}
			return 0;
		case 70:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs70[i];
				Weights[i] = LegWeight70[i];
			}
			return 0;
		case 75:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs75[i];
				Weights[i] = LegWeight75[i];
			}
			return 0;
		case 80:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs80[i];
				Weights[i] = LegWeight80[i];
			}
			return 0;
		case 85:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs85[i];
				Weights[i] = LegWeight85[i];
			}
			return 0;
		case 90:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs90[i];
				Weights[i] = LegWeight90[i];
			}
			return 0;
		case 95:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs95[i];
				Weights[i] = LegWeight95[i];
			}
			return 0;
		case 100:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LegAbs100[i];
				Weights[i] = LegWeight100[i];
			}
			return 0;
	}

	int m = (n+1)/2;  // The roots are symmetric about the origin, so we only have to find half of them.

	for (int i = 1; i <= m; i++) {
		long double xprev = 0.0, x, xprevprev;
		long double pprime;  // Derivative of P_n(x) at x.
		long double p1, p2, p3;  // p1 = P_n, p2 = P_{n-1}, and p3 = temporary placeholder.

		// Approximation to the roots from Numerical Recipes in Fortran, Second Edition:
		//  x = cos(pi*(i-1/4)/(n+1/2)),
		//  where i is the ith root, and n is from P_n(x).
		x = cos(PI*(i-0.25)/(n+.5));
		//		cout << (double)x << endl;

		do {
			p1 = 1.0;  // P_0(x)
			p2 = 0.0;  // Using P_{-1}(x) as 0 (equation 4.5.6 of Numerical Methods, page 142).
			// Loop up the recurrence relation to get the Legendre polynomial evaluated at x.
			for (int j = 0; j < n; j++) {
				p3 = p2;
				p2 = p1;
				// Recurrence relation for the Legendre polynomials:
				//  (n+1)*P_{n+1}(x) = (2n+1)*x*P_n(x) - n*P_{n-1}(x)
				//  from Abramowitz and Stegun page 333, equation 8.5.3 with \mu = 0.
				//  This also appears on http://mathworld.wolfram.com/LegendrePolynomial.html.
				p1 = ((2.0*j+1.0)*x*p2 - j*p3)/(j+1);
			}

			// Recurrence relation for the Legendre polynomials involving the first derivative:
			//  (1-x^2)P_n'(x) = -n*x*P_n(x) + n*P_{n-1}(x) = (n+1)*x*P_n(x) - (n+1)*P_{n+1}(x)
			//  from Abramowitz and Stegun page 334, equation 8.5.4 with \mu = 0
			//  and also on http://mathworld.wolfram.com/LegendrePolynomial.html.
			pprime = n*(x*p1-p2)/(x*x-1.0);

			xprevprev = xprev;
			xprev = x;  // Keep track of the previous value of x for error analysis.
			// Newton's method:
			//  x_{n+1} = x_n - f(x_n)/f'(x_n).
			x = x - p1 / pprime;
			//cout << abs((double)x-(double)xprev) << endl;

			// When we get to the limits of double precision, x and xprev will sometimes swap values continuously.
			if (xprevprev == x)
				break;
		} while(fabs(x-xprev) > Tolerance);

		Abscissas[i-1] = x;  // Root is calculated at this point.

		// Weights are calculated from Abramowitz and Stegun page 887, equation 25.4.29:
		//  w_i = 2/((1-x_i^2) * P'_n(x_i)^2).
		Weights[i-1] = 2.0 / ((1.0 - x*x) * pprime*pprime);
	}

	// The roots are symmetric about x = 0, so we fill in the other half of the arrays.
	for (int i = n-1; i >= m; i--) {
		Abscissas[i] = -Abscissas[n-i-1];
		Weights[i] = Weights[n-i-1];
	}

	//@TODO: Do this earlier in the function while building the array up!
	vector <long double> AbsTemp(n);
	for (int i = 0; i < n; i++) {
		AbsTemp[i] = Abscissas[n-i-1];
	}
	Abscissas = AbsTemp;

	return 0;
}


// Generates Gauss-Laguerre quadrature abscissas and weights
//!//  This function is based on the gaulag function on pages 145-146 of
//  Numerical Recipes in Fortran, Second Edition.
//@TODO: We could also make use of the fact that for odd-degree polynomials, one root is always 0.
int GaussLaguerre(vector <long double> & Abscissas, vector <long double> & Weights, int n)
{
	double Tolerance = 1e-14;  //@TODO: Adjust tolerance accordingly.

	//@TODO: Remove
	Abscissas.resize(n);
	Weights.resize(n);

	switch (n)
	{
		case 2:
			Abscissas[0] = (double)0.585786437626904951198311275;
			Abscissas[1] = (double)3.414213562373095048801688724;
			Weights[0] = (double)0.853553390593273762200422181;
			Weights[1] = (double)0.146446609406726237799577818;
			return 0;
		case 5:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs5[i];
				Weights[i] = LagWeight5[i];
			}
			return 0;
		case 10:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs10[i];
				Weights[i] = LagWeight10[i];
			}
			return 0;
		case 15:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs15[i];
				Weights[i] = LagWeight15[i];
			}
			return 0;
		case 20:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs20[i];
				Weights[i] = LagWeight20[i];
			}
			return 0;
		case 25:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs25[i];
				Weights[i] = LagWeight25[i];
			}
			return 0;
		case 30:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs30[i];
				Weights[i] = LagWeight30[i];
			}
			return 0;
		case 35:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs35[i];
				Weights[i] = LagWeight35[i];
			}
			return 0;
		case 40:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs40[i];
				Weights[i] = LagWeight40[i];
			}
			return 0;
		case 45:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs45[i];
				Weights[i] = LagWeight45[i];
			}
			return 0;
		case 50:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs50[i];
				Weights[i] = LagWeight50[i];
			}
			return 0;
		case 55:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs55[i];
				Weights[i] = LagWeight55[i];
			}
			return 0;
		case 60:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs60[i];
				Weights[i] = LagWeight60[i];
			}
			return 0;
		case 65:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs65[i];
				Weights[i] = LagWeight65[i];
			}
			return 0;
		case 70:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs70[i];
				Weights[i] = LagWeight70[i];
			}
			return 0;
		case 75:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs75[i];
				Weights[i] = LagWeight75[i];
			}
			return 0;
		case 80:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs80[i];
				Weights[i] = LagWeight80[i];
			}
			return 0;
		case 85:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs85[i];
				Weights[i] = LagWeight85[i];
			}
			return 0;
		case 90:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs90[i];
				Weights[i] = LagWeight90[i];
			}
			return 0;
		case 95:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs95[i];
				Weights[i] = LagWeight95[i];
			}
			return 0;
		case 100:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs100[i];
				Weights[i] = LagWeight100[i];
			}
			return 0;
		case 105:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs105[i];
				Weights[i] = LagWeight105[i];
			}
			return 0;
		case 110:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs110[i];
				Weights[i] = LagWeight110[i];
			}
			return 0;
		case 115:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs115[i];
				Weights[i] = LagWeight115[i];
			}
		case 120:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs120[i];
				Weights[i] = LagWeight120[i];
			}
		case 125:
			for (int i = 0; i < n; i++) {
				Abscissas[i] = LagAbs125[i];
				Weights[i] = LagWeight125[i];
			}
		return 0;
	}

	long double x = 0.0;
	for (int i = 1; i <= n; i++) {
		long double xprev = 0.0, xprevprev;
		long double pprime;  // Derivative of P_n(x) at x.
		long double p1, p2, p3;  // p1 = P_n, p2 = P_{n-1}, and p3 = temporary placeholder.
		int ai;

		// Approximation to the roots from Numerical Recipes in Fortran, Second Edition:
		//  x = cos(pi*(i-1/4)/(n+1/2)),
		//  where i is the ith root, and n is from P_n(x).
		//		x = cos(PI*(i-0.25)/(n+.5));

		switch(i)
		{
			case 1:
				x = 3.0 / (1.0+2.4*n);
				break;
			case 2:
				x = x + 15.0 / (1.0 + 2.5*n);
				break;
			default:
				ai = i - 2;
				x = x + ((1.0+2.55*ai)/(1.9*ai))*(x-Abscissas[i-3]);
		}

		do {
			p1 = 1.0;  // P_0(x)
			p2 = 0.0;  // Using P_{-1}(x) as 0 (equation 4.5.6 of Numerical Methods, page 142).
			// Loop up the recurrence relation to get the Legendre polynomial evaluated at x.
			for (int j = 0; j < n; j++) {
				p3 = p2;
				p2 = p1;
				// Recurrence relation for the Laguerre polynomials:
				//  (n+1)*L_{k+1}(x) = (2n+1-x)*L_k(x) - n*L_{n-1}(x)
				//!				//  from Abramowitz and Stegun page 333, equation 8.5.3 with \mu = 0.
				//  This also appears on http://mathworld.wolfram.com/LaguerrePolynomial.html.
				p1 = ((2.0*j+1.0-x)*p2 - j*p3)/(j+1);
			}

			// Recurrence relation for the Legendre polynomials involving the first derivative:
			//  x*P_n'(x) = n*L_n(x) - n*L_{n-1}(x)
			//  from Abramowitz and Stegun page 334, equation 8.5.4 with \mu = 0
			//  and also on http://mathworld.wolfram.com/LegendrePolynomial.html.
			pprime = n*(p1-p2)/x;

			xprevprev = xprev;
			xprev = x;  // Keep track of the previous value of x for error analysis.
			// Newton's method:
			//  x_{n+1} = x_n - f(x_n)/f'(x_n).
			x = x - p1 / pprime;
			//cout << abs((double)x-(double)xprev) << endl;

			// When we get to the limits of double precision, x and xprev will sometimes swap values continuously.
			if (xprevprev == x)
				break;
		} while(fabs(x-xprev) > Tolerance);

		Abscissas[i-1] = x;  // Root is calculated at this point.

		//!		// Weights are calculated from Abramowitz and Stegun page 887, equation 25.4.29:
		//!		//  w_i = 2/((1-x_i^2) * P'_n(x_i)^2).
		//Weights[i-1] = 2.0 / ((1.0 - x*x) * pprime*pprime);
		//Weights[i-1] = x / ((n+1)*(n+1)*
		Weights[i-1] = 1.0 / (x*pprime*pprime);
	}

	return 0;
}


//@TODO: Take out hardcoded values for "BETAH" and "BETAPS"
//@TODO: Actually use a CuspR1.
double GaussIntegrationPhi23Perimetric_SLC(int nX, int nY, int nZ, int nR3Leg, int nR3Lag, int nR13, int nPhi23, double CuspR3, double kappa, double mu, int sf)
{
	double xSum = 0.0, ySum, zSum, r3r13Sum;
	double x, y, z;
	double r1, r2;
	double r12;
	vector <long double> LaguerreAbscissasX, LaguerreWeightsX;
	vector <long double> LaguerreAbscissasY, LaguerreWeightsY;
	vector <long double> LaguerreAbscissasZ, LaguerreWeightsZ;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3;
	vector <long double> LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	int Prog = 0;

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(LaguerreAbscissasX, LaguerreWeightsX, nX);
	GaussLaguerre(LaguerreAbscissasY, LaguerreWeightsY, nY);
	GaussLaguerre(LaguerreAbscissasZ, LaguerreWeightsZ, nZ);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);

	#pragma omp parallel for shared(xSum) private(x,y,z,r1,r2,r12,ySum,zSum,r3r13Sum) schedule(guided,1)
	for (int i = 0; i < nX; i++) {  // x integration
		WriteProgress(string("CLS"), Prog, nX);
		x = LaguerreAbscissasX[i] / 0.5;  //@TODO: Above

		ySum = 0.0;
		for (int j = 0; j < nY; j++) {  // y integration
			y = LaguerreAbscissasY[j] / (0.5 * (1.0 + 0.5));

			zSum = 0.0;
			for (int k = 0; k < nZ; k++) {  // z integration
				z = LaguerreAbscissasZ[k] / (0.5 * 0.5);
				// Relates our inter-particle coordinates to the perimetric coordinates for r1, r2 and r12.
				r1 = 0.5 * (x + z);
				r2 = 0.5 * (x + y);
				r12 = 0.5 * (y + z);

				double TempCoeff = LaguerreWeightsX[i] * LaguerreWeightsY[j] * LaguerreWeightsZ[k];
				//r3r13Sum = r3r13phi23IntegrationPerimetric(TempCoeff, Powers, r1, r2, r12, kappa, mu, sf, LegendreAbscissasR3, LegendreWeightsR3, LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR13, LegendreWeightsR13, nR3Leg, nR3Lag, nR13, nPhi23, CuspR3, f);


				double a13, b13;
				long double r3Sum = 0.0, r13Sum, r3;
				vector <long double> r13, r3Abscissas, r3Weights;
				double r23, Phi23, Cos12, Sin12, Cos13, Sin13, rho, rhop;
				long double Phi23Sum;
				int NumR3Points;

				Cos12 = (r1*r1 + r2*r2 - r12*r12) / (2.0*r1*r2);
				Sin12 = sqrt(1.0 - Cos12*Cos12);
				rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12*r12);

				r3Abscissas.resize(nR3Leg + nR3Lag);
				r3Weights.resize(nR3Leg + nR3Lag);
				r13.resize(nR13);

				if (r1 > CuspR3) {
					NumR3Points = nR3Lag;
					for (int a = 0; a < nR3Lag; a++) {  //@TODO: Can we just do this as a single line ala Fortran?
						r3Abscissas[a] = LaguerreAbscissasR3[a];
						r3Weights[a] = LaguerreWeightsR3[a] / 1.0;
					}
				}
				else {
					ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0, r1);
					NumR3Points = nR3Leg + nR3Lag;
					for (int a = 0; a < nR3Leg; a++) {
						r3Weights[a] = LegendreWeightsR3[a] * exp(-1.0*r3Abscissas[a]) * (r1-0.0)/2.0;
					}
					for (int a = nR3Leg; a < NumR3Points; a++) {
						r3Abscissas[a] = (LaguerreAbscissasR3[a-nR3Leg] + r1*1.0) / 1.0;
						r3Weights[a] = LaguerreWeightsR3[a-nR3Leg] * exp(-1.0*r1) / 1.0;
					}
				}

				for (int a = 0; a < NumR3Points; a++) {  // r3 integration
					r3 = r3Abscissas[a];
					r13Sum = 0.0;
					a13 = fabs(r1-r3);
					b13 = fabs(r1+r3);
					ChangeOfIntervalNoResize(LegendreAbscissasR13, r13, a13, b13);

					for (int b = 0; b < nR13; b++) {  // r13 integration
						Cos13 = (r1*r1 + r3*r3 - r13[b]*r13[b]) / (2.0*r1*r3);
						Sin13 = sqrt(1.0 - Cos13*Cos13);
						Phi23Sum = 0.0;
						rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13[b]*r13[b]);
						for (int c = 1; c <= nPhi23; c++) {  // phi_23 integration
							Phi23 = (2.0*c - 1.0)*PI/(2.0*nPhi23);
							r23 = sqrt(r2*r2 + r3*r3 - 2.0*r2*r3*(Sin12*Sin13*cos(Phi23) + Cos12*Cos13));
							//Phi23Sum += CalcCLS(Powers, r1, r2, r3, r12, r13[b], r23, kappa, rho, rhop, mu, sf) * r2 * r3 * r12 * r13[b];

							//double Ret = CP23(r2, r13[b], kappa, rhop, mu) * expl(0.5*r12+r3+r2) * (2.0/r1 - 2.0/r2 - 2.0/r13[b]) * S(r3, r12, kappa, rho);
							//Ret = sf*Ret;

							double ExpMuRho = exp(-mu*rho);
							double S = sin(kappa*rhop) * sqrt(kappa) / (kappa*rhop);
							double C = cos(kappa*rho) / (kappa*rho) * (1.0 - ExpMuRho*(1.0 + 0.50*mu*rho));
							double ExpR12R3 = exp(-(0.5*r12 + r3));
							double CosFactor = cos(kappa*rho) / (kappa*rho);
							double Part = ExpMuRho * mu*mu*mu * rho / 4.0 * CosFactor + ExpMuRho / 2.0 * kappa * mu * (1 + mu*rho) * sin(kappa*rho) / (kappa*rho);
							double Ret = S * ExpR12R3 * ((2.0/r1 - 2.0/r2 - 2.0/r13[b]) * C + Part);

							Phi23Sum += Ret * r2 * r3 * r12 * r13[b];
						}
						Phi23Sum /= nPhi23;  //@TODO: Move to return statement.
						r13Sum += LegendreWeightsR13[b] * Phi23Sum * (b13-a13)/2.0;
					}
					r3Sum += r3Weights[a] * r13Sum;
				}

				//@TODO: Does the / 4.0 need to be here?
				r3r13Sum = r3Sum / 4.0;


				zSum += LaguerreWeightsZ[k] * r3r13Sum;
			}
			ySum += LaguerreWeightsY[j] * zSum;
		}
		//@TODO: Need to have atomic?
		//#pragma omp atomic
		xSum += LaguerreWeightsX[i] * ySum;
	}

	xSum /= (0.5 * 0.5 * (1.0 + 0.5) * 0.5 * 0.5);

	return xSum;
}


double GaussIntegrationPhi12_SLC_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nPhi12, int nR13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int sf)
{
	double r1Sum, r2Sum, r3Sum, r12r13Sum;
	double r1, r2, r3;
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR23, LegendreWeightsR23;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r3Abscissas, r2Weights, r3Weights;
	int NumR2Points, NumR3Points, Prog = 0;

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(r1Abscissas, r1Weights, nR1);
	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR23, LegendreWeightsR23, nR23);
	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);

	r1Abscissas.resize(nR1);
	r1Weights.resize(nR1);

	r1Sum = 0.0;
	#pragma omp parallel for shared(r1Sum,r1Abscissas,r1Weights) private(r1,r2,r3,r2Sum,r3Sum,r12r13Sum,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		WriteProgress(string("CLS R23"), Prog, nR1);

		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);

		r1 = r1Abscissas[i];
		r2Sum = 0.0;

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
				r3Weights[m] = LegendreWeightsR3[m] * exp(-r3Abscissas[m]) * (r1-0.0)/2.0 /** 4.0*/;
			}
			for (int m = nR3Leg; m < NumR3Points; m++) {
				r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r1);
				r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-r1);
			}
		}

		r3Sum = 0.0;
		for (int j = 0; j < NumR3Points; j++) {  // r3 integration
			r3 = r3Abscissas[j];
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
					r2Weights[m] = LegendreWeightsR2[m] * exp(-r2Abscissas[m]) * (r3-0.0)/2.0 /** 4.0*/;
				}
				for (int m = nR2Leg; m < NumR2Points; m++) {
					r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r3);
					r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-r3);
				}
			}

			r2Sum = 0.0;
			for (int g = 0; g < NumR2Points; g++) {  // r2 integration
				r2 = r2Abscissas[g];
				double TempCoeff = r3Weights[j] * r2Weights[g] * r1Weights[i];
				// Integrations over the 3 interparticle coordinates
				//r12r13Sum = r23r13phi12Integration(TempCoeff, Powers, r1, r2, r3, kappa, mu, sf, LegendreAbscissasR23, LegendreWeightsR23, LegendreAbscissasR13, LegendreWeightsR13, nPhi12, nR13, nR23, f);
				long double ExpR1R2R3 = expl(r1+r2+r3);


				double a23, b23, a13, b13;
				double r23Sum = 0.0, r13Sum;
				vector <long double> r23, r13;
				double r12, Phi12, Cos23, Sin23, Cos13, Sin13, rho, rhop;
				long double Phi12Sum;

				a23 = fabs(r2-r3);
				b23 = fabs(r2+r3);
				a13 = fabs(r1-r3);
				b13 = fabs(r1+r3);
				ChangeOfInterval(LegendreAbscissasR23, r23, a23, b23);
				ChangeOfInterval(LegendreAbscissasR13, r13, a13, b13);

				for (int k = 0; k < nR23; k++) {  // r23 integration
					r13Sum = 0.0;
					Cos23 = (r2*r2 + r3*r3 - r23[k]*r23[k]) / (2.0*r2*r3);
					Sin23 = sqrt(1.0 - Cos23*Cos23);
					for (int p = 0; p < nR13; p++) {  // r13 integration
						Cos13 = (r1*r1 + r3*r3 - r13[p]*r13[p]) / (2.0*r1*r3);
						Sin13 = sqrt(1.0 - Cos13*Cos13);
						Phi12Sum = 0.0;
						rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13[p]*r13[p]);

						for (int m = 1; m <= nPhi12; m++) {  // phi_12 integration
							Phi12 = (2.0*m - 1.0)*PI/(2.0*nPhi12);
							r12 = sqrt(r1*r1 + r2*r2 - 2.0*r1*r2*(Sin13*Sin23*cos(Phi12) + Cos13*Cos23));
							rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12*r12);
							//Phi12Sum += CalcCLS_R23Term(Powers, r1, r2, r3, r12, r13[p], r23[k], kappa, rho, rhop, mu, sf) * r1 * r2 * r23[k] * r13[p];


							//@TODO: Could reduce calculations by combining S and C terms.
							double Ret;
							Ret = CP23(r2, r13[p], kappa, rhop, mu) * S(r3, r12, kappa, rho);
							Ret += C(r3, r12, kappa, rho, mu) * SP23(r2, r13[p], kappa, rhop);
							Ret *= 2.0/r23[k] * ExpR1R2R3 * sf;

							Phi12Sum += Ret * r1 * r2 * r23[k] * r13[p];
						}
						Phi12Sum /= nPhi12;
						r13Sum += LegendreWeightsR13[p] * Phi12Sum * (b13-a13)/2.0;
					}
					r23Sum += LegendreWeightsR23[k] * r13Sum * (b23-a23)/2.0;
				}

				//@TODO: Does the / 2.0 need to be here?
				//return r23Sum / 4.0;
				r12r13Sum = r23Sum / 2.0;


				r2Sum += r2Weights[g] * r12r13Sum;
			}

			r3Sum += r3Weights[j] * r2Sum;
		}

		//@TODO: Do I really need the atomic?
		//#pragma omp atomic
		r1Sum += r1Weights[i] * r3Sum;
	}

	return r1Sum;
}



void GaussIntegrationPhi23Perimetric_LongLong(int nX, int nY, int nZ, int nR3Leg, int nR3Lag, int nR13, double CuspR2, double CuspR3, double kappa, double mu, int sf, double &CLC, double &CLS, double &SLS)
{
	double xSumCLC, xSumCLS, xSumSLS;
	double x, y, z;
	double r1, r2;
	double r12;
	vector <long double> LaguerreAbscissasX, LaguerreWeightsX;
	vector <long double> LaguerreAbscissasY, LaguerreWeightsY;
	vector <long double> LaguerreAbscissasZ, LaguerreWeightsZ;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3;
	vector <long double> LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(LaguerreAbscissasX, LaguerreWeightsX, nX);
	GaussLaguerre(LaguerreAbscissasY, LaguerreWeightsY, nY);
	GaussLaguerre(LaguerreAbscissasZ, LaguerreWeightsZ, nZ);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);

	xSumCLC = 0.0, xSumCLS = 0.0, xSumSLS = 0.0;
	#pragma omp parallel for shared(xSumCLC, xSumCLS, xSumSLS) private(x,y,z,r1,r2,r12) schedule(guided,1)
	for (int i = 0; i < nX; i++) {  // x integration
		//WriteProgress(i, nX);
		x = LaguerreAbscissasX[i] / 0.5;  //@TODO: Above

		double ySumCLC = 0.0, ySumCLS = 0.0, ySumSLS = 0.0;
		for (int j = 0; j < nY; j++) {  // y integration
			y = LaguerreAbscissasY[j] / (0.5 * (1.0 + 0.5));

			double zSumCLC = 0.0, zSumCLS = 0.0, zSumSLS = 0.0;
			for (int k = 0; k < nZ; k++) {  // z integration
				z = LaguerreAbscissasZ[k] / (0.5 * 0.5);
				// Relates our inter-particle coordinates to the perimetric coordinates for r1, r2 and r12.
				r1 = 0.5 * (x + z);
				r2 = 0.5 * (x + y);
				r12 = 0.5 * (y + z);

				double TempCoeff = LaguerreWeightsX[i] * LaguerreWeightsY[j] * LaguerreWeightsZ[k];
				long double r3;
				vector <long double> r13, r3Abscissas, r3Weights;
				int NumR3Points;
				double rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12*r12);

				r3Abscissas.resize(nR3Leg + nR3Lag);
				r3Weights.resize(nR3Leg + nR3Lag);
				r13.resize(nR13);

				if (r1 > CuspR3) {
					NumR3Points = nR3Lag;
					for (int a = 0; a < nR3Lag; a++) {  //@TODO: Can we just do this as a single line ala Fortran?
						r3Abscissas[a] = LaguerreAbscissasR3[a];
						r3Weights[a] = LaguerreWeightsR3[a];
					}
				}
				else {
					ChangeOfIntervalNoResize(LegendreAbscissasR3, r3Abscissas, 0.0, r1);
					NumR3Points = nR3Leg + nR3Lag;
					for (int a = 0; a < nR3Leg; a++) {
						r3Weights[a] = LegendreWeightsR3[a] * exp(-r3Abscissas[a]) * r1 / 2.0;
					}
					for (int a = nR3Leg; a < NumR3Points; a++) {
						r3Abscissas[a] = (LaguerreAbscissasR3[a-nR3Leg] + r1);
						r3Weights[a] = LaguerreWeightsR3[a-nR3Leg] * exp(-r1);
					}
				}

				long double r3SumCLC = 0.0, r3SumCLS = 0.0, r3SumSLC = 0.0, r3SumSLS = 0.0;
				for (int a = 0; a < NumR3Points; a++) {  // r3 integration
					r3 = r3Abscissas[a];
					double a13 = fabs(r1-r3);
					double b13 = fabs(r1+r3);
					ChangeOfIntervalNoResize(LegendreAbscissasR13, r13, a13, b13);
					//double Ret1 = /*exp(-(r12/2.0 + r3)) **/ /*sqrt(kappa)*/ kappa * sin(kappa*rho) / (kappa*rho);

					long double r13SumCLC = 0.0, r13SumCLS = 0.0, r13SumSLC = 0.0, r13SumSLS = 0.0;
					for (int b = 0; b < nR13; b++) {  // r13 integration
						double rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13[b]*r13[b]);
						double SqrtKappa = sqrt(kappa);
						double ExpMuRho = exp(-mu*rho);
						double ExpMuRhop = exp(-mu*rhop);
						double Pot = (2.0/r1 - 2.0/r2 - 2.0/r13[b]);
						double CLC_C, CLC_CP23, CLC_S;
						
						// CLC
						CLC_C = cos(kappa*rho) * SqrtKappa / (kappa*rho) * (1.0 - ExpMuRho*(1.0 + 0.50*mu*rho));
						CLC_CP23 = cos(kappa*rhop) * SqrtKappa / (kappa*rhop) * (1.0 - ExpMuRhop*(1.0 + 0.50*mu*rhop));
						CLC_S = sin(kappa*rho) * SqrtKappa / (kappa*rho);

						double RetCLC1 = CLC_C * exp(1.0*(r2-r3)-0.5*r12) *
							(kappa*mu/2.0*ExpMuRho*(1.0 + mu*rho) * CLC_S + SqrtKappa*ExpMuRho*mu*mu*mu/4.0*cos(kappa*rho)/kappa);
						double RetCLC2 = CLC_CP23 * exp(-0.5*r13[b]) *
							(kappa*mu/2.0*ExpMuRho*(1.0 + mu*rho) * CLC_S + SqrtKappa*ExpMuRho*mu*mu*mu/4.0*cos(kappa*rho)/kappa + (2.0/r1 - 2.0/r2 - 2.0/r13[b]) * CLC_C);
						double Phi23SumCLC = (RetCLC1 + sf*RetCLC2) * r2 * r3 * r12 * r13[b];
						r13SumCLC += LegendreWeightsR13[b] * Phi23SumCLC * (b13-a13)/2.0;

						// CLS
						double RetCLS = CP23(r2, r13[b], kappa, rhop, mu) * expl(0.5*r12+r3+r2) * (2.0/r1 - 2.0/r2 - 2.0/r13[b]) * S(r3, r12, kappa, rho);
						double Phi23SumCLS = sf * RetCLS * r2 * r3 * r12 * r13[b];
						r13SumCLS += LegendreWeightsR13[b] * Phi23SumCLS * (b13-a13)/2.0;

						// SLS
						double RetSLS1 = kappa * sin(kappa*rho) / (kappa*rho);
						double RetSLS2 = exp(-(r13[b]/2.0)) * sin(kappa*rhop) / (kappa*rhop);
						double RetSLS = RetSLS1 * RetSLS2 * Pot;
						double Phi23SumSLS = sf * RetSLS * r2 * r3 * r12 * r13[b];
						r13SumSLS += LegendreWeightsR13[b] * Phi23SumSLS * (b13-a13)/2.0;
					}
					r3SumCLC += r3Weights[a] * r13SumCLC;
					r3SumCLS += r3Weights[a] * r13SumCLS;
					r3SumSLS += r3Weights[a] * r13SumSLS;
				}
				zSumCLC += LaguerreWeightsZ[k] * r3SumCLC / 4.0;
				zSumCLS += LaguerreWeightsZ[k] * r3SumCLS / 4.0;
				zSumSLS += LaguerreWeightsZ[k] * r3SumSLS / 4.0;
			}
			ySumCLC += LaguerreWeightsY[j] * zSumCLC;
			ySumCLS += LaguerreWeightsY[j] * zSumCLS;
			ySumSLS += LaguerreWeightsY[j] * zSumSLS;
		}
		xSumCLC += LaguerreWeightsX[i] * ySumCLC;
		xSumCLS += LaguerreWeightsX[i] * ySumCLS;
		xSumSLS += LaguerreWeightsX[i] * ySumSLS;
	}
	//@TODO: Need to have atomic?
	//#pragma omp atomic
	CLC = xSumCLC / (0.5 * 0.5 * (1.0 + 0.5) * 0.5 * 0.5);
	CLS = xSumCLS / (0.5 * 0.5 * (1.0 + 0.5) * 0.5 * 0.5);
	SLS = xSumSLS/ (0.5 * 0.5 * (1.0 + 0.5) * 0.5 * 0.5);

	return;
}


void GaussIntegrationPhi12_LongLong_R23Term(int nR1, int nR2Leg, int nR2Lag, int nR3Leg, int nR3Lag, int nPhi12, int nR13, int nR23, double CuspR2, double CuspR3, double kappa, double mu, int sf, double &CLC, double &CLS, double &SLS)
{
	double r1, r2, r3;
	vector <long double> LaguerreAbscissasR2, LaguerreWeightsR2, LegendreAbscissasR2, LegendreWeightsR2;
	vector <long double> LaguerreAbscissasR3, LaguerreWeightsR3, LegendreAbscissasR3, LegendreWeightsR3;
	vector <long double> LegendreAbscissasR23, LegendreWeightsR23;
	vector <long double> LegendreAbscissasR13, LegendreWeightsR13;
	vector <long double> r1Abscissas, r1Weights, r2Abscissas, r3Abscissas, r2Weights, r3Weights;
	int NumR2Points, NumR3Points;
	double SqrtKappa = sqrt(kappa);

	// Create the abscissas and weights for the needed number of points.
	GaussLaguerre(r1Abscissas, r1Weights, nR1);
	GaussLegendre(LegendreAbscissasR2, LegendreWeightsR2, nR2Leg);
	GaussLaguerre(LaguerreAbscissasR2, LaguerreWeightsR2, nR2Lag);
	GaussLegendre(LegendreAbscissasR3, LegendreWeightsR3, nR3Leg);
	GaussLaguerre(LaguerreAbscissasR3, LaguerreWeightsR3, nR3Lag);
	GaussLegendre(LegendreAbscissasR23, LegendreWeightsR23, nR23);
	GaussLegendre(LegendreAbscissasR13, LegendreWeightsR13, nR13);

	r1Abscissas.resize(nR1);
	r1Weights.resize(nR1);
	for (int m = 0; m < nR1; m++) {  //@TODO: Can we just do this as a single line ala Fortran?
		r1Abscissas[m] = r1Abscissas[m];
		r1Weights[m] = r1Weights[m];
	}

	long double r1SumCLC = 0.0, r1SumCLS = 0.0, r1SumSLC = 0.0, r1SumSLS = 0.0;
	#pragma omp parallel for shared(r1SumCLC,r1SumCLS,r1SumSLC,r1SumSLS,r1Abscissas,r1Weights) private(r1,r2,r3,r2Abscissas,r2Weights,r3Abscissas,r3Weights,NumR2Points,NumR3Points) schedule(guided,1)
	for (int i = 0; i < nR1; i++) {  // r1 integration
		//WriteProgress(i, nR1);

		r2Abscissas.resize(nR2Leg + nR2Lag);
		r2Weights.resize(nR2Leg + nR2Lag);
		r3Abscissas.resize(nR3Leg + nR3Lag);
		r3Weights.resize(nR3Leg + nR3Lag);

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
				r3Weights[m] = LegendreWeightsR3[m] * exp(-r3Abscissas[m]) * (r1-0.0)/2.0 /** 4.0*/;
			}
			for (int m = nR3Leg; m < NumR3Points; m++) {
				r3Abscissas[m] = (LaguerreAbscissasR3[m-nR3Leg] + r1);
				r3Weights[m] = LaguerreWeightsR3[m-nR3Leg] * exp(-r1);
			}
		}

		long double r3SumCLC = 0.0, r3SumCLS = 0.0, r3SumSLC = 0.0, r3SumSLS = 0.0;
		for (int j = 0; j < NumR3Points; j++) {  // r3 integration
			r3 = r3Abscissas[j];
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
					r2Weights[m] = LegendreWeightsR2[m] * exp(-r2Abscissas[m]) * (r3-0.0)/2.0 /** 4.0*/;
				}
				for (int m = nR2Leg; m < NumR2Points; m++) {
					r2Abscissas[m] = (LaguerreAbscissasR2[m-nR2Leg] + r3);
					r2Weights[m] = LaguerreWeightsR2[m-nR2Leg] * exp(-r3);
				}
			}

			long double r2SumCLC = 0.0, r2SumCLS = 0.0, r2SumSLC = 0.0, r2SumSLS = 0.0;
			for (int g = 0; g < NumR2Points; g++) {  // r2 integration
				r2 = r2Abscissas[g];
				double TempCoeff = r3Weights[j] * r2Weights[g] * r1Weights[i];
				long double ExpR1R2R3 = expl(r1+r2+r3);
				vector <long double> r23, r13;
				double r12, Phi12, Cos23, Sin23, Cos13, Sin13, rho, rhop;

				double a23 = fabs(r2-r3);
				double b23 = fabs(r2+r3);
				double a13 = fabs(r1-r3);
				double b13 = fabs(r1+r3);
				ChangeOfInterval(LegendreAbscissasR23, r23, a23, b23);
				ChangeOfInterval(LegendreAbscissasR13, r13, a13, b13);

				long double r23SumCLC = 0.0, r23SumCLS = 0.0, r23SumSLC = 0.0, r23SumSLS = 0.0;
				for (int k = 0; k < nR23; k++) {  // r23 integration
					Cos23 = (r2*r2 + r3*r3 - r23[k]*r23[k]) / (2.0*r2*r3);
					Sin23 = sqrt(1.0 - Cos23*Cos23);

					long double r13SumCLC = 0.0, r13SumCLS = 0.0, r13SumSLC = 0.0, r13SumSLS = 0.0;
					for (int p = 0; p < nR13; p++) {  // r13 integration
						Cos13 = (r1*r1 + r3*r3 - r13[p]*r13[p]) / (2.0*r1*r3);
						Sin13 = sqrt(1.0 - Cos13*Cos13);
						rhop = 0.5 * sqrt(2.0*(r1*r1 + r3*r3) - r13[p]*r13[p]);

						long double Phi12SumCLC = 0.0, Phi12SumCLS = 0.0, Phi12SumSLS = 0.0;
						for (int m = 1; m <= nPhi12; m++) {  // phi_12 integration
							Phi12 = (2.0*m - 1.0)*PI/(2.0*nPhi12);
							r12 = sqrt(r1*r1 + r2*r2 - 2.0*r1*r2*(Sin13*Sin23*cos(Phi12) + Cos13*Cos23));
							rho = 0.5 * sqrt(2.0*(r1*r1 + r2*r2) - r12*r12);

							// CLC
							double RetCLC1 = C(r3, r12, kappa, rho, mu);
							double RetCLC2 = CP23(r2, r13[p], kappa, rhop, mu);
							double RetCLC = RetCLC1 * RetCLC2;
							RetCLC *= exp(r1) * exp(r2+r3) * 4.0 / r23[k];

							// CLS
							double RetCLS = CP23(r2, r13[p], kappa, rhop, mu) * S(r3, r12, kappa, rho);
							RetCLS += C(r3, r12, kappa, rho, mu) * SP23(r2, r13[p], kappa, rhop);
							RetCLS *= 2.0/r23[k] * ExpR1R2R3;

							// SLS
							double RetSLS1 = exp(-(r13[p]/2.0 + r2)) * SqrtKappa * sin(kappa*rhop) / (kappa*rhop); //SP23(r2, r13[p], kappa, rhop);
							double RetSLS2 = exp(-(r12/2.0 + r3)) * SqrtKappa * sin(kappa*rho) / (kappa*rho); //S(r3, r12, kappa, rho);
							double RetSLS = RetSLS1 * ExpR1R2R3 * RetSLS2 / r23[k] * 4.0;

							Phi12SumCLC += sf*RetCLC * r1 * r2 * r23[k] * r13[p] / nPhi12;
							Phi12SumCLS += sf*RetCLS * r1 * r2 * r23[k] * r13[p] / nPhi12;
							Phi12SumSLS += sf*RetSLS * r1 * r2 * r23[k] * r13[p] / nPhi12;
						}
						r13SumCLC += LegendreWeightsR13[p] * Phi12SumCLC * (b13-a13)/2.0;
						r13SumCLS += LegendreWeightsR13[p] * Phi12SumCLS * (b13-a13)/2.0;
						r13SumSLS += LegendreWeightsR13[p] * Phi12SumSLS * (b13-a13)/2.0;
					}
					r23SumCLC += LegendreWeightsR23[k] * r13SumCLC * (b23-a23)/2.0;
					r23SumCLS += LegendreWeightsR23[k] * r13SumCLS * (b23-a23)/2.0;
					r23SumSLS += LegendreWeightsR23[k] * r13SumSLS * (b23-a23)/2.0;
				}
				r2SumCLC += r2Weights[g] * r23SumCLC / 2.0;
				r2SumCLS += r2Weights[g] * r23SumCLS / 2.0;
				r2SumSLS += r2Weights[g] * r23SumSLS / 2.0;
			}
			r3SumCLC += r3Weights[j] * r2SumCLC;
			r3SumCLS += r3Weights[j] * r2SumCLS;
			r3SumSLS += r3Weights[j] * r2SumSLS;
		}
		//@TODO: Do I really need the atomic?
		//#pragma omp atomic
		r1SumCLC += r1Weights[i] * r3SumCLC;
		r1SumCLS += r1Weights[i] * r3SumCLS;
		r1SumSLS += r1Weights[i] * r3SumSLS;
	}

	CLC += r1SumCLC;
	//SLC += r1SumSLC;
	CLS += r1SumCLS;
	SLS += r1SumSLS;

	return;
}
