//
// Gaussian Integration.cpp: Has Gauss-Legendre and Gauss-Laguerre integration routines
//

#include <vector>
#include <iostream>
#include <iomanip>
//#include <math.h>
#include <float.h>
#include <cstdio>
#include "Ps-H Scattering.h"
#include "Gaussian Integration.h"

extern	long double PI;

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
	int Size = Abscissas.size();
	ChangedAbscissas.resize(Size);
	for (int i = 0; i < Size; i++) {
		ChangedAbscissas[i] = Abscissas[i] * (b-a)/2.0 + (a+b)/2.0;
	}
}

void ChangeOfIntervalNoResize(vector <long double> &Abscissas, vector <long double> &ChangedAbscissas, double a, double b)
{
	int Size = Abscissas.size();
	for (int i = 0; i < Size; i++) {
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
