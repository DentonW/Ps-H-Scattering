//
// Long-Range.cpp: This file contains the asymptotic long-range (S and C) functions.
//

#include <math.h>
#include <iostream>
#include "Ps-H Scattering.h"

//@TODO: In "Ps-H Scattering.h"
using namespace std;
extern double PI;


double S(double r3, double r12, double kappa, double rho)
{
	//double Ret = PsWaveFn(r12) * HWaveFn(r3) * sqrt(kappa) * sin(kappa*rho) / (kappa*rho);
	double Ret = exp(-(r12/2.0 + r3)) * sqrt(kappa) * sin(kappa*rho) / (kappa*rho);
	return Ret;
}


double SP23(double r2, double r13, double kappa, double rhop)
{
	//double Ret = PsWaveFn(r13) * HWaveFn(r2) * sqrt(kappa) * sin(kappa*rhop) / (kappa*rhop);
	double Ret = exp(-(r13/2.0 + r2)) * sqrt(kappa) * sin(kappa*rhop) / (kappa*rhop);
	return Ret;
}


double C(double r3, double r12, double kappa, double rho, double mu)
{
	return sqrt(kappa) * PsWaveFn(r12) * HWaveFn(r3) * cos(kappa*rho) / (kappa*rho)
		* (1.0 - exp(-mu*rho)*(1.0 + mu/2.0*rho));
}


double CP23(double r2, double r13, double kappa, double rhop, double mu)
{
	return sqrt(kappa) * PsWaveFn(r13) * HWaveFn(r2) * cos(kappa*rhop) / (kappa*rhop)
		* (1.0 - exp(-mu*rhop)*(1.0 + mu/2.0*rhop));
}


double LOnCBar(double r1, double r2, double r3, double r12, double r13, double r23, double kappa, double rho, double rhop, double mu, int s)
{
	long double Ret1, Ret2;
	long double ExpMuRho = exp(-(mu*rho)), ExpMuRhop = exp(-(mu*rhop));
	long double SinKRho, CosKRho, SinKRhop, CosKRhop;

	SinKRho = sin(kappa*rho);  CosKRho = cos(kappa*rho);  
	SinKRhop = sin(kappa*rhop);  CosKRhop = cos(kappa*rhop);  
	

	Ret1 = kappa*mu/2.0*ExpMuRho*(1.0+mu*rho)*SinKRho/(kappa*rho);
	Ret1 += mu*mu*mu*rho/4.0*ExpMuRho*CosKRho/(kappa*rho);
	Ret1 += (2.0/r1-2.0/r2-2.0/r13)*CosKRho/(kappa*rho)*(1.0-ExpMuRho*(1.0+mu*rho/2.0));
	//Ret1 *= PsWaveFn(r12) * HWaveFn(r3);
	Ret1 *= exp(-(r12/2.0 + r3));

	Ret2 = kappa*mu/2.0*ExpMuRhop*(1.0+mu*rhop)*SinKRhop/(kappa*rhop);
	Ret2 += mu*mu*mu*rhop/4.0*ExpMuRhop*CosKRhop/(kappa*rhop);
	Ret2 += (2.0/r1-2.0/r3-2.0/r12)*CosKRhop/(kappa*rhop)*(1.0-ExpMuRhop*(1.0+mu*rhop/2.0));
	//Ret2 *= PsWaveFn(r13) * HWaveFn(r2);
	Ret2 *= exp(-(r13/2.0 + r2));

	long double Ret = (Ret1 + s*Ret2) * sqrt(kappa);
	return Ret;
}
