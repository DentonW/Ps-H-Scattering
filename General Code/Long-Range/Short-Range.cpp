//
// Short-Range.cpp: This file contains the short-range (phi) functions.
//

#include <math.h>
#include <iostream>
#include "Ps-H Scattering.h"
#ifndef NO_MPI
	#include <mpi.h>
#endif

//@TODO: In "Ps-H Scattering.h"
using namespace std;
extern double PI;


// From http://www.google.com/codesearch/p?hl=en#6QwnD1uQIUw/Programming/Libraries/Scientific/fxt-2006.05.26.tgz%7CXpA1r3Yktrg/fxt/src/mod/ipow.h&q=function:ipow%20lang:c%2B%2B
double ipow(double a, int ex)
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

// The Phi terms describe the short-range interactions.  I have not included the 1/(4*pi) spherical harmonic in this definition.
double Phi(rPowers &rp, double r1, double r2, double r3, double r12, double r13, double r23)
{
	return /*exp(-(rp.alpha*r1 + rp.beta*r2 + rp.gamma*r3)) **/ ipow(r1,rp.ki) * ipow(r2,rp.li) * ipow(r12,rp.mi)
		* ipow(r3,rp.ni) * ipow(r13,rp.pi) * ipow(r23,rp.qi);
}


// Same as above but with the 2<->3 permutation
double PhiP23(rPowers &rp, double r1, double r2, double r3, double r12, double r13, double r23)
{
	return /*exp(-(rp.alpha*r1 + rp.beta*r3 + rp.gamma*r2)) **/ ipow(r1,rp.ki) * ipow(r3,rp.li) * ipow(r13,rp.mi)
		* ipow(r2,rp.ni) * ipow(r12,rp.pi) * ipow(r23,rp.qi);
}


// Combines Phi and the 2<->3 permutation
double PhiBar(rPowers &Powers, double r1, double r2, double r3, double r12, double r13, double r23, int s)
{
	return Phi(Powers, r1, r2, r3, r12, r13, r23) + s * PhiP23(Powers, r1, r2, r3, r12, r13, r23);
}


// Same as the above Phi function but with a full list of parameters.
double CalcWithPhi(rPowers &rp, double r1, double r2, double r3, double r12, double r13,
				 double r23, double kappa, double rho, double rhop, double mu, int s)
{
	return /*exp(-(rp.alpha*r1 + rp.beta*r2 + rp.gamma*r3)) **/ ipow(r1,rp.ki) * ipow(r2,rp.li) * ipow(r12,rp.mi)
		* ipow(r3,rp.ni) * ipow(r13,rp.pi) * ipow(r23,rp.qi);
}


double CalcWithPhiP23(rPowers &rp, double r1, double r2, double r3, double r12, double r13,
				 double r23, double kappa, double rho, double rhop, double mu, int s)
{
	return /*exp(-(rp.alpha*r1 + rp.beta*r3 + rp.gamma*r2)) **/ ipow(r1,rp.ki) * ipow(r3,rp.li) * ipow(r13,rp.mi)
		* ipow(r2,rp.ni) * ipow(r12,rp.pi) * ipow(r23,rp.qi);
}


double CalcWithPhiBar(rPowers &rp, double r1, double r2, double r3, double r12, double r13,
				 double r23, double kappa, double rho, double rhop, double mu, int s)
{
	return Phi(rp, r1, r2, r3, r12, r13, r23) + s * PhiP23(rp, r1, r2, r3, r12, r13, r23);
}


// Returns the number of terms for a given omega.  This could use the formula for combination with repetition,
//  except then it would be unable to use a restricted set of terms if we needed.
int CalcPowerTableSize(int Omega)
{
	int NumTerms = 0;  // The total number of terms
	int om, ki, li, mi, ni, pi, qi;  // These are the exponents we are determining.

	if (Omega == -1)  // Special case of no short-range terms
		return 0;

	for (om = 0; om <= Omega; om++) {
		for (ki = 0; ki <= Omega; ki++) {
			for (li = 0; li <= Omega; li++) {
				for (mi = 0; mi <= Omega; mi++) {
					for (ni = 0; ni <= Omega; ni++) {
						for (pi = 0; pi <= Omega; pi++) {
							for (qi = 0; qi <= Omega; qi++) {
								if (ki + li + mi + ni + pi + qi == om)
									NumTerms = NumTerms + 1;
							}
						}
					}
				}
			}
		}
	}
	
	return NumTerms;
}


// Same as CalcPowerTableSize, except only returns the number of terms with qi == 0.
int CalcPowerTableSizeQi0(int Omega, int Ordering, int Start, int End)
{
	int NumTerms = 0;  // The total number of terms
	int om, ki, li, mi, ni, pi, qi;  // These are the exponents we are determining.
	int Cur = 0;

	if (Start == -1)
		Start = 0;

	if (Ordering == 0) {  // My ordering
		for (om = 0; om <= Omega; om++) {
			for (ki = 0; ki <= Omega; ki++) {
				for (li = 0; li <= Omega; li++) {
					for (mi = 0; mi <= Omega; mi++) {
						for (ni = 0; ni <= Omega; ni++) {
							for (pi = 0; pi <= Omega; pi++) {
								for (qi = 0; qi <= Omega; qi++) {
									if (ki + li + mi + ni + pi + qi == om) {
										if (qi == 0) {
											if (Cur >= Start)
												NumTerms = NumTerms + 1;
										}
										Cur = Cur + 1;
										if (Cur > End)
											return NumTerms;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {  // Peter Van Reeth's ordering
		int IHDPP1 = Omega + 1;
		int INX = 0;
		for (int I = 1; I <= IHDPP1; I++) {
			for (int I23P1 = 1; I23P1 <= I; I23P1++) {
				int I23 = I23P1 - 1;
				int I12P1M = I - I23;
				for (int I12P1 = 1; I12P1 <= I12P1M; I12P1++) {
					int I12 = I12P1 - 1;
					int I2P1M = I12P1M - I12;
					for (int I2P1 = 1; I2P1 <= I2P1M; I2P1++) {
						int I2 = I2P1-1;
						int I13P1M = I2P1M-I2;
						for (int I13P1 = 1; I13P1 <= I13P1M; I13P1++) {
							int I13 = I13P1 -1;
							int I3P1M = I13P1M -I13;
							for (int I3P1 = 1; I3P1 <= I3P1M; I3P1++) {
								int I3=I3P1-1;
								int I1 = I3P1M -I3P1;

								//cout << Node << " looking for qi == 0 at " << Cur << endl;
								if (Cur >= Start) {
									if (I23 == 0) {
										INX = INX + 1;
										//cout << Node << " found qi == 0 at " << Cur << endl;
									}
								}
								Cur = Cur + 1;
								if (Cur > End)
									return INX;
							}
						}
					}
				}
			}
		}
		NumTerms = INX;
	}

	return NumTerms;
}


// Same as CalcPowerTableSize, except only returns the number of terms with qi > 0.
int CalcPowerTableSizeQiGt0(int Omega, int Ordering, int Start, int End)
{
	int NumTerms = 0;  // The total number of terms
	int om, ki, li, mi, ni, pi, qi;  // These are the exponents we are determining.
	int Cur = 0;

	if (Start == -1)
		Start = 0;

	if (Ordering == 0) {  // My ordering
		for (om = 0; om <= Omega; om++) {
			for (ki = 0; ki <= Omega; ki++) {
				for (li = 0; li <= Omega; li++) {
					for (mi = 0; mi <= Omega; mi++) {
						for (ni = 0; ni <= Omega; ni++) {
							for (pi = 0; pi <= Omega; pi++) {
								for (qi = 0; qi <= Omega; qi++) {
									if (ki + li + mi + ni + pi + qi == om) {
										if (qi > 0) {
											if (Cur >= Start)
												NumTerms = NumTerms + 1;
										}
										Cur = Cur + 1;
										if (Cur > End)
											return NumTerms;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {  // Peter Van Reeth's ordering
		int IHDPP1 = Omega + 1;
		int INX = 0;
		for (int I = 1; I <= IHDPP1; I++) {
			for (int I23P1 = 1; I23P1 <= I; I23P1++) {
				int I23 = I23P1 - 1;
				int I12P1M = I - I23;
				for (int I12P1 = 1; I12P1 <= I12P1M; I12P1++) {
					int I12 = I12P1 - 1;
					int I2P1M = I12P1M - I12;
					for (int I2P1 = 1; I2P1 <= I2P1M; I2P1++) {
						int I2 = I2P1-1;
						int I13P1M = I2P1M-I2;
						for (int I13P1 = 1; I13P1 <= I13P1M; I13P1++) {
							int I13 = I13P1 -1;
							int I3P1M = I13P1M -I13;
							for (int I3P1 = 1; I3P1 <= I3P1M; I3P1++) {
								int I3=I3P1-1;
								int I1 = I3P1M -I3P1;

								//cout << Node << " looking for qi > 0 at " << Cur << endl;
								if (Cur >= Start) {
									if (I23 > 0) {
										INX = INX + 1;
										//cout << Node << " found qi > 0 at " << Cur << endl;
									}
								}
								Cur = Cur + 1;
								if (Cur > End)
									return INX;
							}
						}
					}
				}
			}
		}
		NumTerms = INX;
	}

	return NumTerms;
}


// This subroutine calculates values of k_i, l_i, m_i, n_i, p_i and q_i of
//  (3.14).  The summation of these is given by equation (3.15), but we do not have
//  the restriction on q of being even.
// This is written to always have the terms in order of increasing omega, i.e. the terms
//  for omega = 0 are first, etc.
void GenOmegaPowerTable(int Omega, int l, int Ordering, vector <rPowers> &PowerTable)
{
	int NumTerm = 0;  // The number of the current term
	int om, ki, li, mi, ni, pi, qi;  // These are the exponents we are determining.

	if (Ordering == 0) {  // My ordering
		for (om = 0; om <= Omega; om++) {
			for (ki = 0; ki <= Omega; ki++) {
				for (li = 0; li <= Omega; li++) {
					for (mi = 0; mi <= Omega; mi++) {
						for (ni = 0; ni <= Omega; ni++) {
							for (pi = 0; pi <= Omega; pi++) {
								for (qi = 0; qi <= Omega; qi++) {
									if (ki + li + mi + ni + pi + qi == om) {
										PowerTable[NumTerm].ki = ki;
										PowerTable[NumTerm].li = li;
										PowerTable[NumTerm].mi = mi;
										PowerTable[NumTerm].ni = ni;
										PowerTable[NumTerm].pi = pi;
										PowerTable[NumTerm].qi = qi;
										PowerTable[NumTerm].Index = NumTerm;
										NumTerm = NumTerm + 1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {  // Peter Van Reeth's ordering
		int IHDPP1 = Omega + 1;
		int INX = 0;
		for (int I = 1; I <= IHDPP1; I++) {
			for (int I23P1 = 1; I23P1 <= I; I23P1++) {
				int I23 = I23P1 - 1;
				int I12P1M = I - I23;
				for (int I12P1 = 1; I12P1 <= I12P1M; I12P1++) {
					int I12 = I12P1 - 1;
					int I2P1M = I12P1M - I12;
					for (int I2P1 = 1; I2P1 <= I2P1M; I2P1++) {
						int I2 = I2P1-1;
						int I13P1M = I2P1M-I2;
						for (int I13P1 = 1; I13P1 <= I13P1M; I13P1++) {
							int I13 = I13P1 -1;
							int I3P1M = I13P1M -I13;
							for (int I3P1 = 1; I3P1 <= I3P1M; I3P1++) {
								int I3=I3P1-1;
								int I1 = I3P1M -I3P1;

								PowerTable[INX].ki = I1;
								PowerTable[INX].li = I2;
								PowerTable[INX].mi = I12;
								PowerTable[INX].ni = I3;
								PowerTable[INX].pi = I13;
								PowerTable[INX].qi = I23;
								PowerTable[INX].Index = INX;
								INX = INX + 1;
								NumTerm = INX;
							}
						}
					}
				}
			}
		}
	}

	int size = PowerTable.size();
	for (int i = 0; i < NumTerm; i++) {
		PowerTable[NumTerm+i].ki = PowerTable[i].ki;
		// Increase power of r2 for second symmetry
		PowerTable[NumTerm+i].li = PowerTable[i].li + l;
		PowerTable[NumTerm+i].mi = PowerTable[i].mi;
		PowerTable[NumTerm+i].ni = PowerTable[i].ni;
		PowerTable[NumTerm+i].pi = PowerTable[i].pi;
		PowerTable[NumTerm+i].qi = PowerTable[i].qi;

		// Increase power of r1 for first symmetry
		PowerTable[i].ki = PowerTable[i].ki + l;
	}

	return;
}


// Same as GenOmegaPowerTable, except only returns terms with qi == 0.
void GenOmegaPowerTableQi0(int Omega, int l, int Ordering, vector <rPowers> &PowerTable, int Start, int End)
{
	int NumTerm = 0;  // The number of the current term
	int om, ki, li, mi, ni, pi, qi;  // These are the exponents we are determining.
	int Cur = 0;

	if (Start == -1)
		Start = 0;

	if (Ordering == 0) {  // My ordering
		for (om = 0; om <= Omega; om++) {
			for (ki = 0; ki <= Omega; ki++) {
				for (li = 0; li <= Omega; li++) {
					for (mi = 0; mi <= Omega; mi++) {
						for (ni = 0; ni <= Omega; ni++) {
							for (pi = 0; pi <= Omega; pi++) {
								for (qi = 0; qi <= Omega; qi++) {
									if (ki + li + mi + ni + pi + qi == om) {
										if (Cur >= Start) {
											if (qi == 0) {
												PowerTable[NumTerm].ki = ki;
												PowerTable[NumTerm].li = li;
												PowerTable[NumTerm].mi = mi;
												PowerTable[NumTerm].ni = ni;
												PowerTable[NumTerm].pi = pi;
												PowerTable[NumTerm].qi = qi;
												PowerTable[NumTerm].Index = NumTerm;
												NumTerm = NumTerm + 1;
											}
										}
										Cur = Cur + 1;
										if (Cur > End)
											goto finish;  //@TODO: Do this without a goto statement.
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {  // Peter Van Reeth's ordering
		int IHDPP1 = Omega + 1;
		int INX = 0;
		for (int I = 1; I <= IHDPP1; I++) {
			for (int I23P1 = 1; I23P1 <= I; I23P1++) {
				int I23 = I23P1 - 1;
				int I12P1M = I - I23;
				for (int I12P1 = 1; I12P1 <= I12P1M; I12P1++) {
					int I12 = I12P1 - 1;
					int I2P1M = I12P1M - I12;
					for (int I2P1 = 1; I2P1 <= I2P1M; I2P1++) {
						int I2 = I2P1-1;
						int I13P1M = I2P1M-I2;
						for (int I13P1 = 1; I13P1 <= I13P1M; I13P1++) {
							int I13 = I13P1 -1;
							int I3P1M = I13P1M -I13;
							for (int I3P1 = 1; I3P1 <= I3P1M; I3P1++) {
								int I3=I3P1-1;
								int I1 = I3P1M -I3P1;

								if (Cur >= Start) {
									if (I23 == 0) {  // We only want the qi == 0 terms right now.
										PowerTable[INX].ki = I1;
										PowerTable[INX].li = I2;
										PowerTable[INX].mi = I12;
										PowerTable[INX].ni = I3;
										PowerTable[INX].pi = I13;
										PowerTable[INX].qi = I23;
										PowerTable[INX].Index = INX;
										INX = INX + 1;
									}
								}
								Cur = Cur + 1;
								if (Cur > End)
									goto finish;  //@TODO: Do this without a goto statement.
							}
						}
					}
				}
			}
		}
	}

finish:
	int size = PowerTable.size();
	NumTerm = size / 2;  //@TODO: Use Cur.
	for (int i = 0; i < NumTerm; i++) {
		PowerTable[NumTerm+i].ki = PowerTable[i].ki;
		// Increase power of r2 for second symmetry
		PowerTable[NumTerm+i].li = PowerTable[i].li + l;
		PowerTable[NumTerm+i].mi = PowerTable[i].mi;
		PowerTable[NumTerm+i].ni = PowerTable[i].ni;
		PowerTable[NumTerm+i].pi = PowerTable[i].pi;
		PowerTable[NumTerm+i].qi = PowerTable[i].qi;

		// Increase power of r1 for first symmetry
		PowerTable[i].ki = PowerTable[i].ki + l;
	}

	return;
}


// Same as GenOmegaPowerTable, except only returns terms with qi > 0.
void GenOmegaPowerTableQiGt0(int Omega, int l, int Ordering, vector <rPowers> &PowerTable, int Start, int End)
{
	int NumTerm = 0;  // The number of the current term
	int om, ki, li, mi, ni, pi, qi;  // These are the exponents we are determining.
	int Cur = 0;

	if (Start == -1)
		Start = 0;

	if (Ordering == 0) {  // My ordering
		for (om = 0; om <= Omega; om++) {
			for (ki = 0; ki <= Omega; ki++) {
				for (li = 0; li <= Omega; li++) {
					for (mi = 0; mi <= Omega; mi++) {
						for (ni = 0; ni <= Omega; ni++) {
							for (pi = 0; pi <= Omega; pi++) {
								for (qi = 0; qi <= Omega; qi++) {
									if (ki + li + mi + ni + pi + qi == om) {
										if (Cur >= Start) {
											if (qi > 0) {
												PowerTable[NumTerm].ki = ki;
												PowerTable[NumTerm].li = li;
												PowerTable[NumTerm].mi = mi;
												PowerTable[NumTerm].ni = ni;
												PowerTable[NumTerm].pi = pi;
												PowerTable[NumTerm].qi = qi;
												PowerTable[NumTerm].Index = NumTerm;
												NumTerm = NumTerm + 1;
											}
										}
										Cur = Cur + 1;
										if (Cur > End)
											goto finish;  //@TODO: Do this without a goto statement.
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else {  // Peter Van Reeth's ordering
		int IHDPP1 = Omega + 1;
		int INX = 0;
		for (int I = 1; I <= IHDPP1; I++) {
			for (int I23P1 = 1; I23P1 <= I; I23P1++) {
				int I23 = I23P1 - 1;
				int I12P1M = I - I23;
				for (int I12P1 = 1; I12P1 <= I12P1M; I12P1++) {
					int I12 = I12P1 - 1;
					int I2P1M = I12P1M - I12;
					for (int I2P1 = 1; I2P1 <= I2P1M; I2P1++) {
						int I2 = I2P1-1;
						int I13P1M = I2P1M-I2;
						for (int I13P1 = 1; I13P1 <= I13P1M; I13P1++) {
							int I13 = I13P1 -1;
							int I3P1M = I13P1M -I13;
							for (int I3P1 = 1; I3P1 <= I3P1M; I3P1++) {
								int I3=I3P1-1;
								int I1 = I3P1M -I3P1;

								if (Cur >= Start) {
									if (I23 > 0) {  // We only want the qi > 0 terms right now.
										PowerTable[INX].ki = I1;
										PowerTable[INX].li = I2;
										PowerTable[INX].mi = I12;
										PowerTable[INX].ni = I3;
										PowerTable[INX].pi = I13;
										PowerTable[INX].qi = I23;
										PowerTable[INX].Index = INX;
										INX = INX + 1;
									}
								}
								Cur = Cur + 1;
								if (Cur > End)
									goto finish;  //@TODO: Do this without a goto statement.
							}
						}
					}
				}
			}
		}
	}

finish:
	int size = PowerTable.size();
	NumTerm = size / 2;  //@TODO: Use Cur.
	for (int i = 0; i < NumTerm; i++) {
		PowerTable[NumTerm+i].ki = PowerTable[i].ki;
		// Increase power of r2 for second symmetry
		PowerTable[NumTerm+i].li = PowerTable[i].li + l;
		PowerTable[NumTerm+i].mi = PowerTable[i].mi;
		PowerTable[NumTerm+i].ni = PowerTable[i].ni;
		PowerTable[NumTerm+i].pi = PowerTable[i].pi;
		PowerTable[NumTerm+i].qi = PowerTable[i].qi;

		// Increase power of r1 for first symmetry
		PowerTable[i].ki = PowerTable[i].ki + l;
	}

	return;
}

