//
// Phase Shift.cpp
//

#include <iostream>
#include <iomanip>
#include <string>
#include <cerrno>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <mkl_lapack.h>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <boost/filesystem.hpp>

using namespace std;

typedef complex<double> dcmplx;

int		CalcPowerTableSize(int Omega);
int		ReadMatrixElem(ifstream &FileMatrixElem, int NumShortTerms, vector <double> &ARow, vector <double> &B, double &SLS, int &IsTriplet, int &Ordering, int &LValue, int &Formalism, int &Omega, int &NumSets, double &Alpha1, double &Beta1, double &Gamma1, double &Alpha2, double &Beta2, double &Gamma2, double &Kappa, double &Mu, string &LString, int &Shielding, string &Lambda, double &Epsilon12, double &Epsilon13, bool &ExtraExponential);
int		ReadShortHeader(ifstream &FileShortRange, int &Omega, int &LValue, int &IsTriplet, int &Formalism, int &Ordering, int &NumShortTerms, int &NumSets, int &Integration, double &Alpha1, double &Beta1, double &Gamma1, double &Alpha2, double &Beta2, double &Gamma2, bool &ExtraExponential, double &Epsilon12, double &Epsilon13, vector <int> &ExpLen);
void	WriteHeader(ofstream &OutFile, string &LString, int &LValue, char *FileShortName, char *FileMatrixElemName, char *EnergyFileName, bool &Paired, bool &Resorted,
				int &ShortInt, int &NumTerms, double &Kappa, double &Mu, int &Shielding, string &Lambda, double &Alpha, double &Beta, double &Gamma, string &ProgName);
int		CreateSubset(vector <double> &ARow, vector <double> &B, vector <double> &ShortTerms, vector <double> &ARowSub, vector <double> &BSub, vector <double> &ShortTermsSub, int NumShortTerms, int NSub);
int		FindOrderedToddTerm(string EnergyFilename, int TermToFind);
int		LoadToddTerms(int LValue, vector <double> &ARow, vector <double> &B, vector <double> &ShortTerms, vector <double> &ARowSub, vector <double> &BSub, vector <double> &ShortTermsSub, int NumShortTerms, int NSub, string EnergyFilename, bool Resorted, bool Paired, int ResortedSize);
int		TestToddFile(string EnergyFilename, int NumShortTerms);
void	uGenKohn(dcmplx (&u)[2][2], double Tau);
void	uGenTKohn(dcmplx (&u)[2][2], double Tau);
void	uGenSKohn(dcmplx (&u)[2][2], double Tau);
double	CombinedKohn(dcmplx (&u)[2][2], int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> &ShortTerms, double SLS, int LValue, int IsTriplet);
string	ShortIntString(int &Integration);
void	FixPhase(double &PhaseShift, int &LValue, int &IsTriplet);
string	GetDateTime(void);


#define NUM_TAUARRAY 35
double TauArray[NUM_TAUARRAY] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7853981633974483, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.570796326794897, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.356194490192345, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.141592653589793};
const char *DataHeader = "       n |         Kohn            |       Inverse Kohn      |     Complex Kohn (S)    |     Complex Kohn (T) "
	"   |    Gen Kohn tau = 0.0   |    Gen Kohn tau = 0.1   |    Gen Kohn tau = 0.2   |    Gen Kohn tau = 0.3   |    Gen Kohn tau = 0.4   |    Gen Kohn tau = 0.5"
	"   |    Gen Kohn tau = 0.6   |    Gen Kohn tau = 0.7   |   Gen Kohn tau = pi/4   |    Gen Kohn tau = 0.8   |    Gen Kohn tau = 0.9   |    Gen Kohn tau = 1.0"
	"   |    Gen Kohn tau = 1.1   |    Gen Kohn tau = 1.2   |    Gen Kohn tau = 1.3   |    Gen Kohn tau = 1.4   |    Gen Kohn tau = 1.5   |   Gen Kohn tau = pi/2"
	"   |    Gen Kohn tau = 1.6   |    Gen Kohn tau = 1.7   |    Gen Kohn tau = 1.8   |    Gen Kohn tau = 1.9   |    Gen Kohn tau = 2.0"
	"   |    Gen Kohn tau = 2.1   |    Gen Kohn tau = 2.2   |    Gen Kohn tau = 2.3   |  Gen Kohn tau = 3*pi/4  |    Gen Kohn tau = 2.4   |    Gen Kohn tau = 2.5"
	"   |    Gen Kohn tau = 2.6   |    Gen Kohn tau = 2.7   |    Gen Kohn tau = 2.8   |    Gen Kohn tau = 2.9   |    Gen Kohn tau = 3.0   |    Gen Kohn tau = pi "
	"   |   Gen T Kohn tau = 0.0  |   Gen T Kohn tau = 0.1  |   Gen T Kohn tau = 0.2  |   Gen T Kohn tau = 0.3  |   Gen T Kohn tau = 0.4  |   Gen T Kohn tau = 0.5"
	"  |   Gen T Kohn tau = 0.6  |   Gen T Kohn tau = 0.7  |  Gen T Kohn tau = pi/4  |   Gen T Kohn tau = 0.8  |   Gen T Kohn tau = 0.9  |   Gen T Kohn tau = 1.0"
	"  |   Gen T Kohn tau = 1.1  |   Gen T Kohn tau = 1.2  |   Gen T Kohn tau = 1.3  |   Gen T Kohn tau = 1.4  |   Gen T Kohn tau = 1.5  |  Gen T Kohn tau = pi/2"
	"  |   Gen T Kohn tau = 1.6  |   Gen T Kohn tau = 1.7  |   Gen T Kohn tau = 1.8  |   Gen T Kohn tau = 1.9  |   Gen T Kohn tau = 2.0"
	"  |   Gen T Kohn tau = 2.1  |   Gen T Kohn tau = 2.2  |   Gen T Kohn tau = 2.3  | Gen T Kohn tau = 3*pi/4 |   Gen T Kohn tau = 2.4  |   Gen T Kohn tau = 2.5"
	"  |   Gen T Kohn tau = 2.6  |   Gen T Kohn tau = 2.7  |   Gen T Kohn tau = 2.8  |   Gen T Kohn tau = 2.9  |   Gen T Kohn tau = 3.0  |   Gen T Kohn tau = pi "
	"  |   Gen S Kohn tau = 0.0  |   Gen S Kohn tau = 0.1  |   Gen S Kohn tau = 0.2  |   Gen S Kohn tau = 0.3  |   Gen S Kohn tau = 0.4  |   Gen S Kohn tau = 0.5"
	"  |   Gen S Kohn tau = 0.6  |   Gen S Kohn tau = 0.7  |  Gen S Kohn tau = pi/4  |   Gen S Kohn tau = 0.8  |   Gen S Kohn tau = 0.9  |   Gen S Kohn tau = 1.0"
	"  |   Gen S Kohn tau = 1.1  |   Gen S Kohn tau = 1.2  |   Gen S Kohn tau = 1.3  |   Gen S Kohn tau = 1.4  |   Gen S Kohn tau = 1.5  |  Gen S Kohn tau = pi/2"
	"  |   Gen S Kohn tau = 1.6  |   Gen S Kohn tau = 1.7  |   Gen S Kohn tau = 1.8  |   Gen S Kohn tau = 1.9  |   Gen S Kohn tau = 2.0"
	"  |   Gen S Kohn tau = 2.1  |   Gen S Kohn tau = 2.2  |   Gen S Kohn tau = 2.3  | Gen S Kohn tau = 3*pi/4 |   Gen S Kohn tau = 2.4  |   Gen S Kohn tau = 2.5"
	"  |   Gen S Kohn tau = 2.6  |   Gen S Kohn tau = 2.7  |   Gen S Kohn tau = 2.8  |   Gen S Kohn tau = 2.9  |   Gen S Kohn tau = 3.0  |   Gen S Kohn tau = pi";


dcmplx uKohn[2][2] = {{dcmplx(1,0),dcmplx(0,0)}, {dcmplx(0,0),dcmplx(1,0)}};
dcmplx uInvKohn[2][2] = {{dcmplx(0,0),dcmplx(1,0)}, {dcmplx(-1,0),dcmplx(0,0)}};
dcmplx uCompSKohn[2][2] = {{dcmplx(0,1), dcmplx(-1,0)}, {dcmplx(0,1), dcmplx(1,0)}};
dcmplx uCompTKohn[2][2] = {{dcmplx(1,0), dcmplx(0,0)}, {dcmplx(0,1), dcmplx(1,0)}};


std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
    const int strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const int strEnd = str.find_last_not_of(whitespace);
    const int strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


int main(int argc, char *argv[])
{
	ifstream FileMatrixElem, FileShortRange;
	ofstream OutFile;
	string LString, Lambda;
	double ShortAlpha1, ShortBeta1, ShortGamma1, ShortAlpha2, ShortBeta2, ShortGamma2, ShortEpsilon12, ShortEpsilon13, Kappa, Mu;
	double LongAlpha1, LongBeta1, LongGamma1, LongAlpha2, LongBeta2, LongGamma2, LongEpsilon12, LongEpsilon13;
	int ShortOmega, ShortLValue, ShortIsTriplet, ShortOrdering, /*ShortNumSets,*/ ShortFormalism;
	int LongOmega, LongLValue, LongIsTriplet, LongOrdering, LongNumSets, LongFormalism;
	int /*Ordering,*/ NumSets, NumShort, NumShortTotal, NumShortTermsFile, ShortInt, Shielding;
	vector <int> ExpLen;
	bool ExtraExponential;
	int TotalTerms;

	vector <double> ARow, B, ShortTerms;
	double *PhiPhi, *PhiHPhi;
	double SLS;
	vector <double> ARowSub, BSub, ShortTermsSub;
	double GenKohnPhase;
	bool Paired, Resorted = false;
	int TermStep;
	dcmplx u[2][2];

	string ProgName = boost::filesystem::canonical(argv[0]).string();  // Get the absolute path of this program

	// Initialize the second set of nonlinear parameters for the files that don't use them.
	ShortAlpha2 = 0.0; LongAlpha2 = 0.0; ShortBeta2 = 0.0; LongBeta2 = 0.0; ShortGamma2 = 0.0; LongGamma2 = 0.0;

	char *FileMatrixElemName = argv[2];
	char *FileShortName = argv[3];
	char *OutFileName = argv[4];
	char *EnergyFileName = argv[6];

	if (argc < 7) {
		cerr << "Not enough parameters on the command line." << endl;
		cerr << "Usage: Phase pairing matrixelements.txt shortrangefile.bin results.txt #terms (energyfile.txt) (resorted?)" << endl;
		cerr << "Example: Phase 1 matrixelements.txt shortrangefile.bin results.txt 84 energyfile.txt true" << endl << endl;
		cerr << " The pairing parameter is 0 for no pairing of terms for the two symmetries and" << endl;
		cerr << " 1 for pairing." << endl;
		return 1;
	}

	if (atoi(argv[1]) == 0) {
		Paired = false;
		TermStep = 1;
	}
	else if (atoi(argv[1]) == 1) {
		Paired = true;
		TermStep = 2;
	}
	else {
		cout << "The pairing entry must be either 0 or 1." << endl;
		return 2;
	}

	FileMatrixElem.open(FileMatrixElemName);
	if (FileMatrixElem.fail()) {
		cerr << "Unable to open file " << FileMatrixElemName << " for reading." << endl;
		return 2;
	}

	FileShortRange.open(FileShortName, ios::in | ios::binary);
	if (FileShortRange.fail()) {
		cerr << "Unable to open file " << FileShortName << " for reading." << endl;
		return 3;
	}

	OutFile.open(OutFileName);
	if (!OutFile.is_open()) {
		cout << "Could not open output file...exiting." << endl;
		return 4;
	}
	TotalTerms = atoi(argv[5]);

	if (argc > 6) {
		int ToddTermNum = TestToddFile(EnergyFileName, TotalTerms);
		if (ToddTermNum < TotalTerms) {
			cout << "Using less than the requested number of terms: " << ToddTermNum << " instead of " << TotalTerms << endl;
			TotalTerms = ToddTermNum;
		}
	}

	if (argc > 7) {
		//@TODO: Case-insensitive string compare
		if (string(argv[7]) == "true") {
			Resorted = true;
			cout << "Computations will be performed with the terms resorted." << endl;
		}
		// Any other string just sets Resorted to false.
		else {
			cout << "Computations will be performed with the ordering specified in the energy file." << endl;
		}
	}

	// Include trailing zeros so the columns line up in the output file.
	cout.setf(ios::showpoint);
	OutFile.setf(ios::showpoint);
	cout << setprecision(18);
	OutFile << setprecision(18);

	int err = ReadShortHeader(FileShortRange, ShortOmega, ShortLValue, ShortIsTriplet, ShortFormalism, ShortOrdering, NumShortTermsFile, NumSets, ShortInt, ShortAlpha1, ShortBeta1, ShortGamma1, ShortAlpha2, ShortBeta2, ShortGamma2, ExtraExponential, ShortEpsilon12, ShortEpsilon13, ExpLen);
	if (err == -1) {
		return 8;
	}

	// Calculate number of terms for a given omega and generate the r-powers.
	NumShort = CalcPowerTableSize(ShortOmega);
	if (NumShort != NumShortTermsFile) {
		//cout << "Number of terms does not match in files...exiting." << endl;
		//return 2;
	}

	// The P-wave files have double the number of elements.
	//NumShortTerms = NumShortTerms*(ShortLValue+1);

	if (ShortLValue == 0)
		NumShortTotal = NumShort;  // The S-wave only has one symmetry
	else
		NumShortTotal = NumShort * 2;

	// Allocate PhiPhi and PhiHPhi matrices
	PhiPhi = new double[NumShortTotal*NumShortTotal];
	PhiHPhi = new double[NumShortTotal*NumShortTotal];
	if (PhiPhi == NULL || PhiHPhi == NULL) {
		cout << "Memory allocation error" << endl;
		return 5;
	}

	// Read in the <phi|phi> and <phi|H|phi> matrix elements.
	FileShortRange.read((char*)PhiPhi, NumShortTotal*NumShortTotal*sizeof(double));
	FileShortRange.read((char*)PhiHPhi, NumShortTotal*NumShortTotal*sizeof(double));

	ARow.resize(NumShortTotal+1);
	B.resize(NumShortTotal+1);
	ShortTerms.resize(NumShortTotal*NumShortTotal);
	err = ReadMatrixElem(FileMatrixElem, NumShortTotal, ARow, B, SLS, LongIsTriplet, LongOrdering, LongLValue, LongFormalism, LongOmega, LongNumSets, LongAlpha1, LongBeta1, LongGamma1, LongAlpha2, LongBeta2, LongGamma2, Kappa, Mu, LString, Shielding, Lambda, LongEpsilon12, LongEpsilon13, ExtraExponential);
	if (err == -1)
		return 6;

	// Compare the short-range and long-range files to make sure they are describing the same problem.
	if ((LongLValue != ShortLValue && ShortLValue != 0) || LongIsTriplet != ShortIsTriplet || LongOrdering != ShortOrdering || LongOmega != ShortOmega || LongAlpha1 != ShortAlpha1 || LongBeta1 != ShortBeta1 || LongGamma1 != ShortGamma1 || LongAlpha2 != ShortAlpha2 || LongBeta2 != ShortBeta2 || LongGamma2 != ShortGamma2) {
		cout << "Short-range and long-range files describe different problems...exiting." << endl;
		return 8;
	}

	for (int i = 0; i < NumShortTotal*NumShortTotal; i++) {
		ShortTerms[i] = PhiHPhi[i] - 0.5*Kappa*Kappa * PhiPhi[i] + 1.5*PhiPhi[i];
		//ShortTerms[i] = PhiHPhi[i] - Kappa*Kappa * PhiPhi[i] + 1.5*PhiPhi[i];  // For electron or positron scattering
	}

	delete [] PhiPhi;
	delete [] PhiHPhi;

	WriteHeader(OutFile, LString, ShortLValue, FileShortName, FileMatrixElemName, EnergyFileName, Paired, Resorted,
				ShortInt, TotalTerms, Kappa, Mu, Shielding, Lambda, ShortAlpha1, ShortBeta1, ShortGamma1, ProgName);

	int FieldWidth = 25;

	for (int i = 0; i <= TotalTerms; i++) {
		LoadToddTerms(ShortLValue, ARow, B, ShortTerms, ARowSub, BSub, ShortTermsSub, NumShortTotal, i, EnergyFileName, Resorted, Paired, TotalTerms);
		double KohnPhase = CombinedKohn(uKohn, i*TermStep, ARowSub, BSub, ShortTermsSub, SLS, ShortLValue, ShortIsTriplet);
		double InvKohnPhase = CombinedKohn(uInvKohn, i*TermStep, ARowSub, BSub, ShortTermsSub, SLS, ShortLValue, ShortIsTriplet);
		double CompKohnSPhase = CombinedKohn(uCompSKohn, i*TermStep, ARowSub, BSub, ShortTermsSub, SLS, ShortLValue, ShortIsTriplet);
		double CompKohnTPhase = CombinedKohn(uCompTKohn, i*TermStep, ARowSub, BSub, ShortTermsSub, SLS, ShortLValue, ShortIsTriplet);
		if ((KohnPhase == 0.0) || (InvKohnPhase == 0.0) || (CompKohnSPhase == 0.0) || (CompKohnTPhase == 0.0)) {
			cout << "Terminating loop early due to LAPACK errors." << endl;
			OutFile << "Terminating loop early due to LAPACK errors." << endl;
			break;
		}
		cout << i << " " << KohnPhase << " " << InvKohnPhase << " " << CompKohnSPhase << " " << CompKohnTPhase << endl;
		OutFile << setw(8) << i << setw(1) << " " << setw(FieldWidth) << KohnPhase << setw(1) << " " << setw(FieldWidth) << InvKohnPhase << setw(1) << " " << setw(FieldWidth)
				<< CompKohnSPhase << setw(1) << " " << setw(FieldWidth) << CompKohnTPhase;

		// Generalized Kohn
		for (int t = 0; t < NUM_TAUARRAY; t++) {
			LoadToddTerms(ShortLValue, ARow, B, ShortTerms, ARowSub, BSub, ShortTermsSub, NumShortTotal, i, EnergyFileName, Resorted, Paired, TotalTerms);
			uGenKohn(u, TauArray[t]);
			GenKohnPhase = CombinedKohn(u, i*TermStep, ARowSub, BSub, ShortTermsSub, SLS, ShortLValue, ShortIsTriplet);
			if (GenKohnPhase == 0.0)
				break;
			OutFile << setw(1) << " " << setw(FieldWidth) << GenKohnPhase;
		}

		// Generalized T-matrix
		for (int t = 0; t < NUM_TAUARRAY; t++) {
			LoadToddTerms(ShortLValue, ARow, B, ShortTerms, ARowSub, BSub, ShortTermsSub, NumShortTotal, i, EnergyFileName, Resorted, Paired, TotalTerms);
			uGenTKohn(u, TauArray[t]);
			GenKohnPhase = CombinedKohn(u, i*TermStep, ARowSub, BSub, ShortTermsSub, SLS, ShortLValue, ShortIsTriplet);
			if (GenKohnPhase == 0.0)
				break;
			OutFile << setw(1) << " " << setw(FieldWidth) << GenKohnPhase;
		}

		// Generalized S-matrix
		for (int t = 0; t < NUM_TAUARRAY; t++) {
			LoadToddTerms(ShortLValue, ARow, B, ShortTerms, ARowSub, BSub, ShortTermsSub, NumShortTotal, i, EnergyFileName, Resorted, Paired, TotalTerms);
			uGenSKohn(u, TauArray[t]);
			GenKohnPhase = CombinedKohn(u, i*TermStep, ARowSub, BSub, ShortTermsSub, SLS, ShortLValue, ShortIsTriplet);
			if (GenKohnPhase == 0.0)
				break;
			OutFile << setw(1) << " " << setw(FieldWidth) << GenKohnPhase;
		}
		if (GenKohnPhase == 0.0) {
			cout << "Terminating loop early due to LAPACK errors." << endl;
			OutFile << "Terminating loop early due to LAPACK errors." << endl;
			break;
		}

		OutFile << setw(1) << " " << endl;
	}

	OutFile << "</data>" << endl << "</psh_data>" << endl;

	FileMatrixElem.close();
	FileShortRange.close();
	OutFile.close();

	return 0;
}


int CreateSubset(vector <double> &ARow, vector <double> &B, vector <double> &ShortTerms, vector <double> &ARowSub, vector <double> &BSub, vector <double> &ShortTermsSub, int NumShortTerms, int NSub)
{
	ARowSub.resize(NSub*2+1);
	BSub.resize(NSub*2+1);
	ShortTermsSub.resize(NSub*NSub*4);

	vector <int> UsedTerms, UsedTermsSub;

	UsedTerms.resize(NumShortTerms*2);
	for (int i = 0; i < NumShortTerms*2; i++) {
		UsedTerms[i] = i+1;
	}
	UsedTermsSub.resize(NSub*2);
	for (int i = 0; i < NSub; i++) {
		UsedTermsSub[i*2] = i+1;
		UsedTermsSub[i*2+1] = NumShortTerms+i+1;
	}

	ARowSub[0] = ARow[0];
	BSub[0] = B[0];
	for (int i = 0; i < NSub*2; i++) {
		for (int j = 0; j < NSub*2; j++) {
			// The Fortran output counts from 1, hence the -1 on the RHS.
			ShortTermsSub[i*NSub*2 + j] = ShortTerms[(UsedTermsSub[i]-1)*NumShortTerms*2 + (UsedTermsSub[j]-1)];
		}
		// Skips the 0 entry in ARow, so they line up.
		ARowSub[i+1] = ARow[UsedTermsSub[i]];
	}
	for (int i = 0; i < NSub*2; i++) {
		BSub[i+1] = B[UsedTermsSub[i]];
	}

	return 0;
}


// Sort algorithm example code from http://www.cplusplus.com/reference/algorithm/sort/
struct myclass {
	bool operator() (int i,int j) { return (i<j);}
} myobject;

int FindOrderedToddTerm(string EnergyFilename, int TermToFind, int NumTerms)
{
	ifstream EnergyFile;
	vector <int> UsedTerms, UsedTermsSub;
	string Line;
	int Term, Index;
	double Energy;

	if (TermToFind < 1) {
		cout << "TermToFind must be 1 or greater." << endl;
		return 1;
	}

	EnergyFile.open(EnergyFilename.c_str());
	getline(EnergyFile, Line);
	getline(EnergyFile, Line);  // Skip the first 4 lines
	getline(EnergyFile, Line);  //  (unimportant for this)
	getline(EnergyFile, Line);

	UsedTerms.resize(NumTerms);
	for (int i = 0; i < NumTerms; i++) {
		EnergyFile >> Term >> Index >> Energy;
		UsedTerms[i] = Term;
	}

	sort(UsedTerms.begin(), UsedTerms.end(), myobject);

	EnergyFile.close();

	// Now search for the term (or find where it would go).
	for (int i = 0; i < NumTerms; i++) {
		if (UsedTerms[i] == TermToFind)
			return i+1;
		if (UsedTerms[i] > TermToFind) {
			if (i == 0)
				return 1;  // Don't want to return 0.
			return i;
		}
	}

	return NumTerms+1;  // Term not found (larger than last used term).
}


// We could just open the energy file in the main program instead of reopening it many times, but this is just easier.
int TestToddFile(string EnergyFilename, int NumShortTerms)
{
	ifstream EnergyFile;
	string Line;
	int Term, Index;
	double Energy;

	if (NumShortTerms < 1) {
		cout << "NumShortTerms must be 1 or greater." << endl;
		return 1;
	}

	EnergyFile.open(EnergyFilename.c_str());
	if (EnergyFile.fail())
		return 0;
	getline(EnergyFile, Line);
	getline(EnergyFile, Line);  // Skip the first 4 lines
	getline(EnergyFile, Line);  //  (unimportant for this)
	getline(EnergyFile, Line);

	for (int i = 0; i < NumShortTerms; i++) {
		EnergyFile >> Term >> Index >> Energy;
		if (EnergyFile.fail())  // End of terms to use
			return i;
	}

	EnergyFile.close();

	return NumShortTerms;
}


// We could just open the energy file in the main program instead of reopening it many times, but this is just easier.
int LoadToddTerms(int ShortLValue, vector <double> &ARow, vector <double> &B, vector <double> &ShortTerms, vector <double> &ARowSub, vector <double> &BSub, 
					vector <double> &ShortTermsSub, int NumShortTotal, int NSub, string EnergyFilename, bool Resorted, bool Paired, int ResortedSize)
{
	ifstream EnergyFile;
	vector <int> UsedTerms, UsedTermsSub;
	string Line;
	int Term, Index, Size;
	double Energy;

	if (NSub < 0) {
		cout << "NSub must be 0 or greater." << endl;
		return 1;
	}

	EnergyFile.open(EnergyFilename.c_str());
	if (EnergyFile.fail())
		return 0;
	getline(EnergyFile, Line);
	getline(EnergyFile, Line);  // Skip the first 4 lines
	getline(EnergyFile, Line);  //  (unimportant for this)
	getline(EnergyFile, Line);

	if (Resorted == false)
		Size = NSub;
	else
		Size = ResortedSize;

	UsedTerms.resize(Size);
	for (int i = 0; i < Size; i++) {
		EnergyFile >> Term >> Index >> Energy;
		if (EnergyFile.fail())  // End of terms to use
			return 0;
		UsedTerms[i] = Term;
	}

	if (Paired == true) {
		UsedTermsSub.resize(NSub*2);
	
		if (Resorted) {
			cout << "Reordering terms" << endl;
			sort(UsedTerms.begin(), UsedTerms.end(), myobject);
		}

		for (int i = 0; i < NSub; i++) {
			UsedTermsSub[i*2] = UsedTerms[i];
			UsedTermsSub[i*2+1] = NumShortTotal/2+UsedTerms[i];
		}

		EnergyFile.close();

		ARowSub.resize(NSub*2+1);
		BSub.resize(NSub*2+1);
		ShortTermsSub.resize(NSub*NSub*4);

		ARowSub[0] = ARow[0];
		BSub[0] = B[0];
		for (int i = 0; i < NSub*2; i++) {
			for (int j = 0; j < NSub*2; j++) {
				// The Fortran output counts from 1, hence the -1 on the RHS.
				ShortTermsSub[i*NSub*2 + j] = ShortTerms[(UsedTermsSub[i]-1)*NumShortTotal + (UsedTermsSub[j]-1)];
			}
			// Skips the 0 entry in ARow, so they line up.
			ARowSub[i+1] = ARow[UsedTermsSub[i]];
		}
		for (int i = 0; i < NSub*2; i++) {
			BSub[i+1] = B[UsedTermsSub[i]];
		}
	}
	else {
		UsedTermsSub.resize(NSub);
	
		if (Resorted) {
			cout << "Reordering terms" << endl;
			sort(UsedTerms.begin(), UsedTerms.end(), myobject);
		}

		for (int i = 0; i < NSub; i++) {
			UsedTermsSub[i] = UsedTerms[i];
		}

		EnergyFile.close();

		ARowSub.resize(NSub+1);
		BSub.resize(NSub+1);
		ShortTermsSub.resize(NSub*NSub);

		ARowSub[0] = ARow[0];
		BSub[0] = B[0];
		for (int i = 0; i < NSub; i++) {
			for (int j = 0; j < NSub; j++) {
				// The Fortran output counts from 1, hence the -1 on the RHS.
				ShortTermsSub[i*NSub + j] = ShortTerms[(UsedTermsSub[i]-1)*NumShortTotal + (UsedTermsSub[j]-1)];
			}
			// Skips the 0 entry in ARow, so they line up.
			ARowSub[i+1] = ARow[UsedTermsSub[i]];
		}
		for (int i = 0; i < NSub; i++) {
			BSub[i+1] = B[UsedTermsSub[i]];
		}
	}

	return NSub;
}


// Returns the number of terms for a given omega.  This could use the formula for combination with repetition,
//  except then it would be unable to use a restricted set of terms if we needed.
int CalcPowerTableSize(int Omega)
{
	int NumTerms = 0;  // The total number of terms
	int om, ki, li, mi, ni, pi, qi;  // These are the exponents we are determining.

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


// Reads in the output from the scattering program (short-range - long-range and long-range - long-range terms)
int ReadMatrixElem(ifstream &FileMatrixElem, int NumShortTerms, vector <double> &ARow, vector <double> &B, double &SLS, int &IsTriplet, int &Ordering, int &LValue, int &Formalism, int &Omega, int &NumSets, double &Alpha1, double &Beta1, double &Gamma1, double &Alpha2, double &Beta2, double &Gamma2, double &Kappa, double &Mu, string &LString, int &Shielding, string &Lambda, double &Epsilon12, double &Epsilon13, bool &ExtraExponential)
{
	string Line, Line1, Line2, Line3, Line4, OrderString;
	int NumTerms, Offset;

	getline(FileMatrixElem, Line);
	getline(FileMatrixElem, Line);
	getline(FileMatrixElem, LString);
	getline(FileMatrixElem, OrderString);

	FileMatrixElem >> Line >> Omega;
	FileMatrixElem >> Line1 >> Line2 >> Line3 >> NumTerms;

	cout << LString << endl;
	if (LString == "S-Wave Singlet Ps-H") {
		LValue = 0;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (LString == "S-Wave Triplet Ps-H") {
		LValue = 0;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (LString == "S-Wave Singlet Ps-H - Laplacian Formalism - Exponential") {
		LValue = 0;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (LString == "S-Wave Triplet Ps-H - Laplacian Formalism - Exponential") {
		LValue = 0;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (LString == "P-Wave Singlet Ps-H: 1st formalism" || LString == "P-Wave Singlet Ps-H") {  // Second is the older type
		LValue = 1;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (LString == "P-Wave Triplet Ps-H: 1st formalism" || LString == "P-Wave Triplet Ps-H") {  // Second is the older type
		LValue = 1;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (LString == "P-Wave Singlet Ps-H: 2nd formalism") {
		LValue = 1;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (LString == "P-Wave Triplet Ps-H: 2nd formalism") {
		LValue = 1;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (LString == "P-Wave Singlet Ps-H: 1st formalism / 2 sets") {
		LValue = 1;
		IsTriplet = 0;
		NumSets = 2;
	}
	else if (LString == "P-Wave Triplet Ps-H: 1st formalism / 2 sets") {
		LValue = 1;
		IsTriplet = 1;
		NumSets = 2;
	}
	else if (LString == "P-Wave Singlet Ps-H: 2nd formalism / 2 sets") {
		LValue = 1;
		IsTriplet = 0;
		NumSets = 2;
	}
	else if (LString == "P-Wave Triplet Ps-H: 2nd formalism / 2 sets") {
		LValue = 1;
		IsTriplet = 1;
		NumSets = 2;
	}
	else if (LString == "D-Wave Singlet Ps-H: 1st formalism") {
		LValue = 2;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (LString == "D-Wave Triplet Ps-H: 1st formalism") {
		LValue = 2;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (LString == "F-Wave Singlet Ps-H") {
		LValue = 3;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (LString == "F-Wave Triplet Ps-H") {
		LValue = 3;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (LString == "G-Wave Singlet Ps-H") {
		LValue = 4;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (LString == "G-Wave Triplet Ps-H") {
		LValue = 4;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (LString == "H-Wave Singlet Ps-H") {
		LValue = 5;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (LString == "H-Wave Triplet Ps-H") {
		LValue = 5;
		IsTriplet = 1;
		NumSets = 1;
	}
	else {
		cout << "Problem string in matrix element file has an unknown value...exiting." << endl;
		return -1;
	}

	if (ExtraExponential) {
		FileMatrixElem >> Line >> Alpha1 >> Line1 >> Beta1 >> Line2 >> Gamma1 >> Line3 >> Epsilon12 >> Line4 >> Epsilon13;
		cout << Alpha1 << " " << Beta1 << " " << Gamma1 << endl;
	}
	else {
		if (NumSets == 1) {
			FileMatrixElem >> Line >> Alpha1 >> Line1 >> Beta1 >> Line2 >> Gamma1;
			cout << Alpha1 << " " << Beta1 << " " << Gamma1 << endl;
		}
		else if (NumSets == 2) {
			FileMatrixElem >> Line >> Alpha2 >> Line1 >> Beta2 >> Line2 >> Gamma2;
		}
	}

	//getline(FileMatrixElem, Line);
	//getline(FileMatrixElem, Line);
	FileMatrixElem >> Line >> Mu;
	FileMatrixElem >> Line;
	Shielding = -1;
	if (Line.find("Shielding") != string::npos) {  // Skip this line - not yet to kappa line
		FileMatrixElem >> Line;
		FileMatrixElem >> Shielding;
		//getline(FileMatrixElem, Line);
		FileMatrixElem >> Line;
	}
	FileMatrixElem >> Kappa;
	getline(FileMatrixElem, Line);

	if (OrderString == "Using Denton's ordering") {
		Ordering = 0;
	}
	else if (OrderString == "Using Peter Van Reeth's ordering") {
		Ordering = 1;
	}
	else {
		cout << "Ordering string in matrix element file has an unknown value...exiting." << endl;
		return -1;
	}

	//getline(FileMatrixElem, Line);
	FileMatrixElem >> Line;
	if (Line.find("Lambda") != string::npos) {  // Has the extra lambda line here
		getline(FileMatrixElem, Lambda);
	}
	getline(FileMatrixElem, Line);
	
	for (int i = 0; i < 11; i++) {
		getline(FileMatrixElem, Line);
		if (Line == "A matrix row")  // Some files have one less extra line
			break;
	}

	// Reads in first row (and column) of A
	for (int i = 0; i < NumShortTerms+1; i++) {
		getline(FileMatrixElem, Line);
		istringstream iss(Line);
		iss >> Offset >> ARow[i];
	}

	// Skips extra lines
	getline(FileMatrixElem, Line);
	getline(FileMatrixElem, Line);

	// Reads in B vector
	for (int i = 0; i < NumShortTerms+1; i++) {
		getline(FileMatrixElem, Line);
		istringstream iss(Line);
		iss >> Offset >> B[i];
	}

	getline(FileMatrixElem, Line);
	getline(FileMatrixElem, Line);

	// SLS is not in A or B, so read it in.
	FileMatrixElem >> SLS;

	return 0;
}


// Reads in the short-range file header
int ReadShortHeader(ifstream &FileShortRange, int &Omega, int &LValue, int &IsTriplet, int &Formalism, int &Ordering, int &NumShortTerms, int &NumSets, int &Integration, double &Alpha1, double &Beta1, double &Gamma1, double &Alpha2, double &Beta2, double &Gamma2, bool &ExtraExponential, double &Epsilon12, double &Epsilon13, vector <int> &ExpLen)
{
	int MagicNum, Version, HeaderLen, DataFormat, NumShortTerms1, NumShortTerms2;
	int VarLen;

	FileShortRange.read((char*)&MagicNum, 4);
	FileShortRange.read((char*)&Version, 4);
	FileShortRange.read((char*)&HeaderLen, 4);
	FileShortRange.read((char*)&DataFormat, 4);
	FileShortRange.read((char*)&Omega, 4);
	FileShortRange.read((char*)&NumShortTerms1, 4);
	FileShortRange.read((char*)&NumShortTerms2, 4);
	FileShortRange.read((char*)&LValue, 4);
	FileShortRange.read((char*)&Formalism, 4);
	FileShortRange.read((char*)&IsTriplet, 4);
	FileShortRange.read((char*)&Ordering, 4);
	FileShortRange.read((char*)&Integration, 4);
	FileShortRange.read((char*)&NumSets, 4);

	//@TODO: More descriptive errors for each

	if (MagicNum != 0x31487350) {  // "PsH1" in hexadecimal (with reverse due to endianness)
		cout << "This is not a valid Ps-H file (MagicNum)...exiting." << endl;
		return -1;
	}

	if (Version < 1 || Version > 9) {
		cout << "This is not a valid Ps-H file (Version)...exiting." << endl;
		cout << Version << endl;
		return -1;
	}

	if (HeaderLen != 80 && HeaderLen != 104) {
		cout << "This is not a valid Ps-H file (HeaderLen)...exiting." << endl;
		return -1;
	}

	if (DataFormat != 8) {
		cout << "This is not a valid Ps-H file (Dataformat)...exiting." << endl;
		return -1;
	}

	if ((Formalism != 1 && Formalism != 2) || (IsTriplet != 0 && IsTriplet != 1) || (Ordering != 0 && Ordering != 1)) {
		cout << "This is not a valid Ps-H file (Formalism/IsTriplet/Ordering)...exiting." << endl;
		return -1;
	}

	if (NumSets != 1 && NumSets != 2) {
		cout << "This is not a valid Ps-H file (NumSets)...exiting." << endl;
		return -1;
	}

	/*if (NumShortTerms1 != NumShortTerms2) {
		cout << "This is not a valid Ps-H file...exiting." << endl;  // Cannot handle two different values for this yet.
		return -1;
	}*/
	NumShortTerms = NumShortTerms1;

	// @TODO: Set up to work properly with sectors
	FileShortRange.read((char*)&Alpha1, 8);
	FileShortRange.read((char*)&Beta1, 8);
	FileShortRange.read((char*)&Gamma1, 8);
	ExtraExponential = false;
	if (Version == 9) {  // Extra exponentials
		double BlankDouble;
		int BlankInt;
		ExtraExponential = true;
		FileShortRange.read((char*)&Epsilon12, 8);
		FileShortRange.read((char*)&Epsilon13, 8);
		FileShortRange.read((char*)&BlankDouble, 8);  // To be reserved for Epsilon23 at some point in the future

		ExpLen.resize(2);  // No r23 exponential right now
		FileShortRange.read((char*)&ExpLen[0], 4);
		FileShortRange.read((char*)&ExpLen[1], 4);
		FileShortRange.read((char*)&BlankInt, 4);
	}
	if (NumSets == 2) {
		//read (FileShortRange) Alpha2, Beta2, Gamma2
		FileShortRange.read((char*)&Alpha2, 8);
		FileShortRange.read((char*)&Beta2, 8);
		FileShortRange.read((char*)&Gamma2, 8);
	}

	FileShortRange.read((char*)&VarLen, 4);

	return 0;
}


void WriteHeader(ofstream &OutFile, string &LString, int &LValue, char *FileShortName, char *FileMatrixElemName, char *EnergyFileName, bool &Paired, bool &Resorted,
				int &ShortInt, int &NumTerms, double &Kappa, double &Mu, int &Shielding, string &Lambda, double &Alpha, double &Beta, double &Gamma, string &ProgName)
{
	OutFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?> " << endl;
	OutFile << "<psh_data>" << endl << "<header>" << endl;
	OutFile << "	<problem>" << LString << "</problem>" << endl;
	OutFile << "	<lvalue>" << LValue << "</lvalue>" << endl;
	OutFile << "	<shortfile>" << FileShortName << "</shortfile>" << endl;
	OutFile << "	<longfile>" << FileMatrixElemName << "</longfile>" << endl;
	//if (argc >= 6)  // Energy file present  //@TODO: Are we even accepting runs without now?
		OutFile << "	<energyfile>" << EnergyFileName << "</energyfile>" << endl;
	if (Paired == true)
		OutFile << "	<paired>" << "true" << "</paired>" << endl;
	else
		OutFile << "	<paired>" << "false" << "</paired>" << endl;
	OutFile << "	<ordering>" << "Peter" << "</ordering>" << endl;
	if (Resorted == true)
		OutFile << "	<reorder>" << "true" << "</reorder>" << endl;
	else
		OutFile << "	<reorder>" << "false" << "</reorder>" << endl;
	OutFile << "	<shortint>" << ShortIntString(ShortInt) << "</shortint>" << endl;
	OutFile << "	<numterms>" << NumTerms << "</numterms>" << endl;
	OutFile << "	<numsets>" << 1 << "</numsets>" << endl;
	OutFile << "	<kappa>" << Kappa << "</kappa>" << endl;
	OutFile << "	<mu>" << Mu << "</mu>" << endl;
	if (Shielding == -1)  // No shielding value specified in file - assume default
		Shielding = 2*LValue + 1;
	OutFile << "	<shielding>" << Shielding << "</shielding>" << endl;
	if (trim(Lambda) != "" || Lambda.size() > 0)
		OutFile << "	<lambda>" << trim(Lambda) << "</lambda>" << endl;
	OutFile << "	<nonlinear>" << endl;
	OutFile << "		<alpha>" << Alpha << "</alpha>" << endl;
	OutFile << "		<beta>" << Beta << "</beta>" << endl;
	OutFile << "		<gamma>" << Gamma << "</gamma>" << endl;
	OutFile << "	</nonlinear>" << endl;
	OutFile << "	<program>" << ProgName << "</program>" << endl;
	OutFile << "	<datetime>" << GetDateTime() << "</datetime>" << endl;

	cout << "n          Kohn              Inverse Kohn          Complex Kohn (S)       Complex Kohn (T)" << endl;
	OutFile << "</header>" << endl << "<dataheader>" << endl << DataHeader << endl << "</dataheader>" << endl << "<data>" << endl;
	OutFile << setprecision(16);
	OutFile << scientific;

	return;
}


// Generalized real Kohn
void uGenKohn(dcmplx (&u)[2][2], double Tau)
{
	u[0][0] = dcmplx(cos(Tau),0);
	u[0][1] = dcmplx(sin(Tau),0);
	u[1][0] = dcmplx(-sin(Tau),0);
	u[1][1] = dcmplx(cos(Tau),0);
	return;
}


// Generalized T-matrix Kohn
void uGenTKohn(dcmplx (&u)[2][2], double Tau)
{
	u[0][0] = dcmplx(cos(Tau),0);
	u[0][1] = dcmplx(sin(Tau),0);
	u[1][0] = dcmplx(-sin(Tau),cos(Tau));
	u[1][1] = dcmplx(cos(Tau),sin(Tau));
	return;
}


// Generalized S-matrix Kohn
void uGenSKohn(dcmplx (&u)[2][2], double Tau)
{
	u[0][0] = dcmplx(-sin(Tau),-cos(Tau));
	u[0][1] = dcmplx(cos(Tau),-sin(Tau));
	u[1][0] = dcmplx(-sin(Tau),cos(Tau));
	u[1][1] = dcmplx(cos(Tau),sin(Tau));
	return;
}


// ARow and BVec are sent from the main program as the A and B from the Kohn method.  This function rearranges everything into
//  the matrix equation (7) of their paper and solves.
double CombinedKohn(dcmplx (&u)[2][2], int NumShortTerms, vector <double> &ARow, vector <double> &B, vector <double> &ShortTerms, double SLS, int LValue, int IsTriplet)
{
	MKL_INT n, nrhs, lda, ldb, info;
	//double *A = new double[(NumShortTerms+1)*(NumShortTerms+1)];
	vector <dcmplx> A((NumShortTerms+1)*(NumShortTerms+1));
	vector <dcmplx> X(NumShortTerms+1);
	double CLC = ARow[0], CLS = B[0];
	double SLC = CLS + 1.0;  // Use (S,LC) = (C,LS) + 1
	dcmplx SLSt, SLCt, CLSt, CLCt;

	dcmplx detu = u[0][0]*u[1][1] - u[0][1]*u[1][0];  // Determinant

	SLSt = u[0][0]*u[0][0]*SLS + u[0][0]*u[0][1]*SLC + u[0][1]*u[0][0]*CLS + u[0][1]*u[0][1]*CLC;
	SLCt = u[0][0]*u[1][0]*SLS + u[0][0]*u[1][1]*SLC + u[0][1]*u[1][0]*CLS + u[0][1]*u[1][1]*CLC;
	CLSt = u[1][0]*u[0][0]*SLS + u[1][0]*u[0][1]*SLC + u[1][1]*u[0][0]*CLS + u[1][1]*u[0][1]*CLC;
	CLCt = u[1][0]*u[1][0]*SLS + u[1][0]*u[1][1]*SLC + u[1][1]*u[1][0]*CLS + u[1][1]*u[1][1]*CLC;

	// Copy short-range terms to bottom-right NumShortTerms x NumShortTerms submatrix of A.
	for (int i = 0; i < NumShortTerms; i++) {
		for (int j = 0; j < NumShortTerms; j++) {
			A[(i+1)*(NumShortTerms+1) + (j+1)] = ShortTerms[i*NumShortTerms + j];
		}
	}

	// Fill in the rest of A
	A[0] = CLCt;
	for (int i = 1; i < NumShortTerms+1; i++) {
		//A[i] = dcmplx(ARow[i], B[i]);
		A[i] = u[1][0]*B[i] + u[1][1]*ARow[i];
		A[i*(NumShortTerms+1)] = A[i];
	}

	// Fill in B (or X)
	X[0] = -CLSt;
	for (int i = 1; i < NumShortTerms+1; i++) {
		X[i] = - u[0][0]*B[i] - u[0][1]*ARow[i];
	}

	// LAPACK requires calls by reference, so we have to define all these variables.
	vector <int> ipiv(NumShortTerms+1);
	n = lda = ldb = NumShortTerms+1;
	nrhs = 1;
	zgesv(&n, &nrhs, (MKL_Complex16*)&A[0], &lda, &ipiv[0], (MKL_Complex16*)&X[0], &ldb, &info);
	if (info != 0) {
		cout << "LAPACK Error: " << info << endl;
		return 0.0;
	}

	// Equation () of notes
	dcmplx PsiLS = X[0] * CLSt;
	for (int i = 1; i < (NumShortTerms+1); i++) {
		PsiLS += X[i] * (u[0][0]*B[i] + u[0][1]*ARow[i]);
	}
	dcmplx L = -(PsiLS + SLSt) / detu;
	// Go from general L matrix element to K.
	dcmplx K = (u[0][1] + u[1][1]*L) / (u[0][0] + u[1][0]*L);

	//@TODO: Check K for imaginary part.
	double PhaseShift = atan(K.real());
	FixPhase(PhaseShift, LValue, IsTriplet);
	return PhaseShift;
}


string ShortIntString(int &Integration)
{
	switch (Integration)
	{
		case 1:
			return string("Direct summation");
		case 2:
			return string("Asymptotic expansion");
		case 3:
			return string("Recursion relations");
	}
	return string("Unknown integration");
}


// Since we are finding atan(delta) instead of delta directly, some of the results
//  are in the wrong range.
void FixPhase(double &PhaseShift, int &LValue, int &IsTriplet)
{
	double Pi = 4.0 * atan(1.0);

	if (LValue == 0 && IsTriplet == 0) {  // ^1S
		if (PhaseShift > 0.0) {
			PhaseShift = PhaseShift - Pi;
		}
	}
	else if (LValue == 0 && IsTriplet == 1) {  // ^3S
		if (PhaseShift > 0.0) {
			PhaseShift = PhaseShift - Pi;
		}
	}
	else if (LValue >= 3) {  // ^1F and higher
		if (PhaseShift > Pi) {
			PhaseShift = PhaseShift - Pi;
		}
	}
}


// Modified from http://www.dreamincode.net/code/snippet1102.htm
string GetDateTime(void)
{
	//Find the current time
	time_t curtime = time(0); 

	//convert it to tm
	tm now=*localtime(&curtime); 

	//BUFSIZ is standard macro that expands to a integer constant expression 
	//that is greater then or equal to 256. It is the size of the stream buffer 
	//used by setbuf()
	char dest[BUFSIZ]={0};

	//Format string determines the conversion specification's behaviour
	const char format[]="%x %X"; 

	//strftime - converts date and time to a string
	if (strftime(dest, sizeof(dest)-1, format, &now)>0) {
		return string(dest);
	}
	else 
		cerr << "strftime failed. Errno code: " << errno << endl;
	return string("");
}
