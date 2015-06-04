#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>

using namespace std;

bool DetermineProblem(string &Line, int &Version, int &LValue, int &IsTriplet, int &NumSets);
bool DetermineOrdering(string &Line, int &Ordering);
bool DetermineIntegration(string &Line, int &Integration);
bool GetNonlinearParams(ifstream &InFile, int Version, int NumSets, vector <double> &Alpha, vector <double> &Beta, vector <double> &Gamma, vector <double> &Eps12, vector <double> &Eps13, int &ExpR12Len, int &ExpR13Len, int &NumTerms, vector <int> &NumTermsSets);


// trim and reduce are from http://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
const std::string trim(const std::string& pString,
					   const std::string& pWhitespace = " \t")
{
	const size_t beginStr = pString.find_first_not_of(pWhitespace);
	if (beginStr == std::string::npos)
	{
		// no content
		return "";
	}

	const size_t endStr = pString.find_last_not_of(pWhitespace);
	const size_t range = endStr - beginStr + 1;

	return pString.substr(beginStr, range);
}

const std::string reduce(const std::string& pString,
						 const std::string& pFill = " ",
						 const std::string& pWhitespace = " \t")
{
	// trim first
	std::string result(trim(pString, pWhitespace));

	// replace sub ranges
	size_t beginSpace = result.find_first_of(pWhitespace);
	while (beginSpace != std::string::npos)
	{
		const size_t endSpace =
						result.find_first_not_of(pWhitespace, beginSpace);
		const size_t range = endSpace - beginSpace;

		result.replace(beginSpace, range, pFill);

		const size_t newStart = beginSpace + pFill.length();
		beginSpace = result.find_first_of(pWhitespace, newStart);
	}

	return result;
}


// split is from http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while(std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	return split(s, delim, elems);
}


string replaceStrChar(string str, const string& replace, char ch) {
	// set our locator equal to the first appearance of any character in replace
	size_t found = str.find_first_of(replace);

	while (found != string::npos) { // While our position in the sting is in range.
		str[found] = ch; // Change the character at position.
		found = str.find_first_of(replace, found+1); // Relocate again.
	}

	return str; // return our new string.
}



int main(int argc, char *argv[])
{
	ifstream InFile;
	ofstream OutFile;
	string OutName, Temp1, Temp2, Temp3, Line;
	//double Alpha1, Beta1, Gamma1, Alpha2, Beta2, Gamma2, PhiPhi, PhiHPhi;
	vector <double> Alpha, Beta, Gamma, Eps12, Eps13;
	int Omega, NumTerms, Ordering, IsTriplet, Formalism;
	int LValue, Integration, Version, HeaderLen, DataFormat, NumSets;
	int ExpR12Len, ExpR13Len;
	vector <int> NumTermsSets;
	
	int MagicNum = 0x31487350; // PsH1 in hexadecimal (with reverse due to endianness)

	if (argc < 2) {
		cout << "Usage: Binary inputfile.txt outputfile.psh" << endl;
		exit(1);
	}
	
#ifdef __linux__
	cout << "Calling dos2unix to resolve any EOL problems" << endl;
	string d2ustr = string("dos2unix ") + string(argv[1]);
	int d2uret = system(d2ustr.c_str());
	cout << endl;
#endif//__linux__

	InFile.open(argv[1]);
	if (InFile.fail()) {
		cerr << "Unable to open file " << argv[1] << " for reading." << endl;
		exit(2);
	}

	if (argc == 3) {
		OutName = argv[2];
	}
	else {  //TODO: This is very hackish! Replace with a much more robust solution, such as the solution from the Boost library.
		OutName = argv[1];
		int OutNameLen = OutName.length();
		OutName[OutNameLen-3] = 'p';
		OutName[OutNameLen-2] = 's';
		OutName[OutNameLen-1] = 'h';
		cout << "Outputting on file " << OutName << endl;
	}

	OutFile.open(OutName.c_str(), ios::out | ios::binary);
	if (OutFile.fail()) {
		cerr << "Unable to open file " << argv[2] << " for writing." << endl;
		exit(3);
	}

	//
	// Read input file
	//

	getline(InFile, Temp1);
	Temp1 = trim(Temp1);
	if (!DetermineProblem(Temp1, Version, LValue, IsTriplet, NumSets)) {
		cout << "Could not determine type of file: " << Temp1 << endl;
		exit(4);
	}

	getline(InFile, Temp1);
	Temp1 = trim(Temp1);
	if (!DetermineOrdering(Temp1, Ordering)) {
		cout << "Could not determine ordering: " << Temp1 << endl;
		exit(5);
	}

	getline(InFile, Temp1);
	Temp1 = trim(Temp1);
	if (!DetermineIntegration(Temp1, Integration)) {
		cout << "Could not determine type of file: " << Temp1 << endl;
		exit(4);
	}

	getline(InFile, Temp1);  // We do not need to know the type of eigenvalue routine used, since the energies are not included.

	//@TODO: Do this in a function.
	getline(InFile, Temp1);
	Temp1 = reduce(Temp1);
	vector<string> s = split(Temp1, ' ');
	if (s.size() < 3) {
		cout << "Cannot determine omega value: " << Temp1 << endl;
		exit(5);
	}
	Omega = atoi(s[2].c_str());  //@TODO: Could use method from http://forums.codeguru.com/showthread.php?t=231054

	Alpha.resize(NumSets);
	Beta.resize(NumSets);
	Gamma.resize(NumSets);
	Eps12.resize(NumSets);
	cout << NumSets << endl;
	Eps13.resize(NumSets);

	if (!GetNonlinearParams(InFile, Version, NumSets, Alpha, Beta, Gamma, Eps12, Eps13, ExpR12Len, ExpR13Len, NumTerms, NumTermsSets)) {
		cout << "Could not read nonlinear parameters." << endl;
		exit(6);
	}

	cout << "Number of terms: " << NumTerms << endl;
	if (NumTerms < 0) {
		cout << "Cannot determine number of terms: " << NumTerms << endl;
		exit(5);
	}
	

	//
	// Begin writing output file
	//

	OutFile.write((char*)&MagicNum, 4);  // Assuming 32-bit integer
	//Version = 1;
	OutFile.write((char*)&Version, 4);
	HeaderLen = 80;  // Just trying something right now
	OutFile.write((char*)&HeaderLen, 4);
	DataFormat = 8;  // Only have the option of double precision right now
	OutFile.write((char*)&DataFormat, 4);

	OutFile.write((char*)&Omega, 4);
	OutFile.write((char*)&NumTerms, 4);  // Number of terms for omega value
	if (NumSets == 5) {
		OutFile.write((char*)&NumSets, 4);  // Number of terms in each set
		for (int i = 0; i < 5; i++) 
			OutFile.write((char*)&NumTermsSets[i], 4);
	}
	else {
		int Temp = 0;
		OutFile.write((char*)&Temp, 4);  // Number of terms used
	}
	OutFile.write((char*)&LValue, 4);
	Formalism = 1;  // Only supporting the first formalism right now.
	OutFile.write((char*)&Formalism, 4);
	OutFile.write((char*)&IsTriplet, 4);
	OutFile.write((char*)&Ordering, 4);
	OutFile.write((char*)&Integration, 4);
	OutFile.write((char*)&NumSets, 4);
	for (int i = 0; i < NumSets; i++) {
		OutFile.write((char*)&Alpha[i], 8);
		OutFile.write((char*)&Beta[i], 8);
		OutFile.write((char*)&Gamma[i], 8);
		if (Version == 9 || Version == 10) {  // Exponential terms in r12 and r13
			OutFile.write((char*)&Eps12[i], 8);
			OutFile.write((char*)&Eps13[i], 8);
			double Blank = 0.0;
			OutFile.write((char*)&Blank, 8);  // Possibility for r23 expansion
		}
	}
	if (Version == 9 || Version == 10) {  // Exponential terms in r12 and r13
		OutFile.write((char*)&ExpR12Len, 4);
		OutFile.write((char*)&ExpR13Len, 4);
		int Blank = 1;
		OutFile.write((char*)&Blank, 4);  // Possibility for r23 expansion
	}

	int VarLen = 0;
	OutFile.write((char*)&VarLen, 4);  // Descriptor field length

	//
	// Read matrix elements
	//
	cout << "PhiPhi matrix elements" << endl;
	getline(InFile, Temp1);  // Skip these
	getline(InFile, Temp1);  //  lines
	
	int Pos = (int)InFile.tellg();
	// Read / write PhiPhi
	for (int i = 0; i < NumTerms; i++) {
		cout << i+1 << endl;
		for (int j = 0; j < NumTerms; j++) {
			getline(InFile, Temp1);
			Temp1 = reduce(Temp1);
			//cout << Temp1 << endl;
			s = split(Temp1, ' ');
			s[2] = replaceStrChar(s[2], string("Dd"), 'e');
			double PhiPhi = atof(s[2].c_str());
			OutFile.write((char*)&PhiPhi, 8);
		}
	}

	// Seek back to read / write skipped PhiHPhi entries
#ifdef _WIN32
	//InFile.seekg(Pos-5, ios::beg);  // Why?!
	InFile.seekg(Pos, ios::beg);
#else
	InFile.seekg(Pos, ios::beg);
#endif

//#ifdef _WIN32test
//	// Not sure why I have to skip the following lines.
//	//  tellg/seekg does not take this back to the proper spot in Windows.
//	//  It probably has something to do with getline.
//	getline(InFile, Temp1);
//	getline(InFile, Temp1);
//	getline(InFile, Temp1);
//#endif//_WIN32
	cout << endl << "PhiHPhi matrix elements" << endl;
	for (int i = 0; i < NumTerms; i++) {
		cout << i+1 << endl;
		for (int j = 0; j < NumTerms; j++) {
			getline(InFile, Temp1);
			//cout << Temp1 << endl;
			Temp1 = reduce(Temp1);
			s = split(Temp1, ' ');
			s[3] = replaceStrChar(s[3], string("Dd"), 'e');
			double PhiHPhi = atof(s[3].c_str());
			OutFile.write((char*)&PhiHPhi, 8);
		}
	}

	InFile.close();
	OutFile.close();

	return 0;
}


bool DetermineProblem(string &Line, int &Version, int &LValue, int &IsTriplet, int &NumSets)
{
	if (Line == "S-Wave Singlet Ps-H") {
		Version = 1;
		LValue = 0;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "S-Wave Triplet Ps-H") {
		Version = 1;
		LValue = 0;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "S-Wave Singlet Ps-H - Grad-Grad Formalism") {
		Version = 1;
		LValue = 0;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "S-Wave Triplet Ps-H - Grad-Grad Formalism") {
		Version = 1;
		LValue = 0;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "S-Wave Singlet Ps-H - Laplacian Formalism") {
		Version = 1;
		LValue = 0;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "S-Wave Triplet Ps-H - Laplacian Formalism") {
		Version = 1;
		LValue = 0;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "S-Wave Singlet Ps-H - Grad-Grad Formalism - Sectors") {
		Version = 7;
		LValue = 0;
		IsTriplet = 0;
		NumSets = 5;
	}
	else if (Line == "S-Wave Triplet Ps-H - Grad-Grad Formalism - Sectors") {
		Version = 7;
		LValue = 0;
		IsTriplet = 1;
		NumSets = 5;
	}
	else if (Line == "S-Wave Singlet Ps-H - Laplacian Formalism - Sectors") {
		Version = 7;
		LValue = 0;
		IsTriplet = 0;
		NumSets = 5;
	}
	else if (Line == "S-Wave Triplet Ps-H - Laplacian Formalism - Sectors") {
		Version = 7;
		LValue = 0;
		IsTriplet = 1;
		NumSets = 5;
	}
	else if (Line == "S-Wave Singlet Ps-H - Grad-Grad Formalism - Exponential") {
		Version = 9;
		LValue = 0;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "S-Wave Triplet Ps-H - Grad-Grad Formalism - Exponential") {
		Version = 9;
		LValue = 0;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "S-Wave Singlet Ps-H - Laplacian Formalism - Exponential") {
		Version = 9;
		LValue = 0;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "S-Wave Triplet Ps-H - Laplacian Formalism - Exponential") {
		Version = 9;
		LValue = 0;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "P-Wave Singlet Ps-H: 1st formalism") {
		Version = 2;
		LValue = 1;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "P-Wave Triplet Ps-H: 1st formalism") {
		Version = 2;
		LValue = 1;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "P-Wave Singlet Ps-H: 2nd formalism") {
		Version = 3;
		LValue = 1;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "P-Wave Triplet Ps-H: 2nd formalism") {
		Version = 3;
		LValue = 1;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "P-Wave Singlet Ps-H: 1st formalism / 2 sets") {
		Version = 4;
		LValue = 1;
		IsTriplet = 0;
		NumSets = 2;
	}
	else if (Line == "P-Wave Triplet Ps-H: 1st formalism / 2 sets") {
		Version = 4;
		LValue = 1;
		IsTriplet = 1;
		NumSets = 2;
	}
	else if (Line == "P-Wave Singlet Ps-H: 2nd formalism / 2 sets") {
		Version = 5;
		LValue = 1;
		IsTriplet = 0;
		NumSets = 2;
	}
	else if (Line == "P-Wave Triplet Ps-H: 2nd formalism / 2 sets") {
		Version = 5;
		LValue = 1;
		IsTriplet = 1;
		NumSets = 2;
	}
	else if (Line == "D-Wave Singlet Ps-H: 1st formalism") {
		Version = 6;
		LValue = 2;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "D-Wave Triplet Ps-H: 1st formalism") {
		Version = 6;
		LValue = 2;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "D-Wave Singlet Ps-H - Exponential") {
		Version = 10;
		LValue = 2;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "D-Wave Triplet Ps-H - Exponential") {
		Version = 10;
		LValue = 2;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "F-Wave Singlet Ps-H") {
		Version = 8;
		LValue = 3;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "F-Wave Triplet Ps-H") {
		Version = 8;
		LValue = 3;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "G-Wave Singlet Ps-H") {
		Version = 8;
		LValue = 4;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "G-Wave Triplet Ps-H") {
		Version = 8;
		LValue = 4;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "H-Wave Singlet Ps-H") {
		Version = 8;
		LValue = 5;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "H-Wave Triplet Ps-H") {
		Version = 8;
		LValue = 5;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "I-Wave Singlet Ps-H") {
		Version = 8;
		LValue = 6;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "I-Wave Triplet Ps-H") {
		Version = 8;
		LValue = 6;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "K-Wave Singlet Ps-H") {
		Version = 8;
		LValue = 7;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "K-Wave Triplet Ps-H") {
		Version = 8;
		LValue = 7;
		IsTriplet = 1;
		NumSets = 1;
	}
	else if (Line == "L-Wave Singlet Ps-H") {
		Version = 8;
		LValue = 8;
		IsTriplet = 0;
		NumSets = 1;
	}
	else if (Line == "L-Wave Triplet Ps-H") {
		Version = 8;
		LValue = 8;
		IsTriplet = 1;
		NumSets = 1;
	}
	else {
		return false;
	}

	return true;
}


bool DetermineOrdering(string &Line, int &Ordering)
{
	if (Line == "Using Denton's ordering") {
		Ordering = 0;
	}
	else if (Line == "Using Peter Van Reeth's ordering") {
		Ordering = 1;
	}
	else {
		return false;
	}

	return true;
}


bool DetermineIntegration(string &Line, int &Integration)
{
	if (Line == "Integration Technique: Direct Summation") {
		Integration = 1;
	}
	else if (Line == "Integration Technique: Asymptotic Expansion") {
		Integration = 2;
	}
	else if (Line == "Integration Technique: Recursion Relations") {
		Integration = 3;
	}
	else {
		return false;
	}

	return true;
}


bool GetNonlinearParams(ifstream &InFile, int Version, int NumSets, vector <double> &Alpha, vector <double> &Beta, vector <double> &Gamma, vector <double> &Eps12, vector <double> &Eps13, int &ExpR12Len, int &ExpR13Len, int &NumTerms, vector <int> &NumTermsSets)
{
	string Line;

	if (Version == 9) {  // Exponential falloff terms
		InFile >> Line;
		InFile >> Line;
		InFile >> Alpha[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Beta[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Gamma[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Eps12[0];
		InFile >> Eps13[0];
		for (int i = 0; i < 5; i++)
			InFile >> Line;
		InFile >> ExpR12Len;
		for (int i = 0; i < 5; i++)
			InFile >> Line;
		InFile >> ExpR13Len;
		cout << Alpha[0] << " " << Beta[0] << " " << Gamma[0] << " " << Eps12[0] << " " << Eps13[0] << endl;
	} else if (Version == 10) {  // D-wave exponential falloff terms
		InFile >> Line;
		InFile >> Line;
		InFile >> Alpha[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Beta[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Gamma[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Eps12[0];
		InFile >> Eps13[0];
		for (int i = 0; i < 5; i++)
			InFile >> Line;
		InFile >> ExpR12Len;
		for (int i = 0; i < 5; i++)
			InFile >> Line;
		InFile >> ExpR13Len;
		cout << Alpha[0] << " " << Beta[0] << " " << Gamma[0] << " " << Eps12[0] << " " << Eps13[0] << endl;
	} else if (NumSets == 1) {
		InFile >> Line;
		InFile >> Line;
		InFile >> Alpha[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Beta[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Gamma[0];
	}
	else if (NumSets == 2) {
		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		InFile >> Alpha[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		InFile >> Beta[0];
		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		InFile >> Gamma[0];

		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		InFile >> Alpha[1];
		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		InFile >> Beta[1];
		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		InFile >> Gamma[1];
	}
	else if (NumSets == 5) {
		getline(InFile, Line);
		for (int i = 0; i < NumSets; i++) {
			InFile >> Line;
			InFile >> Line;
			InFile >> Line;
			InFile >> Alpha[i];
			InFile >> Beta[i];
			InFile >> Gamma[i];
		}
	}

	/*getline(InFile, Line);
	Line = reduce(Line);
	vector<string> s = split(Line, ' ');
	if (s.size() < 3)
		return false;
	Alpha = atof(s[2].c_str());  //@TODO: Could use method from http://forums.codeguru.com/showthread.php?t=231054*/


	if (NumSets == 5) {
		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		InFile >> NumTerms;
		NumTermsSets.resize(NumSets);
		for (int i = 0; i < 5; i++) InFile >> Line;
		for (int i = 0; i < 5; i++) InFile >> NumTermsSets[i];
	}
	else {
		InFile >> Line;
		InFile >> Line;
		InFile >> Line;
		//InFile >> Line;
		InFile >> NumTerms;
	}

	return true;
}
