//
// Ps-H P-Wave Scattering Code
//  written by Denton Woods
//


#ifndef NO_MPI
	#define USE_MPI
#endif

#include "Ps-H Scattering.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctime>
#include <cerrno>
#include <sstream>
#include "tinyxml2.h"
#include <algorithm>
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include <omp.h>

long double	PI;
#ifdef USE_MPI
	MPI_File MpiLog;
#endif

using namespace tinyxml2;

//#define VERBOSE


int main(int argc, char *argv[])
{
	//vector<rPowers> PowerTable;
	int Omega, NumShortTerms, Ordering, IsTriplet;
	double Alpha, Beta, Gamma, Mu, Kappa;
	//double **PhiHPhi, **PhiPhi;
	vector <double> B(1, 0.0), ARow(1, 0.0);//, ShortTerms;
	double SLS = 0.0, SLC = 0.0;
	QuadPoints q;
	int sf;
	double r2Cusp, r3Cusp;
	int Node = 0, TotalNodes = 1;
	ifstream ParameterFile, FileShortRange;
	ofstream OutFile;
	//int NodeStart, NodeEnd;
	time_t TimeStart, TimeEnd;
#ifdef USE_MPI
	char ProcessorName[MPI_MAX_PROCESSOR_NAME];
	MPI_Status MpiStatus;
	int MpiError, ProcNameLen;
#endif

	//@TODO: Will MPI time be different?
	TimeStart = time(NULL);
	//omp_set_num_threads(1);

	// Initialize global constants.
	PI = 3.1415926535897932384626433832795029L;

	if (argc < 5) {
		cerr << "Not enough parameters on the command line." << endl;
		cerr << "Usage: Scattering kappa parameterfile.txt shortrangefile.psh results.txt" << endl;
		exit(1);
	}
	cout << setprecision(18);
	cout.setf(ios::showpoint);

	// MPI initialization
	//@TODO: Check MpiError.
#ifdef USE_MPI
	MpiError = MPI_Init(&argc, &argv);  // All MPI programs start with MPI_Init; all 'N' processes exist thereafter.
	MpiError = MPI_Comm_size(MPI_COMM_WORLD, &TotalNodes);  // Find out how big the SPMD world is
	MpiError = MPI_Comm_rank(MPI_COMM_WORLD, &Node);  // and what this process's rank is.
	MpiError = MPI_Get_processor_name(ProcessorName, &ProcNameLen);
	string LogName;
	if (argc > 5)
		LogName = argv[5];
	else
		LogName = "test.log";
	//MpiError = MPI_File_open(MPI_COMM_WORLD, (char*)LogName.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &MpiLog);
	//cout << "Node " << Node << " of " << TotalNodes << " on " << ProcessorName << endl;

	// The following code here just shows what the nodes are and how many threads each has.
	stringstream ThreadStringStream;
	int ThreadStrLen;
	int MaxThreads = omp_get_max_threads();
	ThreadStringStream << MaxThreads << " threads on Node " << Node << " of " << TotalNodes << " on " << ProcessorName;
	if (Node == 0) {
		cout << endl << ThreadStringStream.str() << endl;
		for (int i = 1; i < TotalNodes; i++) {
			MpiError = MPI_Recv(&ThreadStrLen, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &MpiStatus);
			char *StringStream = new char[ThreadStrLen+1];
			MpiError = MPI_Recv(StringStream, ThreadStrLen, MPI_CHAR, i, 0, MPI_COMM_WORLD, &MpiStatus);
			StringStream[ThreadStrLen] = NULL;
			cout << StringStream << endl;
		}
		cout << endl;
	}
	else {
		ThreadStrLen = ThreadStringStream.str().length();
		const char *ThreadString = ThreadStringStream.str().c_str();
		MpiError = MPI_Send(&ThreadStrLen, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MpiError = MPI_Send((void*)ThreadString, ThreadStrLen, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}
#else
	Node = 0;
#endif

	if (Node == 0) {
		ParameterFile.open(argv[2]);
		OutFile.open(argv[4]);
		if (!ParameterFile.is_open()) {
			cout << "Could not open parameter file...exiting." << endl;
			FinishMPI();
			exit(3);
		}
		if (!OutFile.is_open()) {
			cout << "Could not open output file...exiting." << endl;
			FinishMPI();
			exit(4);
		}
		OutFile << setprecision(18);
		OutFile.setf(ios::showpoint);
		ShowDateTime(OutFile);

		ReadParamFile(ParameterFile, q, Mu, r2Cusp, r3Cusp);

		// Read in short-range short-range elements.  These have already been calculated by the PsHBound program.
		//  We are reading in only the binary versions (they were originally text files).
		//
		/*FileShortRange.open(argv[3], ios::in | ios::binary);
		if (FileShortRange.fail()) {
			cerr << "Unable to open file " << argv[3] << " for reading." << endl;
			FinishMPI();
			exit(2);
		}*/

		// Create short-range short-range terms
		//
		// Read in omega and nonlinear parameters.
		if (!ReadShortHeader(string(argv[3]), Omega, IsTriplet, Ordering, NumShortTerms, Alpha, Beta, Gamma)) {
			cerr << argv[3] << " is not a valid Ps-H P-wave short-range file...exiting." << endl;
			FinishMPI();
			exit(3);
		}

		// Calculate number of terms for a given omega. Do we need to compare to the one in the short-range file above?
		NumShortTerms = CalcPowerTableSize(Omega);
		Kappa = atof(argv[1]);

		if (IsTriplet == 0) {
			sf = 1;
			cout << "P-Wave Singlet Ps-H" << endl;
		}
		else if (IsTriplet == 1) {
			sf = -1;
			cout << "P-Wave Triplet Ps-H" << endl;
		}
		else {
			cout << "IsTriplet parameter incorrect in " << argv[3] << endl;
			exit(6);
		}

		cout << "Omega = " << Omega << endl;
		cout << "Mu = " << Mu << endl;
		cout << "Alpha = " << Alpha << endl;
		cout << "Beta  = " << Beta << endl;
		cout << "Gamma = " << Gamma << endl;
		cout << "Kappa = " << Kappa << endl;

		cout << endl << "Number of quadrature points" << endl;
		cout << "Long-long:                     " << q.LongLong_r1 << " " << q.LongLong_r2Leg << " " << q.LongLong_r2Lag << " " << q.LongLong_r3Leg << " " << q.LongLong_r3Lag << " " << q.LongLong_r12 << " " << q.LongLong_r13 << " " << q.LongLong_phi23 << endl;
		cout << "Long-long 2/r23 term:          " << q.LongLongr23_r1 << " " << q.LongLongr23_r2Leg << " " << q.LongLongr23_r2Lag << " " << q.LongLongr23_r3Leg << " " << q.LongLongr23_r3Lag << " " << q.LongLongr23_phi12 << " " << q.LongLongr23_r13 << " " << q.LongLongr23_r23 << endl;
		cout << "Short-long with qi = 0:        " << q.ShortLong_r1 << " " << q.ShortLong_r2Leg << " " << q.ShortLong_r2Lag << " " << q.ShortLong_r3Leg << " " << q.ShortLong_r3Lag << " " << q.ShortLong_r12 << " " << q.ShortLong_r13 << " " << q.ShortLong_phi23 << endl;
		cout << "Short-long 2/r23 with qi = 0:  " << q.ShortLongr23_r1 << " " << q.ShortLongr23_r2Leg << " " << q.ShortLongr23_r2Lag << " " << q.ShortLongr23_r3Leg << " " << q.ShortLongr23_r3Lag << " " << q.ShortLongr23_r12 << " " << q.ShortLongr23_phi13 << " " << q.ShortLongr23_r23 << endl;
		cout << "Short-long (full) with qi > 0: " << q.ShortLongQiGt0_r1 << " " << q.ShortLongQiGt0_r2Leg << " " << q.ShortLongQiGt0_r2Lag << " " << q.ShortLongQiGt0_r3Leg << " " << q.ShortLongQiGt0_r3Lag << " " << q.ShortLongQiGt0_r12 << " " << q.ShortLongQiGt0_r13 << " " << q.ShortLongQiGt0_phi23 << endl;
		cout << endl;
		cout << "Cusp parameters" << endl;
		cout << r2Cusp << " " << r3Cusp << endl;
		cout << endl;

		/*// Allocate memory for the overlap matrix and point PhiPhiP to rows of PhiPhi so we
		//  can access it like a 2D array but have it in contiguous memory for LAPACK.
		PhiPhi = new double*[NumShortTerms*2];
		PhiPhi[0] = new double[NumShortTerms*NumShortTerms*4];
		for (int i = 1; i < NumShortTerms*2; i++) {
			PhiPhi[i] = PhiPhi[i-1] + NumShortTerms*2;
		}
		// Read in the <phi|phi> matrix elements.
		FileShortRange.read((char*)PhiPhi[0], NumShortTerms*NumShortTerms*4*sizeof(double));

		// Allocate memory for the <phi|H|phi> matrix and point PhiHPhiP to rows of PhiHPhi so we
		//  can access it like a 2D array but have it in contiguous memory for LAPACK.
		PhiHPhi = new double*[NumShortTerms*2];
		PhiHPhi[0] = new double[NumShortTerms*NumShortTerms*4];
		for (int i = 1; i < NumShortTerms*2; i++) {
			PhiHPhi[i] = PhiHPhi[i-1] + NumShortTerms*2;
		}
		// Read in the <phi|H|phi> matrix elements.
		FileShortRange.read((char*)PhiHPhi[0], NumShortTerms*NumShortTerms*4*sizeof(double));

		// Allocate matrix of short-range - short-range terms.
		ShortTerms.resize(NumShortTerms*NumShortTerms*4);

		//@TODO: Remove next line.
		//memset(ShortTerms, 0, NumShortTerms*NumShortTerms*4*sizeof(double));  // Initialize to all 0.*/
	}

	//@TODO: Check results of MpiError.
#ifdef USE_MPI
	MpiError = MPI_Bcast(&Omega, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&NumShortTerms, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&Ordering, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&IsTriplet, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&Mu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&Kappa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&Alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&Beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&Gamma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&r2Cusp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&r3Cusp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&sf, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MpiError = MPI_Bcast(&q, sizeof(q), MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
#endif


	ARow.resize(NumShortTerms*2+1);
	B.resize(NumShortTerms*2+1);
	//@TODO: Remove next line.
	//memset(ARow, 0, (NumShortTerms+1)*sizeof(double));  // Initialize to all 0.
	//memset(B, 0, (NumShortTerms+1)*sizeof(double));  // Initialize to all 0.


	vector <rPowers> PowerTableQi0, PowerTableQiGt0;
	vector <vector <rPowers> > PowerTableQi0Array(TotalNodes), PowerTableQiGt0Array(TotalNodes);
	vector <double> AResultsQi0, BResultsQi0, AResultsQiGt0, BResultsQiGt0;
	vector <double> AResultsQi0Final, BResultsQi0Final, AResultsQiGt0Final, BResultsQiGt0Final;
	int NumTerms, NumTermsQi0, NumTermsQiGt0;

	NumTerms = CalcPowerTableSize(Omega);

	if (Node == 0) {
		NumTermsQi0 = CalcPowerTableSizeQi0(Omega, Ordering, 0, NumShortTerms);
		NumTermsQiGt0 = CalcPowerTableSizeQiGt0(Omega, Ordering, 0, NumShortTerms);

		// The *2 comes from the 2 types of symmetry
		AResultsQi0.resize(NumTermsQi0*2, 0.0);
		AResultsQiGt0.resize(NumTermsQiGt0*2, 0.0);
		BResultsQi0.resize(NumTermsQi0*2, 0.0);
		BResultsQiGt0.resize(NumTermsQiGt0*2, 0.0);

		AResultsQi0Final.resize(NumTermsQi0*2, 0.0);
		AResultsQiGt0Final.resize(NumTermsQiGt0*2, 0.0);
		BResultsQi0Final.resize(NumTermsQi0*2, 0.0);
		BResultsQiGt0Final.resize(NumTermsQiGt0*2, 0.0);

		PowerTableQi0.resize(NumTermsQi0*2, rPowers(Alpha, Beta, Gamma));
		PowerTableQiGt0.resize(NumTermsQiGt0*2, rPowers(Alpha, Beta, Gamma));
		GenOmegaPowerTableQi0(Omega, Ordering, PowerTableQi0, 0, NumShortTerms-1);
		GenOmegaPowerTableQiGt0(Omega, Ordering, PowerTableQiGt0, 0, NumShortTerms-1);
	}


#ifdef USE_MPI
	double NumTermsQi0Proc, NumTermsQiGt0Proc;
	vector <int> NumTermsQi0Array(TotalNodes), NumTermsQiGt0Array(TotalNodes);

	MpiError = MPI_Barrier(MPI_COMM_WORLD);
	// Tell all processes what terms they should be evaluating.
	if (Node == 0) {
		NumTermsQi0Proc = (double)NumTermsQi0 / (double)TotalNodes;  //@TODO: Need the typecast?
		NumTermsQiGt0Proc = (double)NumTermsQiGt0 / (double)TotalNodes;  //@TODO: Need the typecast?
		int Qi0Pos = (int)NumTermsQi0Proc, QiGt0Pos = (int)NumTermsQiGt0Proc;
		for (int i = 1; i < TotalNodes; i++) {
			//@TODO: Do we need to save these values in an array?
			NumTermsQi0Array[i] = int(NumTermsQi0Proc * (i+1)) - int(NumTermsQi0Proc * i);
			NumTermsQiGt0Array[i] = int(NumTermsQiGt0Proc * (i+1)) - int(NumTermsQiGt0Proc * i);
			MpiError = MPI_Send(&NumTermsQi0Array[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MpiError = MPI_Send(&NumTermsQiGt0Array[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			cout << "Node " << i << ": " << NumTermsQi0Array[i] << " " << NumTermsQiGt0Array[i] << endl;

			PowerTableQi0Array[i].resize(NumTermsQi0Array[i]);
			PowerTableQiGt0Array[i].resize(NumTermsQiGt0Array[i]);
			for (int j = 0; j < NumTermsQi0Array[i]; j++, Qi0Pos++) {
				PowerTableQi0Array[i][j] = PowerTableQi0[Qi0Pos];
				//cout << i << ": " << PowerTableQi0[Qi0Pos].ki + PowerTableQi0[Qi0Pos].li + PowerTableQi0[Qi0Pos].mi + PowerTableQi0[Qi0Pos].ni + PowerTableQi0[Qi0Pos].pi + PowerTableQi0[Qi0Pos].qi << " - " <<
				//	PowerTableQi0[Qi0Pos].ki << " " << PowerTableQi0[Qi0Pos].li << " " << PowerTableQi0[Qi0Pos].mi << " " << PowerTableQi0[Qi0Pos].ni << " " << PowerTableQi0[Qi0Pos].pi << " " << PowerTableQi0[Qi0Pos].qi << endl;
			}
#ifdef VERBOSE
			cout << endl;
#endif
			for (int j = 0; j < NumTermsQiGt0Array[i]; j++, QiGt0Pos++) {
				PowerTableQiGt0Array[i][j] = PowerTableQiGt0[QiGt0Pos];
#ifdef VERBOSE
				cout << i << ": " << PowerTableQiGt0[QiGt0Pos].ki + PowerTableQiGt0[QiGt0Pos].li + PowerTableQiGt0[QiGt0Pos].mi + PowerTableQiGt0[QiGt0Pos].ni + PowerTableQiGt0[QiGt0Pos].pi + PowerTableQiGt0[QiGt0Pos].qi << " - " <<
					PowerTableQiGt0[QiGt0Pos].ki << " " << PowerTableQiGt0[QiGt0Pos].li << " " << PowerTableQiGt0[QiGt0Pos].mi << " " << PowerTableQiGt0[QiGt0Pos].ni << " " << PowerTableQiGt0[QiGt0Pos].pi << " " << PowerTableQiGt0[QiGt0Pos].qi << endl;
#endif//VERBOSE
			}

			//@TODO: How safe is this?
			MpiError = MPI_Send(&PowerTableQi0Array[i][0], sizeof(rPowers)*NumTermsQi0Array[i], MPI_BYTE, i, 0, MPI_COMM_WORLD);
			MpiError = MPI_Send(&PowerTableQiGt0Array[i][0], sizeof(rPowers)*NumTermsQiGt0Array[i], MPI_BYTE, i, 0, MPI_COMM_WORLD);
		}

		NumTermsQi0Array[0] = int(NumTermsQi0Proc);
		NumTermsQiGt0Array[0] = int(NumTermsQiGt0Proc);
		NumTermsQi0 = NumTermsQi0Array[0];
		NumTermsQiGt0 = NumTermsQiGt0Array[0];

		//AResultsQi0.resize(NumTermsQi0*2, 0.0);
		//AResultsQiGt0.resize(NumTermsQiGt0*2, 0.0);
		//BResultsQi0.resize(NumTermsQi0*2, 0.0);
		//BResultsQiGt0.resize(NumTermsQiGt0*2, 0.0);

		//@TODO: Temporary?
		int NumTermsQi0Temp = CalcPowerTableSizeQi0(Omega, Ordering, 0, NumShortTerms);
		for (int i = 0; i < NumTermsQi0; i++) {
			PowerTableQi0[NumTermsQi0+i].ki = PowerTableQi0[NumTermsQi0Temp+i].ki;
			PowerTableQi0[NumTermsQi0+i].li = PowerTableQi0[NumTermsQi0Temp+i].li;
			PowerTableQi0[NumTermsQi0+i].mi = PowerTableQi0[NumTermsQi0Temp+i].mi;
			PowerTableQi0[NumTermsQi0+i].ni = PowerTableQi0[NumTermsQi0Temp+i].ni;
			PowerTableQi0[NumTermsQi0+i].pi = PowerTableQi0[NumTermsQi0Temp+i].pi;
			PowerTableQi0[NumTermsQi0+i].qi = PowerTableQi0[NumTermsQi0Temp+i].qi;
		}
		int NumTermsQiGt0Temp = CalcPowerTableSizeQiGt0(Omega, Ordering, 0, NumShortTerms);
		for (int i = 0; i < NumTermsQiGt0; i++) {
			PowerTableQiGt0[NumTermsQiGt0+i].ki = PowerTableQiGt0[NumTermsQiGt0Temp+i].ki;
			PowerTableQiGt0[NumTermsQiGt0+i].li = PowerTableQiGt0[NumTermsQiGt0Temp+i].li;
			PowerTableQiGt0[NumTermsQiGt0+i].mi = PowerTableQiGt0[NumTermsQiGt0Temp+i].mi;
			PowerTableQiGt0[NumTermsQiGt0+i].ni = PowerTableQiGt0[NumTermsQiGt0Temp+i].ni;
			PowerTableQiGt0[NumTermsQiGt0+i].pi = PowerTableQiGt0[NumTermsQiGt0Temp+i].pi;
			PowerTableQiGt0[NumTermsQiGt0+i].qi = PowerTableQiGt0[NumTermsQiGt0Temp+i].qi;
		}

		cout << "Node " << 0 << ": " << NumTermsQi0Array[0] << " " << NumTermsQiGt0Array[0] << endl;
		for (int i = 0; i < NumTermsQiGt0; i++) {
			//cout << 0 << ": " << PowerTableQi0[i].ki + PowerTableQi0[i].li + PowerTableQi0[i].mi + PowerTableQi0[i].ni + PowerTableQi0[i].pi + PowerTableQi0[i].qi << " - " <<
			//		PowerTableQi0[i].ki << " " << PowerTableQi0[i].li << " " << PowerTableQi0[i].mi << " " << PowerTableQi0[i].ni << " " << PowerTableQi0[i].pi << " " << PowerTableQi0[i].qi << endl;
#ifdef VERBOSE
			cout << 0 << ": " << PowerTableQiGt0[i].ki + PowerTableQiGt0[i].li + PowerTableQiGt0[i].mi + PowerTableQiGt0[i].ni + PowerTableQiGt0[i].pi + PowerTableQiGt0[i].qi << " - " <<
					PowerTableQiGt0[i].ki << " " << PowerTableQiGt0[i].li << " " << PowerTableQiGt0[i].mi << " " << PowerTableQiGt0[i].ni << " " << PowerTableQiGt0[i].pi << " " << PowerTableQiGt0[i].qi << endl;
#endif//VERBOSE
		}
		cout << endl;
	}
	else {
		MpiError = MPI_Recv(&NumTermsQi0, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &MpiStatus);
		MpiError = MPI_Recv(&NumTermsQiGt0, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &MpiStatus);

		PowerTableQi0.resize(NumTermsQi0*2);
		PowerTableQiGt0.resize(NumTermsQiGt0*2);
		AResultsQi0.resize(NumTermsQi0*2);
		AResultsQiGt0.resize(NumTermsQiGt0*2);
		BResultsQi0.resize(NumTermsQi0*2);
		BResultsQiGt0.resize(NumTermsQiGt0*2);

		MpiError = MPI_Recv(&PowerTableQi0[0], sizeof(rPowers)*NumTermsQi0, MPI_INT, 0, 0, MPI_COMM_WORLD, &MpiStatus);
		MpiError = MPI_Recv(&PowerTableQiGt0[0], sizeof(rPowers)*NumTermsQiGt0, MPI_INT, 0, 0, MPI_COMM_WORLD, &MpiStatus);

		// Create phi2 terms
		for (int i = 0; i < NumTermsQi0; i++) {
			PowerTableQi0[NumTermsQi0+i].ki = PowerTableQi0[i].ki - 1;
			PowerTableQi0[NumTermsQi0+i].li = PowerTableQi0[i].li + 1;
			PowerTableQi0[NumTermsQi0+i].mi = PowerTableQi0[i].mi;
			PowerTableQi0[NumTermsQi0+i].ni = PowerTableQi0[i].ni;
			PowerTableQi0[NumTermsQi0+i].pi = PowerTableQi0[i].pi;
			PowerTableQi0[NumTermsQi0+i].qi = PowerTableQi0[i].qi;
		}
		for (int i = 0; i < NumTermsQiGt0; i++) {
			PowerTableQiGt0[NumTermsQiGt0+i].ki = PowerTableQiGt0[i].ki - 1;
			PowerTableQiGt0[NumTermsQiGt0+i].li = PowerTableQiGt0[i].li + 1;
			PowerTableQiGt0[NumTermsQiGt0+i].mi = PowerTableQiGt0[i].mi;
			PowerTableQiGt0[NumTermsQiGt0+i].ni = PowerTableQiGt0[i].ni;
			PowerTableQiGt0[NumTermsQiGt0+i].pi = PowerTableQiGt0[i].pi;
			PowerTableQiGt0[NumTermsQiGt0+i].qi = PowerTableQiGt0[i].qi;
		}
	}
	//cout << "Node " << Node << ": " << NodeStart << " " << NodeEnd << endl;
	MpiError = MPI_Barrier(MPI_COMM_WORLD);
#endif
//#ifdef USE_MPI
//	MpiError = MPI_Barrier(MPI_COMM_WORLD);
//
//	// Tell all processes what terms they should be evaluating.
//	if (Node == 0) {
//		NumTermsProc = (double)(NumShortTerms+1) / (double)TotalNodes;  //@TODO: Need the typecast?
//		for (int i = 1; i < TotalNodes; i++) {
//			NodeStart = NumTermsProc * i;
//			NodeEnd = NumTermsProc * (i+1) - 1;
//			MpiError = MPI_Send(&NodeStart, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
//			MpiError = MPI_Send(&NodeEnd, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
//		}
//		// This is the set of values for the root node to take.
//		NodeStart = 0;
//		NodeEnd = NumTermsProc * (Node+1) - 1;
//		}
//	else {
//		MpiError = MPI_Recv(&NodeStart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &MpiStatus);
//		MpiError = MPI_Recv(&NodeEnd, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &MpiStatus);
//	}
//	MpiError = MPI_Barrier(MPI_COMM_WORLD);
//	//@TODO: To clean this up some, we could just output this in the loop above for node 0.
//	cout << "Node " << Node << ": " << NodeStart << " " << NodeEnd << endl;
//#else
//	NumTermsProc = (double)(NumShortTerms+1);
//	NodeStart = 0;
//	NodeEnd = NumTermsProc * (Node+1) - 1;
//#endif//USE_MPI

	//TimeStart = time(NULL);
	//CalcARowAndBVector(Node, NumTerms, Omega, PowerTable, AResults, ARow, BResults, B, SLS, q, r2Cusp, r3Cusp, Alpha, Beta, Gamma, Kappa, Mu, sf);
	CalcARowAndBVector(Node, NumTermsQi0, NumTermsQiGt0, Omega, PowerTableQi0, PowerTableQiGt0, AResultsQi0, AResultsQiGt0, ARow, BResultsQi0, BResultsQiGt0, B, SLS, SLC, q, r2Cusp, r3Cusp, Alpha, Beta, Gamma, Kappa, Mu, sf);
	//TimeEnd = time(NULL);
	//cout << "Time elapsed: " << difftime(TimeEnd, TimeStart) << endl;
	//OutFile << "Time elapsed: " << difftime(TimeEnd, TimeStart) << endl;
	

#ifdef USE_MPI
	if (Node == 0) {
		//@TODO: Temporary
		int NumTermsQi0Temp = CalcPowerTableSizeQi0(Omega, Ordering, 0, NumShortTerms);
		int NumTermsQiGt0Temp = CalcPowerTableSizeQiGt0(Omega, Ordering, 0, NumShortTerms);

		for (int i = 0; i < NumTermsQi0; i++) {
			AResultsQi0Final[i] = AResultsQi0[i];
			AResultsQi0Final[NumTermsQi0Temp+i] = AResultsQi0[NumTermsQi0+i];
			BResultsQi0Final[i] = BResultsQi0[i];
			BResultsQi0Final[NumTermsQi0Temp+i] = BResultsQi0[NumTermsQi0+i];
		}
		for (int i = 0; i < NumTermsQiGt0; i++) {
			AResultsQiGt0Final[i] = AResultsQiGt0[i];
			AResultsQiGt0Final[NumTermsQiGt0Temp+i] = AResultsQiGt0[NumTermsQiGt0+i];
			BResultsQiGt0Final[i] = BResultsQiGt0[i];
			BResultsQiGt0Final[NumTermsQiGt0Temp+i] = BResultsQiGt0[NumTermsQiGt0+i];
		}

		cout << " Node " << Node << " (" << ProcessorName << ") finished computation at " << ShowTime() << endl;
		// ResultsQi0 and ResultsQiGt0 already have process 0's results.
		//  Gather the results from the rest of the processes.
		int nQi0 = NumTermsQi0Array[0], nQiGt0 = NumTermsQiGt0Array[0];
		// IMPORTANT NOTE! i++ must be at the end, or else it increments before nQi0 and nQiGt0.
		for (int i = 1; i < TotalNodes; nQi0 += NumTermsQi0Array[i], nQiGt0 += NumTermsQiGt0Array[i], i++) {
			// Phi1 terms
			MpiError = MPI_Recv(&AResultsQi0Final[nQi0], NumTermsQi0Array[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &MpiStatus);
			MpiError = MPI_Recv(&BResultsQi0Final[nQi0], NumTermsQi0Array[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &MpiStatus);
			MpiError = MPI_Recv(&AResultsQiGt0Final[nQiGt0], NumTermsQiGt0Array[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &MpiStatus);
			MpiError = MPI_Recv(&BResultsQiGt0Final[nQiGt0], NumTermsQiGt0Array[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &MpiStatus);
			// Phi2 terms
			MpiError = MPI_Recv(&AResultsQi0Final[NumTermsQi0Temp+nQi0], NumTermsQi0Array[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &MpiStatus);
			MpiError = MPI_Recv(&BResultsQi0Final[NumTermsQi0Temp+nQi0], NumTermsQi0Array[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &MpiStatus);
			MpiError = MPI_Recv(&AResultsQiGt0Final[NumTermsQiGt0Temp+nQiGt0], NumTermsQiGt0Array[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &MpiStatus);
			MpiError = MPI_Recv(&BResultsQiGt0Final[NumTermsQiGt0Temp+nQiGt0], NumTermsQiGt0Array[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &MpiStatus);
		}
	}
	else {
		cout << " Node " << Node << " (" << ProcessorName << ") finished computation at " << ShowTime() << endl;
		// Phi1 terms
		MpiError = MPI_Send(&AResultsQi0[0], NumTermsQi0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MpiError = MPI_Send(&BResultsQi0[0], NumTermsQi0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MpiError = MPI_Send(&AResultsQiGt0[0], NumTermsQiGt0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MpiError = MPI_Send(&BResultsQiGt0[0], NumTermsQiGt0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		// Phi2 terms
		MpiError = MPI_Send(&AResultsQi0[NumTermsQi0], NumTermsQi0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MpiError = MPI_Send(&BResultsQi0[NumTermsQi0], NumTermsQi0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MpiError = MPI_Send(&AResultsQiGt0[NumTermsQiGt0], NumTermsQiGt0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MpiError = MPI_Send(&BResultsQiGt0[NumTermsQiGt0], NumTermsQiGt0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
#else
	AResultsQi0Final = AResultsQi0;
	BResultsQi0Final = BResultsQi0;
	AResultsQiGt0Final = AResultsQiGt0;
	BResultsQiGt0Final = BResultsQiGt0;
#endif

	if (Node == 0) {
		//@TODO: Put directly into A and B
		vector <double> AResults(NumShortTerms*2), BResults(NumShortTerms*2);
		CombineResults(Omega, Ordering, AResultsQi0Final, AResultsQiGt0Final, AResults, 0, NumShortTerms);
		CombineResults(Omega, Ordering, BResultsQi0Final, BResultsQiGt0Final, BResults, 0, NumShortTerms);

		for (int i = 0; i < NumShortTerms*2; i++) {
			ARow[i+1] = AResults[i];
		}
		for (int i = 0; i < NumShortTerms*2; i++) {
			B[i+1] = BResults[i];
		}

		cout << endl << endl;
		for (int i = 0; i < NumShortTerms*2+1; i++) {
			cout << i << " " << ARow[i] << endl;
		}
		cout << endl << endl;

		for (int i = 0; i < NumShortTerms*2+1; i++) {
			cout << i << " " << B[i] << endl;
		}
		cout << endl;

		/*// Construct the rest of the matrix A with the short-range - short-range terms.
		for (int i = 0; i < NumShortTerms*2; i++) {
			for (int j = 0; j < NumShortTerms*2; j++) {
				ShortTerms[i*NumShortTerms*2+j] = PhiHPhi[i][j] - 0.5*Kappa*Kappa * PhiPhi[i][j] + 1.5*PhiPhi[i][j];
			}
		}*/
	}

	if (Node == 0) {
		//@TODO: Move further up.
		// Output results to file
		if (sf == 1)
			OutFile << "P-Wave Singlet Ps-H" << endl;
		else if (sf == -1)
			OutFile << "P-Wave Triplet Ps-H" << endl;

		if (Ordering == 0)
			OutFile << "Using Denton's ordering" << endl;
		else
			OutFile << "Using Peter Van Reeth's ordering" << endl;

		OutFile << "Omega: " << Omega << endl;
		OutFile << "Number of terms: " << NumShortTerms << endl;
		OutFile << "Alpha: " << Alpha << "  Beta: " << Beta << "  Gamma: " << Gamma << endl;
		OutFile << "Mu: " << Mu << endl;
		OutFile << "Kappa: " << Kappa << endl;

		OutFile << endl << "Number of quadrature points" << endl;
		OutFile << "Long-long:                     " << q.LongLong_r1 << " " << q.LongLong_r2Leg << " " << q.LongLong_r2Lag << " " << q.LongLong_r3Leg << " " << q.LongLong_r3Lag << " " << q.LongLong_r12 << " " << q.LongLong_r13 << " " << q.LongLong_phi23 << endl;
		OutFile << "Long-long 2/r23 term:          " << q.LongLongr23_r1 << " " << q.LongLongr23_r2Leg << " " << q.LongLongr23_r2Lag << " " << q.LongLongr23_r3Leg << " " << q.LongLongr23_r3Lag << " " << q.LongLongr23_phi12 << " " << q.LongLongr23_r13 << " " << q.LongLongr23_r23 << endl;
		OutFile << "Short-long with qi = 0:        " << q.ShortLong_r1 << " " << q.ShortLong_r2Leg << " " << q.ShortLong_r2Lag << " " << q.ShortLong_r3Leg << " " << q.ShortLong_r3Lag << " " << q.ShortLong_r12 << " " << q.ShortLong_r13 << " " << q.ShortLong_phi23 << endl;
		OutFile << "Short-long 2/r23 with qi = 0:  " << q.ShortLongr23_r1 << " " << q.ShortLongr23_r2Leg << " " << q.ShortLongr23_r2Lag << " " << q.ShortLongr23_r3Leg << " " << q.ShortLongr23_r3Lag << " " << q.ShortLongr23_r12 << " " << q.ShortLongr23_phi13 << " " << q.ShortLongr23_r23 << endl;
		OutFile << "Short-long (full) with qi > 0: " << q.ShortLongQiGt0_r1 << " " << q.ShortLongQiGt0_r2Leg << " " << q.ShortLongQiGt0_r2Lag << " " << q.ShortLongQiGt0_r3Leg << " " << q.ShortLongQiGt0_r3Lag << " " << q.ShortLongQiGt0_r12 << " " << q.ShortLongQiGt0_r13 << " " << q.ShortLongQiGt0_phi23 << endl;
		OutFile << endl;
		OutFile << "Cusp parameters" << endl;
		OutFile << r2Cusp << " " << r3Cusp << endl;
		OutFile << endl;

		OutFile << "A matrix row" << endl;
		for (int i = 0; i < NumShortTerms*2+1; i++) {
			OutFile << i << " " <<  ARow[i] << endl;
		}
		OutFile << endl << "B vector" << endl;
		for (int i = 0; i < NumShortTerms*2+1; i++) {
			OutFile << i << " " << B[i] << endl;
		}

		cout << "SLS Term: " << SLS << endl;
		OutFile << endl << "SLS Term" << endl << SLS << endl;
		cout << "SLC Term: " << SLC << endl;
		OutFile << endl << "SLC Term" << endl << SLC << endl;
		cout << "SLC - CLS: " << SLC - B[0] << endl << endl;
		OutFile << "SLC - CLS: " << SLC - B[0] << endl << endl;

		/*double KohnPhase, InvKohnPhase, ComplexKohnPhase;
		KohnPhase = Kohn(NumShortTerms, ARow, B, ShortTerms, SLS);
		InvKohnPhase = InverseKohn(NumShortTerms, ARow, B, ShortTerms, SLS);
		ComplexKohnPhase = ComplexKohnT(NumShortTerms, ARow, B, ShortTerms, SLS);
		cout << "Kohn phase shift: " << KohnPhase << endl;
		cout << "Inverse Kohn phase shift: " << InvKohnPhase << endl;
		cout << "Complex Kohn phase shift: " << ComplexKohnPhase << endl << endl;
		OutFile << "Kohn phase shift: " << KohnPhase << endl;
		OutFile << "Inverse Kohn phase shift: " << InvKohnPhase << endl;
		OutFile << "Complex (T-matrix) Kohn phase shift: " << ComplexKohnPhase << endl << endl;

		double CrossSection = 4.0 * PI * sin(KohnPhase) * sin(KohnPhase);  // (2l+1) = 1 with l = 0
		cout << "Kohn partial wave cross section: " << CrossSection << endl;
		OutFile << "Kohn Partial wave cross section: " << CrossSection << endl;
		CrossSection = 4.0 * PI * sin(InvKohnPhase) * sin(InvKohnPhase);
		cout << "Inverse Kohn partial wave cross section: " << CrossSection << endl;
		OutFile << "Inverse Kohn partial wave cross section: " << CrossSection << endl;
		CrossSection = 4.0 * PI * sin(ComplexKohnPhase) * sin(ComplexKohnPhase);
		cout << "Complex (T-matrix) Kohn partial wave cross section: " << CrossSection << endl << endl;
		OutFile << "Complex (T-matrix) Kohn partial wave cross section: " << CrossSection << endl << endl;*/
	}


	// Cleanup
	if (Node == 0) {
		TimeEnd = time(NULL);
		cout << "Time elapsed: " << difftime(TimeEnd, TimeStart) << endl;
		OutFile << "Time elapsed: " << difftime(TimeEnd, TimeStart) << endl;
		ParameterFile.close();
		OutFile.close();
		/*FileShortRange.close();
		delete [] PhiHPhi[0];
		delete [] PhiHPhi;
		delete [] PhiPhi[0];
		delete [] PhiPhi;*/
	}

#ifdef USE_MPI
	//MpiError = MPI_File_close(&MpiLog);
	MpiError = MPI_Finalize();
#endif

	return 0;
}


// Reads in the parameter file, getting the number of integration points, mu, kappa, etc.
//  Comments are in the example parameterfile.txt.
void ReadParamFile(ifstream &ParameterFile, QuadPoints &q, double &Mu, double &r2Cusp, double &r3Cusp)
{
	string Line;

	getline(ParameterFile, Line);
	getline(ParameterFile, Line);

	getline(ParameterFile, Line);
	ParameterFile >> q.LongLong_r1 >> q.LongLong_r2Leg >> q.LongLong_r2Lag >> q.LongLong_r3Leg >> q.LongLong_r3Lag >> q.LongLong_r12 >> q.LongLong_r13 >> q.LongLong_phi23;
	getline(ParameterFile, Line);  getline(ParameterFile, Line);
	ParameterFile >> q.LongLongr23_r1 >> q.LongLongr23_r2Leg >> q.LongLongr23_r2Lag >> q.LongLongr23_r3Leg >> q.LongLongr23_r3Lag >> q.LongLongr23_phi12 >> q.LongLongr23_r13 >> q.LongLongr23_r23;
	getline(ParameterFile, Line);
	getline(ParameterFile, Line);
	getline(ParameterFile, Line);

	ParameterFile >> q.ShortLong_r1 >> q.ShortLong_r2Leg >> q.ShortLong_r2Lag >> q.ShortLong_r3Leg >> q.ShortLong_r3Lag >> q.ShortLong_r12 >> q.ShortLong_r13 >> q.ShortLong_phi23;
	getline(ParameterFile, Line);
	getline(ParameterFile, Line);
	ParameterFile >> q.ShortLongr23_r1 >> q.ShortLongr23_r2Leg >> q.ShortLongr23_r2Lag >> q.ShortLongr23_r3Leg >> q.ShortLongr23_r3Lag >> q.ShortLongr23_r12 >> q.ShortLongr23_phi13 >> q.ShortLongr23_r23;
	getline(ParameterFile, Line);
	getline(ParameterFile, Line);
	ParameterFile >> q.ShortLongQiGt0_r1 >> q.ShortLongQiGt0_r2Leg >> q.ShortLongQiGt0_r2Lag >> q.ShortLongQiGt0_r3Leg >> q.ShortLongQiGt0_r3Lag >> q.ShortLongQiGt0_r12 >> q.ShortLongQiGt0_r13 >> q.ShortLongQiGt0_phi23;

	getline(ParameterFile, Line);
	getline(ParameterFile, Line);
	getline(ParameterFile, Line);
	//@TODO: Probably want more than one set (another for the 1/r23 term integration).
	ParameterFile >> r2Cusp;
	ParameterFile >> r3Cusp;
	getline(ParameterFile, Line);
	getline(ParameterFile, Line);
	ParameterFile >> Mu;

	return;
}


#ifndef XMLCheckResult
	#define XMLCheckResult(a_eResult) if (a_eResult != XML_SUCCESS) { cout << "Error parsing short-range file." << endl; return false; }
#endif

#ifndef XMLGetElem
	#define XMLGetElem(Element, Head, Str) { Element = Head->FirstChildElement(Str);  if (Element == NULL) { cout << "Could not parse the " << Str << " element in short-range file." << endl; return false; } }
#endif

#ifndef XMLCheckResult2
	#define XMLCheckResult2(Element) 	if (Element == NULL) {	cout << "Could not parse element in short-range file." << endl; return false;	}
#endif

bool ReadShortHeader(string FileShort, int &Omega, int &IsTriplet, int &Ordering, int &NumShortTerms, double &Alpha, double &Beta, double &Gamma)
{
	int MagicNum, Version, HeaderLen, DataFormat, NumShortTerms1, NumShortTerms2, LValue, Formalism, IntType, NumSets, VarLen;

	//FileShortRange.read((char*)&MagicNum, 4);
	//FileShortRange.read((char*)&Version, 4);
	//FileShortRange.read((char*)&HeaderLen, 4);
	//FileShortRange.read((char*)&DataFormat, 4);

	//FileShortRange.read((char*)&Omega, 4);
	//FileShortRange.read((char*)&NumShortTerms1, 4);
	//FileShortRange.read((char*)&NumShortTerms2, 4);
	//FileShortRange.read((char*)&LValue, 4);
	//FileShortRange.read((char*)&Formalism, 4);
	//FileShortRange.read((char*)&IsTriplet, 4);
	//FileShortRange.read((char*)&Ordering, 4);
	//FileShortRange.read((char*)&IntType, 4);

	//FileShortRange.read((char*)&NumSets, 4);
	//FileShortRange.read((char*)&Alpha, 8);
	//FileShortRange.read((char*)&Beta, 8);
	//FileShortRange.read((char*)&Gamma, 8);

	//FileShortRange.read((char*)&VarLen, 4);

	//if (MagicNum != 0x31487350 || /*Version != 1 ||*/ HeaderLen != 80 || DataFormat != 8 || Omega < 0 || NumShortTerms1 < 1
	//	|| Alpha <= 0.0 || Beta <= 0.0 || Gamma <= 0.0 /*|| NumShortTerms1 != NumShortTerms2*/ || LValue != 1 || Formalism != 1 || NumSets != 1) {
	//	cout << Version << endl;
	//	return false;
	//}

	XMLDocument doc;
	XMLError eResult;
	XMLElement *Element;
	string Entry;

	eResult = doc.LoadFile(FileShort.c_str());
	XMLCheckResult(eResult);
	XMLElement *titleElement = doc.FirstChildElement("psh_data");
	XMLCheckResult2(titleElement);

	XMLElement *HeaderElement = titleElement->FirstChildElement("header");
	XMLCheckResult2(HeaderElement);

	XMLGetElem(Element, HeaderElement, "problem");
	string LString = Element->GetText();  //@TODO: Do we want to check this?

	XMLGetElem(Element, HeaderElement, "ordering");
	Entry = Element->GetText();
	if (Entry == "Denton")
		Ordering = 0;
	else if (Entry == "Peter")
		Ordering = 1;
	else {
		cout << "Ordering in short-range file has an invalid value...exiting." << endl;
		return false;
	}

	//@TODO: Create macro for reading in an integer
	XMLGetElem(Element, HeaderElement, "lvalue");
	eResult = Element->QueryIntText(&LValue);
	XMLCheckResult(eResult);

	XMLGetElem(Element, HeaderElement, "spin");
	Entry = Element->GetText();
	if (Entry == "Singlet")
		IsTriplet = 0;
	else if (Entry == "Triplet")
		IsTriplet = 1;

	XMLGetElem(Element, HeaderElement, "omega");
	eResult = Element->QueryIntText(&Omega);
	XMLCheckResult(eResult);

	XMLGetElem(Element, HeaderElement, "numterms");
	eResult = Element->QueryIntText(&NumShortTerms);
	XMLCheckResult(eResult);

	XMLGetElem(Element, HeaderElement, "numsets");
	eResult = Element->QueryIntText(&NumSets);
	XMLCheckResult(eResult);
	
	if (NumSets > 1) {
		//@TODO: May want to use NextSibling to read these in
		cout << "This can only use 1 set of nonlinear parameters." << endl;
		return false;
	}

	XMLElement *NonlinearElement = HeaderElement->FirstChildElement("nonlinear");
	XMLCheckResult2(NonlinearElement);

	XMLGetElem(Element, NonlinearElement, "alpha");
	eResult = Element->QueryDoubleText(&Alpha);
	XMLCheckResult(eResult);

	XMLGetElem(Element, NonlinearElement, "beta");
	eResult = Element->QueryDoubleText(&Beta);
	XMLCheckResult(eResult);

	XMLGetElem(Element, NonlinearElement, "gamma");
	eResult = Element->QueryDoubleText(&Gamma);
	XMLCheckResult(eResult);

	XMLElement *ExpElement = HeaderElement->FirstChildElement("extraexp");
	if (ExpElement != NULL) {
		cout << "This cannot use extra exponentials." << endl;
		return false;
	}

	if (Omega < 0 || NumShortTerms < 0 || Alpha <= 0.0 || Beta <= 0.0 || Gamma <= 0.0 || LValue != 1 || NumSets != 1)
		return false;

	return true;
}


// Modified from http://www.dreamincode.net/code/snippet1102.htm
void ShowDateTime(ofstream &OutFile)
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
	const char format[]="Program started on %A, %B %d %Y at %X"; 

	//strftime - converts date and time to a string
	if (strftime(dest, sizeof(dest)-1, format, &now)>0) {
		cout << dest << endl << endl;
		OutFile << dest << endl << endl;
	}
	else 
		cerr << "strftime failed. Errno code: " << errno << endl;
	return;
}


string ShowTime()
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
	const char format[]="%X"; 

	//strftime - converts date and time to a string
	if (strftime(dest, sizeof(dest)-1, format, &now)>0) {
		return string(dest);
	}
	else 
		cerr << "strftime failed. Errno code: " << errno << endl;
	return string();
}


void CalcARowAndBVector(int Node, int NumTermsQi0, int NumTermsQiGt0, int Omega, vector <rPowers> &PowerTableQi0, vector <rPowers> &PowerTableQiGt0, vector <double> &AResultsQi0,
			  vector <double> &AResultsQiGt0, vector <double> &ARow, vector <double> &BResultsQi0, vector <double> &BResultsQiGt0, vector <double> &B, double &SLS, double &SLC, QuadPoints &q, double r2Cusp,
			  double r3Cusp, double alpha, double beta, double gamma, double kappa, double mu, int sf)
{
	//int NumTermsSub = NodeEnd-NodeStart+1;
	vector <rPowers> PowerTableSub;
#ifdef USE_MPI
	int MpiError;
	MPI_Status MpiStatus;
	string Buffer;
#endif

	// Calculate CLC term separately (requires different integration than the PhiLS terms).
	if (Node == 0) {
		ARow[0] = 0.0;
		B[0] = 0.0;
		SLS = 0.0;

		if (Node == 0) cout << "Starting long-long calculations at " << ShowTime() << endl;
		GaussIntegrationPhi23_LongLong(q.LongLong_r1, q.LongLong_r2Leg, q.LongLong_r2Lag, q.LongLong_r3Leg, q.LongLong_r3Lag, q.LongLong_r12, q.LongLong_r13, q.LongLong_phi23, r2Cusp, r3Cusp, kappa, mu, sf, ARow[0], SLC, B[0], SLS);

		cout << "SLS w/o r23 term: " << SLS << endl;
		cout << "SLC w/o r23 term: " << SLC << endl;
		cout << "CLS w/o r23 term: " << B[0] << endl;
		cout << "CLC w/o r23 term: " << ARow[0] << endl;
		cout << endl;

		if (Node == 0) cout << endl << "Starting long-long r23 term calculations at " << ShowTime() << endl;
		double CLCTemp = 0.0, SLSTemp = 0.0, SLCTemp = 0.0, CLSTemp = 0.0;
		//GaussIntegrationPhi12_LongLong_R23Term(q.LongLongr23_r1, q.LongLongr23_r2Leg, q.LongLongr23_r2Lag, q.LongLongr23_r3Leg, q.LongLongr23_r3Lag, q.LongLongr23_phi12, q.LongLongr23_r13, q.LongLongr23_r23, r2Cusp, r3Cusp, kappa, mu, sf, CLCTemp, SLCTemp, CLSTemp, SLSTemp);
		GaussIntegrationPhi13_LongLong_R23Term(q.LongLongr23_r1, q.LongLongr23_r2Leg, q.LongLongr23_r2Lag, q.LongLongr23_r3Leg, q.LongLongr23_r3Lag, q.LongLongr23_phi12, q.LongLongr23_r13, q.LongLongr23_r23, r2Cusp, r3Cusp, kappa, mu, sf, CLCTemp, SLCTemp, CLSTemp, SLSTemp);

		cout << "SLS r23 term: " << SLSTemp << endl;
		cout << "SLC r23 term: " << SLCTemp << endl;
		cout << "CLS r23 term: " << CLSTemp << endl;
		cout << "CLC r23 term: " << CLCTemp << endl << endl;
		ARow[0] += CLCTemp;
		SLS += SLSTemp;
		SLC += SLCTemp;
		B[0] += CLSTemp;

		cout << "SLS Term: " << SLS << endl;
		cout << "SLC Term: " << SLC << endl;
		cout << "CLS Term: " << B[0] << endl;
		cout << "CLC Term: " << ARow[0] << endl << endl;
		cout << "SLC - CLS = " << SLC - B[0] << endl << endl;

		//// The first node is a special case, because the 0th term is CLC or CLS, which is calculated above.
		////NumTermsSub = NodeEnd-NodeStart;
		////NodeStart = 1;  // Original NodeStart is not modified, since it is passed by value (be careful of this in Fortran).
	}

	////////// Create a copy of a portion of the power table (only the portion that this node is computing).
	////////PowerTableSub.resize(NumTermsSub*2);
	////////AResults.resize(NumTermsSub*2);
	////////BResults.resize(NumTermsSub*2);
	//////////cout << endl << "Node " << Node << " set: [" << NumTermsSub << " " << NumTerms << " " << NodeStart << " " << NodeEnd << "]" << endl;
	////////for (int i = 0; i < NumTermsSub; i++) {
	////////	// -1 to compensate for 0 being CLC or CLS (not in PowerTable).
	////////	PowerTableSub[i] = PowerTable[NodeStart+i-1];  // Phi1 terms
	////////	PowerTableSub[NumTermsSub+i] = PowerTable[NumTerms+NodeStart+i-1];  // Phi2 terms
	////////}
	//////////cout << endl << "Node " << Node << " PowerTableSub first entry: [" << PowerTableSub[0].ki << " " << PowerTableSub[0].li << " " << PowerTableSub[0].mi << " " << PowerTableSub[0].ni << " " << PowerTableSub[0].pi << " " << PowerTableSub[0].qi << "]" << endl;

	////////if (Node == 0) cout << endl << "Starting short-long calculations at " << ShowTime() << endl;
	////////VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar(AResults, BResults, q.ShortLong_r1, q.ShortLong_r2Leg, q.ShortLong_r2Lag, q.ShortLong_r3Leg, q.ShortLong_r3Lag, q.ShortLong_r12, q.ShortLong_r13, q.ShortLong_phi23, r2Cusp, r3Cusp, kappa, mu, sf, NumTermsSub, PowerTableSub, Omega);
	//////////#ifdef USE_MPI
	////////////MpiError = MPI_Barrier(MPI_COMM_WORLD);
	//////////Buffer = "Finished PhiLCBar and PhiLSBar on node " + to_string(Node) + "\n";
	//////////MpiError = MPI_File_write_shared(MpiLog, (void*)Buffer.c_str(), Buffer.length(), MPI_CHAR, &MpiStatus);
	//////////#endif

	////////if (Node == 0) cout << endl << "Starting short-long r23 term calculations at " << ShowTime() << endl;
	////////VecGaussIntegrationPhi13_PhiLCBar_PhiLSBar_R23Term(AResults, BResults, q.ShortLongr23_r1, q.ShortLongr23_r2Leg, q.ShortLongr23_r2Lag, q.ShortLongr23_r3Leg, q.ShortLongr23_r3Lag, q.ShortLongr23_r12, q.ShortLongr23_phi13, q.ShortLongr23_r23, r2Cusp, r3Cusp, kappa, mu, sf, NumTermsSub, PowerTableSub, Omega);
	//////////#ifdef USE_MPI
	////////////MpiError = MPI_Barrier(MPI_COMM_WORLD);
	//////////Buffer = "Finished PhiLCBar and PhiLSBar R23 term on node " + to_string(Node) + "\n";
	//////////MpiError = MPI_File_write_shared(MpiLog, (void*)Buffer.c_str(), Buffer.length(), MPI_CHAR, &MpiStatus);
	//////////#endif


	if (NumTermsQi0 > 0) {  // Skips when no terms with qi == 0
		if (Node == 0) cout << "Starting short-long calculations at " << ShowTime() << endl;
		VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar(AResultsQi0, BResultsQi0, q.ShortLong_r1, q.ShortLong_r2Leg, q.ShortLong_r2Lag, q.ShortLong_r3Leg, q.ShortLong_r3Lag, q.ShortLong_r12, q.ShortLong_r13, q.ShortLong_phi23, r2Cusp, r3Cusp, kappa, mu, sf, NumTermsQi0, PowerTableQi0, Omega);
		#ifdef USE_MPI
		//MpiError = MPI_Barrier(MPI_COMM_WORLD);
		Buffer = "Finished short-long on node " + to_string(Node) + "\n";
		//MpiError = MPI_File_write_shared(MpiLog, (void*)Buffer.c_str(), Buffer.length(), MPI_CHAR, &MpiStatus);
		#endif

		if (Node == 0) cout << "Starting short-long r23 term calculations at " << ShowTime() << endl;
		VecGaussIntegrationPhi13_PhiLCBar_PhiLSBar_R23Term(AResultsQi0, BResultsQi0, q.ShortLongr23_r1, q.ShortLongr23_r2Leg, q.ShortLongr23_r2Lag, q.ShortLongr23_r3Leg, q.ShortLongr23_r3Lag, q.ShortLongr23_r12, q.ShortLongr23_phi13, q.ShortLongr23_r23, r2Cusp, r3Cusp, kappa, mu, sf, NumTermsQi0, PowerTableQi0, Omega);
		//VecGaussIntegrationPhi12_PhiLCBar_PhiLSBar_R23Term(AResultsQi0, BResultsQi0, q.ShortLongr23_r1, q.ShortLongr23_r2Leg, q.ShortLongr23_r2Lag, q.ShortLongr23_r3Leg, q.ShortLongr23_r3Lag, q.ShortLongr23_r12, q.ShortLongr23_phi13, q.ShortLongr23_r23, r2Cusp, r3Cusp, kappa, mu, sf, NumTermsQi0, PowerTableQi0, Omega);
		#ifdef USE_MPI
		//MpiError = MPI_Barrier(MPI_COMM_WORLD);
		Buffer = "Finished short-long r23 term on node " + to_string(Node) + "\n";
		//MpiError = MPI_File_write_shared(MpiLog, (void*)Buffer.c_str(), Buffer.length(), MPI_CHAR, &MpiStatus);
		#endif
	}

	if (NumTermsQiGt0 > 0) {  // Skips when no terms with qi > 0
		if (Node == 0) cout << "Starting short-long full calculations at " << ShowTime() << endl;
		VecGaussIntegrationPhi23_PhiLCBar_PhiLSBar_Full(AResultsQiGt0, BResultsQiGt0, q.ShortLongQiGt0_r1, q.ShortLongQiGt0_r2Leg, q.ShortLongQiGt0_r2Lag, q.ShortLongQiGt0_r3Leg, q.ShortLongQiGt0_r3Lag, q.ShortLongQiGt0_r12,q. ShortLongQiGt0_r13, q.ShortLongQiGt0_phi23, r2Cusp, r3Cusp, kappa, mu, sf, NumTermsQiGt0, PowerTableQiGt0, Omega);
		#ifdef USE_MPI
		Buffer = "Finished short-long full on node " + to_string(Node) + "\n";
		//MpiError = MPI_File_write_shared(MpiLog, (void*)Buffer.c_str(), Buffer.length(), MPI_CHAR, &MpiStatus);
		#endif
	}

	return;
}


// We calculate terms with qi == 0 and qi > 0 separately, so this recombines them back into a single set.
void CombineResults(int Omega, int Ordering, vector <double> &ResultsQi0, vector <double> &ResultsQiGt0, vector <double> &Results, int Start, int End)
{
	int NumTerm = 0;  // The number of the current term
	int om, ki, li, mi, ni, pi, qi;  // These are the exponents we are determining.
	int j = 0, k = 0;
	int Cur = 0;

	int NumTermsQi0 = ResultsQi0.size() / 2;
	int NumTermsQiGt0 = ResultsQiGt0.size() / 2;
	int NumTerms = Results.size() / 2;

	if (Ordering == 0) {  // My ordering
		int r = 0;
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
												Results[r] = ResultsQi0[j];
												Results[NumTerms+r] = ResultsQi0[NumTermsQi0+j];
												j++;
											}
											else {  // qi > 0
												Results[r] = ResultsQiGt0[k];
												Results[NumTerms+r] = ResultsQiGt0[NumTermsQiGt0+k];
												k++;
											}

											r = r + 1;
										}
										Cur = Cur + 1;
										if (Cur > End)
											return;
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
									if (I23 == 0) {
										Results[INX] = ResultsQi0[j];
										Results[NumTerms+INX] = ResultsQi0[NumTermsQi0+j];
										j++;
									}
									else {  // I23 > 0
										Results[INX] = ResultsQiGt0[k];
										Results[NumTerms+INX] = ResultsQiGt0[NumTermsQiGt0+k];
										k++;
									}

									INX = INX + 1;
								}
								Cur = Cur + 1;
								if (Cur > End)
									return;
							}
						}
					}
				}
			}
		}
	}

	return;
}


//void WriteProgress(int i, int N)
//{
//	int Node = 0;
//#ifdef USE_MPI
//	int MpiError = MPI_Comm_rank(MPI_COMM_WORLD, &Node); 
//#endif
//	if (Node != 0)
//		return;
//	cout << setprecision(5) << float(i)/float(N)*100.0;
//	fflush(stdout);
//	cout << "\b\b\b\b\b\b";
//	cout << setprecision(18);  // Restore original
//	return;
//}


void WriteProgress(string Desc, int &Prog, int i, int N)
{
	int Node = 0;
#ifdef USE_MPI
	int MpiError = MPI_Comm_rank(MPI_COMM_WORLD, &Node);
	MPI_Status MpiStatus;
#endif
	//if (Node != 0)
	//	return;
#pragma omp critical(output)
	{
		Prog++;
#ifdef USE_MPI
		string Buffer = ShowTime() + " - " + Desc + " on node " + to_string(Node) + ": " + to_string(Prog) + " of " + to_string(N) + "\n";
		//MpiError = MPI_File_write_shared(MpiLog, (void*)Buffer.c_str(), Buffer.length(), MPI_CHAR, &MpiStatus);
		//MpiError = MPI_File_sync(MpiLog);
		cerr << Buffer;
#endif
		if (Node == 0) {
			cout << setprecision(5) << float(Prog)/float(N)*100.0 << "\r";
			fflush(stdout);
			//cout << "\b\b\b\b\b\b";
			cout << setprecision(18);  // Restore original
		}
		else {
			cout << "\r" << endl;
		}
	}
	return;
}


// From http://notfaq.wordpress.com/2006/08/30/c-convert-int-to-string/
template <class T> string to_string(const T& t)
{
	stringstream ss;
	ss << t;
	return ss.str();
}
