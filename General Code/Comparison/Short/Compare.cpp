//
//
// @TODO: Enable my ordering as well
//
//

#include <iostream>
#include <iomanip>
#include <string>
#include <cerrno>
#include <vector>
#include <fstream>
#include <math.h>
//#include <mkl_lapack.h>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <IL/il.h>
#include <IL/ilu.h>

using namespace std;
int ReadShortHeader(ifstream &FileShortRange, int &Omega, int &LValue, int &IsTriplet, int &Formalism, int &Ordering, int &NumShortTerms, int &NumSets, double &Alpha1, double &Beta1, double &Gamma1, double &Alpha2, double &Beta2, double &Gamma2);


int main(int argc, char *argv[])
{
	ifstream FileMatrixElem1, FileMatrixElem2;
	int Omega1, NumTerms1, LValue1, IsTriplet1, Formalism1, Ordering1, NumSets1;
	double Alpha11, Beta11, Gamma11, Alpha12, Beta12, Gamma12;
	int Omega2, NumTerms2, LValue2, IsTriplet2, Formalism2, Ordering2, NumSets2;
	double Alpha21, Beta21, Gamma21, Alpha22, Beta22, Gamma22;
	ILuint	ImgId;
	ILenum	Error;
	ILubyte *ImageData;
	double Elem1, Elem2, Diff, MaxDiff = 0.0;
	//double Multiplier;

	//Multiplier = atof(argv[4]);

	if (argc < 5) {
		cerr << "Not enough parameters on the command line." << endl;
		cerr << "Usage: Compare <matrixelements1> <matrixelements2> <PhiPhiResults> <PhiHPhiResults>" << endl;
		cerr << "Example: Compare matrixelements1.psh matrixelements2.psh phiphiresults.png phihphiresults.png" << endl << endl;
		return 1;
	}
	
	if (string(argv[1]) == string(argv[2])) {
		cerr << "There really is no point in comparing a file to itself, is there?" << endl;
		return 2;
	}

	FileMatrixElem1.open(argv[1], ios::in | ios::binary);
	if (FileMatrixElem1.fail()) {
		cerr << "Unable to open file " << argv[1] << " for reading." << endl;
		return 3;
	}

	FileMatrixElem2.open(argv[2], ios::in | ios::binary);
	if (FileMatrixElem2.fail()) {
		cerr << "Unable to open file " << argv[2] << " for reading." << endl;
		return 4;
	}

	
	// Initialize DevIL.
	ilInit();
	// Generate the main image name to use.
	ilGenImages(1, &ImgId);
	// Bind this image name.
	ilBindImage(ImgId);
	// Enable file overwrites
	ilEnable(IL_FILE_OVERWRITE);
	


	
	//FileMatrixElem1 >> Omega;
	//FileMatrixElem1 >> Alpha;
	//FileMatrixElem1 >> Beta;
	//FileMatrixElem1 >> Gamma;
	//FileMatrixElem1 >> NumShortTerms;

	int err = ReadShortHeader(FileMatrixElem1, Omega1, LValue1, IsTriplet1, Formalism1, Ordering1, NumTerms1, NumSets1, Alpha11, Beta11, Gamma11, Alpha12, Beta12, Gamma12);
	if (err == -1) {
		return 5;
	}

	err = ReadShortHeader(FileMatrixElem2, Omega2, LValue2, IsTriplet2, Formalism2, Ordering2, NumTerms2, NumSets2, Alpha21, Beta21, Gamma21, Alpha22, Beta22, Gamma22);
	if (err == -1) {
		return 6;
	}


	cout << Omega1 << " " << Alpha11 << " " << Beta11 << " " << Gamma11 << " " << NumTerms1 << endl;
	cout << Omega2 << " " << Alpha21 << " " << Beta21 << " " << Gamma21 << " " << NumTerms2 << endl;

	if ((Omega1 != Omega2) || (Alpha11 != Alpha21) || (Beta11 != Beta21) || (Gamma11 != Gamma21) || (NumTerms1 != NumTerms2)) {
		//@TODO: Cleanup!
		cout << "File parameters are not the same...exiting." << endl;
		return 7;
	}




	//@TODO: Check for failures.
	ilOriginFunc(IL_ORIGIN_UPPER_LEFT);
	ilTexImage(NumTerms1, NumTerms1, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);
	ImageData = ilGetData();

	//for (int i = 0; i < NumShortTerms1*NumShortTerms1; i++) {
	//	ImageData[i*3+0] = 0;
	//	ImageData[i*3+1] = 0;
	//	FileMatrixElem1.read((char*)&Elem1, sizeof(double));
	//	FileMatrixElem2.read((char*)&Elem2, sizeof(double));
	//	Diff = fabs(Elem1-Elem2);
	//	if (Diff > MaxDiff)
	//		MaxDiff = Diff;
	//	//cout << i << " " << Diff << " " << min(int(Diff * Multiplier), 255) << endl;
	//	ImageData[i*3+2] = min(int(Diff * Multiplier), 255);
	//}

	// Read PhiPhi parts
	for (int i = 0; i < NumTerms1*NumTerms1; i++) {
		FileMatrixElem1.read((char*)&Elem1, sizeof(double));
		FileMatrixElem2.read((char*)&Elem2, sizeof(double));
		Diff = fabs((Elem1-Elem2)/((Elem1+Elem2)/2.0))*100;

		// We will have division by 0 if this is the case.
		if (Elem1 == 0.0 || Elem2 == 0.0)
			Diff = 0.0;
		if (Diff > MaxDiff)
			MaxDiff = Diff;

		if (Diff <= 1e-13) {
			// Black
			ImageData[i*3+0] = 0;
			ImageData[i*3+1] = 0;
			ImageData[i*3+2] = 0;
		}
		else if (Diff <= 1e-12) {
			// Dark Blue
			ImageData[i*3+0] = 51;
			ImageData[i*3+1] = 102;
			ImageData[i*3+2] = 153;
		}
		else if (Diff <= 1e-11) {
			// Blue
			ImageData[i*3+0] = 0;
			ImageData[i*3+1] = 143;
			ImageData[i*3+2] = 204;
		}
		else if (Diff <= 1e-10) {
			// Light Blue
			ImageData[i*3+0] = 71;
			ImageData[i*3+1] = 200;
			ImageData[i*3+2] = 255;
		}
		else if (Diff <= 1e-9) {
			// Cyan
			ImageData[i*3+0] = 71;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 218;
		}
		else if (Diff <= 1e-8) {
			// Dark Green
			ImageData[i*3+0] = 71;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 126;
		}
		else if (Diff <= 1e-7) {
			// Green
			ImageData[i*3+0] = 108;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 71;
		}
		else if (Diff <= 1e-6) {
			// Olive Green
			ImageData[i*3+0] = 200;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 71;
		}
		else if (Diff <= 1e-5) {
			// Yellow
			ImageData[i*3+0] = 255;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 71;
		}
		else if (Diff <= 1e-4) {
			// Light Orange
			ImageData[i*3+0] = 255;
			ImageData[i*3+1] = 163;
			ImageData[i*3+2] = 71;
		}
		else if (Diff <= 1e-3) {
			// Dark Orange
			ImageData[i*3+0] = 255;
			ImageData[i*3+1] = 133;
			ImageData[i*3+2] = 10;
		}
		else {  // Really large difference!
			// Red
			ImageData[i*3+0] = 255;
			ImageData[i*3+1] = 61;
			ImageData[i*3+2] = 61;
		}
	}
	//for (int i = 0; i < NumShortTerms1*NumShortTerms1; i++) {
	//	FileMatrixElem1.read((char*)&Elem1, sizeof(double));
	//	FileMatrixElem2.read((char*)&Elem2, sizeof(double));
	//	Diff = fabs(Elem1-Elem2);
	//	if (Diff > MaxDiff)
	//		MaxDiff = Diff;

	//	if (Diff <= 1e-13) {
	//		// Black
	//		ImageData[i*3+0] = 0;
	//		ImageData[i*3+1] = 0;
	//		ImageData[i*3+2] = 0;
	//	}
	//	else if (Diff <= 1e-10) {
	//		// Green
	//		ImageData[i*3+0] = 0;
	//		ImageData[i*3+1] = 255;
	//		ImageData[i*3+2] = 0;
	//	}
	//	else if (Diff <= 1e-9) {
	//		// Grey
	//		ImageData[i*3+0] = 128;
	//		ImageData[i*3+1] = 128;
	//		ImageData[i*3+2] = 128;
	//	}
	//	else if (Diff <= 1e-8) {
	//		// Yellow
	//		ImageData[i*3+0] = 255;
	//		ImageData[i*3+1] = 255;
	//		ImageData[i*3+2] = 0;
	//	}
	//	else if (Diff <= 1e-7) {
	//		// Magenta
	//		ImageData[i*3+0] = 255;
	//		ImageData[i*3+1] = 0;
	//		ImageData[i*3+2] = 255;
	//	}
	//	else if (Diff <= 1e-6) {
	//		// Blue
	//		ImageData[i*3+0] = 0;
	//		ImageData[i*3+1] = 0;
	//		ImageData[i*3+2] = 255;
	//	}
	//	else if (Diff <= 1e-5) {
	//		// Cyan
	//		ImageData[i*3+0] = 0;
	//		ImageData[i*3+1] = 255;
	//		ImageData[i*3+2] = 255;
	//	}
	//	else if (Diff <= 1e-4) {
	//		// White
	//		ImageData[i*3+0] = 255;
	//		ImageData[i*3+1] = 255;
	//		ImageData[i*3+2] = 255;
	//	}
	//	else {  // Really large difference!
	//		// Red
	//		ImageData[i*3+0] = 255;
	//		ImageData[i*3+1] = 0;
	//		ImageData[i*3+2] = 0;
	//	}
	//}
	
	cout << "PhiPhi Maximum Difference: " << MaxDiff << endl;

	iluFlipImage();
	if (!ilSaveImage(argv[3])) {
		cout << "Could not write output file...exiting." << endl;
		//cout << ilGetError() << endl;  //@TODO: Use ILU to generate error string.
		//@TODO: Cleanup!
		return 5;
	}


	MaxDiff = 0.0;
	// Read PhiHPhi parts
	for (int i = 0; i < NumTerms1*NumTerms1; i++) {
		FileMatrixElem1.read((char*)&Elem1, sizeof(double));
		FileMatrixElem2.read((char*)&Elem2, sizeof(double));
		Diff = fabs((Elem1-Elem2)/((Elem1+Elem2)/2.0))*100;

		// We will have division by 0 if this is the case.
		if (Elem1 == 0.0 || Elem2 == 0.0)
			Diff = 0.0;
		if (Diff > MaxDiff)
			MaxDiff = Diff;

		if (Diff <= 1e-13) {
			// Black
			ImageData[i*3+0] = 0;
			ImageData[i*3+1] = 0;
			ImageData[i*3+2] = 0;
		}
		else if (Diff <= 1e-12) {
			// Dark Blue
			ImageData[i*3+0] = 51;
			ImageData[i*3+1] = 102;
			ImageData[i*3+2] = 153;
		}
		else if (Diff <= 1e-11) {
			// Blue
			ImageData[i*3+0] = 0;
			ImageData[i*3+1] = 143;
			ImageData[i*3+2] = 204;
		}
		else if (Diff <= 1e-10) {
			// Light Blue
			ImageData[i*3+0] = 71;
			ImageData[i*3+1] = 200;
			ImageData[i*3+2] = 255;
		}
		else if (Diff <= 1e-9) {
			// Cyan
			ImageData[i*3+0] = 71;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 218;
		}
		else if (Diff <= 1e-8) {
			// Dark Green
			ImageData[i*3+0] = 71;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 126;
		}
		else if (Diff <= 1e-7) {
			// Green
			ImageData[i*3+0] = 108;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 71;
		}
		else if (Diff <= 1e-6) {
			// Olive Green
			ImageData[i*3+0] = 200;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 71;
		}
		else if (Diff <= 1e-5) {
			// Yellow
			ImageData[i*3+0] = 255;
			ImageData[i*3+1] = 255;
			ImageData[i*3+2] = 71;
		}
		else if (Diff <= 1e-4) {
			// Light Orange
			ImageData[i*3+0] = 255;
			ImageData[i*3+1] = 163;
			ImageData[i*3+2] = 71;
		}
		else if (Diff <= 1e-3) {
			// Dark Orange
			ImageData[i*3+0] = 255;
			ImageData[i*3+1] = 133;
			ImageData[i*3+2] = 10;
		}
		else {  // Really large difference!
			// Red
			ImageData[i*3+0] = 255;
			ImageData[i*3+1] = 61;
			ImageData[i*3+2] = 61;
		}
	}
	//for (int i = 0; i < NumShortTerms1*NumShortTerms1; i++) {
	//	FileMatrixElem1.read((char*)&Elem1, sizeof(double));
	//	FileMatrixElem2.read((char*)&Elem2, sizeof(double));
	//	Diff = fabs(Elem1-Elem2);
	//	if (Diff > MaxDiff)
	//		MaxDiff = Diff;

	//	if (Diff <= 1e-13) {
	//		// Black
	//		ImageData[i*3+0] = 0;
	//		ImageData[i*3+1] = 0;
	//		ImageData[i*3+2] = 0;
	//	}
	//	else if (Diff <= 1e-10) {
	//		// Green
	//		ImageData[i*3+0] = 0;
	//		ImageData[i*3+1] = 255;
	//		ImageData[i*3+2] = 0;
	//	}
	//	else if (Diff <= 1e-9) {
	//		// Grey
	//		ImageData[i*3+0] = 128;
	//		ImageData[i*3+1] = 128;
	//		ImageData[i*3+2] = 128;
	//	}
	//	else if (Diff <= 1e-8) {
	//		// Yellow
	//		ImageData[i*3+0] = 255;
	//		ImageData[i*3+1] = 255;
	//		ImageData[i*3+2] = 0;
	//	}
	//	else if (Diff <= 1e-7) {
	//		// Magenta
	//		ImageData[i*3+0] = 255;
	//		ImageData[i*3+1] = 0;
	//		ImageData[i*3+2] = 255;
	//	}
	//	else if (Diff <= 1e-6) {
	//		// Blue
	//		ImageData[i*3+0] = 0;
	//		ImageData[i*3+1] = 0;
	//		ImageData[i*3+2] = 255;
	//	}
	//	else if (Diff <= 1e-5) {
	//		// Cyan
	//		ImageData[i*3+0] = 0;
	//		ImageData[i*3+1] = 255;
	//		ImageData[i*3+2] = 255;
	//	}
	//	else if (Diff <= 1e-4) {
	//		// White
	//		ImageData[i*3+0] = 255;
	//		ImageData[i*3+1] = 255;
	//		ImageData[i*3+2] = 255;
	//	}
	//	else {  // Really large difference!
	//		// Red
	//		ImageData[i*3+0] = 255;
	//		ImageData[i*3+1] = 0;
	//		ImageData[i*3+2] = 0;
	//	}
	//}

	iluFlipImage();
	if (!ilSaveImage(argv[4])) {
		cout << "Could not write output file...exiting." << endl;
		//cout << ilGetError() << endl;  //@TODO: Use ILU to generate error string.
		//@TODO: Cleanup!
		return 5;
	}


	ilDeleteImages(1, &ImgId);
	FileMatrixElem1.close();
	FileMatrixElem2.close();

	cout << "PhiHPhi Maximum Difference: " << MaxDiff << endl;

	return 0;
}


// Reads in the short-range file header
int ReadShortHeader(ifstream &FileShortRange, int &Omega, int &LValue, int &IsTriplet, int &Formalism, int &Ordering, int &NumShortTerms, int &NumSets, double &Alpha1, double &Beta1, double &Gamma1, double &Alpha2, double &Beta2, double &Gamma2)
{
	int MagicNum, Version, HeaderLen, DataFormat, NumShortTerms1, NumShortTerms2;
	int Integration, VarLen;

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
		cout << "This is not a valid Ps-H file...exiting." << endl;
		return -1;
	}

	if (Version < 1 || Version > 8) {
		cout << "This is not a valid Ps-H file...exiting." << endl;
		return -1;
	}

	if (HeaderLen != 80 && HeaderLen != 104) {
		cout << "This is not a valid Ps-H file...exiting." << endl;
		return -1;
	}

	if (DataFormat != 8) {
		cout << "This is not a valid Ps-H file...exiting." << endl;
		return -1;
	}

	if ((Formalism != 1 && Formalism != 2) || (IsTriplet != 0 && IsTriplet != 1) || (Ordering != 0 && Ordering != 1)) {
		cout << "This is not a valid Ps-H file...exiting." << endl;
		return -1;
	}

	if (NumSets != 1 && NumSets != 2) {
		cout << "This is not a valid Ps-H file...exiting." << endl;
		return -1;
	}

	/*if (NumShortTerms1 != NumShortTerms2) {
		cout << "This is not a valid Ps-H file...exiting." << endl;  // Cannot handle two different values for this yet.
		return -1;
	}*/
	NumShortTerms = NumShortTerms1;

	FileShortRange.read((char*)&Alpha1, 8);
	FileShortRange.read((char*)&Beta1, 8);
	FileShortRange.read((char*)&Gamma1, 8);
	if (NumSets == 2) {
		//read (FileShortRange) Alpha2, Beta2, Gamma2
		FileShortRange.read((char*)&Alpha2, 8);
		FileShortRange.read((char*)&Beta2, 8);
		FileShortRange.read((char*)&Gamma2, 8);
	}

	FileShortRange.read((char*)&VarLen, 4);

	return 0;
}
