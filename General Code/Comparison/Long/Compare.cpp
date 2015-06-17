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


int main(int argc, char *argv[])
{
	ifstream FileMatrixElem1, FileMatrixElem2;
	int Omega1, NumShortTerms1;
	double Alpha1, Beta1, Gamma1, Eta1, Mu1, Kappa1;
	int Omega2, NumShortTerms2;
	double Alpha2, Beta2, Gamma2, Eta2, Mu2, Kappa2;
	ILuint	ImgId;
	ILenum	Error;
	ILubyte *ImageData;
	double Elem1, Elem2, Diff, MaxDiff = 0.0;
	string Line1, Line2;

	//Multiplier = atof(argv[4]);

	if (argc < 5) {
		cerr << "Not enough parameters on the command line." << endl;
		cerr << "Usage: Compare <matrixelements1> <matrixelements2> <AResults> <BResults>" << endl;
		cerr << "Example: Compare matrixelements1.txt matrixelements2.txt AResults.png BResults.png" << endl << endl;
		return 1;
	}

	if (string(argv[1]) == string(argv[2])) {
		cerr << "There really is no point in comparing a file to itself, is there?" << endl;
		return 2;
	}

	FileMatrixElem1.open(argv[1]);
	if (FileMatrixElem1.fail()) {
		cerr << "Unable to open file " << argv[1] << " for reading." << endl;
		return 3;
	}

	FileMatrixElem2.open(argv[2]);
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

	for (int i = 1; i <= 3; i++) {
		getline(FileMatrixElem1, Line1);
		getline(FileMatrixElem2, Line2);
	}
	if (Line1 != Line2 || ((Line1 != "Singlet Case") && (Line1 != "Triplet Case"))) {
		cout << "Comparing different types of problems: " << Line1 << " " << Line2 << endl;
		return 4;
	}

	// Read in the various parameters for each input file.
	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> Omega1;
	FileMatrixElem2 >> Omega2;

	getline(FileMatrixElem1, Line1);
	getline(FileMatrixElem2, Line2);

	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> NumShortTerms1;
	FileMatrixElem2 >> NumShortTerms2;

	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> Alpha1;
	FileMatrixElem2 >> Alpha2;
	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> Beta1;
	FileMatrixElem2 >> Beta2;
	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> Gamma1;
	FileMatrixElem2 >> Gamma2;
	//FileMatrixElem1 >> Line1;
	//FileMatrixElem2 >> Line2;
	//FileMatrixElem1 >> Eta1;
	//FileMatrixElem2 >> Eta2;
	
	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> Mu1;
	FileMatrixElem2 >> Mu2;
	FileMatrixElem1 >> Line1;
	FileMatrixElem2 >> Line2;
	FileMatrixElem1 >> Kappa1;
	FileMatrixElem2 >> Kappa2;

	if ((Omega1 != Omega2) || (Alpha1 != Alpha2) || (Beta1 != Beta2) || (Gamma1 != Gamma2) || (NumShortTerms1 != NumShortTerms2) || (Mu1 != Mu2) || (Kappa1 != Kappa2)) {
		//@TODO: Cleanup!
		cout << "File parameters are not the same...exiting." << endl;
		return 4;
	}


	// We do not care what the number of integration points are, since they will change from file to file and should be different.
	for (int i = 1; i <= 17; i++) {
		getline(FileMatrixElem1, Line1);
		getline(FileMatrixElem2, Line2);
	}



	//@TODO: Check for failures.
	ilOriginFunc(IL_ORIGIN_UPPER_LEFT);
	ilTexImage(NumShortTerms1+1, 32, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);
	ImageData = ilGetData();


	// Read A-matrix parts
	cout << "A matrix" << endl;
	for (int i = 0; i < NumShortTerms1+1; i++) {
		FileMatrixElem1 >> Elem1;
		FileMatrixElem2 >> Elem2;
		// Computes the relative difference, not the absolute difference.
		Diff = fabs(Elem1-Elem2) / (fabs(Elem1+Elem2)/2.0);
		cout << Diff << " " << Elem1 << " " << Elem2 << endl;
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
	// Duplicate the first row of the image to the next 31.
	for (int i = 0; i < NumShortTerms1+1; i++) {
		for (int j = 1; j < 32; j++) {
			for (int k = 0; k < 3; k++) {
				ImageData[j*(NumShortTerms1+1)*3 + i*3 + k] = ImageData[i*3 + k];
			}
		}
	}


	cout << "A Matrix Maximum Difference: " << MaxDiff << endl;

	iluFlipImage();
	if (!ilSaveImage(argv[3])) {
		cout << "Could not write output file...exiting." << endl;
		Error = ilGetError();
		cout << "DevIL error number: " << Error << endl;
		//cout << iluErrorString(IL_NO_ERROR ) << endl;  //@TODO: Use ILU to generate error string.
		//@TODO: Cleanup!
		return 5;
	}


	getline(FileMatrixElem1, Line1);
	getline(FileMatrixElem2, Line2);
	getline(FileMatrixElem1, Line1);
	getline(FileMatrixElem2, Line2);
	getline(FileMatrixElem1, Line1);
	getline(FileMatrixElem2, Line2);

	MaxDiff = 0.0;
	// Read B vector parts
	cout << endl << "B vector" << endl;
	for (int i = 0; i < NumShortTerms1+1; i++) {
		FileMatrixElem1 >> Elem1;
		FileMatrixElem2 >> Elem2;
		// Computes the relative difference, not the absolute difference.
		Diff = fabs(Elem1-Elem2) / (fabs(Elem1+Elem2)/2.0);
		cout << Diff << " " << Elem1 << " " << Elem2 << endl;
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
	// Duplicate the first row of the image to the next 31.
	for (int i = 0; i < NumShortTerms1+1; i++) {
		for (int j = 1; j < 32; j++) {
			for (int k = 0; k < 3; k++) {
				ImageData[j*(NumShortTerms1+1)*3 + i*3 + k] = ImageData[i*3 + k];
			}
		}
	}

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

	cout << "B Vector Maximum Difference: " << MaxDiff << endl;

	return 0;
}
