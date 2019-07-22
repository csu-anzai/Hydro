/** 
	Multifractal analysis module utilizing boxcounting to determine 
	the multifractal spectrum of a measure in the input text/binary file.

	Prints the data that forms the multifractal spectrum to a text file
	that can later be analyzed in Matlab or another visualization program.

	The algorithm used is described in:
		A. Chhabra and R. V. Jensen, Phys. Rev. Lett. 62, 1330 (1989).

	@author Samuel Brenner
	@version July 11, 2013
**/
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <iostream>
#include <algorithm>
#include "FlashAmrMesh.h"
#include "Driver_main.h"

using namespace std;

double qMin = -10;	//minimum and maximum values of the moment used in 
double qMax = 10;
double qIncrement = 1; //the step used between each value of q tested.

void printArray(int** arrayIn, int height, int width){
	cout << endl;
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			printf("%7d", arrayIn[i][j]);
		}
		printf("\n");
	 }
}

void printArray(double** arrayIn, int arrayheight, int width){
	cout << endl;
 	for(int i = 0; i < arrayheight; i++){
		for(int j = 0; j < width; j++){
			printf("%3.3f ", arrayIn[i][j]);
		}
		printf("\n");
 	}
}

void printToFile(double** arrayIn, int level, int height, char* fileOutName){
	FILE * pFile;
	char* fileOut = new char[100];
	sprintf(fileOut, "%s%s%d%s", fileOutName, "_falpha_level_", level, ".txt");
	pFile = fopen(fileOut, "w");
	for(int i = 0; i < height; i++){
		fprintf(pFile, "%.8f %.8f", arrayIn[i][0], arrayIn[i][1]);
		fprintf(pFile, "\n");
	}
	fclose(pFile);
	delete[] fileOut;
}

/** 
	Remove zeros from an array and normalize it.
	Implemented in multifracBoxCounter3D.cpp:
		Multifractal analysis module utilizing boxcounting to determine 
		the multifractal spectrum of a measure in the input text/binary file.

		Prints the data that forms the multifractal spectrum to a text file
		that can later be analyzed in Matlab or another visualization program.

		The algorithm used is described in:
			A. Chhabra and R. V. Jensen, Phys. Rev. Lett. 62, 1330 (1989).
**/
void dataCorrection (double***& elements, int height, int width, int depth, bool haveZeros, double arraySum){
	//creates new int to deal with normalizing array when there are zeros.
	double arraySumWithCorrection = arraySum;
	//cout << arraySum << endl;
	//cout << haveZeros;
	if(haveZeros){
		arraySumWithCorrection += height * width * depth;
	}

	//check to see if array is normalized
	//if(((arraySum - 1)/10.0 != 0) || haveZeros){
		for(int i = 0; i < height; i++){
			for(int j = 0; j < width; j++){
				for(int k = 0; k < depth; k++){
					//adds one if there's a zero in the array
					if(haveZeros){
						elements[i][j][k]++;
					}
					elements[i][j][k] /= double(arraySumWithCorrection);
					//cout << elements[i][j][k] << endl;
				}
			}
		}
	//}
}

/**	
	Returns in a double** the data points that make up the plot of f(alpha) vs. alpha.
	The analysis is performed using the method of moments described in Chhabra and Jensen, 
	1989 (Phys. Rev. Lett.).
	@param arrayIn is the 2-/3-D array of doubles that contains the data
	@param height
	@param width
	@param DPETH
	@param level is the level of analysis to be performed. Level 0, for example, 
	examines only the individual data points, whereas at level 1 they are merged
	into boxes of side length 2^1, and for level l size 2^l.
	@return alpha vs. f(alpha)
**/
double** boxCounting(double*** arrayIn, int height, int width, int depth, int level){
	int boxDimension = (int) pow(2, level); //side length of the boxes that we'll use to coarse-grain the data, in pixels.
	int boxDimensionZ;

	if(depth != 1){
		boxDimensionZ = boxDimension;
	}
	else{
		boxDimensionZ = 1;
	}

	int nDataPts = height * width * depth / pow(boxDimension, 2) / boxDimensionZ;
	int falphaSize = qMax - qMin + 1;

	//creates array to hold falpha
	double** falpha = new double*[falphaSize];
	for(int i = 0; i < qMax - qMin + 1; i++){
		falpha[i] = new double[2];
	}

	/**	iterates over all values of q to be tested, where q is the "microscope for exploring different regions of the measure",
		an exaggerating exponent that gives us the qth moment of the measure. We parametrize the relationship f(alpha) vs. alpha
		in terms of q to find functions f(q) and alpha(q). These values are then added to the falpha out-array.
	**/
	for(double q = qMin; q <= qMax; q += qIncrement){
		double* probability = new double[nDataPts]; //an array of all the non-normalized measures taken in.
		double sumOfProbabilitiesQthMoment = 0.0;
		int arrayCounter = 0;
		double fOfQ = 0.0;
		double alphaOfQ = 0.0;

		//iterates through all boxes
		for(int i = 0; i < height; i += boxDimension){
			for(int j = 0; j < width; j += boxDimension){
				for(int k = 0; k < depth; k += boxDimensionZ){
					double boxSum = 0.0;

					//sums each box
					for(int boxSumX = 0; boxSumX < boxDimension; boxSumX++){
						for(int boxSumY = 0; boxSumY < boxDimension; boxSumY++){
							for(int boxSumZ = 0; boxSumZ < boxDimensionZ; boxSumZ++){
								
								//This is to deal with dimensions that aren't powers of two.
									//if(boxSumX + i >= height || boxSumY + j >= width || boxSumZ + k >= depth){ 
									//	cout << endl << "Dimensions must be powers of two!" << endl;				
									//	boxSum += 0;
									//}
									//else{
									boxSum += arrayIn[boxSumX + i][boxSumY + j][boxSumZ + k];
									//}
							}
						}
					}
					probability[arrayCounter] = boxSum;
					arrayCounter++;
				}
			}
		}
	
		for(int i = 0; i < nDataPts; i++){
			sumOfProbabilitiesQthMoment += pow(probability[i], q); //the denominator in eqn. 6 of Chhabra and Jensen
		}

		double log2size = log2(double(boxDimension) / height);

		for(int i = 0; i < nDataPts; i++){
			double normalizedMeasure = pow(probability[i], q) / sumOfProbabilitiesQthMoment; //eqn. 6 of Chhabra and Jensen
			//cout << log2(normalizedMeasure) << endl;
			fOfQ += normalizedMeasure * log2(normalizedMeasure); //the numerator of eqn. 7 of Chhabra and Jensen
			alphaOfQ += normalizedMeasure * log2(probability[i]); //the numerator of eqn. 8 of Chhabra and Jensen
		}

		falpha[int(q - qMin)][0] = alphaOfQ / log2size; //the final divisions in each equation
		falpha[int(q - qMin)][1] = fOfQ / log2size;
		delete[] probability;
	}

	return falpha;
}

void allocateVarArrayMain(double****& inData, int numVars, int *nCellsVec) {

	inData = new double***[numVars];

	for (int i = 0; i < numVars; i++) {
		inData[i] = new double**[nCellsVec[IAXIS]];

		for (int j = 0; j < nCellsVec[IAXIS]; j++) {
			inData[i][j] = new double*[nCellsVec[JAXIS]];

			for (int k = 0; k < nCellsVec[JAXIS]; k++) {
				inData[i][j][k] = new double[nCellsVec[KAXIS]];
			}
		}
	}
}

void deallocateVarArrayMain(double****& inData, int numVars, int *nCellsVec) {

	if (inData != NULL) {
		for (int i = 0; i < numVars; i++) {
			for (int j = 0; j < nCellsVec[IAXIS]; j++) {
				for (int k = 0; k < nCellsVec[JAXIS]; k++) {
					delete[] inData[i][j][k];
				}
				delete[] inData[i][j];
			}
			delete[] inData[i];
		}
		delete[] inData;
	}
}

int LAVAFLOW_DRIVER_BASENM (int argc, char* argv[]) {

    /* NOTE:
         Input file must be contain a uniformly defined mesh in an HDF5 file. Use the 'refineToFinest'
         function (preferably in a separate program) to create this file from an AMR mesh. 
         An example implementation of refineToFinest can be found in 'UniformRefinementTest'.
         Note that ExportUniformData does generally contain everything required to create
         a uniform mesh, but does not allow the user to specify a subdomain or new variables.
         -pb
    */
    char* fileInName = "/home/pb12c/Research/LAVAflow/build/drivers/devel/ExportUniformData/sedovHDF5_lref5.hdf5";
    char* fileOutName = "sedovHDF5_fractal";
    char* var = "dens";

    /* LAVAflow-required functions to read in data */

	/* Initialization */
	string inFileNameStr(fileInName);
	string outFileNameStr(fileOutName);

	vector<string> meshVars;
	meshVars.push_back(var);

	FlashAmrMesh mesh(inFileNameStr, outFileNameStr, meshVars);

	int varIndex = mesh.findVarIndex(var);

	//Retrieve the data from the single block.
	double**** blocksMain;
	allocateVarArrayMain(blocksMain, mesh.getNUserVars(), mesh.getNcells());

	int startingPos[MESH_MDIM];
	startingPos[IAXIS] = 0;
	startingPos[JAXIS] = 0;
	startingPos[KAXIS] = 0;

	mesh.getBlkData(0, CENTER, var, INTERIOR, startingPos, blocksMain);

	int height = mesh.getNcells(IAXIS);
	int width  = mesh.getNcells(JAXIS);
	int depth  = mesh.getNcells(KAXIS);

	bool haveZeros = false;
    double arraySum;

	//Check for zeros and sum array
	for(int k = 0; k < depth; k++){
		for(int i = 0; i < height; i++){
		 	for(int j = 0; j < width; j++){
		 		arraySum += blocksMain[0][i][j][k];
		 		if (blocksMain[0][i][j][k] == 0){
		 			haveZeros = true;
		 		}
		 	}
		}
	}

    /* Former, non-LAVAflow reader */
   	// if (string(asciiOrBinary) == string("ascii")) {
   	// 	elements = dataReaderASCII<double>(fileInName, height, width, depth, haveZeros, arraySum); 
   	// }
   	// else {
   	// 	elements = dataReaderBinary<double>(fileInName, height, width, depth, haveZeros, arraySum); //for reading in binary data
   	// }

    /**	
    	Corrects the data by removing zeros (does so by adding one to every value in the array)
    	and normalizing the array so that the sum of all the elements is one.
    **/
    dataCorrection(blocksMain[0], height, width, depth, haveZeros, arraySum);

	int lowestLevel = 0; 	//the lowest level of coarse-graining, where level == 0 examines each individual pixel as its own box.
							//level = log2(box's side length) so that side length = 1 when level == 0.

	//iterates over all levels to be tested, going from lowestLevel to the highest possible level permitted by the arrayIn size.
	for(int k = 0; k < log2(height) - lowestLevel; k++){
		cout << endl << "Level: " << k + lowestLevel << endl;
		double** falpha = boxCounting(blocksMain[0], height, width, depth, k + lowestLevel);
		printArray(falpha, qMax - qMin + 1, 2);
		printToFile(falpha, k + lowestLevel, qMax - qMin + 1, fileOutName); //This spectrum output can then be analyzed in Matlab or another visualization program.
	}

	deallocateVarArrayMain(blocksMain, mesh.getNUserVars(), mesh.getNcells());
	
	return 0;
}
