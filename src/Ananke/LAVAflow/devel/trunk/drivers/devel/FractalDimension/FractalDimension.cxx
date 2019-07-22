/** Fractal analysis module utilizing boxcounting to determine 
	the fractal dimension of a shape in the input text/binary file.

	Prints to terminal:
		- the results at each level of analysis
		- the overall dimension as determined by least-squares regression
		- an R^2 value
		
	@author Samuel Brenner
	@version July 11, 2013

**/



#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <iostream>
#include "FlashAmrMesh.h"
#include "Driver_main.h"
#include "FlashIndexing.h"
#include <iomanip> 
#include <random>

#define OLD_BOXCOUNT

using namespace std;

void printArray(int**, int, int);
void printArray(double**, double, double);
void printToFile(double**, int, int, double, double, const char*);
void printToFile(int***, int, int, const char*);
double mean(double*, int);
double stdev(double*, double, int);
double corrCoeff(double*, double*, double, double, double, double, int);
void slope(double**, int, double*);
int*** interpolate(double***, int, int, int, double);

#ifdef OLD_BOXCOUNT
double boxCounting(int***, int, int, int, int);
#else
double boxCounting(int***, int, int, int, int);
#endif

void allocateVarArrayMain(double***&, int, int);
void deallocateVarArrayMain(double***&, int, int);
double GetIgnitionTime(double , double , double );

void printArray(int** arrayIn, int height, int width){
	cout << endl;
	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			printf("%3d", arrayIn[i][j]);
		}
		printf("\n");
	 }
}

void printArray(double** arrayIn, double height, double width){
	cout << endl;
	printf("%7s%7s\n", "level", "log(n)");
 	for(int i = 0; i < height; i++){
		for(int j = 0; j < width; j++){
			printf("%7.3f ", arrayIn[i][j]);
		}
		printf("\n");
 	}
}

void printToFile(double** arrayIn, int height, int width, double slope, double yint, const char* outFileName){

	std::cout<<"outFileName: "<<outFileName<<std::endl;

	FILE * pFile;
	pFile = fopen(outFileName, "w");
	fprintf(pFile, "slope: %f\nintercept: %f\n", slope, yint);
	for(int i = height - 1; i >= 0; i--){
		for(int j = 0; j < width; j++){
			fprintf(pFile, "%-2.6f ", arrayIn[i][j]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}

void printToFile(int*** arrayIn, int height, int width, const char* outFileName){
	FILE * pFile;
	pFile = fopen(outFileName, "w");
	for(int i = height - 1; i >= 0; i--){
		for(int j = 0; j < width; j++){
			fprintf(pFile, "%d ", arrayIn[i][j][0]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}

/**
	slope, stdev, corrCoeff, and mean: Module to run regression for boxcounting.
	
	Changes a double[3] such that regressionArray[0] = slope,
	regressionArray[1] = R^2, and regressionArray[2] = y-int.

	Implemented in boxCounter3D.cpp:
		Fractal analysis module utilizing boxcounting to determine 
		the fractal dimension of a shape in the input text/binary file.

		Prints to terminal:
			- the results at each level of analysis
			- the overall dimension as determined by least-squares regression
			- an R^2 value
			
		@author Samuel Brenner
		@version July 11, 2013

	@author Samuel Brenner
	@version July 11, 2013

**/

double mean(double* arrayIn, int arrayLength){
	double sum = 0;

	for(int i = 0; i < arrayLength; i++){
		sum += arrayIn[i];
	}

	return sum / arrayLength;
}

double stdev(double* arrayIn, double arrayMean, int arrayLength){

	double squaresSum = 0;

	for(int i = 0; i < arrayLength; i++){
		squaresSum += pow(arrayIn[i] - arrayMean, 2);
	}

	return sqrt(squaresSum / (arrayLength - 1));
}

double corrCoeff(double* arrayX, double* arrayY, double meanX, 
	double meanY, double stdevX, double stdevY, int length){

	double sum = 0;

	for(int i = 0; i < length; i++){
		sum += (arrayX[i] - meanX) * (arrayY[i] - meanY);
	}

	return sum / (stdevX * stdevY * (length - 1));
}

void slope(double** arrayIn, int arrayLength, double* regressionArray){
	
	double* arrayX = new double[3];

	for(int i = 0; i < 3; i++){
		arrayX[i] = arrayIn[i][0];
	}

	double* arrayY = new double[3];
	for(int i = 0; i < 3; i++){
		arrayY[i] = arrayIn[i][1];
	}

	double meanX = mean(arrayX, 3);
	double meanY = mean(arrayY, 3);

	double stdevX = stdev(arrayX, meanX, 3);
	double stdevY = stdev(arrayY, meanY, 3);

	double r = corrCoeff(arrayX, arrayY, meanX, meanY, stdevX, stdevY, 3);

	double slope = r 
		* stdevY / stdevX;


	delete[] arrayX;
 	delete[] arrayY;

 	regressionArray[0] = slope;
 	regressionArray[1] = pow(r, 2);
 	regressionArray[2] = meanY - slope * meanX;
}

/**
	Interpolate the isocontour at some pre-specified value inside a 2- or 3-D array.

	Implemented in boxCounter3D.cpp:
		Fractal analysis module utilizing boxcounting to determine 
		the fractal dimension of a shape in the input text/binary file.

		Prints to terminal:
			- the results at each level of analysis
			- the overall dimension as determined by least-squares regression
			- an R^2 value
			
		@author Samuel Brenner
		@version July 11, 2013

	@author Samuel Brenner
	@version July 11, 2013

**/
int*** interpolate(double*** arrayIn, int height, int width, int depth, double valueToFind){
	//Declare and initialize arrayOut. If any values in the arrayIn happen to be the valueToFind,
	//we'll initialize the corresponding arrayOut to a 1; otherwise it'll be a 0.
	int*** arrayOut = new int**[height];
	for(int i = 0; i < height; i++){
		arrayOut[i] = new int*[width];
		for(int j = 0; j < width; j++){
			arrayOut[i][j] = new int[depth];
			for(int k = 0; k < depth; k++){
				//performs a preliminary check so that we don't miss any values that are exactly == valueToFind
				if(arrayIn[i][j][k] == valueToFind){
					arrayOut[i][j][k] = 1;
				}
				else{
					arrayOut[i][j][k] = 0;
				}
			}
		}
	}

	//Makes the interpolater able to handle planar data
	int startingDepthIndex = 0;
	if(depth != 1){
		startingDepthIndex = 1;
	}
	else{
		depth++;
	}

	/*
		Iterates over all entries in the array that aren't on the edge of the array, unless the array is 2D: for 2D
		arrays the program ignores the edge condition in the third dimension.

		For each entry analyzed, we check to see if the value is lower than valueToFind. If so, the neighbors on all 
		six sides (four if planar data) are checked against it to determine if there is some crossing over valueToFind
		in between the two. Then we determine which cell should contain that crossing and assign it a 1 in the 
		corresponding cell of the outArray.
	*/
	for(int i = 1; i < height - 1; i++){
		for(int j = 1; j < width -1; j++){
			for(int k = startingDepthIndex; k < depth - 1; k++){
				if(arrayIn[i][j][k] < valueToFind){
					if(arrayIn[i + 1][j][k] > valueToFind){
						if(int((valueToFind - arrayIn[i][j][k]) / (arrayIn[i + 1][j][k] - arrayIn[i][j][k])) == 0){
							arrayOut[i][j][k] = 1;
						}
						else{
							arrayOut[i + 1][j][k] = 1;
						}
					}
					if(arrayIn[i][j + 1][k] > valueToFind){
						if(int((valueToFind - arrayIn[i][j][k]) / (arrayIn[i][j + 1][k] - arrayIn[i][j][k])) == 0){
							arrayOut[i][j][k] = 1;
						}
						else{
							arrayOut[i][j + 1][k] = 1;
						}
					}
					if(arrayIn[i][j - 1][k] > valueToFind){
						if(int((valueToFind - arrayIn[i][j][k]) / (arrayIn[i][j - 1][k] - arrayIn[i][j][k])) == 0){
							arrayOut[i][j][k] = 1;
						}
						else{
							arrayOut[i][j - 1][k] = 1;
						}
					}
					if(arrayIn[i - 1][j][k] > valueToFind){
						if(int((valueToFind - arrayIn[i][j][k]) / (arrayIn[i - 1][j][k] - arrayIn[i][j][k])) == 0){
							arrayOut[i][j][k] = 1;
						}
						else{
							arrayOut[i - 1][j][k] = 1;
						}
					}
					if(arrayIn[i][j][k + 1] > valueToFind){
						if(int((valueToFind - arrayIn[i][j][k]) / (arrayIn[i][j][k + 1] - arrayIn[i][j][k])) == 0){
							arrayOut[i][j][k] = 1;
						}
						else{
							arrayOut[i][j][k + 1] = 1;
						}
					}
					if(arrayIn[i][j][k - 1] > valueToFind){
						if(int((valueToFind - arrayIn[i][j][k]) / (arrayIn[i][j][k - 1] - arrayIn[i][j][k])) == 0){
							arrayOut[i][j][k] = 1;
						}
						else{
							arrayOut[i][j][k - 1] = 1;
						}
					}

				}
			}
		}
	}

	return arrayOut;
}

/**
	Returns the log base 2 of the number of cells filled in a given level
	of fractal analysis.
	@param arrayIn is the 2-3-D array of integers that contains the data
	@param height
	@param width
	@param DPETH
	@param level is the level of analysis to be performed. Level 0, for example, 
	examines only the individual data points, whereas at level 1 they are merged
	into boxes of side length 2^1, and for level l size 2^l.
	@return log base 2 of the number of cells filled in a given level.

**/
#ifdef OLD_BOXCOUNT
double boxCounting(int*** arrayIn, int height, int width, int depth, int level){
	int numberFilled = 0; //number of boxes with an element of the figure
	int boxDimension = (int) pow(2, level);
	int boxDimensionZ;

	if(depth != 1){
		boxDimensionZ = boxDimension;
	}
	else{
		boxDimensionZ = 1;
	}
	//iterates through all boxes
	//cout << "here!" << endl;
	for(int i = 0; i < height; i += boxDimension){
		for(int j = 0; j < width; j += boxDimension){
			for(int k = 0; k < depth; k += boxDimensionZ){
				int boxSum = 0;


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

				if(boxSum > 0){
					numberFilled++;
				}

			}
		}
	}
	printf("%s%d\n%s%d\n%s%d\n\n", "level: ", level, 
	"--box dimension: ", boxDimension, "--filled boxes: ", numberFilled);		

	return log2(numberFilled);
}
#else

double boxCounting(int*** arrayIn, int height, int width, int depth, int level){
	const double tol = 1e-3;

	int numberFilled = 0; //number of boxes with an element of the figure
	int boxDimension = (int) pow(2, level);
	int boxDimensionZ;

	int boxVol = pow(boxDimension, 3);

	cout << "Box dimension: " << boxDimension << endl;

	if(depth != 1){
		boxDimensionZ = boxDimension;
	}
	else{
		boxDimensionZ = 1;
	}

	int domainDims[] = {height, width, depth};
	double numBoxesPerDim[] = {double(height) / boxDimension, double(width) / boxDimension, double(depth) / boxDimension};

	std::random_device rd;
    std::mt19937 gen(rd());
    
    std::uniform_int_distribution<> disX(0, domainDims[IAXIS]-1);
    std::uniform_int_distribution<> disY(0, domainDims[JAXIS]-1);
    std::uniform_int_distribution<> disZ(0, domainDims[KAXIS]-1);

    int sampleCoeff = 200;
    double lengthScaleFactor = pow(numBoxesPerDim[IAXIS], 3);
    int numSamples = sampleCoeff * lengthScaleFactor;
    cout << "Num samples: " << numSamples << endl;


	// std::vector<int> boxSums;
	//iterates through all boxes
	int numOverlaps = 0;

	int numOverlapIts = pow(2,numOverlaps);
	for (int overlapIndex=0; overlapIndex<numOverlapIts; overlapIndex++)
	{
		double offsetBoxFrac = 1/double(numOverlapIts)*overlapIndex;
		int offset =  boxDimension * offsetBoxFrac;

		// cout << "Overlap: " << overlapIndex << '\n';
		// cout << "Current offset is " << offset << endl;

		for(int boxIndexX = 0; boxIndexX < numBoxesPerDim[IAXIS]; boxIndexX++)
		{
			for(int boxIndexY = 0; boxIndexY < numBoxesPerDim[JAXIS]; boxIndexY++)
			{
				for(int boxIndexZ = 0; boxIndexZ < numBoxesPerDim[KAXIS]; boxIndexZ++)
				{
				// for (int sampleNum=0; sampleNum<numSamples; sampleNum++)
				// {	

				    // int boxIndexX = disX(gen),
				    // 	boxIndexY = disY(gen),
				    // 	boxIndexZ = disZ(gen);


					int boxSum = 0;

					int lowerBoxBounds[] = {boxIndexX*boxDimension+offset,
											boxIndexY*boxDimension+offset,
											boxIndexZ*boxDimension+offset};

					// int lowerBoxBounds[] = {boxIndexX,
					// 						boxIndexY,
					// 						boxIndexZ};

					// if (boxDimension > 30)
					// {
					// 	cout << "Overlap: " << overlapIndex << '\n';
					// 	cout << "Box lower bounds (x,y,z): " << lowerBoxBounds[IAXIS] << ", " << lowerBoxBounds[JAXIS] << ", " << lowerBoxBounds[KAXIS] << endl;
					// }


					//sums each box
					for(int boxX = 0; boxX < boxDimension; boxX++)
					{
						// if (boxDimension == 32 && boxIndexX == 1 && boxIndexY == 1 && boxIndexZ == 1)
						// {
						// 	int tempX = lowerBoxBounds[IAXIS] + boxX;

						// 	if (tempX >= domainDims[IAXIS])
						// 				tempX -= domainDims[IAXIS];
							
						// 	cout << "x coord: " << tempX << '\n';
						// 	cout << arrayIn[tempX][50][50];
						// 	cout << endl;
						// }
						for(int boxY = 0; boxY < boxDimension; boxY++)
						{
							for(int boxZ = 0; boxZ < boxDimensionZ; boxZ++)
							{
								// int gridX = boxIndexX*boxDimension+offset + boxX,
								// 	gridY = boxIndexY*boxDimension+offset + boxY,
								// 	gridZ = boxIndexZ*boxDimension+offset + boxZ;


								int gridCoord[] = {lowerBoxBounds[IAXIS] + boxX,
												   lowerBoxBounds[JAXIS] + boxY,
												   lowerBoxBounds[KAXIS] + boxZ};

								for (int dimCoord=0; dimCoord<MESH_MDIM; dimCoord++)
								{
									if (gridCoord[dimCoord] >= domainDims[dimCoord])
										gridCoord[dimCoord] -= domainDims[dimCoord];
								}
								
								boxSum += arrayIn[gridCoord[IAXIS]][gridCoord[JAXIS]][gridCoord[KAXIS]];								
							}
						}
					}



					if(boxSum > 0){
						// weight by how much of the box is filled
						// double volWeight = (double) boxVol / boxSum;
						// numberFilled+=volWeight;
						numberFilled++;
					}

					// if (boxDimension == 32)
					// {
					// 	if (boxSum == 0)
					// 		cout << "Empty box indices: " << boxIndexX << ", " << boxIndexY << ", " << boxIndexZ << endl;
					// 	else
					// 		cout << "Filled box indices: " << boxIndexX << ", " << boxIndexY << ", " << boxIndexZ << endl;
					// }

					// boxSums.push_back(boxSum);

					// boxIndex++;
				// }

				}
			}
		}
		cout << "Total number filled for this overlap: " << numberFilled << endl;
		cout << "Scaled number filled for this overlap: " << numberFilled * lengthScaleFactor << endl;
	}

	// for (auto sumThisBox : boxSums)

	// cout << "Average per box: " << double(totalElements) / numberFilled << endl;

	cout << "--box dimension: " << boxDimension << '\n' <<  
			"--filled boxes: " << numberFilled << '\n' << endl;

	// printf("%s%d\n%s%d\n%s%d\n\n", "level: ", level, 
	// "--box dimension: ", boxDimension, "--filled boxes: ", numberFilled);		

	return log2(numberFilled/double(numOverlapIts));
	// return log10(numberFilled) / log10(epsilon);
	// return numberFilled;
}

#endif


void allocateVarArrayMain(double***& inData, int *nCellsVec) 
{

	inData = new double**[nCellsVec[IAXIS]];
	
	for (int i = 0; i < nCellsVec[IAXIS]; i++) 
	{
		inData[i] = new double*[nCellsVec[JAXIS]];

		for (int j = 0; j < nCellsVec[JAXIS]; j++) 
		{
			inData[i][j] = new double[nCellsVec[KAXIS]];
		}
	}
	
}


void deallocateVarArrayMain(double***& inData, int *nCellsVec) 
{

	if (inData != NULL) 
	{
		for (int j = 0; j < nCellsVec[IAXIS]; j++) 
		{
			for (int k = 0; k < nCellsVec[JAXIS]; k++) 
			{
				delete[] inData[j][k];
			}
			delete[] inData[j];
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
    // char* fileInName = "/home/pb12c/Research/LAVAflow/build/drivers/devel/UniformRefinementTest/sedovHDF5_lref5";
    // char* fileOutName = "sedovHDF5_fractal";
    // char* var = "dens";

	double contourVal = 1.e-5;

	/* Initialization */
	// string fileInName = "/data1/df11c/data/data/wdStir/512ppm/512ppmi7_75comp";
	// string fileInName = "/data1/df11c/data/data/sedovHDF5_lref5";
	// string fileInName = "/data1/df11c/Simulations/WDStir/ClusterTest/tburn_hdf5_chk_0000_128";
	string fileInName = "/data1/df11c/data/data/wdStir/512ppm/512ppmi7_75comp_t=0.75";
	string fileOutName = "fracDimComp";

	vector<string> meshVars = {"igtm"};
	// vector<string> meshVars = {"dens", "temp", "c12 "};
	// meshVars.push_back(var);

	FlashAmrMesh mesh(fileInName, meshVars);

	// int densIndex = mesh.findVarIndex("dens");
	// int tempIndex = mesh.findVarIndex("temp");
	// int c12Index = mesh.findVarIndex("c12 ");

	//Retrieve the data from the single block.
	// double**** blocksMain;
	// allocateVarArrayMain(blocksMain, mesh.getNUserVars(), mesh.getNcells());

	// int startingPos[MESH_MDIM];
	// startingPos[IAXIS] = 0;
	// startingPos[JAXIS] = 0;
	// startingPos[KAXIS] = 0;

	// mesh.getBlkData(0, CENTER, meshVars[0].c_str(), INTERIOR, startingPos, blocksMain);

	// int height = mesh.getNcells(IAXIS);
	// int width  = mesh.getNcells(JAXIS);
	// int depth  = mesh.getNcells(KAXIS);

	// bool haveZeros = false;
    // double arraySum;

	// //Check for zeros and sum array
	// for(int k = 0; k < depth; k++){
	// 	for(int i = 0; i < height; i++){
	// 	 	for(int j = 0; j < width; j++){
	// 	 		arraySum += blocksMain[0][i][j][k];
	// 	 		if (blocksMain[i][j][k] == 0){
	// 	 			haveZeros = true;
	// 	 		}
	// 	 	}
	// 	}
	// }

	

  //  	if (string(asciiOrBinary) == string("ascii")) {
		// elements = dataReaderASCII<double>(inFileName, height, width, depth, haveZeros, arraySum);
  //  	}
  //  	else {
  //  		elements = dataReaderBinary<double>(inFileName, height, width, depth, haveZeros, arraySum); //for reading in binary data
  //  	}

	// elementsInterpolated = interpolate(blocksMain[0], height, width, depth, contourVal);
	//printToFile(elementsInterpolated, height, width);


	// new method. works for multiple blocks
	// works as long as the domain is cubic
	double*** uniformMesh;
	int numBlocks = mesh.getNBlocks();

	int nCellsBlockX = mesh.getInt("nxb"),
        nCellsBlockY = mesh.getInt("nyb"),
        nCellsBlockZ = mesh.getInt("nzb");

    int nBlocks[MESH_MDIM] = {mesh.getInt("nblockx"), mesh.getInt("nblocky"), mesh.getInt("nblockz")};
    
	int *nCellsMesh = mesh.getNcells();
	int nCellsUniformMesh[MESH_MDIM];

	for (int i=0; i<MESH_MDIM; i++)
		nCellsUniformMesh[i] = nCellsMesh[i] * nBlocks[i];


	allocateVarArrayMain(uniformMesh, nCellsUniformMesh);

	auto blockLinearIndexMap = BuildBlockLinearIndexMap(mesh);

	
    int numVars = meshVars.size();

    int height = nCellsUniformMesh[IAXIS];
	int width  = nCellsUniformMesh[JAXIS];
	int depth  = nCellsUniformMesh[KAXIS];

	
 //    for (int var=0; var<numVars; var++)
	// {
	for (int i=0; i<height; i++)
	{
		for (int j=0; j<width; j++)
		{
			for (int k=0; k<depth; k++)
			{
				// get the row-major linear block index, assuming an equal number of cells per dimension in the block, and a cubic domain
				int linBlockIndex = GetLinearBlockIndex(i, j, k, nCellsBlockX, nBlocks[0]);
    			int blockID = blockLinearIndexMap[linBlockIndex];

    			int blockI = i % nCellsBlockX,
                	blockJ = j % nCellsBlockY,
                	blockK = k % nCellsBlockZ;

    			//Get pointer to the block's data
    			double ****solnData;
                mesh.getBlkPtr(blockID, solnData, CENTER);

                // double  densThisCell = solnData[densIndex][blockI][blockJ][blockK],
                // 		tempThisCell = solnData[tempIndex][blockI][blockJ][blockK],
                // 		c12ThisCell = solnData[c12Index][blockI][blockJ][blockK];

                // uniformMesh[i][j][k] = GetIgnitionTime(tempThisCell, densThisCell, c12ThisCell);

                uniformMesh[i][j][k] = solnData[0][blockI][blockJ][blockK];

                // // compare to original array
                // if (uniformMesh[var][i][j][k] != blocksMain[var][i][j][k])
                // {
                // 	std::cerr << "Mapped array values do not match at index i, j, k = " << i << ", " << j << ", " << k << '\n';
                // 	std::cerr << "Value in original, new arrays: " << blocksMain[var][i][j][k] << ", " << uniformMesh[var][i][j][k] << std::endl;
                // }
			}		
		}
	}
	// }


	int*** elementsInterpolated;
	elementsInterpolated = interpolate(uniformMesh, height, width, depth, contourVal);

	// int kSlice = 32;
	// for (int i=0; i<height; i++)
	// {
	// 	for (int j=0; j<width; j++)
	// 	{
	// 		// int outVal;
	// 		// if (uniformMesh[i][j][kSlice] < 6e8)
	// 		// 	outVal = 0;
	// 		// else
	// 		// 	outVal = 1;
	// 		// cout << outVal << " ";

	// 		cout << elementsInterpolated[i][j][kSlice] << " ";
	// 	}
	// 	cout << '\n';
	// }
	// cout << endl;


	std::vector<double> filledBoxes;
	int lowestLevel = 0;

	//initialize out-array
	int largestDimension = fmin(height, width);
	if(depth > 1){
		largestDimension = fmin(largestDimension, depth);
	}


	int outArrayLength = int(log2(largestDimension) - lowestLevel + 1);
	int levelDivisor = 2;
	// int outArrayLength = int(log2(largestDimension) - lowestLevel + 1)*levelDivisor;
	double** outDataArray = new double* [outArrayLength];

	//printf("\n\n%d\n\n", outArrayLength);

	for(int i = 0; i < outArrayLength; i++){
		outDataArray[i] = new double[2];
	}

	

	for(int k = 0; k < outArrayLength; k++)
	{
		outDataArray[k][0] = outArrayLength - 1 - (k + lowestLevel);
		outDataArray[k][1] = boxCounting(elementsInterpolated, height, width, depth, k + lowestLevel);

		filledBoxes.push_back(pow(2,outDataArray[k][1]));
	}

	// // // dfenn test ///////////////////////////////////////////////////
	// // // 
	// std::random_device rd;
 //    std::mt19937 gen(rd());
    
 //    std::uniform_int_distribution<> disLengthScales(1, largestDimension);
	

	// // const int numLengthScales = largestDimension;
	
	// int beginBoxLength = 4;
	// int endBoxLength = largestDimension;
	// int step = 2;
	// int numLengthScales = (endBoxLength - beginBoxLength + 1) / step;
	// std::vector<int> lengthScales(numLengthScales, 0);

	// // for (auto& thisLengthScale : lengthScales)
	// // {
	// // 	thisLengthScale = disLengthScales(gen);
	// // }

	// // std::sort(lengthScales.begin(), lengthScales.end());
	// // auto last = std::unique(lengthScales.begin(), lengthScales.end());
 // //    lengthScales.erase(last, lengthScales.end());

    
 //    int i=beginBoxLength;
 //    for (auto& thisLengthScale : lengthScales)
 //    {
 //    	thisLengthScale = i;
 //    	i+=step;
 //    }



	// for (auto thisLengthScale : lengthScales)
	// {
	// 	int numBoxesFilled = boxCounting(elementsInterpolated, height, width, depth, thisLengthScale);
	// 	filledBoxes.push_back(numBoxesFilled);
	// }

	// int outArrayLength = lengthScales.size()-1;
	// int maxLengthScale = lengthScales.back();

	// double** outDataArray = new double* [outArrayLength];

	// //printf("\n\n%d\n\n", outArrayLength);

	// for(int i = 0; i < outArrayLength; i++)
	// {
	// 	outDataArray[i] = new double[2];
	// }

	// double cumFracDim = 0;
	// for (int i=0; i<filledBoxes.size()-1; i++)
	// {
	// 	// double avgScaleChange = (double(lengthScales[i+1])/lengthScales[i] + double(lengthScales[i])/lengthScales[i-1])/2.;
	// 	double scaleChange = double(lengthScales[i+1])/lengthScales[i];
	// 	cout << "Box length: " << lengthScales[i] << endl;
	// 	cout << "Num boxes filled: " << filledBoxes[i] << endl;
	// 	cout << "Scale change: " << scaleChange << '\n';
	// 	// double logBoxChange = log10((double(filledBoxes[i])/filledBoxes[i+1] + double(filledBoxes[i-1])/filledBoxes[i])/2) / log10(scaleChange);
	// 	double logBoxChange = log10(double(filledBoxes[i])/filledBoxes[i+1]) / log10(scaleChange);
	// 	cumFracDim += logBoxChange;
	// 	cout << "Log of box count change: " << logBoxChange << endl << endl;


	// 	outDataArray[i][0] = log10(lengthScales[i]);
	// 	// outDataArray[i][0] = i;
	// 	outDataArray[i][1] = log10(filledBoxes[i]);
	// }

	

	//////////////////////////////////////////////////////////////////////


	// for(int k = 0; k < outArrayLength - levelDivisor+1; k++)
	// {
	// 	cout << "k: " << k << '\n';
	// 	double level = k/double(levelDivisor) + lowestLevel;
	// 	outDataArray[k][0] = outArrayLength/levelDivisor - 1 - level;
	// 	outDataArray[k][1] = boxCounting(elementsInterpolated, height, width, depth, level);		
	// }

	// get the average fractal dimension
	double cumFracDim = 0;
	for (int i=0; i<filledBoxes.size()-1; i++)
	{
		double logBoxChange = log2(double(filledBoxes[i])/filledBoxes[i+1]);
		cumFracDim += logBoxChange;
		cout << "Log of box count change: " << logBoxChange << endl << endl;

		
		// outDataArray[i][0] = log2();
		// outDataArray[i][1] = logBoxChange;
	}

	cout << "Average fractal dimension: " << cumFracDim / (filledBoxes.size()-1) << endl;

	printArray(outDataArray, outArrayLength, 2);

	double regressionArray[3];

	slope(outDataArray, outArrayLength, regressionArray);

	printf("\n\n%s%s\n\n%s%.5f\n%s%1.3f\n\n", "Curve: ", fileInName.c_str(), "Dimension: ", regressionArray[0], "R^2: ", regressionArray[1]);

	printToFile(outDataArray, outArrayLength, 2, regressionArray[0], regressionArray[2], fileOutName.c_str());
	char interpolatedOutName[50];
	sprintf(interpolatedOutName, "%s_interpolated.txt", fileOutName.c_str());
	printToFile(elementsInterpolated, height, width, interpolatedOutName);

	deallocateVarArrayMain(uniformMesh, nCellsUniformMesh);

	delete[] elementsInterpolated;
	delete[] outDataArray;
	
	return 0;
}

double GetIgnitionTime(double temp, double dens, double c12Frac)
{
    double T9 = temp / 1e9;
    double rho8 = dens / 1e8;

    double f = pow(T9 - 0.214, -7.566);
    double igTime = 1.15e-5*pow(c12Frac*rho8,-1.9) * f * (1+1193*f);

    return igTime;
}
