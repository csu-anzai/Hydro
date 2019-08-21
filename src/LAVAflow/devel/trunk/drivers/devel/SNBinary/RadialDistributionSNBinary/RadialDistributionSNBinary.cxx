
#include "CommonSNBinary.h"
using namespace std;

enum MASSSCALAR {COMPANION, EJECTA, BOTH};

/* Parameters */
#define ORIGINCHOICE HIGHESTDENSITY //USERDEF
#define ORIGINX 0 //Used if ORIGINCHOICE is USERDEF (offset a little)
#define ORIGINY 1e-12
#define ORIGINZ 0

// These are values between the objects at the initial time step.
// #define ORIGINX 0 //MS7 L200
// #define ORIGINY 2.72e11
// #define ORIGINZ 0
// #define ORIGINX 0 //MS38 L200
// #define ORIGINY 1.25e11
// #define ORIGINZ 0
// #define ORIGINX 0 //MS54
// #define ORIGINY 9e10
// #define ORIGINZ 0
// #define ORIGINX 0 //MS63
// #define ORIGINY 1.65e11
// #define ORIGINZ 0
// #define ORIGINX 0 //SG
// #define ORIGINY 3.40e11
// #define ORIGINZ 0
// #define ORIGINX 0 //SY319
// #define ORIGINY 8.6e12
// #define ORIGINZ 0
// #define ORIGINX 0 //SY428
// #define ORIGINY 3.5e13
// #define ORIGINZ 0

#define DTHETA  1.f
#define THETAMIN -90.f
#define THETAMAX 90.f

#define MAXREFINE 8

struct sortContainer{
  double t;
  double rho;
  double dr;
  double vol;
  double totEnergy;
  double compMassScalar;
  double ejectaMassScalar;
};

void findOrigin(double* originVec, double***** blocksMain, FlashAmrMesh &grid, std::set<int> leafBlockIDs, 
			    int numberOfBlocks, int* nCellsArr, int varIndex, ORIGINTYPE originType);
bool findIntersection(double* originVec, double* directionVec, double cellBounds[2][MESH_MDIM], 
					  double lowerBoundT, double upperBoundT, double &tMinReturn, double &tMaxReturn);
double findCellVectorDistance(double* originVec, double* directionVec, double tmin, double tmax);
bool pointInBox(double* point, double cellBounds[2][MESH_MDIM]);
void printOut(char* fileName, double* radiusVals, double* densVals, int radiusCount);
bool sortFunction(sortContainer i,sortContainer j);

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv) {

	modelContainer models[NUMMODELS];
	readModels(models);

	/* Declare variables to be read */
	int nVars = 2;
	char** vars;
	int* varIndex;
	vars = new char*[nVars];
	for (int i = 0; i < nVars; i++) {
		vars[i] = new char[MESH_VAR_STRING_SIZE+1];
	}
	varIndex = new int[nVars];
	char inFileName[MAX_FILENAME_LENGTH];
	char outFileName[MAX_FILENAME_LENGTH];
	char inProgramFileName[MAX_FILENAME_LENGTH];
	char outProgramFileName[MAX_FILENAME_LENGTH];
	char outFileExtraName[MAX_FILENAME_LENGTH];
	
	/* --- Specify variables --- */

	/* TODO: Check the method used to find the center of the companion and make sure it's good.
	   Also check out what's going on with the symmetry axis.
	*/
	MODEL currentModel = MS38_L200;
	// MODEL currentModel = MS54;
	// MODEL currentModel = MS63;
	// MODEL currentModel = MS7_L200;
	// MODEL currentModel = SG_L200;
	// MODEL currentModel = SY319_L200;
	// MODEL currentModel = SY428;

	int checkNum = MS38_L200_I;
	// int checkNum = MS54_I;
	// int checkNum = MS63_I;
	// int checkNum = MS7_L200_I;
	// int checkNum = SG_L200_I;
	// int checkNum = SY319_L200_I;
	// int checkNum = SY428_I;

	MASSSCALAR massScalar = COMPANION;

	strcpy(vars[0],"dens");
	strcpy(vars[1],"ms_3");
	strcpy(outFileExtraName, "_radialDistribution.txt");

	/* Useful Variable names:
		dens = Density
		gpot = Gravitational Potential
		ener = Kinetic Energy
		ms_1 = Ambient Medium Mass Scalar
		ms_2 = Ejecta Mass Scalar
		ms_3 = Companion Mass Scalar
		ms_4 = mass...Mass Scalar?
	*/	
	
	strcpy(inFileName, models[currentModel].checkpoint[checkNum]);
	strcpy(outFileName, models[currentModel].checkpoint[checkNum]);
	

	strcat(outFileName, outFileExtraName);

	strcpy(inProgramFileName, inFolderLocation);
	strcat(inProgramFileName, models[currentModel].modelFolder);
	
	strcpy(outProgramFileName, outFolderLocation);
	strcat(outProgramFileName, models[currentModel].modelFolder);
	strcat(outProgramFileName, "radialDistribution/");

	strcat(inProgramFileName,inFileName);
	strcat(outProgramFileName,outFileName);

	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "WARNING: THIS ANALYSIS PROGRAM NOT YET BUILT FOR MPI" << std::endl;
	std::cout << "Input file: " << inProgramFileName << std::endl;
	std::cout << "Output file: " << outProgramFileName << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	/* ------------------------- */

	string inFileNameStr(inProgramFileName);
	string outFileNameStr(outProgramFileName);

	/* Initialization */
	vector<string> meshVars;
	for (int i = 0; i < nVars; i++) {
		meshVars.push_back(vars[i]);
	}

    FlashAmrMesh mesh(inFileNameStr, outFileNameStr, meshVars);

	for (int i = 0; i < nVars; i++) {
		varIndex[i] = mesh.findVarIndex(vars[i]);
	}

	std::cout << std::endl << "--Entering Main Functions--" << std::endl << std::endl;
	int numberOfBlocks = mesh.getNBlocks();

	double***** blocksMain;

	//Allocate data arrays for each block.
	blocksMain = new double****[numberOfBlocks];

	for (int i = 0; i < numberOfBlocks; i++) {
		allocateVarArray(blocksMain[i], mesh.getNUserVars(), mesh.getNcells());
	}

	//Get leaf block IDs
	std::set<int> leafBlockIDs;

	for (int i = 0; i < numberOfBlocks; i++) {
		if (mesh.getNodeType(i) == FR_LEAF_NODE) {
			leafBlockIDs.insert(i);
		}
		// if (mesh.getRefineLevel(i) == MAXREFINE) {
		// 	leafBlockIDs.insert(i);
		// }
	}

	//Retrieve the data from each block.
	int startingPos[MESH_MDIM];
	startingPos[IAXIS] = 0;
	startingPos[JAXIS] = 0;
	startingPos[KAXIS] = 0;

	for (int i = 0; i < numberOfBlocks; i++) {
		for (int j = 0; j < nVars; j++) {
			mesh.getBlkData(i, CENTER, vars[j], INTERIOR, startingPos, blocksMain[i]);
		}
	}

	// ---------------
	// MAIN SETUP
	// ---------------

	double currentTheta = 0;
	double originVec[3];
	double directionVec[3];

	int nCellsX = mesh.getNcells(IAXIS);
	int nCellsY = mesh.getNcells(JAXIS);
	// int nCellsZ = mesh.getNcells(KAXIS);
	int* nCellsArr = mesh.getNcells();


	//Determine the origin.
	ORIGINTYPE originType = ORIGINCHOICE;
	findOrigin(originVec, blocksMain, mesh, 
	 		   leafBlockIDs, numberOfBlocks, nCellsArr, varIndex[0], originType);
	// originVec[IAXIS] = 0;
	// originVec[JAXIS] = -10.00e12;
	// originVec[KAXIS] = 0;

	//Find the maximum range of the ray.
	double domainBounds[2][MESH_MDIM];
	mesh.getDomainBoundBox(domainBounds);

	//Make sure that the maximum radius of the rays is only as far as the shortest ray distance,
	// which will be one of the extrema of the Y-axis.
	double maxRadius;
	if (fabs(domainBounds[UPPER][JAXIS]-originVec[JAXIS]) < 
		fabs(domainBounds[LOWER][JAXIS]-originVec[JAXIS])) {
		maxRadius = fabs(domainBounds[UPPER][JAXIS] - originVec[JAXIS]);
	}
	else {
		maxRadius = fabs(domainBounds[LOWER][JAXIS] - originVec[JAXIS]);
	}

	// ---------
	// MAIN
	// ---------

	std::vector<sortContainer> sortVector;
	sortContainer tempContainer;
	
	//Convert from polar to cartesian.
	directionVec[IAXIS] = maxRadius * cos(currentTheta * PI / 180.f)/* - originVec[IAXIS]*/;
	directionVec[JAXIS] = maxRadius * sin(currentTheta * PI / 180.f)/* - originVec[JAXIS]*/;
	//std::cout << currentTheta << "X: " << directionVec[IAXIS] << " Y: " << directionVec[JAXIS] << std::endl;
	directionVec[KAXIS] = 0;

	std::set<int> intersectedBlkIDs;
	std::set<double> intersectedLowerProp;
	std::set<double> intersectedUpperProp;

	std::set<int>::iterator leaf;
	std::set<int>::iterator interBlk;
	std::set<double>::iterator lowerPropIter;
	std::set<double>::iterator upperPropIter;

	double domainMinT = 0.f;
	double domainMaxT = 1.f;

	//Loop over each leaf block and see if it intersects the line vector.
	for (leaf = leafBlockIDs.begin(); leaf != leafBlockIDs.end(); ++leaf) {	
		//Get physical location of corners of box in grid.
		double blockBounds[2][MESH_MDIM];
		mesh.getBlkBoundBox(*leaf, blockBounds);
		
		double blockMinT;
		double blockMaxT;
		//Find if there's an intersection. If there is, update the propagation value for this block.
		bool foundIntersection = findIntersection(originVec, directionVec, blockBounds, 
												  domainMinT, domainMaxT,
												  blockMinT, blockMaxT);

		if (foundIntersection) {
			intersectedBlkIDs.insert(*leaf);
			intersectedLowerProp.insert(blockMinT);
			intersectedUpperProp.insert(blockMaxT);
		}

	}
	//Loop over cells in each block.
	lowerPropIter = intersectedLowerProp.begin();
	upperPropIter = intersectedUpperProp.begin();

	for (interBlk = intersectedBlkIDs.begin(); interBlk != intersectedBlkIDs.end(); ++interBlk) {

		double blockMinT = domainMinT;
		double blockMaxT = domainMaxT;

		double blockBounds[2][MESH_MDIM];
		mesh.getBlkBoundBox(*interBlk, blockBounds);

		for (int i = 0; i < nCellsX; i++) {
			for (int j = 0; j < nCellsY; j++) {
				double cellBounds[2][MESH_MDIM];
				getCellBoundBoxInMain(i,j,/*k,*/ nCellsArr, blockBounds, cellBounds);

				double cellMinT;
				double cellMaxT;
				bool foundIntersection = findIntersection(originVec, directionVec, cellBounds, 
														  blockMinT, blockMaxT, cellMinT, cellMaxT);
				if (foundIntersection){
					// Get relevant properties of the cell
					tempContainer.t   = cellMinT;

					tempContainer.rho = blocksMain[*interBlk][varIndex[0]][i][j][0];

					tempContainer.compMassScalar = blocksMain[*interBlk][varIndex[1]][i][j][0];

					sortVector.push_back(tempContainer);
				}
			}
		}

		lowerPropIter++;
		upperPropIter++;
	}

	/* 
		Process the vector of density values. Sort it and start counting once a t-value's density,
	*/
	//Sort by t.
	std::sort(sortVector.begin(), sortVector.end(), sortFunction);

	int radiusCount = 0;
	for(std::vector<sortContainer>::iterator it = sortVector.begin(); it != sortVector.end(); ++it) {
		if (it->compMassScalar != 0) { //Only count companion mass scalar material
			radiusCount++;
		}
	}

	//Make array to contain densities.
	double *radiusVals;
	double *densVals;
	
	radiusVals = new double[radiusCount];
	densVals = new double[radiusCount];
	for (int i = 0; i < radiusCount; i++) {
		radiusVals[i] = 0;
		densVals[i] = 0;
	}

	int radiusIter = 0;
	std::cout << "Radius count: " << radiusCount << std::endl;
	for(std::vector<sortContainer>::iterator it = sortVector.begin(); it != sortVector.end(); ++it) {
		if (it->compMassScalar != 0.0) {
			radiusVals[radiusIter] = originVec[IAXIS] + it->t * directionVec[IAXIS];
			densVals[radiusIter] = it->rho;
			// std::cout << radiusVals[radiusIter] << "," << densVals[radiusIter] << std::endl;
		    radiusIter++;
		}
	}

	//Output densities as a function of theta.
	printOut(outProgramFileName, radiusVals, densVals, radiusCount);


	for (int i = 0; i < numberOfBlocks; i++) {
		deallocateVarArray(blocksMain[i], mesh.getNUserVars(), mesh.getNcells());
	}

	for (int i = 0; i < nVars; i++) {
		delete[] vars[i];
	}
	delete[] vars;
	delete[] varIndex;
	delete[] densVals;
	delete[] radiusVals;

	return 0;

}


void findOrigin(double* originVec, double***** blocksMain, FlashAmrMesh &grid, std::set<int> leafBlockIDs, 
			    int numberOfBlocks, int* nCellsArr, int varIndex, ORIGINTYPE originType) {

	if (originType == USERDEF) {
		originVec[IAXIS] = ORIGINX;
		originVec[JAXIS] = ORIGINY;
		originVec[KAXIS] = ORIGINZ;
	}
	else {
		int maxLocation[4]; //[blockID, cellX, cellY, cellZ]
		maxLocation[0] = 0;
		maxLocation[1] = 0; //ORIGINX;
		// maxLocation[2] = 1e-12; //ORIGINY; I don't know why this is 1e-12, but i'm keeping it just in case.
		maxLocation[2] = 0; //ORIGINY;
		maxLocation[3] = 0; //ORIGINZ;

		double maxDensity = 0;
		std::set<int>::iterator leaf;
		
		int densCount = 0;
		//Find the highest density value as well as the location of that cell in the mesh.
		for (leaf = leafBlockIDs.begin(); leaf != leafBlockIDs.end(); ++leaf) {
			for (int i = 0; i < nCellsArr[IAXIS]; i++) {
				for (int j = 0; j < nCellsArr[JAXIS]; j++) {
					//for (int k = 0; k < nCellsArr[KAXIS]; k++) {
						if (blocksMain[*leaf][varIndex][i][j][0] > maxDensity) {
							maxDensity = blocksMain[*leaf][varIndex][i][j][0];
							maxLocation[0] = *leaf;
							maxLocation[1] = i;
							maxLocation[2] = j;
							//maxLocation[3] = k;
							densCount = 0;
						}
						if (blocksMain[*leaf][varIndex][i][j][0] == maxDensity) {
							densCount++;
						}
					//}
				}
			}
		}

		//Set the origin to that location.
		//Find the physical bounds of that block.
		double blockBounds[2][MESH_MDIM];
		grid.getBlkBoundBox(maxLocation[0], blockBounds);
		
		//Find physical bounds of the cell.
		double cellBounds[2][MESH_MDIM];
		getCellBoundBoxInMain(maxLocation[1], maxLocation[2], /*maxLocation[3],*/ 
						      nCellsArr, blockBounds, cellBounds);

		//Need to stay on the symmetry axis, so x will always be at 0. Maybe z will too. No clue at this point.
		originVec[IAXIS] = 0; //(cellBounds[LOWER][IAXIS] + cellBounds[UPPER][IAXIS])/2.f;
		originVec[JAXIS] = (cellBounds[LOWER][JAXIS] + cellBounds[UPPER][JAXIS])/2.f;
		originVec[KAXIS] = 0; //(cellBounds[LOWER][KAXIS] + cellBounds[UPPER][KAXIS])/2.f;

		std::cout << "Density Count: " << densCount << std::endl;
		std::cout << "Highest density point: " << originVec[1] << ". Value: " << maxDensity << std::endl; 

	}
}

bool findIntersection(double* originVec, double* directionVec, double cellBounds[2][MESH_MDIM], 
					  double lowerBoundT, double upperBoundT, double &tMinReturn, double &tMaxReturn) {

	// These are initialized as if they were txmin and txmax, 
	// but are updated throughout for the final output.
	double tmin;
	double tmax;
	double tymin, tymax; /* tzmin, tzmax */
	bool swappedX = false;
	bool swappedY = false;

    tmin = (cellBounds[LOWER][IAXIS] - originVec[IAXIS]) / directionVec[IAXIS];
    tmax = (cellBounds[UPPER][IAXIS] - originVec[IAXIS]) / directionVec[IAXIS];
    if (tmin > tmax) {
    	std::swap(tmin, tmax);
    	swappedX = true;
    }
    
    tymax = (cellBounds[UPPER][JAXIS] - originVec[JAXIS]) / directionVec[JAXIS];
    tymin = (cellBounds[LOWER][JAXIS] - originVec[JAXIS]) / directionVec[JAXIS];
    
    if (tymin > tymax){
    	std::swap(tymin, tymax);
    	swappedY = true;
    }
    if ((tmin > tymax) || (tymin > tmax)){
        return false;
    }
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    // tzmin = (boxBounds[LOWER][KAXIS] - originVec[JAXIS]) / directionVec[JAXIS];
    // tzmax = (boxBounds[UPPER][KAXIS] - originVec[JAXIS]) / directionVec[JAXIS];
    // if (tzmin > tzmax) swap(tzmin, tzmax);
    // if ((tmin > tzmax) || (tzmin > tmax))
    //     return false;
    // if (tzmin > tmin)
    //     tmin = tzmin;
    // if (tzmax < tmax)
    //     tmax = tzmax;

    //Check to see if we've exceeded the bounds. Swap if we swapped either of the t's above.
    // if (swapped) {
    // 	std::swap(lowerBoundT, upperBoundT);
    // }
    // if (tmin > upperBoundT || tmax < lowerBoundT) {
    // 	return false;
    // }
    	
    // //If we intersect the maximum radius, don't count this one.
    // if (tmin > domainMaxT || tmax > domainMaxT) {
    // 	return false;
    // }

    //find point corresponding to the min and the max t.
    
    //Special condition where we're inside the box where the origin of the vector lies,
    //or we're at the end of the ray.
    double pointMin[MESH_MDIM];
    double pointMax[MESH_MDIM];
    bool isInBox;

    pointMin[IAXIS] = originVec[IAXIS] + lowerBoundT * directionVec[IAXIS];
    pointMin[JAXIS] = originVec[JAXIS] + lowerBoundT * directionVec[JAXIS];
    // pointMin[IAXIS] = originVec[IAXIS] + 0 * directionVec[IAXIS];

    pointMax[IAXIS] = originVec[IAXIS] + upperBoundT * directionVec[IAXIS];
    pointMax[JAXIS] = originVec[JAXIS] + upperBoundT * directionVec[JAXIS];
    // pointMin[IAXIS] = originVec[IAXIS] + 0 * directionVec[IAXIS];

    isInBox = pointInBox(pointMin, cellBounds) || pointInBox(pointMax, cellBounds);

    if (isInBox) {
    	tmax = fmin(upperBoundT, tmax);
    	tmin = fmax(lowerBoundT, tmin);
    }

    if (tmin < lowerBoundT) {
    	return false;
    }

    if (tmax > upperBoundT) {
    	return false;
    }
    
    // if (tmax > upperBoundT) {
    // 	tMaxReturn = upperBoundT;
    // }
    // if (r.tmin < tmin) 
    // 	r.tmin = tmin;
    // if (r.tmax > tmax) 
    // 	r.tmax = tmax;
    
    // If we got this far, then there's an intersection.
    tMinReturn = tmin;
    tMaxReturn = tmax;
    return true;
}

double findCellVectorDistance(double* originVec, double* directionVec, double tmin, double tmax) {

	double p1X = originVec[IAXIS] + tmin * directionVec[IAXIS];
	double p1Y = originVec[JAXIS] + tmin * directionVec[JAXIS];
	//double p1Z = originVec[KAXIS] + tmin * directionVec[KAXIS];

	double p2X = originVec[IAXIS] + tmax * directionVec[IAXIS];
	double p2Y = originVec[JAXIS] + tmax * directionVec[JAXIS];
	//double p2Z = originVec[KAXIS] + tmax * directionVec[KAXIS];

	double difference = sqrt((p2X-p1X)*(p2X-p1X) + (p2Y-p1Y)*(p2Y-p1Y)/* + (p2Z-p1Z)*(p2Z-p1Z)*/);

	return difference;
}

bool pointInBox(double* point, double cellBounds[2][MESH_MDIM]) {
	if (point[IAXIS] >= cellBounds[LOWER][IAXIS] && point[IAXIS] <= cellBounds[UPPER][IAXIS] &&
		point[JAXIS] >= cellBounds[LOWER][JAXIS] && point[JAXIS] <= cellBounds[UPPER][JAXIS] /*&&
		point[IAXIS] >= cellBounds[LOWER][IAXIS] && point[IAXIS] <= cellBounds[UPPER][IAXIS] */)
	{
		return true;
	}
	else 
		return false;

}

void printOut(char* fileName, double* radiusVals, double* densVals, int radiusCount) {
	std::ofstream outFile;
	outFile.open(fileName);

	for (int i = 0; i < radiusCount; i++) {
		outFile << radiusVals[i] << "," << densVals[i] << std::endl;
		std::cout << radiusVals[i] << "," << densVals[i] << std::endl;
	}

	outFile.close();
}

bool sortFunction (sortContainer i,sortContainer j) { 
	return (i.t<j.t);
}