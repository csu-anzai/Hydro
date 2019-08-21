
#include "CommonSNBinary.h"
using namespace std;

enum MASSSCALAR {COMPANION, EJECTA, BOTH};

/* Parameters */
#define ORIGINCHOICE HIGHESTDENSITY //USERDEF
// Note that if you use HIGHESTDENSITY, it will start counting the column density 
// at the surface of the star, not from the center.
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

  double rInner;
  double rOuter;
  double sortCellMinT;
};

void findOrigin(double* originVec, double***** blocksMain, FlashAmrMesh &grid, std::set<int> leafBlockIDs, 
			    int numberOfBlocks, int* nCellsArr, int varIndex, ORIGINTYPE originType);
bool findIntersection(double* originVec, double* directionVec, double cellBounds[2][MESH_MDIM], 
					  double lowerBoundT, double upperBoundT, double &tMinReturn, double &tMaxReturn);
double findCellVectorDistance(double* originVec, double* directionVec, double tmin, double tmax);
bool pointInBox(double* point, double cellBounds[2][MESH_MDIM]);
void printOut(char* fileName, double* varValues, double thetaMin, double dTheta, int thetaCount);
bool sortFunction (sortContainer i,sortContainer j);

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv) {

	modelContainer models[NUMMODELS];
	readModels(models);

	/* Declare variables to be read */
	int nVars = 5;
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

	int checkNum = MS38_L200_F;
	// int checkNum = MS54_F;
	// int checkNum = MS63_F;
	// int checkNum = MS7_L200_F;
	// int checkNum = SG_L200_F;
	// int checkNum = SY319_L200_F;
	// int checkNum = SY428_F;

	MASSSCALAR massScalar = COMPANION;

	strcpy(vars[0],"dens");
	strcpy(vars[1],"gpot");
	strcpy(vars[2],"ener");
	strcpy(vars[3],"ms_2");
	strcpy(vars[4],"ms_3");
	strcpy(outFileExtraName, "_colDens_comp_vol.txt");

	/* Useful Variable names:
		dens = Density
		gpot = Gravitational Potential
		ener = Total energy
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
	strcat(outProgramFileName, "columnDensities/");

	strcat(inProgramFileName,inFileName);
	strcat(outProgramFileName,outFileName);

	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "WARNING: THIS ANALYSIS PROGRAM NOT YET BUILT FOR MPI" << std::endl;
	std::cout << "Input file: " << inProgramFileName << std::endl;
	std::cout << "Output file: " << outProgramFileName << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;


	//char inFileName[MAX_FILENAME_LENGTH] = "snb_hdf5_chk_0205";
	//char outFileName[MAX_FILENAME_LENGTH] = "snb_hdf5_chk_0205_colDens_ejecta.txt";
	// char outFileName[MAX_FILENAME_LENGTH] = "colDensOutPaths.txt";
	// char inProgramFileName[MAX_FILENAME_LENGTH] = ".";
	// char outProgramFilename[MAX_FILENAME_LENGTH] = "/"

	//char inProgramFileName[MAX_FILENAME_LENGTH] = "/home/pb12c/Data/SNB_2007_Data/w7_sy319_l200-20070509/";
	//char outProgramFileName[MAX_FILENAME_LENGTH] = "/home/pb12c/Data/SNB_2007_Data/w7_sy319_l200-20070509/columnDensities/";
	// char outFileType[MAX_STRING_LENGTH] = "ascii";
	
	/* ------------------------- */


	/* Parameters */
	double dTheta = DTHETA;
	double thetaMin = THETAMIN;
	double thetaMax = THETAMAX;
	double thetaCount = ceil((thetaMax - thetaMin) / dTheta) + 1; // +1 is for the zero. 0 to 180 is 181 values.
	// double densMaxAllow = DENSMAXALLOW;

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

	// // Temporarily change densities to a constant value for testing purposes.
	// for (int i = 0; i < numberOfBlocks; i++) {
	// 	for (int ii = 0; ii < mesh.getNcells(IAXIS); ii++) {
	// 		for (int jj = 0; jj < mesh.getNcells(JAXIS); jj++) {
	// 			for (int kk = 0; kk < mesh.getNcells(KAXIS); kk++) {
	// 				blocksMain[i][mesh.findVarIndex("dens")][ii][jj][kk] = 1.f;
	// 			}
	// 		}
	// 	}
	// }

	//Make array to contain densities.
	double *densVals;
	int densSize = thetaCount;
	//std::cout << "Size of density array: " << densSize << std::endl;
	densVals = new double[densSize];
	for (int i = 0; i < densSize; i++) {
		densVals[i] = 0.f;
	}
	int densValsIter = 0;

	// ---------------
	// MAIN LOOP SETUP
	// ---------------

	//Loop from -90 to 90 degrees with specified increments.
	double currentTheta = thetaMin;
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
	
	std::cout << "vec: " << originVec[JAXIS] << std::endl;
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

	std::cout << "Max Radius: " << maxRadius << std::endl;

	// ---------
	// MAIN LOOP
	// ---------


	while (currentTheta <= thetaMax) {
		std::vector<sortContainer> sortVector;
		sortContainer tempContainer;
		
		//Convert from polar to cartesian.
		directionVec[IAXIS] = maxRadius * cos(currentTheta * PI / 180.f) - originVec[IAXIS];
		directionVec[JAXIS] = maxRadius * sin(currentTheta * PI / 180.f) - originVec[JAXIS];
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
				//std::cout << "Found intersection. Block: " << *leaf << std::endl;
				//densVals[densValsIter] += 1;
			}

		}
		// std::cout << "Number of blocks intersected: " << intersectedBlkIDs.size() << std::endl;

		//Loop over cells in each block.
		lowerPropIter = intersectedLowerProp.begin();
		upperPropIter = intersectedUpperProp.begin();

		for (interBlk = intersectedBlkIDs.begin(); interBlk != intersectedBlkIDs.end(); ++interBlk) {

			//TODO: Still a problem here when using lower and upper prop iter.
			double blockMinT = domainMinT; //*lowerPropIter;
			double blockMaxT = domainMaxT; //*upperPropIter;
			//std::cout << "block: " << *interBlk << " tmin: " << tmin << " tmax: " << tmax << std::endl;

			double blockBounds[2][MESH_MDIM];
			mesh.getBlkBoundBox(*interBlk, blockBounds);

			for (int i = 0; i < nCellsX; i++) {
				for (int j = 0; j < nCellsY; j++) {
					// for (int k = 0; k < nCellsZ; k++) {
					double cellBounds[2][MESH_MDIM];
					getCellBoundBoxInMain(i,j,/*k,*/ nCellsArr, blockBounds, cellBounds);

					double cellMinT;
					double cellMaxT;
					bool foundIntersection = findIntersection(originVec, directionVec, cellBounds, 
															  blockMinT, blockMaxT, cellMinT, cellMaxT);
					// std::cout << "cellMinT " << cellMinT << std::endl;

					if (foundIntersection){
						tempContainer.sortCellMinT = cellMinT;
						// Get relevant properties of the cell
						double dist = findCellVectorDistance(originVec, directionVec, cellMinT, cellMaxT); 
						//blocksMain[*interBlk][varIndex][i][j][0] = currentTheta;
						tempContainer.t   = cellMinT;

						tempContainer.rho = blocksMain[*interBlk][varIndex[0]][i][j][0];

						double minR = /*originVec[IAXIS] + */cellMinT * directionVec[IAXIS];
						double minZ = /*originVec[JAXIS] + */cellMinT * directionVec[JAXIS];
	
						double maxR = /*originVec[IAXIS] + */cellMaxT * directionVec[IAXIS];
						double maxZ = /*originVec[JAXIS] + */cellMaxT * directionVec[JAXIS];
						
						tempContainer.rInner = sqrt(minR*minR + minZ*minZ);
						tempContainer.rOuter = sqrt(maxR*maxR + maxZ*maxZ);

						double volumeOuter = 4./3. * PI * pow(tempContainer.rOuter,3);

						// cellBounds[UPPER][IAXIS] * cellBounds[UPPER][IAXIS] * 
						// 				cellBounds[UPPER][IAXIS] * 
						// 				(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);
						
						double volumeInner = 4./3. * PI * pow(tempContainer.rInner,3);

						// cellBounds[LOWER][IAXIS] * cellBounds[LOWER][IAXIS] * 
						// 				cellBounds[UPPER][JAXIS] * 
						// 				(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);

						tempContainer.vol = volumeOuter - volumeInner;

						tempContainer.dr  = dist;


						tempContainer.totEnergy = blocksMain[*interBlk][varIndex[1]][i][j][0] + 
												  blocksMain[*interBlk][varIndex[2]][i][j][0];

						tempContainer.ejectaMassScalar = blocksMain[*interBlk][varIndex[3]][i][j][0];

						tempContainer.compMassScalar = blocksMain[*interBlk][varIndex[4]][i][j][0];

						sortVector.push_back(tempContainer);
					}
					//}
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

		// double endX, endY;
		//bool startCounting = false;
		double sumTotalMass = 0;
		bool startCounting = false;
		for(std::vector<sortContainer>::iterator it = sortVector.begin(); it != sortVector.end(); ++it) {
			if (!startCounting) {
				if (originType == USERDEF) { //Start counting automatically if not using companion as origin.
					startCounting = true;
				}
				else if (it->totEnergy > 0) { //Start counting at the surface of companion if using it as the origin.
				//if (it->ejectaMassScalar <= 0) {
					startCounting = true;
					
					// endX = originVec[IAXIS] + it->t * directionVec[IAXIS];
					// endY = originVec[JAXIS] + it->t * directionVec[JAXIS];
				}
			}
			if (startCounting && massScalar == COMPANION) {
				densVals[densValsIter] += it->compMassScalar * it->rho * it->vol;//it->dr;
				// densVals[densValsIter] += it->compMassScalar * it->rho * it->dr;
			}
			if (startCounting && massScalar == EJECTA) {
				densVals[densValsIter] += it->ejectaMassScalar * it->rho * it->vol;//it->dr;
				// densVals[densValsIter] += it->ejectaMassScalar * it->rho * it->dr;
			}
			if (startCounting && massScalar == BOTH) {
			// if (startCounting && massScalar == BOTH && it->rInner >= 3e12) {
				densVals[densValsIter] += it->ejectaMassScalar * it->rho * it->vol;//it->dr;
				// densVals[densValsIter] += it->vol;
				// densVals[densValsIter] += it->vol;//it->dr;
				densVals[densValsIter] += it->compMassScalar * it->rho * it->vol;//it->dr;
				// densVals[densValsIter] += it->ejectaMassScalar * it->rho * it->dr;
				// densVals[densValsIter] += it->compMassScalar * it->rho * it->dr;
			}
			// std::cout << "CellMinT: " << it->sortCellMinT << std::endl;
			// std::cout << "Inner: " << it->rInner << " : " << "Outer: " << it->rOuter << std::endl;

			// if (it == sortVector.end()-1) {
			// 	// std::cout << "Dividing" << std::endl;
			// 	densVals[densValsIter] /= (maxRadius * maxRadius * 4 * PI);
			// }
		}

		//Multiply by the area to allow for comparison across models.
		//double snejRadius = sqrt(endX*endX + endY*endY);
		//double tempDens = densVals[densValsIter];

		//densVals[densValsIter] = densVals[densValsIter] * (4 * PI * snejRadius * snejRadius);

		//std::cout << "End: (" << endX << "," << endY << ") | Mass Val: " << densVals[densValsIter] << " | Mass/Area: " << tempDens << 
		//" | Radius used: " << snejRadius << std::endl;
		// Increment theta and the position in the density values array.
		std::cout << "Theta: " << currentTheta << " | Mass Val: " << densVals[densValsIter] << std::endl;

		densValsIter++;
		currentTheta += dTheta;

	}

	std::cout << "vec: " << originVec[IAXIS] << std::endl;

	// std::cout << "Radius Used: " << maxRadius << std::endl;

	//Output densities as a function of theta.
	printOut(outProgramFileName, densVals, thetaMin, dTheta, thetaCount);

	// /* Test stuff for checking vector trajectories */
	// //Put trajectory information from the blocks in main to the one in the class.
	// for (int i = 0; i < numberOfBlocks; i++) {
	// 	mesh.putBlkData(i,CENTER,var[0],INTERIOR,startingPos,blocksMain[i]);
	// }

	// double**** tempData;
	// mesh.refineToFinest(tempData, MAXREFINE);

	// mesh.printRefinedData(tempData, outProgramFileName, outFileType);

	// mesh.deallocateFinest(tempData);

	// delete[] densVals;

	for (int i = 0; i < numberOfBlocks; i++) {
		deallocateVarArray(blocksMain[i], mesh.getNUserVars(), mesh.getNcells());
	}

	for (int i = 0; i < nVars; i++) {
		delete[] vars[i];
	}
	delete[] vars;
	delete[] varIndex;
	delete[] densVals;

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
		// maxLocation[2] = 1e-12; //ORIGINY;
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

void printOut(char* fileName, double* varValues, double thetaMin, double dTheta, int thetaCount) {
	std::ofstream outFile;
	outFile.open(fileName);

	for (int i = 0; i < thetaCount; i++) {
		outFile << thetaMin + (i*dTheta) << "," << varValues[i] << std::endl;
		//std::cout << thetaMin + (i*dTheta) << "," << varValues[i] << std::endl;
	}

	outFile.close();
}

bool sortFunction (sortContainer i,sortContainer j) { 
	return (i.t<j.t);
}