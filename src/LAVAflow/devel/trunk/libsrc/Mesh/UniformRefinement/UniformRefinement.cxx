#include "FlashAmrMesh.h"
// #include "Driver_main.h"

int compare (const void * a, const void * b);
void findFactors(int *primeFactors, int procs, int &numFactors);
bool intersectsSubdomain(int localBlockCornerIDs[MESH_MDIM], int blockMaxIndexExtent[MESH_MDIM], 
				 int subdomainBlockCornerIDs[2][MESH_MDIM], bool flag);
void findDestinationProcessor(int procsX, int procsY, int procsZ, 
							  int offBlkX, int offBlkY, int offBlkZ,
							  int ii, int jj, int kk,
							  int dim, int *procIndexMin[MESH_MDIM],
							  bool subdomainFlag, bool containedInSubdomainFlag, int subdomainBlockCornerIDs[2][MESH_MDIM],
							  int &destinationProcessor);
bool inSubdomain(int i, int j, int k, int subdomainBlockCornerIDs[2][MESH_MDIM]);
void allocateVarArray(double****&, int, int*);
void deallocateVarArray(double****&, int, int*);

void FlashAmrMesh::refineToFinest(double****& inData, double subdomainCoords[2][MESH_MDIM], int refineLevel) {
	/* Pass a NON-ALLOCATED 4-d pointer to this function.
	 * Deallocate that pointer using deallocateFinest().
	 *
	 * Takes an AMR mesh with ALL variables specified in Main (using addUserVariable or addAllVariables)
	 * and returns an array of cells indexed with [variable index][x][y][z] refined to the value specified 
	 * in the argument 'refineLevel'. 
	 * Ex: Using this function on a mesh with block refinement levels ranging from 1 to 9 and a 
	 *     'refineLevel' value set to 5 will return a mesh of cells at refinement level 5. 
	 *     To refine the entire mesh to a refinement level of 9, enter either 9 or -1. uniform refinement. 
	 * This data array can be output in HDF5 format using the printRefinedDataHDF5 function, and
	 * the image for the hdf5 file can be immediately viewed in HDFView, but not VisIt. 
	 * More details regarding what is output in this function are available in its documentation. 
	 * The printRefinedData function can be modified to output variable data in either ascii or binary form, 
	 * but this function is there for testing purposes and is not robust enough for regular use. 
	 *
	 * A subdomain (as physical coordinates) can be specified if only a small area of the mesh is of interest. This will be
	 * declared in main and passed in as an argument. Note that due to the discretization and the
	 * way that a cell is determined to intersect the specified subdomain, the output subdomain coordinates
	 * will be slightly offset from the input coordinates. Also note that the subdomain can't be too small 
	 * (smaller than one cell, or smaller than one cell per processor) or too large (larger than the original domain).
	 * The program will exit with relevant info if either of these cases are found.
	 *
	 * This function can be run in serial or in parallel with any number of processors. For parallel processing, 
	 * the processors are split up using a prime factorization method, where processors are iteratively distributed 
	 * using these factors to the dimensions with the largest amount of remaining cells. An example of this might be a domain with 
	 * 120 cells in the y-dimension and 60 in the x-dimension. Split across 12 processors you'll have 3 assigned to 
	 * the y-dimension (3 total, 40 per processor), then 2 to the x-dimension (2 total, 30 per processor), 
	 * then 3*2 to the y-dimension (6 total, 20 per processor). 2 * (2*3) = 12.
	 *
	 * This function will refine 1-, 2-, and 3-dimensional datasets. 
	 * 
	 * Example program is in the src directory: testRefinement.cpp.
	 * 
	 * Installation notes are in src/INSTALL.txt, setup for NERSC is in src/README_NERSC.txt.
	 * Example scripts for the interactive queue and for submitting a job are also available 
	 * in the src directory.
	 *
	 * Inputs: 
	 *     Non-allocated 4d pointer
	 *     Refinement level (use -1 if the most fine uniform mesh is wanted)
	 *     Subdomain physical coordinates
	 *
	 * Outputs:
	 *     inData: 4-d non-allocated data array that will contain the refined data.
	 */

	typedef struct mpiVarStruct {
		double value;
		int x,y,z;
	} mpiVarStruct;

	
	/* Crude checking mechanism to see if we need to focus on a subdomain rather than
	   the whole domain. If no subdomain is desired, set all the values to 0.
	*/
	bool subdomainFlag = false;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < MESH_MDIM; j++) {
			if (subdomainCoords[i][j] != 0) {
				subdomainFlag = true;
			}
		}
	}

	if (subdomainFlag) {
		/* Make sure the subdomain is actually within the original domain. */
		if (getDim() >= 0) {
			if (   subdomainCoords[LOWER][IAXIS] < realScalarsMap["xmin"]
				|| subdomainCoords[UPPER][IAXIS] > realScalarsMap["xmax"]
			   ){
				std::cout << "Subdomain X-coordinates exceed domain coordinates: " << std::endl;
				std::cout << "\tChosen subdomain X-coordinate range: (" << subdomainCoords[LOWER][IAXIS] << " - "
																   		<< subdomainCoords[UPPER][IAXIS] << ").";
				std::cout << std::endl;
				std::cout << "\tDomain X-coordinate range: (" << realScalarsMap["xmin"] << " - " 
															  << realScalarsMap["xmax"] << ").";
				std::cout << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		if (getDim() >= 1) {
			if (   subdomainCoords[LOWER][JAXIS] < realScalarsMap["ymin"]
				|| subdomainCoords[UPPER][JAXIS] > realScalarsMap["ymax"]
			   ){
				std::cout << "Subdomain Y-coordinates exceed domain coordinates: " << std::endl;
				std::cout << "\tChosen subdomain Y-coordinate range: (" << subdomainCoords[LOWER][JAXIS] << " - " 
																   		<< subdomainCoords[UPPER][JAXIS] << ").";
				std::cout << std::endl;
				std::cout << "\tDomain Y-coordinate range: (" << realScalarsMap["ymin"] << " - " 
															  << realScalarsMap["ymax"] << ").";
				std::cout << std::endl;
				exit(EXIT_FAILURE);	
			}
		}
		if (getDim() >= 2) {
			if (   subdomainCoords[LOWER][KAXIS] < realScalarsMap["zmin"]
				|| subdomainCoords[UPPER][KAXIS] > realScalarsMap["zmax"]
			   ){
				std::cout << "Subdomain Z-coordinates exceed domain coordinates: " << std::endl;
				std::cout << "\tChosen subdomain Z-coordinate range: (" << subdomainCoords[LOWER][KAXIS] << " - " 
																   		<< subdomainCoords[UPPER][KAXIS] << ").";
				std::cout << std::endl;
				std::cout << "\tDomain Z-coordinate range: (" << realScalarsMap["zmin"] << " - " 
															  << realScalarsMap["zmax"] << ").";
				std::cout << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// std::cout << realScalarsMap["xmin"] << std::endl;
		// std::cout << realScalarsMap["xmax"] << std::endl;
		// std::cout << realScalarsMap["ymin"] << std::endl;
		// std::cout << realScalarsMap["ymax"] << std::endl;
		// std::cout << realScalarsMap["zmin"] << std::endl;
		// std::cout << realScalarsMap["zmax"] << std::endl;





		//  If we choose a subdomain that's a dimension smaller than our original domain
		//    (say a 2D subdomain within the original 3D domain), change the dimension to 
		//    that subdomain's dimension. (Not yet implemented, no current plans to. 
		//    May be implemented later.)
		// bool fixedIndexFlag[MESH_MDIM];
		// int fixedIndex[MESH_MDIM];

		// fixedIndexFlag[IAXIS] = false;
		// fixedIndexFlag[JAXIS] = false;
		// fixedIndexFlag[KAXIS] = false;

		// int fixedY, fixedZ;
		// int refinementDim = getDim();
		// for (int i = MESH_MDIM-1; i >= 0; i--) {
		// 	if (subdomainCoords[UPPER][i] == subdomainCoords[LOWER][i]) {
		// 		refinementDim--;
		// 		fixedIndexFlag[i] = true;
		//		fixedIndex[i] = /* the index corresponding to the coordinate chosen. */
		// 	}
		// }
	}

	
	if (mpiGetID() == 0) {
		std::cout << "Starting refinement..." << std::endl;
	}

	int nblockx, nblocky, nblockz;
	nblockx = intScalarsMap["nblockx"];
	nblocky = intScalarsMap["nblocky"];
	nblockz = intScalarsMap["nblockz"];

	/* Find how refined the mesh gets. If it's not specified, refine to the finest level. */
	int maxRefineLevel;
	maxRefineLevel = blocks[0].refineLevel;
	for (int i = 1; i < getLocalNumBlocks() /* getNBlocks() */; i++) {
		maxRefineLevel = std::max(maxRefineLevel,blocks[i].refineLevel);
	}

	// MPI: Get the maximum refinement level across all processors
	//Get the maxiumum across all processors, then broadcast it. Allreduce does all this.

	int globalMaxRefineLevel;
	MPI_Allreduce(&maxRefineLevel, &globalMaxRefineLevel, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
	// std::cout << "Global max refine level: " << globalMaxRefineLevel << std::endl;

	if (refineLevel > globalMaxRefineLevel) {
		std::cout << "WARNING: Specified refinement level is too high. Defaulting to maximum refine level " << globalMaxRefineLevel << std::endl;
		refineLevel = globalMaxRefineLevel;
	}

	if (refineLevel != -1){
		globalMaxRefineLevel = refineLevel;
	}

	/* Find the minimum and maximum index limits of the grid, get block index limits. */
	/* Based off Flash functions Grid_getBlkCornerId.F90 and IO_readCheckpoint.F90 */

	double blkBoundBox[2][MESH_MDIM];
	double gridBoundBox[2][MESH_MDIM];

	//Initialize grid index search
	getBlkBoundBox(0, blkBoundBox);
	for (int j = 0; j < MESH_MDIM; j++) {
		gridBoundBox[LOWER][j] = blkBoundBox[LOWER][j];
		gridBoundBox[UPPER][j] = blkBoundBox[UPPER][j];
	}
	for (int i = 1; i < getLocalNumBlocks(); i++) {
		getBlkBoundBox(i, blkBoundBox);
		
		for (int j = 0; j < MESH_MDIM; j++) {
			gridBoundBox[LOWER][j] = fmin(blkBoundBox[LOWER][j],gridBoundBox[LOWER][j]);
			gridBoundBox[UPPER][j] = fmax(blkBoundBox[UPPER][j],gridBoundBox[UPPER][j]);
		}
	}

	double globalGridBoundBox[2][MESH_MDIM];
	//Distribute the grid bounds to the processors.
	for (int j = 0; j < MESH_MDIM; j++) {
		MPI_Allreduce(&gridBoundBox[LOWER][j], &globalGridBoundBox[LOWER][j], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&gridBoundBox[UPPER][j], &globalGridBoundBox[UPPER][j], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}

	gridDelta[IAXIS] = (globalGridBoundBox[UPPER][IAXIS] - globalGridBoundBox[LOWER][IAXIS])
				     / (1.0 * nCellsVec[IAXIS] * nblockx * pow(2,globalMaxRefineLevel-1));
	gridDelta[JAXIS] = (globalGridBoundBox[UPPER][JAXIS] - globalGridBoundBox[LOWER][JAXIS])
				     / (1.0 * nCellsVec[JAXIS] * nblocky * pow(2,globalMaxRefineLevel-1));
	gridDelta[KAXIS] = (globalGridBoundBox[UPPER][KAXIS] - globalGridBoundBox[LOWER][KAXIS])
				     / (1.0 * nCellsVec[KAXIS] * nblockz * pow(2,globalMaxRefineLevel-1));

	double gridHalfDelta[MESH_MDIM];
	for (int j = 0; j < MESH_MDIM; j++) {
		gridHalfDelta[j] = gridDelta[j] / 2.0;
	}

	int **localBlockCornerIDs;
	int subdomainBlockCornerIDs[2][MESH_MDIM];

	/* Figure out the local block IDs first */
	localBlockCornerIDs = new int*[getLocalNumBlocks()]; //Allocate space for the pointers
	for (int i = 0; i < getLocalNumBlocks(); i++) {
		localBlockCornerIDs[i] = new int[MESH_MDIM];
	}

	/* Finally calculate the corner IDs for each block. */
	for (int i = 0; i < getLocalNumBlocks(); i++) {
		getBlkBoundBox(i, blkBoundBox);
		for (int j = 0; j < MESH_MDIM; j++) {
			localBlockCornerIDs[i][j] = int(((blkBoundBox[LOWER][j] - globalGridBoundBox[LOWER][j] + gridHalfDelta[j])
										/ gridDelta[j]));
		}
	}
	if (subdomainFlag) {
		for (int j = 0; j < MESH_MDIM; j++) {
			if (j < getDim()) 
			{
				subdomainBlockCornerIDs[LOWER][j] = int( (subdomainCoords[LOWER][j] - globalGridBoundBox[LOWER][j]) / gridDelta[j] + 0.5 );
				subdomainBlockCornerIDs[UPPER][j] = int( (subdomainCoords[UPPER][j] - globalGridBoundBox[LOWER][j]) / gridDelta[j] + 0.5 );
				
				// subdomainBlockCornerIDs[LOWER][j] = int(((subdomainCoords[LOWER][j] - globalGridBoundBox[LOWER][j] + gridHalfDelta[j])
				// 								/ gridDelta[j]));
				// subdomainBlockCornerIDs[UPPER][j] = int(((subdomainCoords[UPPER][j] - globalGridBoundBox[LOWER][j] + gridHalfDelta[j])
				// 								/ gridDelta[j]));


				// std::cout << subdomainBlockCornerIDs[LOWER][j] << " " << subdomainBlockCornerIDs[UPPER][j] << std::endl;
			}
			else {
				subdomainBlockCornerIDs[LOWER][j] = 0;
				subdomainBlockCornerIDs[UPPER][j] = 0;
			}
		}
	}

	
	int maxScale = (int)pow(2,globalMaxRefineLevel-1);

	//Calculate how many blocks are in each direction for each processor.
	int fineBlocksX = nblockx * maxScale;
	int fineBlocksY = nblocky * maxScale;
	int fineBlocksZ = nblockz * maxScale;

	int subdomainCellsX, subdomainCellsY, subdomainCellsZ;
	if (subdomainFlag) {
		if (getDim() == 1) {
			subdomainCellsX = subdomainBlockCornerIDs[UPPER][IAXIS] - subdomainBlockCornerIDs[LOWER][IAXIS] + 1;
			subdomainCellsY = 1;
			subdomainCellsZ = 1;	
		}
		else if (getDim() == 2) {
			subdomainCellsX = subdomainBlockCornerIDs[UPPER][IAXIS] - subdomainBlockCornerIDs[LOWER][IAXIS] + 1;
			subdomainCellsY = subdomainBlockCornerIDs[UPPER][JAXIS] - subdomainBlockCornerIDs[LOWER][JAXIS] + 1;
			subdomainCellsZ = 1;
		}
		else {
			subdomainCellsX = subdomainBlockCornerIDs[UPPER][IAXIS] - subdomainBlockCornerIDs[LOWER][IAXIS] + 1;
			subdomainCellsY = subdomainBlockCornerIDs[UPPER][JAXIS] - subdomainBlockCornerIDs[LOWER][JAXIS] + 1;
			subdomainCellsZ = subdomainBlockCornerIDs[UPPER][KAXIS] - subdomainBlockCornerIDs[LOWER][KAXIS] + 1;
		}		
	}

	//Calculate how many cells are in each direction for each processor.
	int fineCellsX, fineCellsY, fineCellsZ;
	int remainCellsX, remainCellsY, remainCellsZ;
	int validCellsX, validCellsY, validCellsZ;
	
	int procsX = 0;
	int procsY = 0;
	int procsZ = 0;

	if (getDim() == 1) {
		
		//No need to factor when there's only one processor.
		
		procsX = mpiGetProcs();

		if (subdomainFlag) {
			fineCellsX = subdomainCellsX / procsX;
			remainCellsX = subdomainCellsX - fineCellsX * procsX;
		}
		else {
			fineCellsX = fineBlocksX * getNcells(IAXIS) / procsX;
			remainCellsX = fineBlocksX * getNcells(IAXIS) - fineCellsX * procsX;
		}

		if (subdomainFlag){
			if (fineCellsX <= 0) {
				//There aren't enough cells to distribute across processors. Exit program.
				std::cout << "Subdomain too small in X-direction. Not enough cells. " << std::endl;
				std::cout << "Number of processors in X-direction: " << procsX << std::endl;
				std::cout << "Number of subdomain cells in X-direction: " << subdomainCellsX << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		fineCellsY = 1;
		fineCellsZ = 1;

		remainCellsY = 0;
		remainCellsZ = 0;

		// if (remainCellsX % procsX != 0) {
		// 	std::cout << "Uneven split in X-direction." << std::endl;
		// }

		fineCellsX += (int)ceil(1.0 * remainCellsX / procsX);

		if (remainCellsX % procsX != 0) {
			if (mpiGetID() % procsX < remainCellsX % procsX) {
				validCellsX = fineCellsX;
			}
			else {
				validCellsX = fineCellsX - 1;
			}
		}
		else {
			validCellsX = fineCellsX;
		}
		validCellsY = 1;
		validCellsZ = 1;

		// std::cout << "Fine cells x: " << fineCellsX << std::endl;
		// std::cout << "Valid cells x: " << validCellsX << std::endl;
	}
	if (getDim() == 2) {

		int maxProcs = mpiGetProcs();
		
		int *primeFactors;
		primeFactors = new int[maxProcs];
		
		int numFactors;

		findFactors(primeFactors, maxProcs, numFactors);

		procsX = 1;
		procsY = 1;
		procsZ = 0;

		for (int i = 0; i < numFactors; i++) {
			if (fineBlocksX * getNcells(IAXIS) / procsX >= fineBlocksY * getNcells(JAXIS) / procsY) {
				procsX = procsX * primeFactors[i];
			}
			else {
				procsY = procsY * primeFactors[i];
			}
		}

		delete[] primeFactors;

		if (subdomainFlag) {
			fineCellsX = subdomainCellsX / procsX;
			fineCellsY = subdomainCellsY / procsY;
			remainCellsX = subdomainCellsX - fineCellsX * procsX;
			remainCellsY = subdomainCellsY - fineCellsY * procsY;
		}
		else {
			fineCellsX = fineBlocksX * getNcells(IAXIS) / procsX;
			fineCellsY = fineBlocksY * getNcells(JAXIS) / procsY;
			remainCellsX = fineBlocksX * getNcells(IAXIS) - fineCellsX * procsX;
			remainCellsY = fineBlocksY * getNcells(JAXIS) - fineCellsY * procsY;
		}

		fineCellsZ = 1;
		remainCellsZ = 0;

		if (subdomainFlag) {
			if (fineCellsX <= 0) {
				//There aren't enough cells to distribute across processors. Exit program.
				std::cout << "Subdomain too small in X-direction. Not enough cells. " << std::endl;
				std::cout << "Number of processors in X-direction: " << procsX << std::endl;
				std::cout << "Number of subdomain cells in X-direction: " << subdomainCellsX << std::endl;
				std::cout << "(The dimensions are checked in a serial manner. The other dimensions may also be too small)" << std::endl;				
				exit(EXIT_FAILURE);
			}
			if (fineCellsY <= 0) {
				//There aren't enough cells to distribute across processors. Exit program.
				std::cout << "Subdomain too small in Y-direction. Not enough cells. " << std::endl;
				std::cout << "Number of processors in Y-direction: " << procsY << std::endl;
				std::cout << "Number of subdomain cells in Y-direction: " << subdomainCellsY << std::endl;
				std::cout << "(The dimensions are checked in a serial manner. The other dimensions may also be too small)" << std::endl;				
				exit(EXIT_FAILURE);
			}
		}

		// std::cout << "Procs in x: " << procsX << std::endl;
		// std::cout << "Procs in y: " << procsY << std::endl;
		// std::cout << "Fine cells x: " << fineCellsX << std::endl;
		// std::cout << "Remain cells x: " << remainCellsX << std::endl;
		// std::cout << "Fine cells y: " << fineCellsY << std::endl;
		// std::cout << "Remain cells y: " << remainCellsY << std::endl;
		// std::cout << "Total cells across procs x: " << fineBlocksX * getNcells(IAXIS) << std::endl;
		// std::cout << "Total cells across procs y: " << fineBlocksY * getNcells(JAXIS) << std::endl;

		// if (remainCellsX % procsX != 0) {
		// 	std::cout << "Uneven split in X-direction." << std::endl;
		// }
		// if (remainCellsY % procsY != 0) {
		// 	std::cout << "Uneven split in Y-direction." << std::endl;
		//}
		// std::cout << "Allocated processors: " << procsX * procsY << std::endl;
		// std::cout << "Remaining X: " << remainCellsX << " Remaining Y: " << remainCellsY << std::endl;
		// std::cout << "Distributed across processors X: " << ceil(remainCellsX / procsX) << std::endl;
		// std::cout << "Distributed across processors Y: " << ceil(remainCellsY / procsY) << std::endl;

		fineCellsX += (int)ceil(1.0 * remainCellsX / procsX);
		fineCellsY += (int)ceil(1.0 * remainCellsY / procsY);

		if (remainCellsX % procsX != 0) {
			if (mpiGetID() % procsX < remainCellsX % procsX) {
				validCellsX = fineCellsX;
			}
			else {
				validCellsX = fineCellsX - 1;
			}
		}
		else {
			validCellsX = fineCellsX;
		}		
		if (remainCellsY % procsY != 0){
			if ((mpiGetID() - (mpiGetID() % procsX)) / procsX < remainCellsY % procsY) {
				validCellsY = fineCellsY;
			}
			else {
				validCellsY = fineCellsY - 1;
			}
		}
		else {
			validCellsY = fineCellsY;
		}
		validCellsZ = 1;
	}
	if (getDim() == 3) {

		int maxProcs = mpiGetProcs();
		
		int *primeFactors;
		primeFactors = new int[maxProcs];
		
		int numFactors;

		findFactors(primeFactors, maxProcs, numFactors);

		procsX = 1;
		procsY = 1;
		procsZ = 1;

		for (int i = 0; i < numFactors; i++) {
			if ((fineBlocksX * getNcells(IAXIS) / procsX >= fineBlocksY * getNcells(JAXIS) / procsY)
			 && (fineBlocksX * getNcells(IAXIS) / procsX >= fineBlocksZ * getNcells(KAXIS) / procsZ)) {
				procsX = procsX * primeFactors[i];
			}
			else if ((fineBlocksY * getNcells(JAXIS) / procsY > fineBlocksX * getNcells(IAXIS) / procsX)
				  && (fineBlocksY * getNcells(JAXIS) / procsY >= fineBlocksZ * getNcells(KAXIS) / procsZ)) {
				procsY = procsY * primeFactors[i];
			}
			else { // z > x, z > y
				procsZ = procsZ * primeFactors[i];
			}
		}

		delete[] primeFactors;

		if (subdomainFlag) {
			fineCellsX = subdomainCellsX / procsX;
			fineCellsY = subdomainCellsY / procsY;
			fineCellsZ = subdomainCellsZ / procsZ;
			remainCellsX = subdomainCellsX - fineCellsX * procsX;
			remainCellsY = subdomainCellsY - fineCellsY * procsY;
			remainCellsZ = subdomainCellsZ - fineCellsZ * procsZ;
		}
		else {
			fineCellsX = fineBlocksX * getNcells(IAXIS) / procsX;
			fineCellsY = fineBlocksY * getNcells(JAXIS) / procsY;
			fineCellsZ = fineBlocksZ * getNcells(KAXIS) / procsZ;
			remainCellsX = fineBlocksX * getNcells(IAXIS) - fineCellsX * procsX;
			remainCellsY = fineBlocksY * getNcells(JAXIS) - fineCellsY * procsY;
			remainCellsZ = fineBlocksZ * getNcells(KAXIS) - fineCellsZ * procsZ;
		}		

		if (subdomainFlag) {
			if (fineCellsX <= 0) {
				//There aren't enough cells to distribute across processors. Exit program.
				std::cout << "Subdomain too small in X-direction. Not enough cells. " << std::endl;
				std::cout << "Number of processors in X-direction: " << procsX << std::endl;
				std::cout << "Number of subdomain cells in X-direction: " << subdomainCellsX << std::endl;
				std::cout << "(The dimensions are checked in a serial manner. The other dimensions may also be too small)" << std::endl;				
				exit(EXIT_FAILURE);
			}
			if (fineCellsY <= 0) {
				//There aren't enough cells to distribute across processors. Exit program.
				std::cout << "Subdomain too small in Y-direction. Not enough cells. " << std::endl;
				std::cout << "Number of processors in Y-direction: " << procsY << std::endl;
				std::cout << "Number of subdomain cells in Y-direction: " << subdomainCellsY << std::endl;
				std::cout << "(The dimensions are checked in a serial manner. The other dimensions may also be too small)" << std::endl;				
				exit(EXIT_FAILURE);
			}
			if (fineCellsZ <= 0) {
				//There aren't enough cells to distribute across processors. Exit program.
				std::cout << "Subdomain too small in Z-direction. Not enough cells. " << std::endl;
				std::cout << "Number of processors in Z-direction: " << procsZ << std::endl;
				std::cout << "Number of subdomian cells in Z-direction: " << subdomainCellsZ << std::endl;
				std::cout << "(The dimensions are checked in a serial manner. The other dimensions may also be too small)" << std::endl;				
				exit(EXIT_FAILURE);
			}			
		}

		// std::cout << "Procs in x: " << procsX << std::endl;
		// std::cout << "Procs in y: " << procsY << std::endl;
		// std::cout << "Procs in z: " << procsZ << std::endl;
		// std::cout << "Fine cells x: " << fineCellsX << std::endl;
		// std::cout << "Remain cells x: " << remainCellsX << std::endl;
		// std::cout << "Fine cells y: " << fineCellsY << std::endl;
		// std::cout << "Remain cells y: " << remainCellsY << std::endl;
		// std::cout << "Fine cells z: " << fineCellsZ << std::endl;
		// std::cout << "Remain cells z: " << remainCellsZ << std::endl;
		// std::cout << "Total cells across procs x: " << fineBlocksX * getNcells(IAXIS) << std::endl;
		// std::cout << "Total cells across procs y: " << fineBlocksY * getNcells(JAXIS) << std::endl;
		// std::cout << "Total cells across procs z: " << fineBlocksZ * getNcells(KAXIS) << std::endl;


		// if (remainCellsX % procsX != 0) {
		// 	std::cout << "Uneven split in X-direction." << std::endl;
		// }
		// if (remainCellsY % procsY != 0) {
		// 	std::cout << "Uneven split in Y-direction." << std::endl;
		// }
		// if (remainCellsZ % procsZ != 0) {
		// 	std::cout << "Uneven split in Z-direction." << std::endl;	
		// }
		// std::cout << "Allocated processors: " << procsX * procsY * procsZ << std::endl;
		// std::cout << "Remaining X: " << remainCellsX << " Remaining Y: " << remainCellsY << " Remaining Z: " << remainCellsZ << std::endl;
		// std::cout << "Distributed across processors X: " << ceil(1.0 * remainCellsX / procsX) << std::endl;
		// std::cout << "Distributed across processors Y: " << ceil(1.0 * remainCellsY / procsY) << std::endl;
		// std::cout << "Distributed across processors Z: " << ceil(1.0 * remainCellsZ / procsZ) << std::endl;

		fineCellsX += (int)ceil(1.0 * remainCellsX / procsX);
		fineCellsY += (int)ceil(1.0 * remainCellsY / procsY);
		fineCellsZ += (int)ceil(1.0 * remainCellsZ / procsZ);

		if (remainCellsX % procsX != 0) {
			if (mpiGetID() % procsX < remainCellsX % procsX) {
				validCellsX = fineCellsX;
			}
			else {
				validCellsX = fineCellsX - 1;
			}
		}
		else {
			validCellsX = fineCellsX;
		}
		if (remainCellsY % procsY != 0){
			if ((mpiGetID() - (mpiGetID() % procsX)) / procsX < remainCellsY % procsY) {
				validCellsY = fineCellsY;
			}
			else {
				validCellsY = fineCellsY - 1;
			}
		}
		else {
			validCellsY = fineCellsY;
		}
		if (remainCellsZ % procsZ != 0) {
			if ((mpiGetID() - (mpiGetID() % (procsX * procsY))) / (procsX * procsY) < remainCellsZ % procsZ) {
				validCellsZ = fineCellsZ;
			}
			else {
				validCellsZ = fineCellsZ - 1;
			}
		}
		else {
			validCellsZ = fineCellsZ;
		}
	}

	//Get the leaf blocks on each processor.
	std::set<int> leafBlockIDs;
	std::set<int>::iterator leaf;
	int blockMaxIndexExtent[MESH_MDIM];

	int leaves = 0;
	//Used if a maximum refinement level is specified.
	if (refineLevel != -1) {
		for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) {
			getBlkBoundBox(i, blkBoundBox);
			for (int j = 0; j < getDim(); j++) {
				blockMaxIndexExtent[j] = int(((blkBoundBox[UPPER][j] - globalGridBoundBox[LOWER][j] + gridHalfDelta[j])
										/ gridDelta[j]));
			}
			for (int j = getDim(); j < MESH_MDIM; j++) {
				blockMaxIndexExtent[j] = 0;
			}
			if (   getNodeType(i) == 1 && getRefineLevel(i) < refineLevel 
			    && intersectsSubdomain(localBlockCornerIDs[i], blockMaxIndexExtent, subdomainBlockCornerIDs, subdomainFlag)) {
				leafBlockIDs.insert(i);
				leaves++;
				// std::cout << i << " " << localBlockCornerIDs[i][KAXIS] << " " << blockMaxIndexExtent[KAXIS] << std::endl;
			}
			if (   getRefineLevel(i) == refineLevel 
				&& intersectsSubdomain(localBlockCornerIDs[i], blockMaxIndexExtent, subdomainBlockCornerIDs, subdomainFlag)) {
				leafBlockIDs.insert(i);
				leaves++;
				// std::cout << i << " " << localBlockCornerIDs[i][KAXIS] << " " << blockMaxIndexExtent[KAXIS] << std::endl;
			}
		}
	}
	//Used if no maximum refinement level is specified.
	else {
		//Get the leaf blocks from the list of blocks.
		for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) 
		{
			getBlkBoundBox(i, blkBoundBox);
			for (int j = 0; j < getDim(); j++) 
			{
				blockMaxIndexExtent[j] = int(((blkBoundBox[UPPER][j] - globalGridBoundBox[LOWER][j] + gridHalfDelta[j])
										/ gridDelta[j]));
			}
			for (int j = getDim(); j < MESH_MDIM; j++) 
			{
				blockMaxIndexExtent[j] = 0;
			}
			if (   getNodeType(i) == 1
				&& intersectsSubdomain(localBlockCornerIDs[i], blockMaxIndexExtent, subdomainBlockCornerIDs, subdomainFlag)) 
			{
				leafBlockIDs.insert(i);
				leaves++;
			}
		}
	}

	//Put the fineCells values in the class to be referenced later by the print and deallocation functions.
	fineCells[IAXIS] = fineCellsX;
	fineCells[JAXIS] = fineCellsY;
	fineCells[KAXIS] = fineCellsZ;

	validCells[IAXIS] = validCellsX;
	validCells[JAXIS] = validCellsY;
	validCells[KAXIS] = validCellsZ;

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < getDim(); j++) {
			if (subdomainFlag) {
				refinedDomainBoundBox[i][j] = globalGridBoundBox[LOWER][j] + (subdomainBlockCornerIDs[i][j] * gridDelta[j]);
			}
			else {
				refinedDomainBoundBox[i][j] = globalGridBoundBox[i][j];
			}
		}
	}	

	if (subdomainFlag) {
		totalCells[IAXIS] = subdomainCellsX;
		totalCells[JAXIS] = subdomainCellsY;
		totalCells[KAXIS] = subdomainCellsZ;
	}
	else {
		if (getDim() >= 1) {
			totalCells[IAXIS] = fineBlocksX * getNcells(IAXIS);
			totalCells[JAXIS] = 1;
			totalCells[KAXIS] = 1;
		}
		if (getDim() >= 2) {
			totalCells[JAXIS] = fineBlocksY * getNcells(JAXIS);
			totalCells[KAXIS] = 1;
		}
		if (getDim() >= 3) {
			totalCells[KAXIS] = fineBlocksZ * getNcells(KAXIS);
		}
	}
	
	// std::cout << "Total cells: " << totalCells[IAXIS] << " " << totalCells[JAXIS] << " " << totalCells[KAXIS] << std::endl;
	
	numProcs[IAXIS] = procsX;
	numProcs[JAXIS] = procsY;
	numProcs[KAXIS] = procsZ;

	//Allocate the 4-d data array per processor. [variable][x][y][z]
	// double**** inData;
	inData = new double***[getNUserVars()];	         		 //Number of variables

	for (int i = 0; i < getNUserVars(); i++) {
		inData[i] = new double**[fineCellsX];                //Max cells in x-direction

		for (int j = 0; j < fineCellsX; j++) {
			inData[i][j] = new double*[fineCellsY];          //Max cells in y-direction

			for (int k = 0; k < fineCellsY; k++) {
				inData[i][j][k] = new double[fineCellsZ];    //Max cells in z-direction	
			}
		}
	}

	// Set the data to 0 so that any extra cells not written to
	// are not given garbage values when they are output.
	for (int var = 0; var < getNUserVars(); var++) {
		for (int i = 0; i < fineCellsX; i++) {
			for (int j = 0; j < fineCellsY; j++) {
				for (int k = 0; k < fineCellsZ; k++) {
					//inData[0][i][j][k] = 1e99;
					inData[var][i][j][k] = 0;
				}
			}
		}
	}

	// std::cout << "Proc: " << mpiGetID() << " Leaves: " << leaves << std::endl;

	/* 
		Determine processor index bounds.
	*/

	int *indexMax;
	indexMax = new int[MESH_MDIM];
	procIndexMin = new int*[mpiGetProcs()];
	for (int i = 0; i < mpiGetProcs(); i++) {
		procIndexMin[i] = new int[MESH_MDIM];
	}

	if (getDim() == 1) {
		
		int dxMin = fineCellsX;
		
		//Can't use validCells here because that's only known locally on each processor.
		int procsWithRemaindersX = remainCellsX % procsX;

		procIndexMin[0][IAXIS] = 0;
		// std::cout << "Proc index mins: " << std::endl;
		for (int i = 1; i < procsX; i++) {
			if (i-1 < procsWithRemaindersX || procsWithRemaindersX == 0) {
				procIndexMin[i][IAXIS] = procIndexMin[i-1][IAXIS] + dxMin;
			}
			else {
				procIndexMin[i][IAXIS] = procIndexMin[i-1][IAXIS] + dxMin - 1;
			}
			// if (i-1 < remainCellsX % procsX) {
			// 	procIndexMin[i][IAXIS] += remainCellsX;
			// }
			// std::cout << procIndexMin[i][IAXIS] << std::endl;
		}

		for (int j = 0; j < mpiGetProcs(); j++) {
			procIndexMin[j][JAXIS] = 0;
			procIndexMin[j][KAXIS] = 0;
		}

		if (subdomainFlag) {
			indexMax[IAXIS] = subdomainCellsX;
		}
		else {
			indexMax[IAXIS] = int(((globalGridBoundBox[UPPER][IAXIS] - globalGridBoundBox[LOWER][IAXIS] + gridHalfDelta[IAXIS])
											/ gridDelta[IAXIS]));
		}
	}
	else if (getDim() == 2) {
		
		int dxMin = fineCellsX; //Original fineCellsX with the remainder added.
		int dyMin = fineCellsY;

		int procsWithRemaindersX = remainCellsX % procsX;
		int procsWithRemaindersY = remainCellsY % procsY;

		//0-index corner
		procIndexMin[0][IAXIS] = 0;
		procIndexMin[0][JAXIS] = 0;

		//0-index corners of box
		for (int i = 1; i < procsX; i++) {
			if (i-1 < procsWithRemaindersX || procsWithRemaindersX == 0) {
				procIndexMin[i][IAXIS] = procIndexMin[i-1][IAXIS] + dxMin;
			}
			else {
				procIndexMin[i][IAXIS] = procIndexMin[i-1][IAXIS] + dxMin - 1;
			}
			procIndexMin[i][JAXIS] = procIndexMin[0][JAXIS];
		}
		for (int j = 1; j < procsY; j++) {
			// std::cout << procsX*j << std::endl;
			procIndexMin[procsX*j][IAXIS] = procIndexMin[0][IAXIS];
			if (j-1 < procsWithRemaindersY || procsWithRemaindersY == 0) {
				procIndexMin[procsX*j][JAXIS] = procIndexMin[procsX*(j-1)][JAXIS] + dyMin;
			}
			else {
				procIndexMin[procsX*j][JAXIS] = procIndexMin[procsX*(j-1)][JAXIS] + dyMin - 1;
			}
			// if (j-1 < remainCellsY) {
			// 	procIndexMin[procsX*j][JAXIS] += remainCellsY;
			// }
		}

		//Remainder of box.
		for (int i = 1; i < procsX; i++) {
			for (int j = 1; j < procsY; j++) {
				if (i-1 < procsWithRemaindersX || procsWithRemaindersX == 0) {
					procIndexMin[procsX*j + i][IAXIS] = procIndexMin[procsX*j + (i-1)][IAXIS] + dxMin;
				}
				else {
					procIndexMin[procsX*j + i][IAXIS] = procIndexMin[procsX*j + (i-1)][IAXIS] + dxMin - 1;	
				}
				if (j-1 < procsWithRemaindersY || procsWithRemaindersY == 0) {
					procIndexMin[procsX*j + i][JAXIS] = procIndexMin[procsX*(j-1) + i][JAXIS] + dyMin;
				}
				else {
					procIndexMin[procsX*j + i][JAXIS] = procIndexMin[procsX*(j-1) + i][JAXIS] + dyMin - 1;
				}
				//std::cout << procsX*i + j << ": " << procIndexMin[procsX*i + j][IAXIS] << "," << procIndexMin[procsX*i + j][JAXIS] << std::endl;
			}
		}
		for (int k = 0; k < mpiGetProcs(); k++) {
			procIndexMin[k][KAXIS] = 0;
		}
		if (subdomainFlag) {
			indexMax[IAXIS] = subdomainCellsX;
			indexMax[JAXIS] = subdomainCellsY;
		}
		else {
			indexMax[IAXIS] = int(((globalGridBoundBox[UPPER][IAXIS] - globalGridBoundBox[LOWER][IAXIS] + gridHalfDelta[IAXIS])
											/ gridDelta[IAXIS]));
			indexMax[JAXIS] = int(((globalGridBoundBox[UPPER][JAXIS] - globalGridBoundBox[LOWER][JAXIS] + gridHalfDelta[JAXIS])
											/ gridDelta[JAXIS]));
		}

		// std::cout << "PRINTING ALL INDICES" << std::endl;
		// for (int j = 0; j < procsY; j++) {
		// 	for (int i = 0; i < procsX; i++) {
		// 		std::cout << procIndexMin[procsX*j + i][IAXIS] << "," << procIndexMin[procsX*j + i][JAXIS] << std::endl;
		// 	}
		// }
		// std::cout << indexMax[IAXIS] << "," << indexMax[JAXIS] << std::endl;

	}

	else if (getDim() == 3) {

		int dxMin = fineCellsX; //Original fineCellsX with the remainder added.
		int dyMin = fineCellsY;
		int dzMin = fineCellsZ;

		int procsWithRemaindersX = remainCellsX % procsX;
		int procsWithRemaindersY = remainCellsY % procsY;
		int procsWithRemaindersZ = remainCellsZ % procsZ;

		// std::cout << "Remainders: " << procsWithRemaindersX << " " << procsWithRemaindersY << " " << procsWithRemaindersZ << std::endl;

		//Single 0-index corner
		procIndexMin[0][IAXIS] = 0;
		procIndexMin[0][JAXIS] = 0;
		procIndexMin[0][KAXIS] = 0;

		//0-index corners of box
		// x-dir, y and z are 0
		for (int i = 1; i < procsX; i++) {
			if (i-1 < procsWithRemaindersX || procsWithRemaindersX == 0) {
				procIndexMin[i][IAXIS] = procIndexMin[i-1][IAXIS] + dxMin;
			}
			else {
				procIndexMin[i][IAXIS] = procIndexMin[i-1][IAXIS] + dxMin - 1;
			}
			procIndexMin[i][JAXIS] = procIndexMin[0][JAXIS];
			procIndexMin[i][KAXIS] = procIndexMin[0][KAXIS];
		}
		// y-dir, x and z are 0
		for (int j = 1; j < procsY; j++) {
			procIndexMin[procsX*j][IAXIS] = procIndexMin[0][IAXIS];
			if (j-1 < procsWithRemaindersY || procsWithRemaindersY == 0) {
				procIndexMin[procsX*j][JAXIS] = procIndexMin[procsX*(j-1)][JAXIS] + dyMin;
			}
			else {
				procIndexMin[procsX*j][JAXIS] = procIndexMin[procsX*(j-1)][JAXIS] + dyMin - 1;
			}
			procIndexMin[procsX*j][KAXIS] = procIndexMin[0][KAXIS];
		}
		// z-dir, x and y are 0
		for (int k = 1; k < procsZ; k++) {
			procIndexMin[(procsX*procsY)*k][IAXIS] = procIndexMin[0][IAXIS];
			procIndexMin[(procsX*procsY)*k][JAXIS] = procIndexMin[0][JAXIS];
			if (k-1 < procsWithRemaindersZ || procsWithRemaindersZ == 0) {
				procIndexMin[(procsX*procsY)*k][KAXIS] = procIndexMin[(procsX * procsY)*(k-1)][KAXIS] + dzMin;
			}
			else {
				procIndexMin[(procsX*procsY)*k][KAXIS] = procIndexMin[(procsX * procsY)*(k-1)][KAXIS] + dzMin - 1;
			}
		}

		//0-index faces of box
		// x-y face
		for (int i = 1; i < procsX; i++) {
			for (int j = 1; j < procsY; j++) {
				if (i-1 < procsWithRemaindersX || procsWithRemaindersX == 0) {
					procIndexMin[procsX*j + i][IAXIS] = procIndexMin[procsX*j + (i-1)][IAXIS] + dxMin;
				}
				else {
					procIndexMin[procsX*j + i][IAXIS] = procIndexMin[procsX*j + (i-1)][IAXIS] + dxMin - 1;
				}
				if (j-1 < procsWithRemaindersY || procsWithRemaindersY == 0) {
					procIndexMin[procsX*j + i][JAXIS] = procIndexMin[procsX*(j-1) + i][JAXIS] + dyMin;
				}
				else {
					procIndexMin[procsX*j + i][JAXIS] = procIndexMin[procsX*(j-1) + i][JAXIS] + dyMin - 1;
				}
				procIndexMin[procsX*j + i][KAXIS] = procIndexMin[0][KAXIS];
			}
		}

		// y-z face
		for (int j = 1; j < procsY; j++) {
			for (int k = 1; k < procsZ; k++) {
				procIndexMin[(procsX*procsY)*k + procsX*j][IAXIS] = procIndexMin[0][IAXIS];
				if (j-1 < procsWithRemaindersY || procsWithRemaindersY == 0) {
					procIndexMin[(procsX*procsY)*k + procsX*j][JAXIS] = procIndexMin[(procsX*procsY)*k + procsX*(j-1)][JAXIS] + dyMin;
				}
				else {
					procIndexMin[(procsX*procsY)*k + procsX*j][JAXIS] = procIndexMin[(procsX*procsY)*k + procsX*(j-1)][JAXIS] + dyMin - 1;
				}
				if (k-1 < procsWithRemaindersZ || procsWithRemaindersZ == 0){
					procIndexMin[(procsX*procsY)*k + procsX*j][KAXIS] = procIndexMin[(procsX*procsY)*(k-1) + procsX*j][KAXIS] + dzMin;
				}
				else {
					procIndexMin[(procsX*procsY)*k + procsX*j][KAXIS] = procIndexMin[(procsX*procsY)*(k-1) + procsX*j][KAXIS] + dzMin - 1;
				}
			}
		}
		
		//x-z face
		for (int i = 1; i < procsX; i++) {
			for (int k = 1; k < procsZ; k++) {
				if (i-1 < procsWithRemaindersX || procsWithRemaindersX == 0) {
					procIndexMin[(procsX*procsY)*k + i][IAXIS] = procIndexMin[(procsX*procsY)*k + (i-1)][IAXIS] + dxMin;
				}
				else {
					procIndexMin[(procsX*procsY)*k + i][IAXIS] = procIndexMin[(procsX*procsY)*k + (i-1)][IAXIS] + dxMin - 1;
				}
				procIndexMin[(procsX*procsY)*k + i][JAXIS] = procIndexMin[0][JAXIS];
				if (k-1 < procsWithRemaindersZ || procsWithRemaindersZ == 0) {
					procIndexMin[(procsX*procsY)*k + i][KAXIS] = procIndexMin[(procsX*procsY)*(k-1) + i][KAXIS] + dzMin;
				}
				else {
					procIndexMin[(procsX*procsY)*k + i][KAXIS] = procIndexMin[(procsX*procsY)*(k-1) + i][KAXIS] + dzMin - 1;
				}
			}
		}

		//Remainder of box
		for (int i = 1; i < procsX; i++) {
			for (int j = 1; j < procsY; j++) {
				for (int k = 1; k < procsZ; k++) {
					if (i-1 < procsWithRemaindersX || procsWithRemaindersX == 0) {
						procIndexMin[(procsX*procsY)*k + procsX*j + i][IAXIS] = procIndexMin[(procsX*procsY)*k + procsX*j + (i-1)][IAXIS] + dxMin;
					}
					else {
						procIndexMin[(procsX*procsY)*k + procsX*j + i][IAXIS] = procIndexMin[(procsX*procsY)*k + procsX*j + (i-1)][IAXIS] + dxMin - 1;
					}
					if (j-1 < procsWithRemaindersY || procsWithRemaindersY == 0) {
						procIndexMin[(procsX*procsY)*k + procsX*j + i][JAXIS] = procIndexMin[(procsX*procsY)*k + procsX*(j-1) + i][JAXIS] + dyMin;
					}
					else {
						procIndexMin[(procsX*procsY)*k + procsX*j + i][JAXIS] = procIndexMin[(procsX*procsY)*k + procsX*(j-1) + i][JAXIS] + dyMin - 1;
					}
					if (k-1 < procsWithRemaindersZ || procsWithRemaindersZ == 0) {
						procIndexMin[(procsX*procsY)*k + procsX*j + i][KAXIS] = procIndexMin[(procsX*procsY)*(k-1) + procsX*j + i][KAXIS] + dzMin;
					}
					else {
						procIndexMin[(procsX*procsY)*k + procsX*j + i][KAXIS] = procIndexMin[(procsX*procsY)*(k-1) + procsX*j + i][KAXIS] + dzMin - 1;
					}
				}
			}
		}

		if (subdomainFlag) {
			indexMax[IAXIS] = subdomainCellsX;
			indexMax[JAXIS] = subdomainCellsY;
			indexMax[KAXIS] = subdomainCellsZ;
		}
		else {
			indexMax[IAXIS] = int(((globalGridBoundBox[UPPER][IAXIS] - globalGridBoundBox[LOWER][IAXIS] + gridHalfDelta[IAXIS])
											/ gridDelta[IAXIS]));
			indexMax[JAXIS] = int(((globalGridBoundBox[UPPER][JAXIS] - globalGridBoundBox[LOWER][JAXIS] + gridHalfDelta[JAXIS])
											/ gridDelta[JAXIS]));
			indexMax[KAXIS] = int(((globalGridBoundBox[UPPER][KAXIS] - globalGridBoundBox[LOWER][KAXIS] + gridHalfDelta[KAXIS])
											/ gridDelta[KAXIS]));
		}

		// std::cout << "PRINTING ALL INDICES" << std::endl;
		// for (int k = 0; k < procsZ; k++) {
		// 	for (int j = 0; j < procsY; j++) {
		// 		for (int i = 0; i < procsX; i++) {
		// 			std::cout << procIndexMin[(procsX*procsY)*k + procsX*j + i][IAXIS] << "," << procIndexMin[(procsX*procsY)*k + procsX*j + i][JAXIS] << "," << procIndexMin[(procsX*procsY)*k + procsX*j + i][KAXIS] << std::endl;
		// 		}
		// 	}
		// }
		// std::cout << indexMax[IAXIS] << "," << indexMax[JAXIS] << "," << indexMax[KAXIS] << std::endl;

	}

	MPI_Datatype mpiMessageStruct, oldtypes[2];
	int          blockcounts[2];
	MPI_Aint     offsets[2], extent;

	offsets[0] = 0;
	oldtypes[0] = MPI_DOUBLE;
	blockcounts[0] = 1;

	MPI_Type_extent(MPI_DOUBLE, &extent);
	offsets[1] = 1 * extent;
	oldtypes[1] = MPI_INT;
	blockcounts[1] = 3;

	/* Now define structured type and commit it */
	MPI_Type_struct(2, blockcounts, offsets, oldtypes, &mpiMessageStruct);
	MPI_Type_commit(&mpiMessageStruct);

	//Iterate over leaf block
	mpiVarStruct **mpiSendVars; //Procs by max refine
	mpiVarStruct *mpiReceiveVars;

	//int transferLevelDifference = pow(2,globalMaxRefineLevel-1);
	//int transferMaxSize = pow(transferLevelDifference,MESH_MDIM);
	int transferMaxSize = fineCellsX*fineCellsY*fineCellsZ;

	//Take the level difference to find the maximum possible refinement for send/receive arrays.
	mpiSendVars = new mpiVarStruct*[mpiGetProcs()];
	for (int mpiProc = 0; mpiProc < mpiGetProcs(); mpiProc++) {
		if (mpiProc != mpiGetID()) { //No need to allocate for own processor.
			mpiSendVars[mpiProc] = new mpiVarStruct[transferMaxSize];
		}
	}
	mpiReceiveVars = new mpiVarStruct[transferMaxSize];

	int *mpiSendIterator;
	mpiSendIterator = new int[mpiGetProcs()];

	int *completionArray;
	completionArray = new int[mpiGetProcs()];

	MPI_Request *mpiRequests;
	mpiRequests = new MPI_Request[mpiGetProcs()];

	
	int resetLeaves = leaves;
	int doneRefining;
	int leafStatus;

	int dataTag = 0;
	int communicationTag = 1;

	for (int var = 0; var < getNUserVars(); var++) {
		if (mpiGetID() == 0) {
			std::cout << "Refining variable: " << findVarName(var) << std::endl;
		}
		
		doneRefining = 0;
		
		leaves = resetLeaves;

		leafStatus = 0;
		if (leaves == 0) {
			leafStatus = 1;
		}

		MPI_Allgather(&leafStatus, 1, MPI_INT, completionArray, 1, MPI_INT, MPI_COMM_WORLD);		
		
		// for (int i = 0; i < mpiGetProcs(); i++) {
		// 	std::cout << completionArray[i] << " ";
		// }
		// std::cout << std::endl;

		while (doneRefining == 0) {
			if (leaves > 0) {
				for (leaf = leafBlockIDs.begin(); leaf != leafBlockIDs.end(); ++leaf) {
					leaves--;
					// std::cout << "Var: " << findVarName(var) << " Remaining leaves: " << leaves << std::endl;
					//Find the lower left finest block index.
					double domainBounds[2][MESH_MDIM];
					getDomainBoundBox(domainBounds);

					double blockBounds[2][MESH_MDIM];
					getBlkBoundBox(*leaf, blockBounds);

					int indX, indY, indZ;
					indX = localBlockCornerIDs[*leaf][IAXIS];
					indY = 0;
					indZ = 0;
					if (getDim() > 1) {
						indY = localBlockCornerIDs[*leaf][JAXIS];
					}
					if (getDim() > 2) {
						indZ = localBlockCornerIDs[*leaf][KAXIS];
					}

					int offBlkX = indX;
					int offBlkY = indY;
					int offBlkZ = indZ;

					for (int mpiProc = 0; mpiProc < mpiGetProcs(); mpiProc++) {
						mpiSendIterator[mpiProc] = 0;
					}

					//Find the difference in refine levels.
					int levelDiff = globalMaxRefineLevel - getRefineLevel(*leaf);
					int scale = (int)pow(2,levelDiff);

					for (int i = 0; i < getNcells(IAXIS); i++) {
					 	for (int j = 0; j < getNcells(JAXIS); j++) {
					 		for (int k = 0; k < getNcells(KAXIS); k++) {

								int iiMin, iiMax, jjMin, jjMax, kkMin, kkMax;
								iiMin = i*scale;
								iiMax = (i+1)*scale;
								jjMin = 0;
								jjMax = 1;
								kkMin = 0;
								kkMax = 1;
								if (getDim() > 1) {
									jjMin = j*scale;
									jjMax = (j+1)*scale;
								}
								if (getDim() > 2) {
									kkMin = k*scale;
									kkMax = (k+1)*scale;
								}

								for (int ii = iiMin; ii < iiMax; ii++) {
									for (int jj = jjMin; jj < jjMax; jj++) {
										for (int kk = kkMin; kk < kkMax; kk++) {
											// if (offBlkX+ii >= fineCellsX) {
												// cout << mpiGetID() << std::endl;
												// cout << leafBlockIDs.size() << endl;
												// cout << *leaf << endl;
												// cout << i << " " << j << " " << k << " " << ii << " " << jj << endl;
												// cout << offBlkX+ii << " " << offBlkY+jj << " " << offBlkZ+kk << endl;	
												// cout << offBlkZ << " " << kk << std::endl;
											//}
											

											int cellProc;
											bool containedInSubdomainFlag = false;

											if (subdomainFlag) {
												containedInSubdomainFlag = inSubdomain(offBlkX+ii, offBlkY+jj, offBlkZ+kk, subdomainBlockCornerIDs);
											}

											findDestinationProcessor(procsX, procsY, procsZ, 
																	 offBlkX, offBlkY, offBlkZ,
																	 ii, jj, kk,
																	 getDim(), procIndexMin,
																	 subdomainFlag, containedInSubdomainFlag, subdomainBlockCornerIDs,
																	 cellProc);

											if (cellProc == mpiGetID()) {
												// std::cout << "Keeping this one" << std::endl;
												// std::cout << offBlkZ + kk - subdomainBlockCornerIDs[LOWER][KAXIS] << std::endl;
												if (subdomainFlag) {
													// if (offBlkZ + kk - subdomainBlockCornerIDs[LOWER][KAXIS] == 72) {
														// std::cout << getVariable(*leaf, var, i, j, k) << std::endl;
													// }
													if (containedInSubdomainFlag) {

														inData[var][offBlkX + ii - subdomainBlockCornerIDs[LOWER][IAXIS] - procIndexMin[cellProc][IAXIS]]
																   [offBlkY + jj - subdomainBlockCornerIDs[LOWER][JAXIS] - procIndexMin[cellProc][JAXIS]]
																   [offBlkZ + kk - subdomainBlockCornerIDs[LOWER][KAXIS] - procIndexMin[cellProc][KAXIS]]
														= getVariable(*leaf, var, i, j, k);
													}
												}
												else {
													inData[var]
														  [offBlkX + ii - procIndexMin[cellProc][IAXIS]]
														  [offBlkY + jj - procIndexMin[cellProc][JAXIS]]
														  [offBlkZ + kk - procIndexMin[cellProc][KAXIS]]
													= getVariable(*leaf, var, i, j, k);
												}
											}
											else {
												//MPI stuff
												//Send the indices, send the variable.
												// std::cout << offBlkZ + kk - subdomainBlockCornerIDs[LOWER][KAXIS] << std::endl;
												if (subdomainFlag) {
													// if (offBlkZ + kk - subdomainBlockCornerIDs[LOWER][KAXIS] == 72) {
														// std::cout << getVariable(*leaf, var, i, j, k) << std::endl;
													// }
													if (containedInSubdomainFlag) {

														mpiSendVars[cellProc][mpiSendIterator[cellProc]].value = getVariable(*leaf, var, i, j, k);
														mpiSendVars[cellProc][mpiSendIterator[cellProc]].x = offBlkX + ii - subdomainBlockCornerIDs[LOWER][IAXIS] - procIndexMin[cellProc][IAXIS];
														mpiSendVars[cellProc][mpiSendIterator[cellProc]].y = offBlkY + jj - subdomainBlockCornerIDs[LOWER][JAXIS] - procIndexMin[cellProc][JAXIS];
														mpiSendVars[cellProc][mpiSendIterator[cellProc]].z = offBlkZ + kk - subdomainBlockCornerIDs[LOWER][KAXIS] - procIndexMin[cellProc][KAXIS];

														mpiSendIterator[cellProc]++;
													}
												}
												else {
													mpiSendVars[cellProc][mpiSendIterator[cellProc]].value = getVariable(*leaf, var, i, j, k);
													mpiSendVars[cellProc][mpiSendIterator[cellProc]].x = offBlkX + ii - procIndexMin[cellProc][IAXIS];
													mpiSendVars[cellProc][mpiSendIterator[cellProc]].y = offBlkY + jj - procIndexMin[cellProc][JAXIS];
													mpiSendVars[cellProc][mpiSendIterator[cellProc]].z = offBlkZ + kk - procIndexMin[cellProc][KAXIS];

													mpiSendIterator[cellProc]++;
												}
											}
										} //kk
									} //jj
								}//End of refined cells for loop ii
							} //k
						} //j
					} //i End of leaf var block

					for (int mpiProc = 0; mpiProc < mpiGetProcs(); mpiProc++) {
						if (mpiProc != mpiGetID() && mpiSendIterator[mpiProc] != 0) {
							//Send off data to the correct processor.
							MPI_Isend(&(mpiSendVars[mpiProc][0]), mpiSendIterator[mpiProc], mpiMessageStruct, mpiProc, dataTag, MPI_COMM_WORLD, &mpiRequests[mpiProc]);
							// std::cout << "Sending from " << mpiGetID() << " to " << mpiProc << " Count: " << mpiSendIterator[mpiProc] << std::endl;
						}
					}

					int flag;
					MPI_Status status;
					int messageLength;
					int sendCompleted = 0;
					int sendCompletedProc;
					while (sendCompleted == 0){ //Receive while waiting for the send to go through.
						for (int mpiProc = 0; mpiProc < mpiGetProcs(); mpiProc++) {
							if (mpiProc != mpiGetID()){
								MPI_Iprobe(mpiProc, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
								if (flag == true && status.MPI_TAG == dataTag) {
									MPI_Get_count(&status, mpiMessageStruct, &messageLength);
									MPI_Recv(&(mpiReceiveVars[0]), messageLength, mpiMessageStruct, mpiProc, dataTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
									for (int i = 0; i < messageLength; i++) {
										inData[var]
											  [mpiReceiveVars[i].x]
											  [mpiReceiveVars[i].y]
											  [mpiReceiveVars[i].z]
										= mpiReceiveVars[i].value;
									}
									// std::cout << "Processor " << mpiGetID() << " Receiving from " << mpiProc << " Count: " << messageLength << std::endl;
								}
								if (flag == true && status.MPI_TAG == communicationTag) {
									MPI_Recv(&completionArray[mpiProc], 1, MPI_INT, mpiProc, communicationTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								}
							}
						}
						sendCompleted = 1;
						for (int mpiProc = 0; mpiProc < mpiGetProcs(); mpiProc++) {
							if (mpiProc != mpiGetID() && mpiSendIterator[mpiProc] != 0) {
								// std::cout << "Proc: " << mpiProc << "Send iterator: " << mpiSendIterator[mpiProc] << std::endl;
								MPI_Test(&(mpiRequests[mpiProc]), &sendCompletedProc, MPI_STATUS_IGNORE);
								sendCompleted = std::min(sendCompleted, sendCompletedProc);
							}
						}
					}
					
					if (leaves == 0) {
						for (int mpiProc = 0; mpiProc < mpiGetProcs(); mpiProc++) {
							completionArray[mpiGetID()] = 1;
							if (mpiProc != mpiGetID()) {
								leafStatus = 1;
								MPI_Request commRequest;
								MPI_Isend(&leafStatus, 1, MPI_INT, mpiProc, communicationTag, MPI_COMM_WORLD, &commRequest);
								MPI_Wait(&commRequest, MPI_STATUS_IGNORE);
							}
						}
					}
				}
			}

			int flag;
			MPI_Status status;
			int messageLength;
			for (int mpiProc = 0; mpiProc < mpiGetProcs(); mpiProc++) {
				if (mpiProc != mpiGetID()) {
					MPI_Iprobe(mpiProc, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
					if (flag == true && status.MPI_TAG == dataTag) {
						MPI_Get_count(&status, mpiMessageStruct, &messageLength);
						MPI_Recv(&(mpiReceiveVars[0]), messageLength, mpiMessageStruct, mpiProc, dataTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						for (int i = 0; i < messageLength; i++) {
							inData[var]
								  [mpiReceiveVars[i].x]
								  [mpiReceiveVars[i].y]
								  [mpiReceiveVars[i].z] 
							= mpiReceiveVars[i].value;
						}
						// std::cout << "Processor " << mpiGetID() << " Receiving from " << mpiProc << " Count: " << messageLength << std::endl;
					}
					if (flag == true && status.MPI_TAG == communicationTag) {
						MPI_Recv(&completionArray[mpiProc], 1, MPI_INT, mpiProc, communicationTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						// std::cout << "Refinement status: " << doneRefining << std::endl;
					}
				}
			}

			doneRefining = 1;
			for (int i = 0; i < mpiGetProcs(); i++) {
				doneRefining = std::min(doneRefining, completionArray[i]);
				// std::cout << completionArray[i] << " ";
			}
			// std::cout << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);

	}

	MPI_Barrier(MPI_COMM_WORLD);

	for (int i = 0; i < getLocalNumBlocks(); i++) {
		delete[] localBlockCornerIDs[i];
	}
	delete[] localBlockCornerIDs;

	delete[] mpiRequests;
	for (int i = 0; i < mpiGetProcs(); i++) {
		if (i != mpiGetID()) {
			delete[] mpiSendVars[i];
		}
	}
	delete[] mpiSendVars;
	delete[] mpiReceiveVars;
	delete[] mpiSendIterator;

	delete[] indexMax;
	delete[] completionArray;
	
	if (mpiGetID() == 0) {
		std::cout << "Refinement completed." << std::endl;
	}
}

int compare (const void * a, const void * b){
  return ( *(int*)b - *(int*)a );
}

void findFactors(int *primeFactors, int procs, int &numFactors) {

	for (int i = 0; i < procs; i++) {
		primeFactors[i] = 0;
	}
	primeFactors[0] = procs;

	int factoring = 0;
	int factorCount = 1;
	int currentFactor;
	int i = 2;

	bool factorized = false;
	bool foundFactor = false;

	if (procs == 1) {
		factorized = true;
	}
	
	while (factorized == false) {
		// Factor and sweep until factorized
		currentFactor = primeFactors[factoring];
		while (i < currentFactor && !foundFactor) {
			if (currentFactor % i == 0){
				//Replace the value we're factoring with the divided amount
				primeFactors[factoring] = currentFactor / i;

				//Put the other factor at the end of the array.
				primeFactors[factorCount] = i;
				factorCount++;
				foundFactor = true;
			}
			i++;
		}
		i = 2;
		foundFactor = false;
		if (primeFactors[factoring] == currentFactor) {
			//Didn't change. Go to the next value.
			factoring++;
		}
		if (primeFactors[factoring] == 0) {
			//Reached the end of the list. Assumably it's factorized now.
			factorized = true;
		}
		// std::cout << factorized << std::endl;
	}

	qsort (primeFactors, factoring, sizeof(int), compare);
	numFactors = factoring;
}

bool intersectsSubdomain(int localBlockCornerIDs[MESH_MDIM], int blockMaxIndexExtent[MESH_MDIM], int subdomainBlockCornerIDs[2][MESH_MDIM], bool flag) {
	// Checks to see if a leaf block intersects the subdomain.

	// dfenn note: I have completely rewritten this routine
	// the old routine checked if the subdomain crossed the block's boundaries. I don't think this is conceptually correct.
	// what the routine should really do is check if any part of the subdomain is present on the block

	// the old method didn't work if the subdomain was contained exclusively by one block
	// the new method works in this case

	// I also don't understand the need for separate arrays for the block's upper and lower bounds, so I've combined them into one
	// a better solution would be to combine them outside this routine, but that would require more widespread changes I'm not willing to undergo

	// If the flag is set to false, there is no subdomain. Automatically return true in this case.
	if (!flag) {
		return true;
	}

	int blockIndexBounds[2][MESH_MDIM];
	for (int i=0; i<2; i++)
	{
		for (int j=0; j<MESH_MDIM; j++)
		{
			if (i == LOWER)
				blockIndexBounds[i][j] = localBlockCornerIDs[j];
			else
				blockIndexBounds[i][j] = blockMaxIndexExtent[j];

		}
	}

	// int nDim = getDim();
	
	// I haven't tested this in 2D -- dfenn
	bool intersects = false;
	for (int dim=0; dim < MESH_MDIM; dim++)
	{
		bool intersects = false;
		if (subdomainBlockCornerIDs[LOWER][dim] <= blockIndexBounds[UPPER][dim] &&
			subdomainBlockCornerIDs[UPPER][dim] >= blockIndexBounds[LOWER][dim])
		{
			intersects = true;
		}
		if (!intersects)
			return false;
	}

	return true;

	// Otherwise, see if any part of the local block intersects with the subdomain. If it does, return true.



	// /* 3D: "Near" side of box */
	// //2D: Top left corner of subdomain box
	// if (   blockMaxIndexExtent[IAXIS] >= subdomainBlockCornerIDs[LOWER][IAXIS]
	// 	&& blockMaxIndexExtent[IAXIS] <= subdomainBlockCornerIDs[UPPER][IAXIS]
	// 	&& localBlockCornerIDs[JAXIS] >= subdomainBlockCornerIDs[LOWER][JAXIS]
	// 	&& localBlockCornerIDs[JAXIS] <= subdomainBlockCornerIDs[UPPER][JAXIS]
	// 	&& blockMaxIndexExtent[KAXIS] >= subdomainBlockCornerIDs[LOWER][KAXIS]
	// 	&& blockMaxIndexExtent[KAXIS] <= subdomainBlockCornerIDs[UPPER][KAXIS]
	// 	){
	// 	return true;
	// }
	// //2D: Top right corner of subdomain box
	// else if (   localBlockCornerIDs[IAXIS] >= subdomainBlockCornerIDs[LOWER][IAXIS]
	// 		 && localBlockCornerIDs[IAXIS] <= subdomainBlockCornerIDs[UPPER][IAXIS]
	// 		 && localBlockCornerIDs[JAXIS] >= subdomainBlockCornerIDs[LOWER][JAXIS]
	// 		 && localBlockCornerIDs[JAXIS] <= subdomainBlockCornerIDs[UPPER][JAXIS]
	// 		 && blockMaxIndexExtent[KAXIS] >= subdomainBlockCornerIDs[LOWER][KAXIS]
	// 		 && blockMaxIndexExtent[KAXIS] <= subdomainBlockCornerIDs[UPPER][KAXIS]
	// 		){
	// 	return true;
	// }
	// //2D: Bottom left corner of subdomain box
	// else if (   blockMaxIndexExtent[IAXIS] >= subdomainBlockCornerIDs[LOWER][IAXIS]
	// 		 && blockMaxIndexExtent[IAXIS] <= subdomainBlockCornerIDs[UPPER][IAXIS]
	// 		 && blockMaxIndexExtent[JAXIS] >= subdomainBlockCornerIDs[LOWER][JAXIS]
	// 		 && blockMaxIndexExtent[JAXIS] <= subdomainBlockCornerIDs[UPPER][JAXIS]
	// 		 && blockMaxIndexExtent[KAXIS] >= subdomainBlockCornerIDs[LOWER][KAXIS]
	// 		 && blockMaxIndexExtent[KAXIS] <= subdomainBlockCornerIDs[UPPER][KAXIS]
	// 		){
	// 	return true;
	// }
	// //2D: Bottom right corner of subdomain box
	// else if (   localBlockCornerIDs[IAXIS] >= subdomainBlockCornerIDs[LOWER][IAXIS]
	// 		 && localBlockCornerIDs[IAXIS] <= subdomainBlockCornerIDs[UPPER][IAXIS]
	// 		 && blockMaxIndexExtent[JAXIS] >= subdomainBlockCornerIDs[LOWER][JAXIS]
	// 		 && blockMaxIndexExtent[JAXIS] <= subdomainBlockCornerIDs[UPPER][JAXIS]
	// 		 && blockMaxIndexExtent[KAXIS] >= subdomainBlockCornerIDs[LOWER][KAXIS]
	// 		 && blockMaxIndexExtent[KAXIS] <= subdomainBlockCornerIDs[UPPER][KAXIS]
	// 		){
	// 	return true;
	// }
	// //3D: "Far" side of box 
	// //2D: Top left corner of subdomain box
	// else if (   blockMaxIndexExtent[IAXIS] >= subdomainBlockCornerIDs[LOWER][IAXIS]
	// 		&& blockMaxIndexExtent[IAXIS] <= subdomainBlockCornerIDs[UPPER][IAXIS]
	// 		&& localBlockCornerIDs[JAXIS] >= subdomainBlockCornerIDs[LOWER][JAXIS]
	// 		&& localBlockCornerIDs[JAXIS] <= subdomainBlockCornerIDs[UPPER][JAXIS]
	// 		&& localBlockCornerIDs[KAXIS] >= subdomainBlockCornerIDs[LOWER][KAXIS]
	// 		&& localBlockCornerIDs[KAXIS] <= subdomainBlockCornerIDs[UPPER][KAXIS]
	// 		){
	// 	return true;
	// }
	// //2D: Top right corner of subdomain box
	// else if (   localBlockCornerIDs[IAXIS] >= subdomainBlockCornerIDs[LOWER][IAXIS]
	// 		 && localBlockCornerIDs[IAXIS] <= subdomainBlockCornerIDs[UPPER][IAXIS]
	// 		 && localBlockCornerIDs[JAXIS] >= subdomainBlockCornerIDs[LOWER][JAXIS]
	// 		 && localBlockCornerIDs[JAXIS] <= subdomainBlockCornerIDs[UPPER][JAXIS]
	// 		 && localBlockCornerIDs[KAXIS] >= subdomainBlockCornerIDs[LOWER][KAXIS]
	// 		 && localBlockCornerIDs[KAXIS] <= subdomainBlockCornerIDs[UPPER][KAXIS]
	// 		 ){
	// 	return true;
	// }
	// //2D: Bottom left corner of subdomain box
	// else if (   blockMaxIndexExtent[IAXIS] >= subdomainBlockCornerIDs[LOWER][IAXIS]
	// 		 && blockMaxIndexExtent[IAXIS] <= subdomainBlockCornerIDs[UPPER][IAXIS]
	// 		 && blockMaxIndexExtent[JAXIS] >= subdomainBlockCornerIDs[LOWER][JAXIS]
	// 		 && blockMaxIndexExtent[JAXIS] <= subdomainBlockCornerIDs[UPPER][JAXIS]
	// 		 && localBlockCornerIDs[KAXIS] >= subdomainBlockCornerIDs[LOWER][KAXIS]
	// 		 && localBlockCornerIDs[KAXIS] <= subdomainBlockCornerIDs[UPPER][KAXIS]
	// 		 ){
	// 	return true;
	// }
	// //2D: Bottom right corner of subdomain box
	// else if (   localBlockCornerIDs[IAXIS] >= subdomainBlockCornerIDs[LOWER][IAXIS]
	// 		 && localBlockCornerIDs[IAXIS] <= subdomainBlockCornerIDs[UPPER][IAXIS]
	// 		 && blockMaxIndexExtent[JAXIS] >= subdomainBlockCornerIDs[LOWER][JAXIS]
	// 		 && blockMaxIndexExtent[JAXIS] <= subdomainBlockCornerIDs[UPPER][JAXIS]
	// 		 && localBlockCornerIDs[KAXIS] >= subdomainBlockCornerIDs[LOWER][KAXIS]
	// 		 && localBlockCornerIDs[KAXIS] <= subdomainBlockCornerIDs[UPPER][KAXIS]){
	// 	return true;
	// }
	// else {
	// 	return false;
	// }
}

void findDestinationProcessor(int procsX, int procsY, int procsZ, 
							  int offBlkX, int offBlkY, int offBlkZ,
							  int ii, int jj, int kk,
							  int dim, int *procIndexMin[MESH_MDIM],
							  bool subdomainFlag, bool containedInSubdomainFlag, int subdomainBlockCornerIDs[2][MESH_MDIM],
							  int &destinationProcessor) {
	/* Finds which processor a cell needs to be sent to */

	int cellProcX = procsX-1;
	int cellProcY = 0;
	int cellProcZ = 0;

	for (int cellsX = procsX-1; cellsX >= 0; cellsX--) {
		if (subdomainFlag && containedInSubdomainFlag) {
			if (offBlkX+ii-subdomainBlockCornerIDs[LOWER][IAXIS] < procIndexMin[cellsX][IAXIS]) {
				cellProcX--;
			}
		}
		else {
			if (offBlkX+ii < procIndexMin[cellsX][IAXIS]) {
				cellProcX--;
			}	
		}	
	}
	if (dim > 1) {
		cellProcY = procsY-1;
		for (int cellsY = procsY-1; cellsY >= 0; cellsY--) {
			if (subdomainFlag && containedInSubdomainFlag) {
				if (offBlkY+jj-subdomainBlockCornerIDs[LOWER][JAXIS] < procIndexMin[procsX * cellsY + cellProcX][JAXIS]){
					cellProcY--;
				}
			}
			else {
				if (offBlkY+jj < procIndexMin[procsX * cellsY + cellProcX][JAXIS]){
					cellProcY--;
				}
			}
			
		}
	}
	if (dim > 2) {
		cellProcZ = procsZ-1;
		for (int cellsZ = procsZ-1; cellsZ >= 0; cellsZ--) {
			if (subdomainFlag && containedInSubdomainFlag) {
				if (offBlkZ+kk-subdomainBlockCornerIDs[LOWER][KAXIS] < procIndexMin[(procsX*procsY)*cellsZ + procsX * cellProcY + cellProcX][KAXIS]){
					cellProcZ--;
				}	
			}
			else {
				if (offBlkZ+kk < procIndexMin[(procsX*procsY)*cellsZ + procsX * cellProcY + cellProcX][KAXIS]){
					cellProcZ--;
				}	
			}
			
		}
	}
	destinationProcessor = (procsX*procsY) * cellProcZ + procsX * cellProcY + cellProcX;
}

bool inSubdomain(int i, int j, int k, int subdomainBlockCornerIDs[2][MESH_MDIM]) {
	//Checks to see if an index is within the subdomain indices.
	
	// std::cout << "Indices: " << std::endl;
	// std::cout << "  " << subdomainBlockCornerIDs[LOWER][IAXIS] << " " << subdomainBlockCornerIDs[UPPER][IAXIS] << std::endl;
	// std::cout << "  " << subdomainBlockCornerIDs[LOWER][JAXIS] << " " << subdomainBlockCornerIDs[UPPER][JAXIS] << std::endl;
	// std::cout << "  " << subdomainBlockCornerIDs[LOWER][KAXIS] << " " << subdomainBlockCornerIDs[UPPER][KAXIS] << std::endl;
	if (   i >= subdomainBlockCornerIDs[LOWER][IAXIS] && i <= subdomainBlockCornerIDs[UPPER][IAXIS]
		&& j >= subdomainBlockCornerIDs[LOWER][JAXIS] && j <= subdomainBlockCornerIDs[UPPER][JAXIS]
		&& k >= subdomainBlockCornerIDs[LOWER][KAXIS] && k <= subdomainBlockCornerIDs[UPPER][KAXIS]){
		return true;
	}
	else {
		return false;
	}
}

void allocateVarArray(double****& inData, int numVars, int *nCellsVec) {

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

void deallocateVarArray(double****& inData, int numVars, int *nCellsVec) {

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
	else {
		if (FORCEDEALLOCATECRASH) {
			std::cout << "[flash_reader] Failure in deallocateVarArray function outside of class." << std::endl;
			exit(EXIT_FAILURE);	
		}
	}
}

void FlashAmrMesh::printRefinedData(double**** outData, const char *outFileName, const char *fileType, const char *var) {
	//Temporary print function used to print the first variable stored
	// in a data array generated by refineToFinest() to a file named fileName.
	// Also only prints 2-d. 

	std::ofstream tempFile;
	std::string binaryStr = "binary";
	std::string asciiStr  = "ascii";
	std::string hdf5Str   = "hdf5";
	std::string inString  = fileType;

	// std::cout << "Dimension: " << getDim() << std::endl;

	int varIndex = findVarIndex(var);

  	// if (getNUserVars() > 1) {
  	// 	cout << "Currently cannot write more than one variable at a time. " << endl;
  	// 	cout << "Nothing has been written to the file. " << endl;
  	// 	exit(EXIT_FAILURE);
  	// }

  	//Transform input string to lowercase.
	std::transform(inString.begin(), inString.end(), inString.begin(), ::tolower);

	if (inString == binaryStr) {
		std::cout << "Writing data to binary file." << std::endl;
		tempFile.open(outFileName, std::ios::out | std::ios::binary );

		// No idea what this is - TAH
        // std::ofstream fs(outFileName, std::ios::out | std::ios::binary | std::ios::app);
		
		if (tempFile.is_open()) {

			// Write the header lines
			// Format is Height (vertical), Width (horizontal), Depth (z-slices)
			tempFile.write(reinterpret_cast<char*>(&fineCells[JAXIS]), sizeof(fineCells[JAXIS]));
			tempFile.write(reinterpret_cast<char*>(&fineCells[IAXIS]), sizeof(fineCells[IAXIS]));
			tempFile.write(reinterpret_cast<char*>(&fineCells[KAXIS]), sizeof(fineCells[KAXIS]));

			for (int k = 0; k < fineCells[KAXIS]; k++) {
				for (int j = 0; j < fineCells[JAXIS]; j++) {
					for (int i = 0; i < fineCells[IAXIS]; i++) {
						tempFile.write(reinterpret_cast<char*>(&outData[varIndex][i][j][k]),sizeof(outData[varIndex][i][j][k]));
					}
				}
	  		}		
		} 
			

        tempFile.close();


	}
	else if (inString == asciiStr) {
		std::cout << "Writing data to ascii file." << std::endl;
		tempFile.open(outFileName, std::ios::out);

		//Make sure file opened properly
		if (tempFile.is_open()) {
			// Write the header lines
			// Format is Height (vertical), Width (horizontal), Depth (z-slices)
			//tempFile<<fineCells[JAXIS]<<endl;
			//tempFile<<fineCells[IAXIS]<<endl;
			//tempFile<<fineCells[KAXIS]<<endl;

			//for (int k = 0; k < fineCells[KAXIS]; k++) {
				for (int j = 0; j < fineCells[JAXIS]; j++) {
					for (int i = 0; i < fineCells[IAXIS]; i++) {
						tempFile << outData[varIndex][i][j][100] << " ";
					}
					tempFile << std::endl;
				}
				tempFile << std::endl;
	  		//}		
		}
	}
	else {
		std::cout << "File type not recognized. Must be 'ascii' or 'binary'." << std::endl;
		std::cout << "Nothing has been written to the file. " << std::endl;
		exit(EXIT_FAILURE);
	}
  	tempFile.close();
	
	std::cout << "Writing complete." << std::endl;
}

void FlashAmrMesh::printRefinedDataHDF5(double**** outData, const char *outFileName) {
	/* As of 3/7/2016 this function is depreciated. Use writeOutFileUniform instead. 
	   This was written to support a specific fractal analysis tool 'boxCounter3D'
       and 'multifracBoxCounter3D', only to be used for that. Note that the driver call,
       'ExportUniformData' does not actually use this function, so though it was developed,
       it may not have actually been used.
	*/

	/* Adapted from the hyperslab dataset HDF5 write example: 
	 * http://www.hdfgroup.org/ftp/HDF5/examples/parallel/Hyperslab_by_row.c
	 * 
	 * Indexing for FLASH output is (Z,Y,X), as opposed to the usual (X,Y,Z). There's some weird
	 * indexing in this function to compensate for this difference and match the FLASH output style.
	 * 
	 * Note: Automatically prints out each variable that was added in the main function. No need
	 *       to call the function for each variable like the other print function.
	 *
	 * Also when printing out data, it's likely that the image will be flipped because most image
	 * viewers start with the 0,0 index at the top left corner rather than the bottom left corner.
	 * There's a note in the data buffer for loop explaining how to change the indexing to unflip
	 * the image.
	 */

    /*
     * HDF5 APIs definitions
     */ 	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[MESH_MDIM];                 /* dataset dimensions */
    hsize_t     chunk_dims[MESH_MDIM];            /* chunk dimensions */
    hsize_t		count[MESH_MDIM];	          /* hyperslab selection parameters */
    hsize_t		offset[MESH_MDIM];
    hid_t		plist_id;                 /* property list identifier */
 
 	double 		*contiguousOutData;

	typedef struct intScalarStruct {
		char name[MAX_STRING_LENGTH+1];
		int  value;
	} intScalarStruct;

    /* 
     * Set up file access property list with parallel I/O access
     */
     plist_id = H5Pcreate(H5P_FILE_ACCESS);
     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    /*
     * Create a new file collectively and release property list identifier.
     */
    file_id = H5Fcreate(outFileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    /*                           */
 	/* First output the metadata */
 	/*                           */
	intScalarStruct outScalars[3];
	strcpy(outScalars[IAXIS].name,"nx");
	strcpy(outScalars[JAXIS].name,"ny");
	strcpy(outScalars[KAXIS].name,"nz");

	outScalars[IAXIS].value = totalCells[IAXIS];
	outScalars[JAXIS].value = totalCells[JAXIS];
	outScalars[KAXIS].value = totalCells[KAXIS];

	/* This is done in parallel, since it won't work in serial. Output isn't affected. */
	hid_t strtype, memtype, filetype, space, dset;
	hsize_t dimsScalar[1] = {3};


    strtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (strtype, MAX_STRING_LENGTH);

    memtype = H5Tcreate (H5T_COMPOUND, sizeof (intScalarStruct));
    H5Tinsert (memtype, "name", HOFFSET (intScalarStruct, name), strtype);
    H5Tinsert (memtype, "value", HOFFSET (intScalarStruct, value), H5T_STD_I32LE);
    
    filetype = H5Tcreate (H5T_COMPOUND, MAX_STRING_LENGTH + 4); // H5T_STD_I32LE is 4 bytes.
    H5Tinsert (filetype, "name", 0, strtype);
    H5Tinsert (filetype, "value", MAX_STRING_LENGTH, H5T_STD_I32LE); //3rd arg is the offset

	space = H5Screate_simple (1, dimsScalar, NULL);

    dset = H5Dcreate1 (file_id, "integer metadata", filetype, space, H5P_DEFAULT);
    H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, outScalars);

    H5Dclose (dset);
    H5Sclose (space);
    H5Tclose (filetype);

    /*                                   */
    /* Output the grid bounding box data */
    /*                                   */
    hsize_t dimsBounding[3];
    dimsBounding[0] = 1;
    dimsBounding[1] = getDim();
    dimsBounding[2] = 2;

    double *boundingBoxOut = new double[getDim()*2];
	
	boundingBoxOut[0] = refinedDomainBoundBox[LOWER][IAXIS];
	boundingBoxOut[1] = refinedDomainBoundBox[UPPER][IAXIS];
    if (getDim() > 1) {
    	boundingBoxOut[2] = refinedDomainBoundBox[LOWER][JAXIS];
    	boundingBoxOut[3] = refinedDomainBoundBox[UPPER][JAXIS];
    }
    if (getDim() > 2) {
    	boundingBoxOut[4] = refinedDomainBoundBox[LOWER][KAXIS];
    	boundingBoxOut[5] = refinedDomainBoundBox[UPPER][KAXIS];
	}
    
    space = H5Screate_simple (3, dimsBounding, NULL);

    dset = H5Dcreate1 (file_id, "bounding box", H5T_IEEE_F64LE, space, H5P_DEFAULT);

    H5Dwrite (dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, boundingBoxOut);

    H5Dclose (dset);
    H5Sclose (space);

    delete[] boundingBoxOut;


   	/*                                       */
	/* Output the uniform mesh variable data */
    /*                                       */

    /*
     * Create the dataspace for the dataset.
     */

    // int procsI, procsJ, procsK; //For loop iterators -- If one dimension is not present, set the proc to 1.
    // procsI = 1;
    // procsJ = 1;
    // procsK = 1;
    // if (getDim() >= 1) {
    	// procsI = numProcs[IAXIS];
    // }
    // if (getDim() >= 2) {
    	// procsJ = numProcs[JAXIS];
    // }
    // if (getDim() >= 3) {
    	// procsK = numProcs[KAXIS];
    // }
    
    /* Z-Y-X ordering */
    dimsf[0] = totalCells[KAXIS];
	dimsf[1] = totalCells[JAXIS];
	dimsf[2] = totalCells[IAXIS];

    chunk_dims[0] = validCells[KAXIS];
    chunk_dims[1] = validCells[JAXIS];
    chunk_dims[2] = validCells[IAXIS];

    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = chunk_dims[0];
    count[1] = chunk_dims[1];
    count[2] = chunk_dims[2];

    offset[0] = procIndexMin[mpiGetID()][KAXIS];
    offset[1] = procIndexMin[mpiGetID()][JAXIS];
    offset[2] = procIndexMin[mpiGetID()][IAXIS];

    // std::cout << "Procs: " << std::endl;
    // std::cout << "\t" << numProcs[KAXIS] << std::endl;
    // std::cout << "\t" << numProcs[JAXIS] << std::endl;
    // std::cout << "\t" << numProcs[IAXIS] << std::endl;

    // std::cout << "chunk_dims: " << std::endl;
    // std::cout << "\t" << chunk_dims[0] << std::endl;
    // std::cout << "\t" << chunk_dims[1] << std::endl;
    // std::cout << "\t" << chunk_dims[2] << std::endl;

    // std::cout << "dims_f: " << std::endl;
    // std::cout << "\t" << dimsf[0] << std::endl;
    // std::cout << "\t" << dimsf[1] << std::endl;
    // std::cout << "\t" << dimsf[2] << std::endl;

    // std::cout << "Proc: " << mpiGetID() << std::endl;
    // std::cout << "Offset Z: " << offset[0] << std::endl;
    // std::cout << "Offset Y: " << offset[1] << std::endl;
    // std::cout << "Offset X: " << offset[2] << std::endl;

	contiguousOutData = new double[chunk_dims[IAXIS]*chunk_dims[JAXIS]*chunk_dims[KAXIS]];
    
    // for (int i = 0; i < (int)chunk_dims[IAXIS]*(int)chunk_dims[JAXIS]*(int)chunk_dims[KAXIS]; i++) {
    // 	contiguousOutData[i] = -1e99;
    // }

	for (int var = 0; var < getNUserVars(); var++) {

	    filespace = H5Screate_simple(MESH_MDIM, dimsf, NULL); 

	    /*
	     * Create the dataset with default properties and close filespace.
	     */	    
	    if (!SUPPRESSOUTPUT) std::cout << "Writing variable: " << findVarName(var).c_str() << std::endl;
	    dset_id = H5Dcreate(file_id, findVarName(var).c_str(), H5T_IEEE_F64LE, filespace,
				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    H5Sclose(filespace);
    	
    	memspace  = H5Screate_simple(MESH_MDIM, count, NULL); 

		/*
		 * Select hyperslab in the file.
		 */
		filespace = H5Dget_space(dset_id);
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

		// std::cout << "Status: " << status << std::endl;

	    /*
	     * Initialize data buffer 
	     */

		//NOTE: FLASH has the 0,0 index defined at the bottom left corner. However,
		//      the 0,0 index in most image viewers is at the top left corner.
		//      To switch, use outData[var][k][(int)chunk_dims[1]-1-j][i] instead.
	    int ind;
	    for (int i = 0; i < (int)chunk_dims[0]; i++) {
	    	for (int j = 0; j < (int)chunk_dims[1]; j++) {
	    		for (int k = 0; k < (int)chunk_dims[2]; k++) {
	    			ind =   i * (int)chunk_dims[2]*(int)chunk_dims[1]
	    				  + j * (int)chunk_dims[2]
	    				  + k;
	    			contiguousOutData[ind] = outData[var][k][j][i];
	    			// std::cout << contiguousOutData[ind] << " ";
	    			// std::cout << outData[var][k][j][i] << std::endl;
	    			// if (outData[var][k][j][i] == 0) {
	    			// 	std::cout << k << " " << j << " " << i << " " << ind << std::endl;
	    			// }
	    		}
	    	}
	    }
	    
	    /*
	     * Create property list for collective dataset write.
	     */

	    plist_id = H5Pcreate(H5P_DATASET_XFER);
	    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	    H5Dwrite(dset_id, H5T_IEEE_F64LE, memspace, filespace,
			      plist_id, contiguousOutData);
	    /*
	     * Close/release resources.
	     */
	    H5Dclose(dset_id);
	    H5Sclose(filespace);
	    H5Sclose(memspace);
	    H5Pclose(plist_id);
	}

    delete[] contiguousOutData;

   	H5Fclose(file_id);
}

void FlashAmrMesh::writeIntsUniform(const char *which){

	typedef struct intScalarStruct {
		char name[MAX_STRING_LENGTH];
		int  value;
	} intScalarStruct;

	char datasetName[50] = {};
	if (std::string(which) == "scalars") {
		strcpy(datasetName,"integer scalars");
	}
	else if (std::string(which) == "parameters") {
		strcpy(datasetName,"integer runtime parameters");
	}

	if (!datasetAvailable(datasetName,"required")) {
		return;
	}

	//Get scalars list.
	hid_t dset_in = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );

	//Make specific type for the 80 character string.
	hid_t strtype_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype_in, MAX_STRING_LENGTH);

	//"integer scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(intScalarStruct) );
	H5Tinsert( memtype_in, "name", HOFFSET(intScalarStruct, name), strtype_in);
	H5Tinsert( memtype_in, "value", HOFFSET(intScalarStruct, value), H5T_STD_I32LE);

	intScalarStruct* intVals_in = new intScalarStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, intVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    /* A couple of values need to be changed before being written. */
    for (int i = 0; i < (int)dims_in[0]; i++) {
    	std::string compareStr = std::string(intVals_in[i].name);
    	if (compareStr.compare(0,3,"nxb") == 0) {
    		intVals_in[i].value = totalCells[IAXIS];
    	}
    	else if (compareStr.compare(0,3,"nyb") == 0) {
    		intVals_in[i].value = totalCells[JAXIS];
    	}
    	else if (compareStr.compare(0,3,"nzb") == 0) {
    		intVals_in[i].value = totalCells[KAXIS];
    	}
    	else if (compareStr.compare(0,15,"globalnumblocks") == 0) {
    		intVals_in[i].value = 1;
    	}
    }

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(intScalarStruct) );
	H5Tinsert( memtype_out, "name", HOFFSET(intScalarStruct, name), strtype_in);
	H5Tinsert( memtype_out, "value", HOFFSET(intScalarStruct, value), H5T_STD_I32LE);

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(intScalarStruct) );
	H5Tinsert( filetype_out, "name", HOFFSET(intScalarStruct, name), strtype_in);
	H5Tinsert( filetype_out, "value", HOFFSET(intScalarStruct, value), H5T_STD_I32LE);

    hid_t dset_out = H5Dcreate(outFile, datasetName, filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, intVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, intVals_in);
	H5Tclose(memtype_in);
	H5Sclose(space_in);
	H5Dclose(dset_in);	
	H5Tclose(strtype_in);

	delete[] intVals_in;
}

void FlashAmrMesh::writeRealsUniform(const char *which){
	//This routine returns a double corresponding to the 
	// requested scalar char array argument from either the 
	// real scalar list or the real runtime parameters list. 
	// These values go into the map object realScalarsMap.

	typedef struct realScalarStruct {
		char 	name[MAX_STRING_LENGTH];
		double  value;
	} realScalarStruct;


	char datasetName[50] = {};
	if (std::string(which) == "scalars") {
		strcpy(datasetName,"real scalars");
	}
	else if (std::string(which) == "parameters") {
		strcpy(datasetName,"real runtime parameters");
	}

	if (!datasetAvailable(datasetName,"required")) {
		return;
	}

	//Get scalars list.
	hid_t dset_in = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );

	//Make specific type for the 80 character string.
	hid_t strtype_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype_in, MAX_STRING_LENGTH);

	//"real scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(realScalarStruct) );
	H5Tinsert( memtype_in, "name", HOFFSET(realScalarStruct, name), strtype_in);
	H5Tinsert( memtype_in, "value", HOFFSET(realScalarStruct, value), H5T_IEEE_F64LE);

	realScalarStruct* realVals_in = new realScalarStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, realVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    /* A couple of values need to be changed before being written. */
    for (int i = 0; i < (int)dims_in[0]; i++) {
    	std::string compareStr = std::string(realVals_in[i].name);	
    	if (compareStr.compare(0,4,"xmin") == 0) {   
    		realVals_in[i].value = refinedDomainBoundBox[LOWER][IAXIS];
    	}
    	else if (compareStr.compare(0,4,"xmax") == 0) {
    		realVals_in[i].value = refinedDomainBoundBox[UPPER][IAXIS];
    	}
    	else if (compareStr.compare(0,4,"ymin") == 0) {
    		realVals_in[i].value = refinedDomainBoundBox[LOWER][JAXIS];
    	}
    	else if (compareStr.compare(0,4,"ymax") == 0) {
    		realVals_in[i].value = refinedDomainBoundBox[UPPER][JAXIS];
    	}
    	else if (compareStr.compare(0,4,"zmin") == 0) {
    		realVals_in[i].value = refinedDomainBoundBox[LOWER][KAXIS];
    	}
    	else if (compareStr.compare(0,4,"zmax") == 0) {
    		realVals_in[i].value = refinedDomainBoundBox[UPPER][KAXIS];
    	}
    }

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(realScalarStruct) );
	H5Tinsert( memtype_out, "name", HOFFSET(realScalarStruct, name), strtype_in);
	H5Tinsert( memtype_out, "value", HOFFSET(realScalarStruct, value), H5T_IEEE_F64LE);

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(realScalarStruct) );
	H5Tinsert( filetype_out, "name", HOFFSET(realScalarStruct, name), strtype_in);
	H5Tinsert( filetype_out, "value", HOFFSET(realScalarStruct, value), H5T_IEEE_F64LE);

    hid_t dset_out = H5Dcreate(outFile, datasetName, filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, realVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, realVals_in);
	H5Dclose(dset_in);
	H5Tclose(memtype_in);
	H5Sclose(space_in);
	H5Tclose(strtype_in);

	delete[] realVals_in;
}

void FlashAmrMesh::writeSimParametersUniform() {
	// This function is for outputting file type 7 information

	typedef struct simParamStruct {
		//Instead of this file being formatted like the rest where there is a name
		// column and a value column, file format 7 has 8 columns with no names.
		// The values have to be read directly into the variables here.
		int totalBlocks;
		float time;
		float timeStep;
		float redShift;
		int numberOfSteps;
		int nxb;
		int nyb;
		int nzb;
	} simParamStruct;

	if (!datasetAvailable("simulation parameters","required")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "simulation parameters", H5P_DEFAULT);
	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );

	//"simulation parameters" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(simParamStruct));

	//NOTE: IF THERE ARE ANY EXTRA PARAMETERS HERE, THEY WILL NOT BE READ.
	//      THIS IS DUE TO THE WEIRD WAY THAT THE PARAMETERS WERE WRITTEN TO THE HDF5 FILE.
	H5Tinsert( memtype_in, "total blocks", 	HOFFSET(simParamStruct, totalBlocks), 	H5T_STD_I32LE);
	H5Tinsert( memtype_in, "time", 			HOFFSET(simParamStruct, time), 			H5T_IEEE_F32LE);
	H5Tinsert( memtype_in, "timestep", 		HOFFSET(simParamStruct, timeStep), 		H5T_IEEE_F32LE);
	H5Tinsert( memtype_in, "redshift", 		HOFFSET(simParamStruct, redShift), 		H5T_IEEE_F32LE);
	H5Tinsert( memtype_in, "number of steps", 	HOFFSET(simParamStruct, numberOfSteps), H5T_STD_I32LE);
	H5Tinsert( memtype_in, "nxb", 				HOFFSET(simParamStruct, nxb), 			H5T_STD_I32LE);
	H5Tinsert( memtype_in, "nyb", 				HOFFSET(simParamStruct, nyb), 			H5T_STD_I32LE);
	H5Tinsert( memtype_in, "nzb", 				HOFFSET(simParamStruct, nzb), 			H5T_STD_I32LE);
	
	simParamStruct* simParamVals_in = new simParamStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, simParamVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    /* A couple of values need to be changed before being written. */
    simParamVals_in[0].totalBlocks = 1;
    simParamVals_in[0].nxb = totalCells[IAXIS];
    simParamVals_in[0].nyb = totalCells[JAXIS];
    simParamVals_in[0].nzb = totalCells[KAXIS];

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(simParamStruct) );
	H5Tinsert( memtype_out, "total blocks", 	HOFFSET(simParamStruct, totalBlocks), 	H5T_STD_I32LE);
	H5Tinsert( memtype_out, "time", 			HOFFSET(simParamStruct, time), 			H5T_IEEE_F32LE);
	H5Tinsert( memtype_out, "timestep", 		HOFFSET(simParamStruct, timeStep), 		H5T_IEEE_F32LE);
	H5Tinsert( memtype_out, "redshift", 		HOFFSET(simParamStruct, redShift), 		H5T_IEEE_F32LE);
	H5Tinsert( memtype_out, "number of steps", 	HOFFSET(simParamStruct, numberOfSteps), H5T_STD_I32LE);
	H5Tinsert( memtype_out, "nxb", 				HOFFSET(simParamStruct, nxb), 			H5T_STD_I32LE);
	H5Tinsert( memtype_out, "nyb", 				HOFFSET(simParamStruct, nyb), 			H5T_STD_I32LE);
	H5Tinsert( memtype_out, "nzb", 				HOFFSET(simParamStruct, nzb), 			H5T_STD_I32LE);

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(simParamStruct) );
	H5Tinsert( filetype_out, "total blocks", 	HOFFSET(simParamStruct, totalBlocks), 	H5T_STD_I32LE);
	H5Tinsert( filetype_out, "time", 			HOFFSET(simParamStruct, time), 			H5T_IEEE_F32LE);
	H5Tinsert( filetype_out, "timestep", 		HOFFSET(simParamStruct, timeStep), 		H5T_IEEE_F32LE);
	H5Tinsert( filetype_out, "redshift", 		HOFFSET(simParamStruct, redShift), 		H5T_IEEE_F32LE);
	H5Tinsert( filetype_out, "number of steps", 	HOFFSET(simParamStruct, numberOfSteps), H5T_STD_I32LE);
	H5Tinsert( filetype_out, "nxb", 				HOFFSET(simParamStruct, nxb), 			H5T_STD_I32LE);
	H5Tinsert( filetype_out, "nyb", 				HOFFSET(simParamStruct, nyb), 			H5T_STD_I32LE);
	H5Tinsert( filetype_out, "nzb", 				HOFFSET(simParamStruct, nzb), 			H5T_STD_I32LE);

    hid_t dset_out = H5Dcreate(outFile, "simulation parameters", filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, simParamVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, simParamVals_in);
	H5Dclose(dset_in);
	H5Sclose(space_in);
	H5Tclose(memtype_in);

	delete[] simParamVals_in;
}

void FlashAmrMesh::writeCoordinatesUniform() {

    /*                                   */
    /* Output the grid coordinate data   */
    /*                                   */
    hsize_t dimsCoordinates[2];
    dimsCoordinates[0] = 1;
    dimsCoordinates[1] = getDim();

    double *coordinatesOut = new double[getDim()*2];
	
	coordinatesOut[0] = (refinedDomainBoundBox[UPPER][IAXIS] + refinedDomainBoundBox[LOWER][IAXIS]) / 2.0;
    if (getDim() > 1) {
    	coordinatesOut[1] = (refinedDomainBoundBox[UPPER][JAXIS] + refinedDomainBoundBox[LOWER][JAXIS]) / 2.0;
    }
    if (getDim() > 2) {
    	coordinatesOut[2] = (refinedDomainBoundBox[UPPER][KAXIS] + refinedDomainBoundBox[LOWER][KAXIS]) / 2.0;
	}
    
    hid_t filespace = H5Screate_simple (2, dimsCoordinates, NULL);
    hid_t dset_id = H5Dcreate (outFile, "coordinates", H5T_IEEE_F64LE, filespace, 
    	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
	    H5Dwrite (dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, plist_id, coordinatesOut);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    H5Dclose (dset_id);

    delete[] coordinatesOut;
}

void FlashAmrMesh::writeBoundsUniform() {

    /*                                   */
    /* Output the grid bounding box data */
    /*                                   */
    hsize_t dimsBounding[3];
    dimsBounding[0] = 1;
    dimsBounding[1] = getDim();
    dimsBounding[2] = 2;

    double *boundingBoxOut = new double[getDim()*2];
	
	boundingBoxOut[0] = refinedDomainBoundBox[LOWER][IAXIS];
	boundingBoxOut[1] = refinedDomainBoundBox[UPPER][IAXIS];
    if (getDim() > 1) {
    	boundingBoxOut[2] = refinedDomainBoundBox[LOWER][JAXIS];
    	boundingBoxOut[3] = refinedDomainBoundBox[UPPER][JAXIS];
    }
    if (getDim() > 2) {
    	boundingBoxOut[4] = refinedDomainBoundBox[LOWER][KAXIS];
    	boundingBoxOut[5] = refinedDomainBoundBox[UPPER][KAXIS];
	}
    
    hid_t filespace = H5Screate_simple (3, dimsBounding, NULL);
    hid_t dset_id = H5Dcreate (outFile, "bounding box", H5T_IEEE_F64LE, filespace, 
    	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite (dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, plist_id, boundingBoxOut);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    H5Dclose (dset_id);

    delete[] boundingBoxOut;
}

void FlashAmrMesh::writeNodeTypeUniform() {

    /*                                   */
    /* Output the node type data         */
    /*                                   */
    hsize_t dimsNodeType[1];
    dimsNodeType[0] = 1;

    // One block, has to be a leaf block.
    int nodeTypeOut = 1;
	
    hid_t filespace = H5Screate_simple (1, dimsNodeType, NULL);
    hid_t dset_id = H5Dcreate (outFile, "node type", H5T_STD_I32LE, filespace, 
    	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite (dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, plist_id, &nodeTypeOut);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    H5Dclose (dset_id);
}

void FlashAmrMesh::writeRefineLevelUniform() {

    /*                                   */
    /* Output the refinement level data  */
    /*                                   */
    hsize_t dimsRefineLevel[1];
    dimsRefineLevel[0] = 1;

    // One block, has to be a leaf block.
    int refineLevelOut = 1;
	
    hid_t filespace = H5Screate_simple (1, dimsRefineLevel, NULL);
    hid_t dset_id = H5Dcreate (outFile, "refine level", H5T_STD_I32LE, filespace, 
    	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite (dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, plist_id, &refineLevelOut);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    H5Dclose (dset_id);
}

void FlashAmrMesh::writeGidUniform() {


    /*                                   */
    /* Output the gid array              */
    /*                                   */
    hsize_t dimsGid[2];
    dimsGid[0] = 1;
    dimsGid[1] = 2*getDim() + 1 + (int)pow(2,getDim());
    			 // number of neighbors + parent + number of children

    int numElements = 2*getDim() + 1 + (int)pow(2,getDim());
    int *gidOut = new int[numElements];
	
    /* Went ahead and described each element in the array here in case
       one in particular needs to be changed later on.
     */
	if (getDim() == 1) {
		/* 1d neighbors (x-direction) */
		gidOut[0] = -1;
		gidOut[1] = -1;

		/* 1d parent */
		gidOut[2] = -1;

		/* 1d children */
		gidOut[3] = -1;
		gidOut[4] = -1;
	}
    else if (getDim() == 2) {
    	/* 2d neighbors (x- and y-direction) */
    	gidOut[0] = -1;
		gidOut[1] = -1;
		gidOut[2] = -1;
		gidOut[3] = -1;

		/* 2d parent */
		gidOut[4] = -1;

		/* 2d children */
		gidOut[5] = -1;
		gidOut[6] = -1;
		gidOut[7] = -1;
		gidOut[8] = -1;
    }
    else if (getDim() == 3) {
    	/* 3d neighbors (x-, y-, and z-direction) */
    	gidOut[0] = -1;
		gidOut[1] = -1;
		gidOut[2] = -1;
		gidOut[3] = -1;
		gidOut[4] = -1;
		gidOut[5] = -1;

		/* 3d parent */
		gidOut[6] = -1;

		/* 3d children */
		gidOut[7] = -1;
		gidOut[8] = -1;
		gidOut[9] = -1;
		gidOut[10] = -1;
		gidOut[11] = -1;
		gidOut[12] = -1;
		gidOut[13] = -1;
		gidOut[14] = -1;		
	}
    
    hid_t filespace = H5Screate_simple (2, dimsGid, NULL);
    hid_t dset_id = H5Dcreate (outFile, "gid", H5T_STD_I32LE, filespace, 
    	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite (dset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, plist_id, gidOut);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    H5Dclose (dset_id);

    delete[] gidOut;
}

void FlashAmrMesh::writeVarDataUniform(const char *outputVariable, double****& outData) {
    
    /* Z-Y-X ordering */
    hsize_t dims_full[4];
    dims_full[0] = 1; // Acting as if there's only one block of data for uniform output
    dims_full[1] = totalCells[KAXIS];
	dims_full[2] = totalCells[JAXIS];
	dims_full[3] = totalCells[IAXIS];

	hsize_t chunk_dims[4];
    chunk_dims[0] = 1;
    chunk_dims[1] = validCells[KAXIS];
    chunk_dims[2] = validCells[JAXIS];
    chunk_dims[3] = validCells[IAXIS];

    hsize_t count_out[4];
	count_out[0] = chunk_dims[0];
	count_out[1] = chunk_dims[1];
	count_out[2] = chunk_dims[2];
	count_out[3] = chunk_dims[3];

    hsize_t offset_out[4];
    offset_out[0] = 0;
    offset_out[1] = procIndexMin[mpiGetID()][KAXIS];
    offset_out[2] = procIndexMin[mpiGetID()][JAXIS];
    offset_out[3] = procIndexMin[mpiGetID()][IAXIS];

    hid_t filespace = H5Screate_simple(4, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, outputVariable, H5T_IEEE_F64LE, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	double* contiguousOutData = new double[chunk_dims[0]*chunk_dims[1]*chunk_dims[2]*chunk_dims[3]];

	double maxVarVal = -1e308;
	double minVarVal = 1e308;
    int ind;
    int var = uVarsMap_sToI[std::string(outputVariable)];

	hid_t memspace_out  = H5Screate_simple(4, chunk_dims, NULL);
	hid_t filespace_out = H5Dget_space(dset_id);

    int counter = 0;
    for (int i = 0; i < (int)chunk_dims[3]; i++) {
    	for (int j = 0; j < (int)chunk_dims[2]; j++) {
    		for (int k = 0; k < (int)chunk_dims[1]; k++) {
    			// Convert from X-Y-Z to Z-Y-X ordering in final array 
    			ind =   k * (int)chunk_dims[3]*(int)chunk_dims[2] 
    			      + j * (int)chunk_dims[3]
    			      + i;
    			contiguousOutData[ind] = outData[var][i][j][k];
				maxVarVal = std::max(maxVarVal,contiguousOutData[ind]);
				minVarVal = std::min(minVarVal,contiguousOutData[ind]);    			
    		}
		}
	}

	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, NULL, 
		                count_out, NULL);

	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	H5Dwrite(dset_id, H5T_IEEE_F64LE, memspace_out, filespace_out,
      		 plist_id, contiguousOutData);

	H5Pclose(plist_id);	


	/* Attributes must be written collectively, so this is done on all
	   processors.
	*/
	double maxVarValAll;
	double minVarValAll;

	MPI_Allreduce(&maxVarVal, &maxVarValAll, 1, MPI_DOUBLE, MPI_MAX, comm);
	MPI_Allreduce(&minVarVal, &minVarValAll, 1, MPI_DOUBLE, MPI_MIN, comm);


	/* Write the min/max attributes */
	hsize_t dims_attr = 1;

	hid_t attrspace_id = H5Screate_simple(1, &dims_attr, NULL);
	hid_t attribute_id = H5Acreate(dset_id, "maximum", H5T_IEEE_F64LE, attrspace_id,
		                      H5P_DEFAULT, H5P_DEFAULT);
	
	

	H5Awrite(attribute_id, H5T_IEEE_F64LE, &maxVarValAll);

	H5Aclose(attribute_id);
	H5Sclose(attrspace_id);

	attrspace_id = H5Screate_simple(1, &dims_attr, NULL);
	attribute_id = H5Acreate(dset_id, "minimum", H5T_IEEE_F64LE, attrspace_id,
		                      H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_IEEE_F64LE, &minVarValAll);

	H5Aclose(attribute_id);
	H5Sclose(attrspace_id);


	H5Sclose(memspace_out);
    H5Sclose(filespace_out);
    H5Dclose(dset_id);

    delete[] contiguousOutData;
}

void FlashAmrMesh::writeVariablesUniform(double****& outData) {
	// Only write the variables the user wants to have written.
	// Default behavior is to write all variables, but the user
	// can use the noWrite() function in main to switch from writing
	// a variable to not writing a variable.

	for (int i = 0; i < nUserVars; i++) {
		if (writeVarMap[variableNames[i]]) {
			if (!SUPPRESSOUTPUT) std::cout << "Writing variable " << variableNames[i] << "...";
			writeVarDataUniform(variableNames[i],outData);
			if (!SUPPRESSOUTPUT) std::cout << " Variable " << variableNames[i] << " written." << std::endl;
		}
	}
}

void FlashAmrMesh::generalWriteUniform() {
	// Writes variables similar to how they are read in the generalRead function.

	if (!SUPPRESSOUTPUT) std::cout << "Writing scalar lists. " << std::endl;
	if (fileType == 7) {
		if (!SUPPRESSOUTPUT) std::cout << "File type 7: Writing file type." << std::endl;
		this->writeFileType();
		if (!SUPPRESSOUTPUT) std::cout << "File type 7: Writing Simulation Parameters." << std::endl;
		this->writeSimParametersUniform();
	}
	else {
		this->writeIntsUniform("scalars");
		this->writeStrings("scalars");
		this->writeRealsUniform("scalars");
		this->writeLogicals("scalars");
	}

	if (!SUPPRESSOUTPUT) std::cout << "Reading parameter lists. " << std::endl;
	this->writeIntsUniform("parameters");
	this->writeStrings("parameters");
	this->writeRealsUniform("parameters");
	this->writeLogicals("parameters");

	this->writeNFileVarsList();
	if (!SUPPRESSOUTPUT) std::cout << "List of variable names written. " << std::endl;
	if (!SUPPRESSOUTPUT) std::cout << "Number of variables in file: " << this->getNUserVars() << std::endl;

	this->writeSimInfo();
	if (!SUPPRESSOUTPUT) std::cout << "Simulation info written. " << std::endl;

	/* These dataset are required in order for VisIt to be able to
	   read the HDF5 file. If they cannot be written, the code WILL crash
	   (see noWrite function, which is called at the top of each
	   of the below functions).
	*/

	this->writeCoordinatesUniform();
	if (!SUPPRESSOUTPUT) std::cout << "Coordinates written. " << std::endl;

	this->writeBoundsUniform();
	if (!SUPPRESSOUTPUT) std::cout << "Block bounds written. " << std::endl;

	this->writeNodeTypeUniform();
	if (!SUPPRESSOUTPUT) std::cout << "Node types written. " << std::endl;

	this->writeRefineLevelUniform();
	if (!SUPPRESSOUTPUT) std::cout << "Refine levels written. " << std::endl;

	this->writeGidUniform();
	if (!SUPPRESSOUTPUT) std::cout << "Gid array written. " << std::endl;

	this->writeWhichChildUniform();
	if (!SUPPRESSOUTPUT) std::cout << "Which child array written. " << std::endl;

	this->writeBFlagsUniform();
	if (!SUPPRESSOUTPUT) std::cout << "bFlags array written. " << std::endl;

	this->writeSizesUniform();
	if (!SUPPRESSOUTPUT) std::cout << "block sizes array written. " << std::endl;

	
	
}

void FlashAmrMesh::writeOutFileUniform(double****& outData) {
	if (mpiGetID() == 0) {
		std::cout << "Writing file" << std::endl;
	}
	createOutFile();
	generalWriteUniform();
	writeVariablesUniform(outData);
	closeOutFile();
	if (mpiGetID() == 0) {
		std::cout << "Writing complete" << std::endl;
	}
}

void FlashAmrMesh::deallocateFinest(double****& outData) {
	//Deallocate the 4-d data array.
	for (int i = 0; i < mpiGetProcs(); i++) {
		delete[] procIndexMin[i];
	}
	delete[] procIndexMin;

	for (int i = 0; i < getNUserVars(); i++) {
		for (int j = 0; j < fineCells[IAXIS]; j++) {
			for (int k = 0; k < fineCells[JAXIS]; k++) {
				delete[] outData [i][j][k];
			}
			delete[]outData [i][j];
		}
		delete[] outData[i];
	}
	delete[] outData;
}

void FlashAmrMesh::writeWhichChildUniform() {
                              
    hsize_t dimsWhichChild[1];
    dimsWhichChild[0] = 1;
    int whichChildOut[1] = {-1};
    
    hid_t filespace = H5Screate_simple (1, dimsWhichChild, NULL);
    hid_t dset_id = H5Dcreate (outFile, "which child", H5T_NATIVE_INT, filespace, 
    	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite (dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, whichChildOut);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    H5Dclose (dset_id);
}

void FlashAmrMesh::writeBFlagsUniform() {
                              
    hsize_t dimsBFlags[1];
    dimsBFlags[0] = 1;
    int bFlagsOut[1] = {-1};
    
    hid_t filespace = H5Screate_simple (1, dimsBFlags, NULL);
    hid_t dset_id = H5Dcreate (outFile, "bflags", H5T_NATIVE_INT, filespace, 
    	    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite (dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, bFlagsOut);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    H5Dclose (dset_id);
}

void FlashAmrMesh::writeSizesUniform() {

	if (!datasetAvailable("block size","optional")) {
		return;
	}

	// there's only one block in this case
	// its size is the total size of the subdomain
	double blockSize[MESH_MDIM];
	for (int i=0; i<MESH_MDIM; i++)
		blockSize[i] = totalCells[i] * gridDelta[i];
	

	/* Output */
	hsize_t dims_full[2];
	dims_full[0] = 1;
	dims_full[1] = 3;

	// Create a data space that expects the total amount of data across all processors.
	hid_t dset_id;
	hid_t filespace;
	if (fileType == 7) {
		filespace = H5Screate_simple(2, dims_full, NULL);
    	dset_id = H5Dcreate(outFile, "block size", H5T_NATIVE_FLOAT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else {
		filespace = H5Screate_simple(2, dims_full, NULL);
   		dset_id = H5Dcreate(outFile, "block size", H5T_NATIVE_DOUBLE, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	H5Sclose(filespace);
    
	// Prepare the hyperslab containing this processor's data
	// hsize_t offset_out[2]; //hyperslab offset in new file
	// hsize_t count_out[2];  //size of hyperslab
	// offset_out[0] = mpiGetGlobalBlockID(0);
	// offset_out[1] = 0;
	// count_out[0] = nBlocksAssigned;
	// count_out[1] = dims_in[1];

	// hid_t memspace_out = H5Screate_simple(2, count_out, NULL);
    
	// hid_t filespace_out = H5Dget_space(dset_id);
	// H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
	// 	                                NULL, count_out, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write hyperslab to dataset
    if (mpiGetID() == 0) 
    {
	    if (fileType == 7) {
	    	H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		  	        	      plist_id, blockSize);    	
	    }
	    else {
	    	H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		  		        	  plist_id, blockSize);
	    }
	}

    H5Pclose(plist_id);
	// H5Sclose(filespace_out);
	// H5Sclose(memspace_out);
    H5Dclose(dset_id);
    // H5Sclose(memspace_in);	
	// H5Dclose(dset_in);	

	// delete[] tempSizesOutD;
	// delete[] tempSizesOutF;

}

int* FlashAmrMesh::getTotalCellsUniform()
{
	return totalCells;
}
