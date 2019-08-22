
#include <cstdlib>
#include <iostream>
#include "FlashAmrMesh.h"
#include "Driver_main.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

int LAVAFLOW_DRIVER_BASENM(int argc, char **argv) {

	/* Parameters */
//	std::string inFileName("/home/pb12c/Research/LAVAflow/drivers/files/testData/cartesianShells/cartesianShells_hdf5_chk_0000");
	std::string inFileName("");

	/* A uniform mesh will be output at this refinement level. 
     * If the refinement level is too high, the highest refinement
     * level will automatically be used.
     */
	int  refineLevel = 5;

	// std::vector<std::string> meshVars = {"vari","rand"};
	std::vector<std::string> meshVars = {"dens"};
	std::vector<std::string> auxVars  = {"new"};

	// char outAsciiName[MAX_FILENAME_LENGTH] = "testAscii";
	char outAsciiName[MAX_FILENAME_LENGTH] = "sedovAscii";
	// char outHDF5Name[MAX_FILENAME_LENGTH] = "testHDF5_multi";
	char outHDF5Name[MAX_FILENAME_LENGTH] = "sedovHDF5_lref5";
	/* End parameters */

	// char inProgramFileName[MAX_FILENAME_LENGTH] = "./";
	// char outProgramAsciiName[MAX_FILENAME_LENGTH] = "./";
	char outProgramHDF5Name[MAX_FILENAME_LENGTH] = "./";

	// strcat(inProgramFileName,inFileName);
	// strcat(outProgramAsciiName,outAsciiName);
	strcat(outProgramHDF5Name,outHDF5Name);

	FlashAmrMesh mesh(inFileName, outProgramHDF5Name, meshVars, auxVars);

	/* Test to verify that new variables can also be refined */
	int newVarIndex = mesh.findVarIndex("new"); //new variable

	int numLocalBlocks = mesh.getLocalNumBlocks();
	int numLeafBlocks = 0;
	int *blockList = new int[numLocalBlocks];
    double* coords;
    coords = new double [8]; 
	mesh.getListOfBlocks(LEAF, blockList, numLeafBlocks);

	double cumVari = 0;
	double cumRand = 0;
	int cellCount = 0;

	int startingPos[MESH_MDIM];
	startingPos[0] = 0;
	startingPos[1] = 0;
	startingPos[2] = 0;
	for (int lb=0; lb<numLeafBlocks; lb++)
	{
		int currentBlock = blockList[lb];
		double ****solnData;
		//Get pointer to the block's data
		mesh.getBlkPtr(currentBlock, solnData, CENTER);

		int blkLimits[2][MESH_MDIM];
		int blkLimitsGC[2][MESH_MDIM];

		mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

		for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++){
			for (int j=blkLimits[LOWER][IAXIS]; j<= blkLimits[UPPER][IAXIS]; j++){
				for (int k=blkLimits[LOWER][IAXIS]; k<= blkLimits[UPPER][IAXIS]; k++){
					// Note that only the leaf blocks are being updated here.
					solnData[newVarIndex][i][j][k] = 1.1;
	        	}
	    	}
		}

		mesh.putBlkData(currentBlock, CENTER, "new ", INTERIOR, startingPos, solnData);
	}

	/* Only need to declare the 4d pointer here. 
 	 * Memory allocation is done within the refineToFinest function. 
 	 */
	double**** fineArray;
	
	/* Subdomain specification -- Must be within the original mesh
	 * and must be large enough to contain at least one cell
	 * per processor. 
	 *
	 * If no subdomain is required, set all values to 0.
	 *
	 * Need to specify all dimensions -- Ex: If you're using a 2D mesh,
	 * Leave the KAXIS values 0 and set the IAXIS and JAXIS values.
	 *
	 * Due to discretization, the output coordinates will be slightly offset
	 * from the original input coordinates.
	 */
 
	double subdomainCoords[2][MESH_MDIM];

	subdomainCoords[LOWER][IAXIS] = 0;
	subdomainCoords[UPPER][IAXIS] = 0;
	subdomainCoords[LOWER][JAXIS] = 0;
	subdomainCoords[UPPER][JAXIS] = 0;
	subdomainCoords[LOWER][KAXIS] = 0;
	subdomainCoords[UPPER][KAXIS] = 0;

	// subdomainCoords[LOWER][IAXIS] = -0.5;
	// subdomainCoords[UPPER][IAXIS] = 1;
	// subdomainCoords[LOWER][JAXIS] = -0.5;
	// subdomainCoords[UPPER][JAXIS] = 0;
	// subdomainCoords[LOWER][KAXIS] = -0.5;
	// subdomainCoords[UPPER][KAXIS] = -0.25;


	mesh.refineToFinest(fineArray, subdomainCoords, refineLevel);

	/* ASCII Output (For 3D mesh, the z-coordinate will be set to 0) 
         * Does not work on 1D mesh.
         * 
         * Will have one output per processor, with the name "testAscii_#.txt"
         * where # is the processor id.
         *
         */

	// char procID[3];
	// sprintf(procID,"%d",mesh.mpiGetID());
	// strcat(outProgramAsciiName, "_");
	// strcat(outProgramAsciiName, procID);
	// strcat(outProgramAsciiName, ".txt");
	// mesh.printRefinedData(fineArray, outProgramAsciiName, "ascii", "vari");	
	
	/* HDF5 output -- Works with all 3 dimensions */
	// mesh.printRefinedDataHDF5(fineArray, "testHDF5_multi_orig");
	
	mesh.writeOutFileUniform(fineArray);

	mesh.deallocateFinest(fineArray);

	std::cout << "Processor " << mesh.mpiGetID() << " at end of main." << std::endl;

	return 0;
}



