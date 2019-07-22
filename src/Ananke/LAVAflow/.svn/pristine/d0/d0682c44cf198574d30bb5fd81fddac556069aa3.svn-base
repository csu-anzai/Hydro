#include <cstdlib>
#include <iostream>
#include "FlashAmrMesh.h"
#include "Driver_main.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

int LAVAFLOW_DRIVER_BASENM(int argc, char **argv) {

  // input file name
	std::string inFileName("super2Dturbulence5000");

  // spatial dimension of data
  int numDim = 1;

  // number of levels of transform
	int numLevels = 5;

  // primitive variables to perform transform on
	std::vector<std::string> meshVars = {"dens","pres"};

  // specify output file names
	char outAsciiName[MAX_FILENAME_LENGTH] = "sedovAscii";
	char outHDF5Name[MAX_FILENAME_LENGTH] = "sedovHDF5_transfrm";
	char outProgramHDF5Name[MAX_FILENAME_LENGTH] = "./";
	strcat(outProgramHDF5Name,outHDF5Name);

  // create AMR mesh object
	FlashAmrMesh mesh(inFileName, outProgramHDF5Name, meshVars, auxVars);

  // get number of local blocks
	int numLocalBlocks = mesh.getLocalNumBlocks();

  // create list of leaf blocks
	int numLeafBlocks = 0;
	int *blockList = new int[numLocalBlocks];
  double* coords;
  coords = new double [8]; 
	mesh.getListOfBlocks(LEAF, blockList, numLeafBlocks);

	double cumVari = 0;
	double cumRand = 0;
	int cellCount = 0;

  // loop through cells in block
	int startingPos[MESH_MDIM];
	startingPos[0] = 0;
	startingPos[1] = 0;
	startingPos[2] = 0;
	for (int lb=0; lb<numLeafBlocks; lb++)
	{
    // get block id
		int currentBlock = blockList[lb];

    // declare quadruple pointer...
		double ****solnData;

		// get pointer to the block's data
		mesh.getBlkPtr(currentBlock, solnData, CENTER);

    // get block limits
		int blkLimits[2][MESH_MDIM];
		int blkLimitsGC[2][MESH_MDIM];
		mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

    // compute total number of cells needed in dyadic storage array
    int nxb = blkLimits[UPPER][IAXIS] - blkLimits[LOWER][IAXIS] + 1;
    int nyb = blkLimits[UPPER][JAXIS] - blkLimits[LOWER][JAXIS] + 1;
    int nzb = blkLimits[UPPER][KAXIS] - blkLimits[LOWER][KAXIS] + 1;

    // store number of cells at each level
    int *iCellsInt = new int[numLevels];
    int *jCellsInt = new int[numLevels];
    int *kCellsInt = new int[numLevels];

    // compute number of cells at each level
    iCellsInt[numLevels-1] = nxb;
    jCellsInt[numLevels-1] = nyb;
    kCellsInt[numLevels-1] = nzb;
    for (int l = 0; i < numLevels; l++)
    {
      iCellsInt[l] = ceil( (double)iCellsInt[numLevels] / (double)pow(2,numLevels-l-1) )
      jCellsInt[l] = ceil( (double)iCellsInt[numLevels] / (double)pow(2,numLevels-l-1) )
      kCellsInt[l] = ceil( (double)iCellsInt[numLevels] / (double)pow(2,numLevels-l-1) )
    }

    // compute total number of cells needed
    int totCells = 0;
    for (int l = 0; i < numLevels; l++)
    {
      // increment counter
      totCells += iCellsInt[l] * jCellsInt[l] * kCellsInt[l];
    }
  
    // declare dyadic grid storage array
    double *cellData = new double[totCells];

		for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++)
    {
      for (int j=blkLimits[LOWER][IAXIS]; j<= blkLimits[UPPER][IAXIS]; j++)
      {
        for (int k=blkLimits[LOWER][IAXIS]; k<= blkLimits[UPPER][IAXIS]; k++)
        {
          // decompose data into dyadic grid data structure
          //solnData[newVarIndex][i][j][k] = 1.1;

          // compute the Q's (Bihari et. al. 1997 eq's 3.15a-c)
          double Qx = ;
          double Qy = ;
          double Qxy = ;

          // compute detail coefficient
	      }
	    }
		}

          //mesh.putBlkData(currentBlock, CENTER, "new ", INTERIOR, startingPos, solnData);
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
