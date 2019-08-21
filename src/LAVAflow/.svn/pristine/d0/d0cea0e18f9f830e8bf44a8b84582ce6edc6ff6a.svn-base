#include <cstdlib>
#include <iostream>
#include "FlashAmrMesh.h"
#include "Driver_main.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

//spatial dimensions of data
const int numDim = 2; 

//number of levels of transform
const int numLevels = 3; 

  // input file name
const std::string inFileName("/home/wkshop24/LAVAflow/data/wdStir_sol_data/wdStir_hdf5_chk_0007");

//calculates linear index
int imap(int l, int i, int j, int k, int *nx,int *ny, int *nz) {

	int map = 0;

	for(int c = 0; c<l; c++)
	{
		map += nx[c]* ny[c] * nz[c];
	}
	map += k + nz[l] * (j + ny[l] * i);
	//map += i;
	return map;
}
int imap(int l, int i, int j, int *nx, int *ny){
	int map = 0;

	for(int c = 0; c<l; c++)
	{
		map += nx[c]* ny[c];
	}
	map += (j + ny[l] * i);
	//map += i;
	return map;
}
int imap(int l, int i, int *nx){
	int map = 0;
	for(int c = 0; c<l; c++)
	{
		map += nx[c];
	}
	map += (ny[l] * i);
	//map += i;
	return map;
}


int LAVAFLOW_DRIVER_BASENM(int argc, char **argv) {

  // primitive variables to perform transform on
	std::vector<std::string> meshVars = {"dens"};
	std::vector<std::string> auxVars  = {};
  // specify output file names
	char outAsciiName[MAX_FILENAME_LENGTH] = "sedovAscii";
	char outHDF5Name[MAX_FILENAME_LENGTH] = "sedovHDF5_transfrm";
	char outProgramHDF5Name[MAX_FILENAME_LENGTH] = "./";
	strcat(outProgramHDF5Name,outHDF5Name);

  	// create AMR mesh object
	FlashAmrMesh mesh(inFileName, outProgramHDF5Name, meshVars, auxVars);

	// compute mapping for dyadic data structure

  	// get number of local blocks
	int numLocalBlocks = mesh.getLocalNumBlocks();

  	// create list of leaf blocks
	int numLeafBlocks = 0;
	int *blockList = new int[numLocalBlocks];
    double* coords;
    coords = new double [8]; 
	mesh.getListOfBlocks(LEAF, blockList, numLeafBlocks);

  	// loop through cells in block
	int startingPos[MESH_MDIM];
	startingPos[0] = 0;
	startingPos[1] = 0;
	startingPos[2] = 0;
	for (int lb=0; lb<numLeafBlocks; lb++)
	{
		//std::cout << lb << std::endl;
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

	    //store number of ghost cells
	    int ngx = blkLimits[LOWER][IAXIS] - blkLimitsGC[LOWER][IAXIS];
	    int ngy = blkLimits[LOWER][JAXIS] - blkLimitsGC[LOWER][JAXIS];
	    int ngz = blkLimits[LOWER][KAXIS] - blkLimitsGC[LOWER][KAXIS];
	    
	    // store number of cells at each level
	    int *iCellsInt = new int[numLevels];
	    int *jCellsInt = new int[numLevels];
	    int *kCellsInt = new int[numLevels];

	    // compute number of cells at each level
	    iCellsInt[numLevels-1] = nxb;
	    jCellsInt[numLevels-1] = nyb;
	    kCellsInt[numLevels-1] = nzb;

	    for (int l = 0; l < numLevels; l++)
	    {
	      iCellsInt[l] = ceil( (double)iCellsInt[numLevels-1] / (double)pow(2,numLevels-l-1) );
	      jCellsInt[l] = ceil( (double)iCellsInt[numLevels-1] / (double)pow(2,numLevels-l-1) );
	      kCellsInt[l] = ceil( (double)iCellsInt[numLevels-1] / (double)pow(2,numLevels-l-1) );
	    }

	    // compute total number of cells needed
	    int totCells = 0;

	    for (int l = 0; l < numLevels; l++)
	    {
	    	totCells += iCellsInt[l] * jCellsInt[l] * kCellsInt[l];
	    }
	  	

	    // declare dyadic grid storage array
	    double *cellData = new double[totCells];
	    double *detailCoeff;

	    switch (numDim){
	    	case 1:
	    		detailCoeff = new double[totCells - iCellsInt[numLevels-1]];
		    	for(int v = 0; v < meshVars.size(); v++)
		    	{
		    	//the variable index in mesh 
			    	int var = mesh.findVarIndex(meshVars[v].c_str()); 
					for (int i = 0; i < iCellsInt[numLevels-1]; i++)
				    {
				      int ii = imap(numLevels-1, i, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      cellData[ii] = solnData[var][i+ngx][blkLimits[LOWER][JAXIS]][blkLimits[LOWER][KAXIS]];
					}
					// average the fine-grid data into our dyadic-grid data structure (1D only)
					for(int l = numLevels-2; l >= 0; l--)
					{	
						for (int i = 0; i < iCellsInt[l]; i++)
				    	{	
				      		int ii = imap(l, i, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      		int ii1 = imap(l+1, 2*i, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      		int ii2 = imap(l+1, 2*i+1, 0, 0, iCellsInt, jCellsInt, kCellsInt);

				      		cellData[ii] = .5 * (cellData[ii1] + cellData[ii2]);
				      		//std::cout << cellData[ii] << std::endl;
				      		//std::cout << "l = " << l << " " << "i = " << i << std::endl;			}
						}
					}
					//computing detail coefficents
					for(int l = numLevels-2; l >= 0; l--)
					{	
						for (int i = 0; i < iCellsInt[l]; i++)
				    	{	
				      		int ii = imap(l+1, 2*i+1, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      		//int ii1 = imap(l+1, 2*i, iCellsInt);
				      		//int ii2 = imap(l+1, 2*i+1, iCellsInt);
				      		int ii1, ii2, ii3;
				      		double ans;
				      		if(i == 0){
				      		//left bias interpolation
				      			ii1 = imap(l, i, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ii2 = imap(l, i+1, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ii3 = imap(l, i+2, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ans = 5*cellData[ii1]/8.0 + .5*cellData[ii2] - cellData[ii3]/8.0;
				      		}
				      		else if(i == iCellsInt[l]-1){
				      		//right bias interpolation
				      			ii1 = imap(l, i, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ii2 = imap(l, i-1, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ii3 = imap(l, i-2, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ans = 11*cellData[ii1]/8.0 - .5*cellData[ii2] + cellData[ii3]/8.0;
				      		}
				      		else{
				      			ii1 = imap(l, i, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ii2 = imap(l, i-1, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ii3 = imap(l, i+1, 0, 0, iCellsInt, jCellsInt, kCellsInt);
				      			ans = cellData[ii1] - (cellData[ii2]-cellData[ii3])/8.0;

				      		}
				      		detailCoeff[ii1] = cellData[ii] - ans;
				      		//detailCoeff[ii] = cellData[ii1] - cellData[ii2];
				      		//std::cout << cellData[ii1] << "-" << cellData[ii2] << " " << detailCoeff[ii] << std::endl;
				      		//std::cout << detailCoeff[ii1] << " " << cellData[ii] << " " << ans << std::endl;
				      		//std::cout << "l = " << l << " " << "i = " << i << std::endl;			}
						}
					}
				}
	    	

	    	case 2:
				// decompose the fine-grid data into our dyadic-grid data structure (2D only)
	    		detailCoeff = new double[totCells];
			    for(int v = 0; v < meshVars.size(); v++){
			    	//the variable index in mesh 
			    	int var = mesh.findVarIndex(meshVars[v].c_str()); 
					for (int i = 0; i < iCellsInt[numLevels-1]; i++)
				    {
				    	for(int j = 0; j< jCellsInt[numLevels-1]; j++)
				    	{

							int ii = imap(numLevels-1, i, j, 0, iCellsInt, jCellsInt, kCellsInt);
							cellData[ii] = solnData[var][i+ngx][j+ngy][blkLimits[LOWER][KAXIS]];
				    	}
				      
					}
					// average the fine-grid data into our dyadic-grid data structure (2D only)
					for(int l = numLevels-2; l >= 0; l--)
					{	
						for (int i = 0; i < iCellsInt[l]; i++)
				    	{	
				    		for(int j = 0; j<jCellsInt[l]; j++)
				    		{
					      		int ii = imap(l, i, j, 0, iCellsInt, jCellsInt, kCellsInt);
					      		int ii1 = imap(l+1, 2*i, 2*j, 0, iCellsInt, jCellsInt, kCellsInt);
					      		int ii2 = imap(l+1, 2*i+1, 2*j, 0, iCellsInt, jCellsInt, kCellsInt);
					      		int ii3 = imap(l+1, 2*i, 2*j+1, 0, iCellsInt, jCellsInt, kCellsInt);
					      		int ii4 = imap(l+1, 2*i+1, 2*j+1, 0, iCellsInt, jCellsInt, kCellsInt);

					      		cellData[ii] =  (cellData[ii1] + cellData[ii2] + cellData[ii3] + cellData[ii4])/4.0;
					      		std::cout << cellData[ii] << std::endl;
					      		//std::cout << "l = " << l << " " << "i = " << i << std::endl;			}
				      		}
						}
					}
			

						//computing detail coefficents
					for(int l = numLevels-2; l >= 0; l--)
					{	

						for (int i = 0; i < iCellsInt[l]; i++)
				    	{	
				    		for(int j = 0; j<jCellsInt[l]; j++)
				    		{
				    			double qx, qy, qxy;

								int ii = imap(l+1, 2*i+1, iCellsInt, iCellsInt, jCellsInt,kCellsInt);
					      		//int ii1 = imap(l+1, 2*i, iCellsInt);
					      		//int ii2 = imap(l+1, 2*i+1, iCellsInt);
					      		int ii1, ii2, ii3, ii4, ii5, ii6, ii7, ii8;
					      		double ans;
					      		if(i == 0){
					      		//left bias interpolation
					      			ii1 = imap(l, i, j, 0, iCellsInt, jCellsInt, kCellsInt);
					      			ii2 = imap(l, i+1, j, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ii3 = imap(l, i+2, j, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ans = 5*cellData[ii1]/8.0 + .5*cellData[ii2] - cellData[ii3]/8.0;
					      		}
					      		else if(j ==0){

					      		}
					      		else if(i == iCellsInt[l]-1){
					      		//right bias interpolation
					      			ii1 = imap(l, i, j, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ii2 = imap(l, i-1, j, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ii3 = imap(l, i-2, j, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ans = 11*cellData[ii1]/8.0 - .5*cellData[ii2] + cellData[ii3]/8.0;
					      		}
					      		else if(j == jCellsInt[l]-1){

					      		}
					      		else{
					      			//calculation of qx
					      			ii1 = imap(l, i+1, j, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ii2 = imap(l, i-1, j, 0, iCellsInt, jCellsInt, kCellsInt);
					      			qx = (cellData[ii1] - cellData[ii2])/8.0;
					      			//calculation of qy
					      			ii3 = imap(l, i, j+1, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ii4 = imap(l, i, j-1, 0, iCellsInt, jCellsInt, kCellsInt);
					      			qy = (cellData[ii3] - cellData[ii4]) / 8.0;
					      			//calculation of qxy
					      			ii5 = imap(l, i+1, j+1, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ii6 = imap(l, i+1, j-1, 0, iCellsInt, jCellsInt, kCellsInt);
					      			ii7 = imap(l, i-1, j+1, 0, iCellsInt, jCellsInt,kCellsInt);
					      			ii8 = imap(l, i-1, j-1, 0, iCellsInt, jCellsInt, kCellsInt);
					      			qxy = (cellData[ii5] - cellData[ii6] - cellData[ii7] + cellData[ii8])/(64.0);

					      		}
					      		detailCoeff[ii1] = cellData[ii] - ans;
					      		//detailCoeff[ii] = cellData[ii1] - cellData[ii2];
					      		//std::cout << cellData[ii1] << "-" << cellData[ii2] << " " << detailCoeff[ii] << std::endl;
					      		std::cout << detailCoeff[ii1] << " " << cellData[ii] << " " << ans << std::endl;
					      		//std::cout << "l = " << l << " " << "i = " << i << std::endl;			}
				    		}
				      		
						}
					}


				}
		}
	    // cleanup memory
	    delete iCellsInt;
	    delete jCellsInt;
	    delete kCellsInt;
	    delete cellData;
	    delete detailCoeff;
	}

          //mesh.putBlkData(currentBlock, CENTER, "new ", INTERIOR, startingPos, solnData);
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
	
	return 0;
}

