#include <iostream>
#include <string>
#include <vector>

#include "FlashAmrMesh.h"
#include "Driver_main.h"
using namespace std;

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{

	//LavaMPI mpiObj(argc, argv);

    /* Replace this with the location of cartesianShells on your machine */
	string inFileName("/home/pb12c/Research/LAVAflow/drivers/files/testData/cartesianShells/cartesianShells_hdf5_chk_0000");
    string outFileName("/home/pb12c/Research/LAVAflow/build/drivers/devel/FlashAMRTest/test_0000");

    /* If all variables present in the file are to be added, use only 'allVars' */
    //vector<string> meshVars = {"allVars"};

	vector<string> meshVars = {"vari", "rand"};
	// vector<string> meshVars = {"vari"};
	// vector<string> auxVars  = {"new" , "n2"};
	vector<string> auxVars  = {"new"};

	//FlashAmrMesh mesh(inFileName, outFileName, meshVars, auxVars);
	FlashAmrMesh mesh(inFileName, outFileName, meshVars, auxVars);
    
    // None of these variables make sense except for rand...
	int variIndex = mesh.findVarIndex("vari");   //radial distance from center?
	int randIndex = mesh.findVarIndex("rand");   //random values between -1 and 1 for each cell
	int newVarIndex = mesh.findVarIndex("new"); //new variable
	// int newVarIndex2 = mesh.findVarIndex("n2"); 

	// mesh.noWrite("rand");
	// mesh.noWrite("vari");
	// mesh.noWrite("n2");

	//cout << variIndex << endl;
	//cout << randIndex << endl;
	//cout << newVarIndex << endl;

	// int numLocalBlocks = mesh.getNBlocks();
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
				for (int k=blkLimits[LOWER][IAXIS]; k<= blkLimits[UPPER][IAXIS]; k++)
				{
					cumVari += solnData[variIndex][i][j][k];
					cumRand += solnData[randIndex][i][j][k];

					// Note that only the leaf blocks are being updated here.
					solnData[newVarIndex][i][j][k] = 1.1;
	        	}
	    	}
		}

		mesh.putBlkData(currentBlock, CENTER, "new ", INTERIOR, startingPos, solnData);
	}

	double avgNewVar = 0;

	for (int lb = 0; lb < numLeafBlocks; lb++) {
		int currentBlock = blockList[lb];
		double ****solnData;
		mesh.getBlkPtr(currentBlock, solnData, CENTER);

		int blkLimits[2][MESH_MDIM];
		int blkLimitsGC[2][MESH_MDIM];
		mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

		for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++){
			for (int j=blkLimits[LOWER][IAXIS]; j<= blkLimits[UPPER][IAXIS]; j++){
				for (int k=blkLimits[LOWER][IAXIS]; k<= blkLimits[UPPER][IAXIS]; k++)
				{
					avgNewVar += solnData[newVarIndex][i][j][k];
					cellCount++;
					// std::cout << solnData[newVarIndex][i][j][k];
	        	}
	        	// std::cout << endl;
	    	}
	    	// std::cout << endl;
		}
		// std::cout << endl;
	}

	avgNewVar = avgNewVar / cellCount;

    //This metric is as good as anything else...
	std::cout << "Cumulative Distances: " << cumVari << std::endl;
	std::cout << "Cumulative Random Values: " << cumRand << std::endl;

	std::cout << "Cell count: " << cellCount << std::endl;
	std::cout << "Average New Variable: " << avgNewVar << std::endl;

	mesh.writeOutFile();

	delete[] blockList;
}
