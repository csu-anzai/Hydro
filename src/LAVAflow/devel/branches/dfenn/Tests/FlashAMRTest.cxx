#include <iostream>
#include <string>
#include <vector>
#include "../libsrc/MPI/LavaFlowMPI.h"
#include "../libsrc/Mesh/FlashAmrMesh.h"

using namespace std;

int main(int argc, char** argv)
{
	
	LavaFlowMPI mpiObj(argc, argv);

	string filename("/home/df11c/yt/data/LocalData/tburn_hdf5_plt_cnt_0006");

	vector<string> meshVars = {"dens", "temp"};
	FlashAmrMesh mesh(filename, meshVars);

	int densIndex = mesh.findVarIndex("dens");
	int tempIndex = mesh.findVarIndex("temp");

	int numLocalBlocks = mesh.getNBlocks();
	int numLeafBlocks = 0;
	int *blockList = new int[numLocalBlocks];

	mesh.getListOfBlocks(ALL_BLKS, blockList, numLeafBlocks);

	double cumDens = 0;
	int cellCount = 0;
	double totalMass = 0;
	
	for (int lb=0; lb<numLeafBlocks; lb++)
	{
		int currentBlock = blockList[lb];
		double ****solnData;
		//Get pointer to the block's data
		mesh.getBlkPtr(currentBlock, solnData, CENTER);

		double cellVol = mesh.getSingleCellVol(currentBlock);

		int blkLimits[2][MESH_MDIM];
		int blkLimitsGC[2][MESH_MDIM];

		mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

		for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++)
			for (int j=blkLimits[LOWER][IAXIS]; j<= blkLimits[UPPER][IAXIS]; j++)
				for (int k=blkLimits[LOWER][IAXIS]; k<= blkLimits[UPPER][IAXIS]; k++)
				{
					double cellMass = solnData[densIndex][i][j][k] * cellVol;
					totalMass += cellMass;
					cumDens += solnData[tempIndex][i][j][k] * cellMass;
					cellCount++;
				}

	}
	
	double avgDens = cumDens / totalMass;
	cout << "Average temperature in domain is: " << avgDens << endl;

	delete blockList;

} 
