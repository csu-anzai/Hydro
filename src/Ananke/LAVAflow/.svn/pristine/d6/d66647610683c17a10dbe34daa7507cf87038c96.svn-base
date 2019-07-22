#include "TrackClusters.h"

void GetOrigMeshInfo(std::string meshFilename, int &nClusters, double &simTime, std::vector<double> &clusterMasses)
{
    std::vector<std::string> origClusterMeshVar = {"clst", "dens"};
    FlashAmrMesh mesh(meshFilename, origClusterMeshVar);

    simTime = mesh.getTime();

    // only works on a uniform mesh!
    int firstCell[] = {0,0,0};
    double cellVol = mesh.getSingleCellVol(0, INTERIOR, firstCell);

    int clusterIDMeshIndex = mesh.findVarIndex("clst");
    int densMeshIndex = mesh.findVarIndex("dens");

    int numLocalBlocks = mesh.getLocalNumBlocks();
    int numLeafBlocks = 0;
    int *blockList = new int[numLocalBlocks];
    mesh.getListOfBlocks(LEAF, blockList, numLeafBlocks);

    std::vector<int> uniqueClusterIDs;

    clusterMasses.resize(2);

    std::vector<int> numCells(2);

    // find the unique cluster IDs
    for (int lb=0; lb<numLeafBlocks; lb++)
    {
        int currentBlock = blockList[lb];
        double ****solnData;
        //Get pointer to the block's data
        mesh.getBlkPtr(currentBlock, solnData, CENTER);

        int blkLimits[2][MESH_MDIM];
        int blkLimitsGC[2][MESH_MDIM];

        mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

        for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++)
        {
            for (int j=blkLimits[LOWER][JAXIS]; j<= blkLimits[UPPER][JAXIS]; j++)
            {
                for (int k=blkLimits[LOWER][KAXIS]; k<= blkLimits[UPPER][KAXIS]; k++)
                {
                    int thisClusterID = int(round(solnData[clusterIDMeshIndex][i][j][k]));
                    double cellMass = solnData[densMeshIndex][i][j][k] * cellVol;
                    
                    if (thisClusterID > 1)
                    {
                        // only add to the list if not already there
                        if ( std::find(uniqueClusterIDs.begin(), uniqueClusterIDs.end(), thisClusterID) ==  uniqueClusterIDs.end())
                        {
                            uniqueClusterIDs.push_back(thisClusterID);
                            
                            // only resize if this ID makes the array bigger
                            if (thisClusterID > clusterMasses.size()-1)
                            {
                                clusterMasses.resize(thisClusterID+1, 0.);
                                numCells.resize(thisClusterID+1, 0.);
                            }

                            // 

                        }

                        clusterMasses[thisClusterID] += cellMass;
                        numCells[thisClusterID]++;
                    }
                }
            }
        }
    }

    delete[] blockList;

    // std::cout << "Unique cluster IDs: " << '\n';
    // for (auto clusterID : uniqueClusterIDs)
    //     std::cout << clusterID << '\n';

    // std::cout << std::endl;

    // std::cout << "Number of cells per cluster: " << '\n';
    // for (auto numCellsThisCluster : numCells)
    //     std::cout << numCellsThisCluster << '\n';

    // std::cout << std::endl;

    nClusters = uniqueClusterIDs.size();
}
