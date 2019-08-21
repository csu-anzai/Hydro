#include "FlashAmrMesh.h"
#include "IndexOperations.h"

// returns a vector mapping a row-major linear block index to a block ID 
// this assumes uniform refinement
std::vector<int> BuildBlockLinearIndexMap(FlashAmrMesh& mesh)
{
    int numLocalBlocks = mesh.getLocalNumBlocks();
    int numLeafBlocks = 0;
    int *blockList = new int[numLocalBlocks];
    
    mesh.getListOfBlocks(LEAF, blockList, numLeafBlocks);

    double domainBoundBox[2][MESH_MDIM];
    mesh.getDomainBoundBox(domainBoundBox);

    
    double *blockSideLengths = mesh.getSizes(blockList[0]);

    double blockBoundBox[2][MESH_MDIM];

    int nBlocksI = mesh.getInt("nblockx"),
        nBlocksJ = mesh.getInt("nblocky"),
        nBlocksK = mesh.getInt("nblockz");



    std::vector<int> blockIndexMap(numLeafBlocks);

    for (int lb=0; lb<numLeafBlocks; lb++)
    {
        int currentBlock = blockList[lb];
        int blkLimits[2][MESH_MDIM];
        int blkLimitsGC[2][MESH_MDIM];

        mesh.getBlkBoundBox(currentBlock, blockBoundBox);

        int iBlock = int(blockBoundBox[LOWER][IAXIS] - domainBoundBox[LOWER][IAXIS]) / (blockSideLengths[IAXIS]),
            jBlock = int(blockBoundBox[LOWER][JAXIS] - domainBoundBox[LOWER][JAXIS]) / (blockSideLengths[JAXIS]),
            kBlock = int(blockBoundBox[LOWER][KAXIS] - domainBoundBox[LOWER][KAXIS]) / (blockSideLengths[KAXIS]);

        int blockLinearInd = Util::sub2ind(iBlock, jBlock, kBlock, nBlocksJ, nBlocksK);
        blockIndexMap[blockLinearInd] = currentBlock;
    }

    delete[] blockList;

    return blockIndexMap;
}

int GetLinearBlockIndex(int i, int j, int k, int nCellBlock, int nBlocksDim)
{
    return nBlocksDim * (int(i/nCellBlock) * nBlocksDim + int(j/nCellBlock)) + int(k/nCellBlock);
}