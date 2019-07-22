#include "EulTimeAutoCorr.h"


using namespace std;

EulTimeAutoCorr::EulTimeAutoCorr(std::string baseFilename, const char** varNames, int nVars, int fileIndexBounds[2], int numPoints, MPI::Intracomm meshComm)
    : TimeAutoCorr(baseFilename, varNames, nVars, fileIndexBounds, meshComm)
{    
    
    // figure out how many points to generate based on the MPI rank
    numSamples = numPoints / meshComm.Get_size() + 
        int(numPoints % meshComm.Get_size() > meshComm.Get_rank());

    // set aside space at sodosopa
    points.resize(numSamples);
    pointBlockIDs = new int[numSamples];
    pointZoneIDs.reserve(numSamples);
}


EulTimeAutoCorr::~EulTimeAutoCorr()
{
    delete[] pointBlockIDs;
}




int EulTimeAutoCorr::GetQuantityAtPoints(int meshVar, std::shared_ptr<FlashIO> mesh, double pointVals[], int maskID, std::string maskVar)
{
    std::shared_ptr<FlashAmrMesh> amrMesh = std::dynamic_pointer_cast<FlashAmrMesh>(mesh);
    int numLocalBlocks = amrMesh->getLocalNumBlocks();
    int numLeafBlocks = 0;
    int *blockList = new int[numLocalBlocks];

    amrMesh->getListOfBlocks(ALL_BLKS, blockList, numLeafBlocks);

    int numPts = 0;

    for (int pointIndex=0; pointIndex<numSamples; pointIndex++)
    {
        int currentBlock = blockList[pointBlockIDs[pointIndex]];
        int currentPoint[MESH_MDIM];

        for (int i=0; i<MESH_MDIM; i++)
            currentPoint[i] = pointZoneIDs[pointIndex][i];

        double ****solnData;
        //Get pointer to the block's data
        amrMesh->getBlkPtr(currentBlock, solnData, CENTER);

        pointVals[pointIndex] = solnData[meshVar][currentPoint[IAXIS]][currentPoint[JAXIS]][currentPoint[KAXIS]];

        numPts++;
    }

    delete[] blockList;

    return numPts;
}


void EulTimeAutoCorr::GeneratePoints(std::shared_ptr<FlashIO> mesh)
{
    // std::cout << "Generating pts" << std::endl;
    std::shared_ptr<FlashAmrMesh> amrMesh = std::dynamic_pointer_cast<FlashAmrMesh>(mesh);

    // Pick a random block
    int numLocalBlocks = amrMesh->getLocalNumBlocks();

    // RNG
    std::random_device rd;

    // use random number to seed the prng
    std::mt19937 gen(rd());

    // define the properties of the distribution
    std::uniform_int_distribution<> dis(0, numLocalBlocks-1);

    int numLeafBlocks = 0;
    int *blockList = new int[numLocalBlocks];
    amrMesh->getListOfBlocks(ALL_BLKS, blockList, numLeafBlocks);

    for (int pointIndex=0; pointIndex<numSamples; pointIndex++)
    {
        int randNum = dis(gen);
        int randBlock = blockList[randNum];

        // save the block ID
        pointBlockIDs[pointIndex] = randBlock;

        // Generate random point in block
        // Get the block's bounds
        double blockBounds[2][MESH_MDIM];
        amrMesh->getBlkBoundBox(randBlock, blockBounds);

        // generate random points between the bounds
        std::uniform_real_distribution<> realDistX(blockBounds[LOWER][IAXIS], blockBounds[UPPER][IAXIS]),
        realDistY(blockBounds[LOWER][JAXIS], blockBounds[UPPER][JAXIS]),
        realDistZ(blockBounds[LOWER][KAXIS], blockBounds[UPPER][KAXIS]);

        points[pointIndex][IAXIS] = realDistX(gen);
        points[pointIndex][JAXIS] = realDistY(gen);
        points[pointIndex][KAXIS] = realDistZ(gen);

        // get the cell IDs for the points
        double ptOffset[MESH_MDIM];
        for (int i=0; i<MESH_MDIM; i++)
            ptOffset[i] = points[pointIndex][i] - blockBounds[LOWER][i];

        // get the beginning index for interior cells
        int blkLimits[2][MESH_MDIM];
        int blkLimitsGC[2][MESH_MDIM];
        amrMesh->getBlkIndexLimits(randBlock, blkLimits, blkLimitsGC);

        double cellSideLengths[MESH_MDIM];
        amrMesh->getCellSideLengths(randBlock, cellSideLengths);

        for (int i=0; i<MESH_MDIM; i++)
            pointZoneIDs[pointIndex][i] = int(ptOffset[i] / cellSideLengths[i]) + blkLimits[LOWER][i];

        int intCellID[3];
        intCellID[IAXIS] = pointZoneIDs[pointIndex][IAXIS] - blkLimits[LOWER][IAXIS];
        intCellID[JAXIS] = pointZoneIDs[pointIndex][JAXIS] - blkLimits[LOWER][JAXIS];
        intCellID[KAXIS] = pointZoneIDs[pointIndex][KAXIS] - blkLimits[LOWER][KAXIS];
        double boundBox[2][MESH_MDIM];
        amrMesh->getCellBoundBox(randBlock, INTERIOR, intCellID[IAXIS], intCellID[JAXIS], intCellID[KAXIS], boundBox);

        for (int j=0; j<MESH_MDIM; j++)
            if (points[pointIndex][j] < boundBox[LOWER][j] || points[pointIndex][j] > boundBox[UPPER][j])
            {
                std::cout << "Point does not fall within cell bounds for dimension " << j << std::endl;
                std::cout << "Point: " << points[pointIndex][j] << std::endl;
                std::cout << "Block bounds: " << blockBounds[LOWER][j] << ", " << blockBounds[UPPER][j] << std::endl;
            }

    }
    delete[] blockList;

}

std::shared_ptr<FlashIO> EulTimeAutoCorr::LoadData(std::string &filename, std::vector<std::string> &meshVars, std::string maskVar, int maskID)
{
    // std::cout << "Loading mesh" << endl;

    std::shared_ptr<FlashIO> meshPtr(new FlashAmrMesh(filename, meshVars));
    double simTime = meshPtr->getTime();

    // cout << "Sim time in LoadData: " << simTime << endl;
    return meshPtr;
}
