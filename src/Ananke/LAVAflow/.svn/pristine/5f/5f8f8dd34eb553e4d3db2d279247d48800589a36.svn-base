#include "TimeAutoCorr.h" 
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>



TimeAutoCorr::TimeAutoCorr(std::string baseFilename, int fileIndexBounds[2], int numPoints, MPI::Intracomm meshComm)
{
	this->meshComm = meshComm;
	this->baseFilename = baseFilename;
	this->fileIndexBounds[0] = fileIndexBounds[0];
	this->fileIndexBounds[1] = fileIndexBounds[1];
	
	this->numPoints = numPoints / meshComm.Get_size() + 
		int(numPoints % meshComm.Get_size() > meshComm.Get_rank());
	
	numTimeSeps = fileIndexBounds[1] - fileIndexBounds[0];
	filenames = new std::string[numTimeSeps+1];

	for (int i=fileIndexBounds[0], index=0; i <= fileIndexBounds[1]; i++, index++)
	{
		// create the filenames
		std::ostringstream fileNumStr;
    	fileNumStr << std::setw( 4 ) << std::setfill( '0' ) << i;
    	std::string filename = baseFilename + fileNumStr.str();
    	
    	filenames[index] = filename;
	}
	
	timeSeps.resize(numTimeSeps);
	autoCorrelations.resize(numTimeSeps);

	points.resize(this->numPoints);
	pointBlockIDs = new int[this->numPoints];
	pointZoneIDs.reserve(this->numPoints);
}

TimeAutoCorr::~TimeAutoCorr()
{
	delete[] filenames;
	delete[] pointBlockIDs;
}

void TimeAutoCorr::GetAutoCorrelations()
{
	// load initial mesh
	std::vector<std::string> meshVars = {"velx", "vely", "velz"};
	FlashAmrMesh initMesh(filenames[0], meshVars);
	double initSimTime = initMesh.getTime();
	GeneratePoints(initMesh);
	
	// do each dimension separately
	for (int dim=0; dim<MESH_MDIM; dim++)
	{
		int meshVarIndex = initMesh.findVarIndex(meshVars[dim].c_str());
	
		// save the values at each point for the original mesh
		double *initPointVals = new double[numPoints];
		GetQuantityAtPoints(meshVarIndex, initMesh, initPointVals);
		
		// for each separation
		for (int i = 0; i < numTimeSeps; i++)
		{
			// load second mesh
			std::vector<std::string> varName = {meshVars[dim]};
			FlashAmrMesh currentMesh(filenames[i+1], varName);

			// get velocities at points
			double *currentPointVals = new double[numPoints];
			GetQuantityAtPoints(0, currentMesh, currentPointVals);

			// Average of above normalized to vel^2
			double cumProduct = 0;
			double cumVelSquared = 0;
			for (int pointIndex=0; pointIndex<numPoints; pointIndex++)
			{
				cumProduct += initPointVals[pointIndex]*currentPointVals[pointIndex];
				cumVelSquared += initPointVals[pointIndex]*initPointVals[pointIndex];
			}

			autoCorrelations[i][dim] = cumProduct / cumVelSquared;
			timeSeps[i][dim] = currentMesh.getTime() - initSimTime;

			delete[] currentPointVals;
		}

		delete[] initPointVals;
	}
}

void TimeAutoCorr::GetQuantityAtPoints(int meshVar, FlashAmrMesh &mesh, double pointVals[])
{
	int numLocalBlocks = mesh.getLocalNumBlocks();
	int numLeafBlocks = 0;
	int *blockList = new int[numLocalBlocks];

	mesh.getListOfBlocks(ALL_BLKS, blockList, numLeafBlocks);

	for (int pointIndex=0; pointIndex<numPoints; pointIndex++)
	{
		int currentBlock = blockList[pointBlockIDs[pointIndex]];
		int currentPoint[MESH_MDIM];

		for (int i=0; i<MESH_MDIM; i++)
			currentPoint[i] = pointZoneIDs[pointIndex][i];

		double ****solnData;
		//Get pointer to the block's data
		mesh.getBlkPtr(currentBlock, solnData, CENTER);

		pointVals[pointIndex] = solnData[meshVar][currentPoint[IAXIS]][currentPoint[JAXIS]][currentPoint[KAXIS]];
	}

	delete[] blockList;
}


void TimeAutoCorr::GeneratePoints(FlashAmrMesh &mesh)
{
	// Pick a random block
	int numLocalBlocks = mesh.getLocalNumBlocks();

	// RNG
 	std::random_device rd;

 	// use random number to seed the prng
    std::mt19937 gen(rd());

    // define the properties of the distribution
    std::uniform_int_distribution<> dis(0, numLocalBlocks-1);

    int numLeafBlocks = 0;
	int *blockList = new int[numLocalBlocks];
	mesh.getListOfBlocks(ALL_BLKS, blockList, numLeafBlocks);

    for (int pointIndex=0; pointIndex<numPoints; pointIndex++)
    {
	    int randNum = dis(gen);
	    int randBlock = blockList[randNum];

	    // save the block ID
	    pointBlockIDs[pointIndex] = randBlock;

		// Generate random point in block
		// Get the block's bounds
		double blockBounds[2][MESH_MDIM];
	    mesh.getBlkBoundBox(randBlock, blockBounds);

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
		mesh.getBlkIndexLimits(randBlock, blkLimits, blkLimitsGC);

		double cellSideLengths[MESH_MDIM];
		mesh.getCellSideLengths(randBlock, cellSideLengths);

		for (int i=0; i<MESH_MDIM; i++)
			pointZoneIDs[pointIndex][i] = int(ptOffset[i] / cellSideLengths[i]) + blkLimits[LOWER][i];

		int intCellID[3];
		intCellID[IAXIS] = pointZoneIDs[pointIndex][IAXIS] - blkLimits[LOWER][IAXIS];
		intCellID[JAXIS] = pointZoneIDs[pointIndex][JAXIS] - blkLimits[LOWER][JAXIS];
		intCellID[KAXIS] = pointZoneIDs[pointIndex][KAXIS] - blkLimits[LOWER][KAXIS];
		double boundBox[2][MESH_MDIM];
		mesh.getCellBoundBox(randBlock, INTERIOR, intCellID[IAXIS], intCellID[JAXIS], intCellID[KAXIS], boundBox);

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
