#include "FloodFill.h"

#include <iostream>
#include <list>
#include <vector>
#include <deque>

#include "../../Utilities/LAVAUtil.h"
#include "../../includes/LAVAconstants.h"


namespace Clustering
{

void FloodFill::process()
{
	// minor error checking
	if(!this->isDataAttached)
	{
		std::cerr<<"[FloodFill::process] ERROR: No data attached! Please call setData or the proper constructor!"<<std::endl;
		return;
	}

	/* =================================================================
		Loop over all pixels in the image
		If the pixel value is equal to TARGET, we've found the start of a new cluster
		If the pixel value is not equal to TARGET, move to the next pixel

		When we've found a new cluster:
		(1) Create a new list<FloodFill::cellIndices>, cellList
			This will store the cell indices for every cell in the cluster
		(2) Call FloodFill::flood(...) with the identified starting index. 
			cellList should be passed by reference to this function, 
			in order to update during the recursive flood method
		(3)	Store cellList in the FloodFill::clusters vector

	   ================================================================= */
	
	// Reserve 16K of memory that can be deleted just in case we run out of memory
    char* _emergencyMemory = new char[16384];

	int replacementValue = FloodFill::TARGET + 1;
	for(int i=minInds[0]; i<maxInds[0]; i++)
	{
		for(int j=minInds[1]; j<maxInds[1]; j++)
		{
			for(int k=minInds[2]; k<maxInds[2]; k++)
			{
				// Determine if we should continue to the next pixel or begin flooding
				int pixelValue = data[i][j][k];
				if(pixelValue != FloodFill::TARGET)
				{
					continue;
				}

				// std::cout << "Found new cluster" << std::endl;

				// Create new list of indices
				std::deque<FloodFill::cellIndices> cellList;

				// std::cout << "Flooding using replacement value " << replacementValue << std::endl;
                // Call flooding procedure on the current pixel
				this->flood(replacementValue, i, j, k, cellList);

				if (cellList.size() >= minClusterSize || minClusterSize < 0)
                {
                    try
    				{
                        // std::cout << "Adding cluster " << replacementValue << " to cluster array" << std::endl;
    					// Push cellList onto the cluster stack
    					this->clusters.push_back(cellList);
    				}
    				catch (std::bad_alloc ex)
    		        {
    		            // Delete the reserved memory so we can print an error message before exiting
    		            delete[] _emergencyMemory;
    		    
    		            std::cerr << "[FloodFill] Out of memory. Exiting..." << std::endl;
    		            exit(EXIT_FAILURE);
    		        }

    				// Increase the replacement value (lets us keep track of things later if we want to)
    				replacementValue++;

                    // std::cout << "New replacement value is " << replacementValue << std::endl;
                }
                else
                {
                    // set all cells in this cell list to the unused value
                    for (auto thisCell : cellList)
                    {
                        data[thisCell.i][thisCell.j][thisCell.k] = FloodFill::UNUSED;
                    }
                }
			}
		}
	}

	// Delete the reserved memory so we can print an error message before exiting
    delete[] _emergencyMemory;

}

void FloodFill::flood(int replacementValue, int i, int j, int k, std::deque<FloodFill::cellIndices>& cellList)
{
	// Reserve 16K of memory that can be deleted just in case we run out of memory
    char* _emergencyMemory = new char[16384];

	// list of pixels that still need to be checked
	std::deque<int> indicesToCheckI, indicesToCheckJ, indicesToCheckK;

    // Call flood on the neighboring indices
	int currentIndex[] = {i,j,k};
	int neighborIndex[3];
	int offset[3];

    int pixelValue = this->data[i][j][k];
    while (pixelValue == FloodFill::TARGET || indicesToCheckI.size() > 0)
    {
		// Store the cell index in the cellList
		FloodFill::cellIndices index;
		index.i = i;
		index.j = j;
		index.k = k;
		index.linearIndex = Util::sub2ind(i,j,k,maxInds[0],maxInds[1]);
		
        if (pixelValue != replacementValue)
        {
            try
    		{
    			cellList.push_back(index);
    		}
    		catch (std::bad_alloc ex)
            {
                // Delete the reserved memory so we can print an error message before exiting
                delete[] _emergencyMemory;
        
                std::cerr << "[FloodFill] Out of memory. Exiting..." << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        // Update the pixel value to replacementValue
        this->data[i][j][k] = replacementValue;

		// check pixel neighbors
		for (int neighbor=0; neighbor<nDim*2; neighbor++)
		{
			offset[0] = 0; offset[1] = 0; offset[2] = 0;
			if (neighbor < 2)
				offset[0] = neighbor*2 - 1;
			else if (neighbor < 4)
				offset[1] = (neighbor-2)*2 - 1;
			else if (neighbor < 6)
				offset[2] = (neighbor-4)*2 - 1;			

			
			if(this->getNeighborIndex(currentIndex, offset, neighborIndex))
			{
				int neighborI = neighborIndex[0],
				neighborJ = neighborIndex[1],
				neighborK = neighborIndex[2];

				// if the neighbor is part of the cluster, add its indices to the list
				if (data[neighborI][neighborJ][neighborK] == FloodFill::TARGET)
				{
					try
					{
						indicesToCheckI.push_back(neighborI);
						indicesToCheckJ.push_back(neighborJ);
						indicesToCheckK.push_back(neighborK);
					}
					catch (std::bad_alloc ex)
		            {
		                // Delete the reserved memory so we can print an error message before exiting
		                delete[] _emergencyMemory;
		        
		                std::cerr << "[FloodFill] Out of memory. Exiting..." << std::endl;
		                exit(EXIT_FAILURE);
		            }
				}
			}
		}

		// move on to the next pixel in the list
		if (indicesToCheckI.size() > 0)
		{
			i = indicesToCheckI.back();
			j = indicesToCheckJ.back();
			k = indicesToCheckK.back();

			indicesToCheckI.pop_back();
			indicesToCheckJ.pop_back();
			indicesToCheckK.pop_back();

			pixelValue = this->data[i][j][k];
			
			currentIndex[0] = i;
			currentIndex[1] = j;
			currentIndex[2] = k;			
		}
		else
		{
			pixelValue = replacementValue;
			continue;
		}

    }

    // Delete the reserved memory so we can print an error message before exiting
    delete[] _emergencyMemory;
}

bool FloodFill::getNeighborIndex(int curr[3], int offset[3], int neigh[3])
{
	// Set the neighbor value the current index + offset
	neigh[0] = curr[0] + offset[0];
	neigh[1] = curr[1] + offset[1];
	neigh[2] = curr[2] + offset[2];

	// std::cout<<"Initial neighbor index: "<<neigh[0]<<"\t"<<neigh[1]<<"\t"<<neigh[2]<<std::endl;

	// Check the boundary conditions for each direction and either
	// 	(1)	return false if the boundary is fixed and the index exceeds the range
	//			This indicates that the neighbor does not exist and should not be used
	//	(2)	Take the modulus of the exceeded index with the maximum when the BC is periodic

	for(int n=0; n<nDim; n++)
	{
		if(neigh[n] < 0)
		{
			if(this->bcType[n][LOW] == FloodFill::PERIODIC)
			{
				// NOTE: This will break if the initial neighbor index is more than one period away
				neigh[n] = maxInds[n]+neigh[n];
			}
			else
			{
				return false;
			}
		}

		if(neigh[n] >= maxInds[n])
		{
			if(this->bcType[n][HIGH] == FloodFill::PERIODIC)
			{
				neigh[n] = neigh[n]%maxInds[n];
			}
			else
			{
				return false;
			}
		}
	}

	return true;
	// std::cout<<"Final neighbor index: "<<neigh[0]<<"\t"<<neigh[1]<<"\t"<<neigh[2]<<std::endl;

}


int*** FloodFill::BuildMaskArray(std::string filename, const double thresholdVarCutoff, std::vector<std::string> meshVars)
{
    // std::vector<std::string> meshVars = {"dens", "temp", "c12 "};

    const int numVars = meshVars.size();
    std::vector<const char*> meshVarConverted = VecStringToVecChar(meshVars);

    FlashAmrMesh mesh(filename, meshVars);

    // const int densIndex = mesh.findVarIndex("dens");
    // const int tempIndex = mesh.findVarIndex("temp");
    // const int c12Index = mesh.findVarIndex("c12 ");

    int varIndices[MVAR_INDICES];
    for (int i=0; i<numVars; i++)
    	varIndices[i] = mesh.findVarIndex(meshVars[i].c_str());


    nDim = mesh.getDim();

    // 3D array representing the entire domain
    int*** image;

    int *numCells = mesh.getNcells();
    
    // get the total number of blocks along each direction
    
    // number of base-level blocks in each direction and refinement level
    const int numBaseLevelBlocks[] = {mesh.getInt("nblockx"), mesh.getInt("nblocky"), mesh.getInt("nblockz")};
    const int lRefine = mesh.getInt("lrefine_max");

    // number of cells per block
    const int nx[] = {mesh.getInt("nxb"), mesh.getInt("nyb"), mesh.getInt("nzb")};

    // std::cout << "Num base level blocks: ";
    // for (int i=0; i<2; i++)
    //     std::cout << numBaseLevelBlocks[i] << ", ";

    // std::cout << numBaseLevelBlocks[2] << std::endl;
    // std::cout << "Refinement level: " << lRefine << std::endl;

    // number of leaf blocks per base-level block and number of cells in each direction
    const int numLeafBlocksPerBL = pow(2, lRefine - 1);
    // std::cout << "Number of leaf blocks per BL block: " << numLeafBlocksPerBL << std::endl;

    int numLeafBlocks[nDim];
    // int totalNumCells[nDim];

    std::array<int, MESH_MDIM> totalNumCells;

    for (int i=0; i<nDim; i++)
    {
        numLeafBlocks[i] = numBaseLevelBlocks[i] * numLeafBlocksPerBL;
        totalNumCells[i] = numLeafBlocks[i] * nx[i];
        maxInds[i] = totalNumCells[i];
        // std::cout << "Total number of cells in dimension: " << totalNumCells[i] << std::endl;
    }

    // Reserve 16K of memory that can be deleted just in case we run out of memory
    char* _emergencyMemory = new char[16384];
    
    // allocate the image array and set mask value to false
    try
    {
        image = new int**[totalNumCells[IAXIS]];
    }
    catch (std::bad_alloc ex)
    {
        // Delete the reserved memory so we can print an error message before exiting
        delete[] _emergencyMemory;

        std::cerr << "[FloodFillTest] Out of memory while allocating masking array. Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i=0; i<totalNumCells[IAXIS]; i++)
    {
        try
        {
            image[i] = new int*[totalNumCells[JAXIS]];
        }
        catch (std::bad_alloc ex)
        {
            // Delete the reserved memory so we can print an error message before exiting
            delete[] _emergencyMemory;
        
            std::cerr << "[FloodFillTest] Out of memory while allocating masking array. Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
        for (int j=0; j<totalNumCells[JAXIS]; j++)
        {
            try
            {
                image[i][j] = new int[totalNumCells[KAXIS]];
            }
            catch (std::bad_alloc ex)
            {
                // Delete the reserved memory so we can print an error message before exiting
                delete[] _emergencyMemory;
        
                std::cerr << "[FloodFillTest] Out of memory while allocating masking array. Exiting..." << std::endl;
                exit(EXIT_FAILURE);
            }
            for (int k=0; k<totalNumCells[KAXIS]; k++)
                image[i][j][k] = Clustering::FloodFill::UNUSED;
            
        }
    }

    delete[] _emergencyMemory;

    // loop through all the blocks
    // for each block, loop through cells
    // check if cell's ig time is below the cutoff threshold
    // if so, mark it. if not, do nothing
    int numLocalBlocks = mesh.getNBlocks();
    int totalNumLeafBlocks = 0;
    int *blockList = new int[numLocalBlocks];

    mesh.getListOfBlocks(ALL_BLKS, blockList, totalNumLeafBlocks);

    // physical properties: domain bounds and cell side lengths
    // this only works for a uniform mesh (it assumes block 0 is representative of all blocks)
    double domainBoundBox[2][MESH_MDIM];
    mesh.getDomainBoundBox(domainBoundBox);

    double cellSideLengths[MESH_MDIM];
    mesh.getCellSideLengths(blockList[0], cellSideLengths);

    for (int lb=0; lb<totalNumLeafBlocks; lb++)
    {
        int currentBlock = blockList[lb];

        // only do if it's a leaf block
        if (mesh.getNodeType(currentBlock) == 1)
        {
            double ****solnData;
            
            //Get pointer to the block's data
            mesh.getBlkPtr(currentBlock, solnData, CENTER);

            int blkLimits[2][MESH_MDIM];
            int blkLimitsGC[2][MESH_MDIM];
            mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

            // physical coordinates of the block's lower bounds
            // double *blockBounds = mesh.getBounds(currentBlock, LOWER);

            for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++)
                for (int j=blkLimits[LOWER][IAXIS]; j<= blkLimits[UPPER][IAXIS]; j++)
                    for (int k=blkLimits[LOWER][IAXIS]; k<= blkLimits[UPPER][IAXIS]; k++)
                    {
                        // get the indices of the corresponding cell in the image array
                        double cellBoundBox[2][MESH_MDIM];
                        mesh.getCellBoundBox(currentBlock, INTERIOR, i, j, k, cellBoundBox);

                        int indicesImageArray[nDim];

                        for (int l=0; l<nDim; l++)
                        {
                            indicesImageArray[l] = lrint((cellBoundBox[LOWER][l] - domainBoundBox[LOWER][l]) / cellSideLengths[l]);
                        }

                        bool markCell = criterionFunc(solnData, thresholdVarCutoff, varIndices, i, j, k);

                        // if the ignition time is less than the cutoff value, mark the image array
                        if (markCell)
                        {
                            int imageI = indicesImageArray[0], 
                                imageJ = indicesImageArray[1],
                                imageK = indicesImageArray[2];
                            
                            image[imageI][imageJ][imageK] = Clustering::FloodFill::TARGET;

                            // cout << "Marking cell as target" << endl;
                        }

                    }
        }

    }

    // cleanup
    delete[] blockList;

    return image;
}


}; // end namespace