#include "FloodFill.h"

#include <iostream>
#include <list>
#include <vector>

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

				// Create new list of indices
				std::list<FloodFill::cellIndices> cellList;

				// Call flooding procedure on the current pixel
				this->flood(replacementValue, i, j, k, cellList);

				// Push cellList onto the cluster stack
				this->clusters.push_back(cellList);

				// Increase the replacement value (let's us keep track of things later if we want to)
				replacementValue++;
			}
		}
	}

}

void FloodFill::flood(int replacementValue, int i, int j, int k, std::list<FloodFill::cellIndices>& cellList)
{
	// Determine if we should end this branch of the flood
	// This happens when the currently evaluated pixel is not equal to the target value
	int pixelValue = this->data[i][j][k];
	if(pixelValue != FloodFill::TARGET)
	{
		return;
	}

	// Update the pixel value to replacementValue
	this->data[i][j][k] = replacementValue;

	// Store the cell index in the cellList
	FloodFill::cellIndices index;
	index.i = i;
	index.j = j;
	index.k = k;
	index.linearIndex = Util::sub2ind(i,j,k,maxInds[0],maxInds[1]);
	cellList.push_back(index);

	// Call flood on the neighboring indices
	int currentIndex[] = {i,j,k};
	int neighborIndex[3];
	int offset[3];
	
	// West
	offset[0] = -1;	offset[1] = 0;	offset[2] = 0;
	if(this->getNeighborIndex(currentIndex, offset, neighborIndex))
	{
		this->flood(replacementValue, neighborIndex[0], neighborIndex[1], neighborIndex[2], cellList);
	}	

	// East
	offset[0] = 1;	offset[1] = 0;	offset[2] = 0;
	if(this->getNeighborIndex(currentIndex, offset, neighborIndex))
	{
		this->flood(replacementValue, neighborIndex[0], neighborIndex[1], neighborIndex[2], cellList);
	}	

	if(nDim>1)
	{
		// South
		offset[0] = 0;	offset[1] = -1;	offset[2] = 0;
		if(this->getNeighborIndex(currentIndex, offset, neighborIndex))
		{
			this->flood(replacementValue, neighborIndex[0], neighborIndex[1], neighborIndex[2], cellList);
		}	

		// North
		offset[0] = 0;	offset[1] = 1;	offset[2] = 0;
		if(this->getNeighborIndex(currentIndex, offset, neighborIndex))
		{
			this->flood(replacementValue, neighborIndex[0], neighborIndex[1], neighborIndex[2], cellList);
		}	
	}

	if(nDim==3)
	{
		// Bottom
		offset[0] = 0;	offset[1] = 0;	offset[2] = -1;
		if(this->getNeighborIndex(currentIndex, offset, neighborIndex))
		{
			this->flood(replacementValue, neighborIndex[0], neighborIndex[1], neighborIndex[2], cellList);
		}	

		// Top
		offset[0] = 0;	offset[1] = 0;	offset[2] = 1;
		if(this->getNeighborIndex(currentIndex, offset, neighborIndex))
		{
			this->flood(replacementValue, neighborIndex[0], neighborIndex[1], neighborIndex[2], cellList);
		}	
	}

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
				neigh[n] = maxInds[n]-1+neigh[n];
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




}; // end namespace