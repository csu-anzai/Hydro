#ifndef LAVA_CLUSTERING_FLOODFILL_H
#define LAVA_CLUSTERING_FLOODFILL_H

#include <list>
#include <ostream>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkDataArray.h>

#include "../../includes/LAVAconstants.h"

namespace Clustering
{

class FloodFill
{
public:

	enum pixelValue 
	{
		UNUSED = 0, 
		TARGET = 1
	};
	enum boundaryCondition
	{
		FIXED = 0,
		PERIODIC = 1
	};


	// Storage for the lists of cell indices comprising clusters
	struct cellIndices
	{
		int i,j,k,linearIndex;
	};


	FloodFill()
	{
		this->init_clean();
	}

	FloodFill(int** data, int nI, int nJ)
	{
		this->init_clean();
		this->setData(data,nI,nJ);
	}

	FloodFill(int*** data, int nI, int nJ, int nK)
	{
		this->init_clean();
		this->setData(data,nI,nJ,nK);
	}

	~FloodFill()
	{
		if(!(this->data==NULL))
		{
			delete[] data;
		}
	}

	// void reset();
	void process();
	// void print(std::ostream& out, int indent = 0);
	// void printBlocks(std::ostream& out);
	// void printConnectivity(std::ostream& out);

	void setData(int** data, int nI, int nJ)
	{
		// Create local data
		maxInds[0] = nI;
		maxInds[1] = nJ;
		maxInds[2] = 1;
		nDim = 2;
		if(!(this->data==NULL))
		{
			delete[] this->data;
		}

		this->data = new int**[nI];
		for(int i=0; i<maxInds[0]; i++)
		{
			this->data[i] = new int*[nJ];
			for(int j=0; j<maxInds[1]; j++)
			{
				this->data[i][j] = new int[1];
				this->data[i][j][0] = data[i][j];	
			}
		}

		isDataAttached = true;

		// std::cout<<"nI: "<<nI<<std::endl;
		// std::cout<<"nJ: "<<nJ<<std::endl;
		// std::cout<<"nK: "<<1<<std::endl;
		// std::cout<<"minInds: "<<minInds[0]<<"\t"<<minInds[1]<<"\t"<<minInds[2]<<"\t"<<std::endl;
		// std::cout<<"maxInds: "<<maxInds[0]<<"\t"<<maxInds[1]<<"\t"<<maxInds[2]<<"\t"<<std::endl;

	}

	void setData(int*** data, int nI, int nJ, int nK)
	{
		// Create local data
		maxInds[0] = nI;
		maxInds[1] = nJ;
		maxInds[2] = nK;
		nDim = 3;
		if(!(this->data==NULL))
		{
			delete[] this->data;
		}

		this->data = new int**[nI];
		for(int i=0; i<maxInds[0]; i++)
		{
			this->data[i] = new int*[nJ];
			for(int j=0; j<maxInds[1]; j++)
			{
				this->data[i][j] = new int[nK];

				for(int k=0; k<maxInds[2]; k++)
				{	
					this->data[i][j][k] = data[i][j][k];	
				}
			}
		}

		isDataAttached = true;
	}

	bool setBoundaryType(int axis, int hilow, int type)
	{
		if( axis<0 || axis>2 )
		{
			return false;
		}

		if( hilow<0 || hilow>1 )
		{
			return false;
		}

		this->bcType[axis][hilow] = type;
	}

	int getPixelValue(int i, int j = 0, int k = 0)
	{
		return data[i][j][k];
	}

	int getNumberOfClusters()
	{
		return clusters.size();
	}

	std::list<FloodFill::cellIndices>* getClusterCellIndices(int clusterNumber)
	{
		if(clusterNumber > (clusters.size()-1))
		{
			return NULL;
		}
		return &clusters[clusterNumber];
	}

	/* data */

private:

	int nDim;
	int numClusters;
	int minInds[3];
	int maxInds[3];
	int bcType[3][2];


	// Attached dataset
	bool isDataAttached;
	int*** data;

	std::vector<std::list<FloodFill::cellIndices> > clusters;

	void init_clean()
	{
		nDim = 2;
		numClusters = 0;
		data = NULL;
		minInds[0] = 0;
		minInds[1] = 0;
		minInds[2] = 0;
		maxInds[0] = 0;
		maxInds[1] = 0;
		maxInds[2] = 0;
		isDataAttached = false;

		for(int i=0; i<3; i++)
		{
			for(int j=0; j<2; j++)
			{
				bcType[i][j] = FloodFill::FIXED;
			}
		}

	}

	void flood(int replacementValue, int i, int j, int k, std::list<FloodFill::cellIndices>& cellList);
	bool getNeighborIndex(int curr[3], int offset[3], int neigh[3]);




}; // End class

}; // End namespace



#endif