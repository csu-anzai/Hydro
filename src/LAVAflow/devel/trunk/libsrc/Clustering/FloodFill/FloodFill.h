#ifndef LAVA_CLUSTERING_FLOODFILL_H
#define LAVA_CLUSTERING_FLOODFILL_H

#define MVAR_INDICES 10

#include <list>
#include <iostream>
#include <vector>
#include <deque>

// #include "../../includes/LAVAconstants.h"
#include "ArrayOperations.h"
#include "FlashAmrMesh.h"

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

	// this constructor will build the masking array by reading in the 
	// mesh specified in filename. currently only implemented for 3D
	FloodFill(std::string filename, const double thresholdValue, bool(*criterionFunc)(double****, double, int[], int, int, int), std::vector<std::string> meshVars)
	{
		this->init_clean();
		SetMaskCriterion(criterionFunc);
		data = BuildMaskArray(filename, thresholdValue, meshVars);
		isDataAttached = true;
	}

	
	// this destructor has not been tested for 2D
	~FloodFill()
	{
		if(!(this->data==NULL))
		{
			for (int i=0; i<maxInds[IAXIS]; i++)
	        {
	            if (nDim == 3)
		            for (int j=0; j<maxInds[JAXIS]; j++)
		                delete[] data[i][j];

	            delete[] data[i];
	        }
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

		this->data = data;

		// this->data = new int**[nI];
		// for(int i=0; i<maxInds[0]; i++)
		// {
		// 	this->data[i] = new int*[nJ];
		// 	for(int j=0; j<maxInds[1]; j++)
		// 	{
		// 		this->data[i][j] = new int[nK];

		// 		for(int k=0; k<maxInds[2]; k++)
		// 		{	
		// 			this->data[i][j][k] = data[i][j][k];	
		// 		}
		// 	}
		// }

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

	void SetMaskCriterion(bool(*criterionFunc)(double****, double, int[], int, int, int))
	{
		this->criterionFunc = criterionFunc;
	}

	int getPixelValue(int i, int j = 0, int k = 0)
	{
		return data[i][j][k];
	}

	int getNumberOfClusters()
	{
		return clusters.size();
	}

	int getClusterSize(int index)
	{
		if (index < clusters.size())
			return clusters[index].size();
		else
		{
			std::cerr << "Requested index exceeds total number of clusters." << std::endl;
			return -1;
		}
	}

	std::deque<FloodFill::cellIndices>* getClusterCellIndices(int clusterNumber)
	{
		if(clusterNumber > (clusters.size()-1))
		{
			return NULL;
		}
		return &clusters[clusterNumber];
	}

	int* GetMaxIndices()
	{
		return maxInds;
	}

	int GetNDim()
	{
		return nDim;
	}

	void SetMinClusterSize(int minClusterSize)
	{
		this->minClusterSize = minClusterSize;
	}

	/* data */
	std::deque<std::deque<FloodFill::cellIndices> > clusters;

private:

	int nDim;
	int numClusters;
	int minInds[3];
	int maxInds[3];
	int bcType[3][2];
	int minClusterSize;


	// Attached dataset
	bool isDataAttached;
	int*** data;

	

	bool (*criterionFunc)(double****, double, int[], int, int, int);

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
		minClusterSize = -1;

		for(int i=0; i<3; i++)
		{
			for(int j=0; j<2; j++)
			{
				bcType[i][j] = FloodFill::FIXED;
			}
		}

	}

	int*** BuildMaskArray(std::string , const double, std::vector<std::string>);

	void flood(int replacementValue, int i, int j, int k, std::deque<FloodFill::cellIndices>& cellList);
	bool getNeighborIndex(int curr[3], int offset[3], int neigh[3]);

	

}; // End class

}; // End namespace



#endif
