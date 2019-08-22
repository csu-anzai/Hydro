#ifndef LAVA_CLUSTERING_SIGBERG_H
#define LAVA_CLUSTERING_SIGBERG_H

#include <list>
#include <ostream>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkDataArray.h>

#include "../../includes/LAVAconstants.h"
#include "../Block.h"
#include "../../Graph/Graph.h"


namespace Clustering
{



class SignatureBerger
{
public:
	SignatureBerger()
	{
		this->init_clean();
	}

	SignatureBerger(int nDim, int nCells[3])
	{
		this->init_clean();

		this->nDim = nDim;
		this->nCells[0] = nCells[0];
		this->nCells[1] = nCells[1];
		this->nCells[2] = nCells[2];

		std::cout<<"[SignatureBerger::constructor] nCells: "<<this->nCells[0]<<"\t"<<this->nCells[1]<<"\t"<<this->nCells[2]<<std::endl;
	}

	~SignatureBerger(){;}

	void reset();
	void process();
	void print(std::ostream& out, int indent = 0);
	void printBlocks(std::ostream& out);
	void printConnectivity(std::ostream& out);



	// Accessors
	void setMinCellsPerDim(int minI, int minJ, int minK)
	{
		minCellsPerDir[0] = minI;
		minCellsPerDir[1] = minJ;
		minCellsPerDir[2] = minK;
	}	

	void setEfficiency(double effthresh)
	{
		this->efficiencyThreshold = effthresh;
	}

	void setFinalPass(bool tf)
	{
		this->doFinalPass = tf;
	}
	bool getFinalPass()
	{
		return this->doFinalPass;
	}

	void setData(vtkSmartPointer<vtkDataArray> data)
	{
		this->data = data;
		isDataAttached = true;
	}

	std::list<Block*>* getBlockList() { return &allBlocks; }

	Graph::Graph<Block*,double>* getConnectivityGraph(){ return &(this->blockConnectivity); }

	int getNumberOfLeafBlocks()
	{
		int nLeafs = 0;
		std::list<Block*>::iterator it;
		for(it = this->allBlocks.begin(); it != this->allBlocks.end(); ++it)
		{
			 if((*it)->isLeaf())
			 	nLeafs++;
		}

		return nLeafs;
	}

	void computeClusters()
	{
		this->blockConnectivity.computeSubgraphRoots();
	}

	int getNumberOfClusters()
	{
		return this->blockConnectivity.getNumberOfSubgraphs();
	}



	/* data */

private:

	int nDim;
	int numIter;
	int numBlocks;
	int minCellsPerDir[3];
	int nCells[3];
	double efficiencyThreshold;
	double efficiencyAlpha;
	bool ignoreMinimumSize;
	bool doFinalPass;

	// Attached dataset
	bool isDataAttached;
	vtkSmartPointer<vtkDataArray> data;

	Graph::Graph<Block*,double> blockConnectivity;
	std::list<Block*> clusterRoots;
	std::list<Block*> allBlocks;


	void shrinkBlock(Block* blk);
	bool scanPossibleSplits(Block* blk, const std::vector<int>& indices, int examinedDir, int& splitDir, int& splitInd);
	double computeBlockEfficiency(Block* blk);
	void computeSignatures(Block* blk, int* sigmaI, int* sigmaJ, int* sigmaK, int* laplaI, int* laplaJ, int* laplaK);
	void refineBlock(Block* blk);
	bool checkBlockSanity(Block* blk, double efficiencyOld);
	bool doBlocksTouch(Block* blk1, Block* blk2, int tolerance = 1);

	bool registerBlock(Block* blk)
	{
		this->allBlocks.push_back(blk);

		// Set the ID of the block to be the number of blocks currently 
		// registered
		blk->setID(numBlocks);

		// increase the number of registered blocks
		numBlocks++;

		return true;
	}






	void init_clean()
	{
		numBlocks = 0;
		numIter = 0;
		nDim = 3;
		minCellsPerDir[0] = 3;
		minCellsPerDir[1] = 3;
		minCellsPerDir[2] = 3;
		efficiencyThreshold = 0.75;
		ignoreMinimumSize = false;
		isDataAttached = false;
		efficiencyAlpha = 0.9;
		doFinalPass = false;
	}







}; // End class

}; // End namespace



#endif