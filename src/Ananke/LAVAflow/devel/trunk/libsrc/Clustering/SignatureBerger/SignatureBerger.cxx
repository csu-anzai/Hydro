#if 0
#include "../../Utilities/LAVAUtil.h"
#include "../../Graph/Graph.h"
#include "SignatureBerger.h"

#include <iostream>
#include <vector>
#include <ostream>
#include <sstream>

namespace Clustering
{

void SignatureBerger::process()
{

	std::clog<<"[SignatureBerger::process] Entering function"<<std::endl;

	// Check that we have data to work with
	if(!this->isDataAttached)
	{
		std::cerr<<"[SignatureBerger::process] ERROR: No data attached!"<<std::endl;
		return;
	}

	// Create initial root block
	int indExtents[3][2];
	indExtents[IAXIS][LOW] = 0;
	indExtents[JAXIS][LOW] = 0;
	indExtents[KAXIS][LOW] = 0;
	indExtents[IAXIS][HIGH] = this->nCells[IAXIS]-1;
	indExtents[JAXIS][HIGH] = this->nCells[JAXIS]-1;
	indExtents[KAXIS][HIGH] = this->nCells[KAXIS]-1;

	if(this->nDim<2)
	{
		indExtents[JAXIS][HIGH] = 0;
	}
	if(this->nDim<3)
	{
		indExtents[KAXIS][HIGH] = 0;
	}

	Clustering::Block* blk = new Clustering::Block(NULL, indExtents);
	this->shrinkBlock(blk);

	// Register this block
	this->registerBlock(blk);

	// Create the original node in the connectivity graph
	this->blockConnectivity.addNode(blk);

	// Begin the depth-first refinement
	std::cout<<"[SignatureBerger::process] Proceeding to refine blocks"<<std::endl;
	this->refineBlock(blk);


	std::clog<<"[SignatureBerger::process] Leaving function"<<std::endl;

}

/*=================================================

	Refine's a block based on directional
	signatures. This is the primary overarching
	function in the Berger algorithm.

	refineBlock takes in a block, 
	computes its signatures, finds possible splits,
	and creates two new Blocks if a suitable
	split is found.



=================================================*/
void SignatureBerger::refineBlock(Block* blk)
{

	// Shrink the block to its minimum size that
	// encompasses all flagged points inside of it
	this->shrinkBlock(blk);

	// Compute the efficiency of the block
	// If the efficiency is greater than the 
	// specified tolerance, leave refineBlocks
	// because there's nothing to be done
	double blkEff = this->computeBlockEfficiency(blk);

	if(blkEff > this->efficiencyThreshold)
	{
		return;
	}

	// Compute the signatures for the block
	// This includes the histograms and Laplacian



	/*=================================================

		Find possible split points.
		(1)	Check for zeros in the histogram
			in each direction
		(2)	Check for inflection points in each direction

		In the event that a check fails, attempt the next check.
		In theory, it is possible to attempt all checks and never
		split.

	=================================================*/
	// Information about a found split
	// Note that splitIndex and possibleSplits are in terms
	// of the Block indices (0->nCells-1)
	// and should be offset appropriately when actually
	// splitting the block to get them in global index space
	bool foundUsablePoint = false;
	int splitDir = -1;
	int splitIndex = -1;
	int examinedDir = -1;
	std::vector<int> possibleSplits;

	// Allocate room for the sum and laplacian arrays
	// The Lapacian arrays only need to be nCells-2 in size
	// but to make this conceptually easier I've let them be the
	// whole size. The 
	int nCellsI = blk->getNumberOfCells(IAXIS);
	int nCellsJ = blk->getNumberOfCells(JAXIS);
	int nCellsK = blk->getNumberOfCells(KAXIS);

	// std::cout<<"[refineBlock] nCells: "<<nCellsI<<"\t"<<nCellsJ<<"\t"<<nCellsK<<std::endl;

	int* sigmaI = new int[nCellsI];
	int* sigmaJ = new int[nCellsJ];
	int* sigmaK = new int[nCellsK];
	int* laplaI = new int[nCellsI];
	int* laplaJ = new int[nCellsJ];
	int* laplaK = new int[nCellsK];

	// Compute signatures
	this->computeSignatures(	blk,
							sigmaI, sigmaJ, sigmaK,
							laplaI, laplaJ, laplaK);

	// std::cout<<"sigmaI: "<<std::endl;
	// for(int i=0; i<nCellsI; i++)
	// 	std::cout<<"\t"<<sigmaI[i];
	// std::cout<<std::endl;


	// std::cout<<"sigmaJ: "<<std::endl;
	// for(int i=0; i<nCellsJ; i++)
	// 	std::cout<<"\t"<<sigmaJ[i];
	// std::cout<<std::endl;


	// std::cout<<"sigmaK: "<<std::endl;
	// for(int i=0; i<nCellsK; i++)
	// 	std::cout<<"\t"<<sigmaK[i];
	// std::cout<<std::endl;


	// Compute when the Laplacian arrays change sign
	bool* signChangedI = new bool[nCellsI];
	bool* signChangedJ = new bool[nCellsJ];
	bool* signChangedK = new bool[nCellsK];

	for(int i=1; i<nCellsI-1; i++)
	{
		signChangedI[i] = false;
		if(	(laplaI[i]>0 && laplaI[i+1]<0) ||
			(laplaI[i]<0 && laplaI[i+1]>0) )
		{
			signChangedI[i] = true;
		}
	}

	for(int i=1; i<nCellsJ-1; i++)
	{
		signChangedJ[i] = false;
		if(	(laplaJ[i]>0 && laplaJ[i+1]<0) ||
			(laplaJ[i]<0 && laplaJ[i+1]>0) )
		{
			signChangedJ[i] = true;
		}
	}

	for(int i=1; i<nCellsK-1; i++)
	{
		signChangedK[i] = false;
		if(	(laplaK[i]>0 && laplaK[i+1]<0) ||
			(laplaK[i]<0 && laplaK[i+1]>0) )
		{
			signChangedK[i] = true;
		}
	}


	// Check for zero-splits before attempting inflection splits

	// Check I-zeros
	possibleSplits.clear(); // Clear the vector
	for(int i=0; i<nCellsI; i++)
	{
		if(sigmaI[i]==0)
		{
			possibleSplits.push_back(i);
		}
	}

	if(!possibleSplits.empty() && !foundUsablePoint)
	{	
		// std::cout<<"Attemping a split using IAXIS zeros"<<std::endl;
		examinedDir = IAXIS;
		foundUsablePoint = this->scanPossibleSplits(blk,possibleSplits, examinedDir, splitDir, splitIndex);
	}

	// Check J-zeros
	if(this->nDim>1)
	{
		possibleSplits.clear(); // Clear the vector
		for(int i=0; i<nCellsJ; i++)
		{
			if(sigmaJ[i]==0)
			{
				possibleSplits.push_back(i);
			}
		}

		if(!possibleSplits.empty() && !foundUsablePoint)
		{
			examinedDir = JAXIS;
			foundUsablePoint = this->scanPossibleSplits(blk,possibleSplits, examinedDir, splitDir, splitIndex);
		}
	}

	// Check K-zeros
	if(this->nDim==3)
	{
		possibleSplits.clear(); // Clear the vector
		for(int i=0; i<nCellsK; i++)
		{
			if(sigmaK[i]==0)
			{
				possibleSplits.push_back(i);
			}
		}

		if(!possibleSplits.empty() && !foundUsablePoint)
		{
			examinedDir = KAXIS;
			foundUsablePoint = this->scanPossibleSplits(blk,possibleSplits, examinedDir, splitDir, splitIndex);
		}
	}

	// Check inflections in I
	possibleSplits.clear(); // Clear the vector
	for(int i=1; i<nCellsI-1; i++)
	{
		if(signChangedI[i]==0)
		{
			possibleSplits.push_back(i);
		}
	}

	if(!possibleSplits.empty() && !foundUsablePoint)
	{
		examinedDir = IAXIS;
		foundUsablePoint = this->scanPossibleSplits(blk,possibleSplits, examinedDir, splitDir, splitIndex);
	}


	// Check inflections in J
	if(this->nDim>1)
	{
		possibleSplits.clear(); // Clear the vector
		for(int i=1; i<nCellsJ-1; i++)
		{
			if(signChangedJ[i]==0)
			{
				possibleSplits.push_back(i);
			}
		}

		if(!possibleSplits.empty() && !foundUsablePoint)
		{
			examinedDir = JAXIS;
			foundUsablePoint = this->scanPossibleSplits(blk,possibleSplits, examinedDir, splitDir, splitIndex);
		}
	}


	// Check inflections in K
	if(this->nDim==3)
	{
		possibleSplits.clear(); // Clear the vector
		for(int i=1; i<nCellsK-1; i++)
		{
			if(signChangedK[i]==0)
			{
				possibleSplits.push_back(i);
			}
		}

		if(!possibleSplits.empty() && !foundUsablePoint)
		{
			examinedDir = KAXIS;
			foundUsablePoint = this->scanPossibleSplits(blk,possibleSplits, examinedDir, splitDir, splitIndex);
		}
	}

	// Free allocated memory
	// Failure to do this before recursing 
	// could lead to crashes, as the memory wont
	// be freed until the depth of the tree is traversed
	delete[] sigmaI;
	delete[] sigmaJ;
	delete[] sigmaK;
	delete[] laplaI;
	delete[] laplaJ;
	delete[] laplaK;
	delete[] signChangedI;
	delete[] signChangedJ;
	delete[] signChangedK;

	// Perform the split operation if we found a usable point
	// This should create two new Blocks whose parent is the current block
	// The appropriate child list in the current block should be updated
	//
	// If we don't find a usable point, do a final pass through ignoring the 
	// minimum block size
	if(foundUsablePoint)
	{
		// std::cout<<"Found a usable split!"<<std::endl;
		// std::cout<<"\tDirection: "<<splitDir<<std::endl;
		// std::cout<<"\t     Index: "<<splitIndex<<std::endl;

		// Create the two new blocks
		int newBoundsLeft[3][2], newBoundsRight[3][2];
		blk->getSplitBounds(splitDir, splitIndex, newBoundsLeft, newBoundsRight);

		Block* blkLeft = new Block(blk, newBoundsLeft);
		Block* blkRight = new Block(blk, newBoundsRight);

		// Add these children to the current block
		blk->addChild(blkLeft);
		blk->addChild(blkRight);

		// Add the children to the list of all blocks
		this->registerBlock(blkLeft);
		this->registerBlock(blkRight);

		// Shrink the blocks so that the graph connectivity will be updated accurately
		this->shrinkBlock(blkLeft);
		this->shrinkBlock(blkRight);
		
		/*=================================================

			Update connectivity graph information.

			(1)	Add blkLeft and blkRight as new nodes in the graph
			(2)	Find the node corresponding to the original Block
			(3)	Join the new left and right nodes together with an edge
			(4)	Add edges between the original neighbors and the new nodes
				if the blocks they represent touch
			(5) Delete the original node

		=================================================*/
		// Get all relevant nodes
		Graph::Node<Block*,double>* nodeLeft	= this->blockConnectivity.addNode(blkLeft);
		Graph::Node<Block*,double>* nodeRight= this->blockConnectivity.addNode(blkRight);
		Graph::Node<Block*,double>* nodeOrig 	= this->blockConnectivity.findNode(blk);

		// Connect the new nodes together
		this->blockConnectivity.addEdge(nodeLeft, nodeRight,0.0);

		// Connect nodes to the original neighbor if they touch
		std::list<Graph::Edge<Block*,double>*>::iterator origEdgeIt;
		for(origEdgeIt = nodeOrig->begin(); origEdgeIt != nodeOrig->end(); ++origEdgeIt)
		{
			Graph::Node<Block*,double>* nodeOrigNeighbor = (*origEdgeIt)->getOtherNode(nodeOrig);

			// Determine if the new node and original node touch
			// If so, add an edge between them
			if(this->doBlocksTouch(nodeLeft->getValue(), nodeOrigNeighbor->getValue()))
			{
				this->blockConnectivity.addEdge(nodeLeft, nodeOrigNeighbor, 0.0);
			}

			if(this->doBlocksTouch(nodeRight->getValue(), nodeOrigNeighbor->getValue()))
			{
				this->blockConnectivity.addEdge(nodeRight, nodeOrigNeighbor, 0.0);
			}

		}

		// Delete original node
		this->blockConnectivity.removeNode(nodeOrig);


		//Call refineBlock on the new blocks
		this->refineBlock(blkLeft);
		this->refineBlock(blkRight);
	}
	else
	{
		// Do a final pass, ignoring the minimum block size
		if(this->getFinalPass() && !ignoreMinimumSize)
		{
			// std::cout<<"[SignatureBerger::refineBlock] Doing a final pass with no minimum size"<<std::endl;
			
			this->ignoreMinimumSize = true;
			this->refineBlock(blk);
			this->ignoreMinimumSize = false;
		}
	}

	

}



/*=================================================

	Determine if two blocks touch

	Taken from:
	http://stackoverflow.com/questions/6664281/detect-if-two-rectangles-can-be-combined-into-a-single-rectangle

     *	tolerance = 0: require at least one pixel overlap
     *	tolerance = 1: accepts 'flush' adjacent neighbours
     *	Negative tolerance require more overlap
     *	tolerance > 1 allows gaps between rects

=================================================*/

bool SignatureBerger::doBlocksTouch(Block* A, Block* B, int tolerance)
{

	int bndsA[3][2], bndsB[3][2];
	A->getBounds(bndsA);
	B->getBounds(bndsB);

	bool areAdjacent = false;

	for(int i=0; i<nDim; i++)
	{
		areAdjacent = areAdjacent || 
					 (bndsA[i][LOW]-bndsB[i][HIGH]) > tolerance ||
					 (bndsB[i][LOW]-bndsA[i][HIGH]) > tolerance;
	}

	return !areAdjacent;
}


void SignatureBerger::shrinkBlock(Block* blk)
{
	int iMin, iMax;
	int jMin, jMax;
	int kMin, kMax;

	// Determine if we need to shrink the block
	// If not, leave immediately
	if(blk->isShrunk())
		return;

	// Set the min/max to be the opposite extents of the block
	blk->getBounds(IAXIS, iMax, iMin);
	blk->getBounds(JAXIS, jMax, jMin);
	blk->getBounds(KAXIS, kMax, kMin);

	// Get the actual bounds of the block
	int bnds[3][2];
	blk->getBounds(bnds);

	// std::cout<<"[SignatureBerger::shrinkBlock] Block bounds: "<<std::endl;
	// std::cout<<"\t"<<bnds[IAXIS][LOW]<<"\t"<<bnds[IAXIS][HIGH]<<std::endl;
	// std::cout<<"\t"<<bnds[JAXIS][LOW]<<"\t"<<bnds[JAXIS][HIGH]<<std::endl;
	// std::cout<<"\t"<<bnds[KAXIS][LOW]<<"\t"<<bnds[KAXIS][HIGH]<<std::endl;

	bool anyFound = false;
	int numCellsVisited = 0;

	for(int k=bnds[KAXIS][LOW]; k<=bnds[KAXIS][HIGH]; k++)
	{
		for(int j=bnds[JAXIS][LOW]; j<=bnds[JAXIS][HIGH]; j++)
		{
			for(int i=bnds[IAXIS][LOW]; i<=bnds[IAXIS][HIGH]; i++)
			{

				// int index = i + this->nCells[IAXIS]*(j + this->nCells[JAXIS]*k);
				int index = Util::sub2ind(i,j,k,this->nCells[IAXIS],this->nCells[JAXIS]);
					
				if(data->GetTuple1(index) > 0.5)
				{

					anyFound = true;

					if(k<kMin)
						kMin = k;

					if(k>kMax)
						kMax = k;

					if(j<jMin)
						jMin = j;

					if(j>jMax)
						jMax = j;

					if(i<iMin)
						iMin = i;

					if(i>iMax)
						iMax = i;

				}

				numCellsVisited++;
	
			}
		}
	}


	// std::cout<<"[SignatureBerger::shrinkBlock] Number of cells visited: "<<numCellsVisited<<std::endl;

	if(anyFound)
	{
	// std::cout<<"[SignatureBerger::shrinkBlock] Old bounds: "<<std::endl;
	// std::cout<<"\t"<<bnds[IAXIS][LOW]<<"\t"<<bnds[IAXIS][HIGH]<<std::endl;
	// std::cout<<"\t"<<bnds[JAXIS][LOW]<<"\t"<<bnds[JAXIS][HIGH]<<std::endl;
	// std::cout<<"\t"<<bnds[KAXIS][LOW]<<"\t"<<bnds[KAXIS][HIGH]<<std::endl;
	// std::cout<<"[SignatureBerger::shrinkBlock] New bounds: "<<std::endl;
	// std::cout<<"\t"<<iMin<<"\t"<<iMax<<std::endl;
	// std::cout<<"\t"<<jMin<<"\t"<<jMax<<std::endl;
	// std::cout<<"\t"<<kMin<<"\t"<<kMax<<std::endl;
	// std::cout<<"[SignatureBerger::shrinkBlock] Any values in range: "<<anyFound<<std::endl;
	

		blk->setBounds(IAXIS, iMin, iMax);
		blk->setBounds(JAXIS, jMin, jMax);
		blk->setBounds(KAXIS, kMin, kMax);
	}

	// Set that the block was shrunk
	blk->setShrunk(true);
}


double SignatureBerger::computeBlockEfficiency(Block* blk)
{
	int numCells = 0;
	int numMarkedCells = 0;


	// Compute the number of flagged cells in the block
	int bnds[3][2];
	blk->getBounds(bnds);
	for(int k=bnds[KAXIS][LOW]; k<=bnds[KAXIS][HIGH]; k++)
	{
		for(int j=bnds[JAXIS][LOW]; j<=bnds[JAXIS][HIGH]; j++)
		{
			for(int i=bnds[IAXIS][LOW]; i<=bnds[IAXIS][HIGH]; i++)
			{
				int index = i + this->nCells[IAXIS]*(j + this->nCells[JAXIS]*k);

				numCells++;
				if(data->GetTuple1(index) > 0.5)
				{
					numMarkedCells++;
				}	
			}
		}
	}


	// Compute the efficiency as the number of flagged cells divided by the number of total cells
	return double(numMarkedCells)/double(numCells);
}


void SignatureBerger::computeSignatures(Block* blk, int* sigmaI, int* sigmaJ, int* sigmaK, int* laplaI, int* laplaJ, int* laplaK)
{

	// Get the number of cells in each direction of the block
	int nCells[3];
	blk->getNumberOfCells(nCells);
	int bounds[3][2];
	blk->getBounds(bounds);

	// Zero the elements in the arrays
	memset(sigmaI, 0, nCells[IAXIS]);
	memset(sigmaJ, 0, nCells[JAXIS]);
	memset(sigmaK, 0, nCells[KAXIS]);
	memset(laplaI, 0, nCells[IAXIS]);
	memset(laplaJ, 0, nCells[JAXIS]);
	memset(laplaK, 0, nCells[KAXIS]);

	// Compute the sums
	for(int i=0; i<nCells[IAXIS]; i++)
	{
		int iGlobal = i + bounds[IAXIS][LOW];

		for(int j=0; j<nCells[JAXIS]; j++)
		{

			int jGlobal = j + bounds[JAXIS][LOW];
			for(int k=0; k<nCells[KAXIS]; k++)
			{
				int kGlobal = k + bounds[KAXIS][LOW];
				int index = Util::sub2ind(iGlobal,jGlobal,kGlobal,this->nCells[IAXIS], this->nCells[JAXIS]);

				if(this->data->GetTuple1(index) > 0.5)
				{

					sigmaI[i]++;
					sigmaJ[j]++;
					sigmaK[k]++;

				}

			}
		}
	}

	// Compute the Laplacian of the sigma distributions
	if(nCells[IAXIS]>=3)
	{
		for(int ind=1; ind<nCells[IAXIS]-1; ind++)
		{
			laplaI[ind] = sigmaI[ind+1] - 2*sigmaI[ind] + sigmaI[ind-1];
		}
	}

	if(nCells[JAXIS]>=3)
	{
		for(int ind=1; ind<nCells[JAXIS]-1; ind++)
		{
			laplaJ[ind] = sigmaJ[ind+1] - 2*sigmaJ[ind] + sigmaJ[ind-1];
		}
	}

	if(nCells[KAXIS]>=3)
	{
		for(int ind=1; ind<nCells[KAXIS]-1; ind++)
		{
			laplaK[ind] = sigmaK[ind+1] - 2*sigmaK[ind] + sigmaK[ind-1];
		}
	}

}


bool SignatureBerger::scanPossibleSplits(Block* blk, const std::vector<int>& possibleSplits, int examinedDir, int& splitDir, int& splitInd)
{

	int nPossibleSplits = possibleSplits.size();

	// std::cout<<"Possible split locations: ";
	// for(int ind=0; ind<nPossibleSplits; ind++)
	// {
	// 	std::cout<<"\t"<<possibleSplits[ind];
	// }
	// std::cout<<std::endl;

	// NOTE: Interesting problem!
	// If you treat the distances as doubles, the algorithm appears to become non-deterministic
	// I presume this is because the distances get a weebit of error in them,
	// which throws the sort out of whack
	// The first time it happens, the new path is laid and one gets increasingly disparate
	// results as the minimum block size decreases

	// Compute the mean split point
	int meanSplit = 0;
	for(int ind=0; ind<nPossibleSplits; ind++)
		meanSplit += possibleSplits[ind];
	// meanSplit /= double(nPossibleSplits);
	meanSplit = floor(meanSplit/nPossibleSplits);

	// Determine the absolute value of the distance from each possible split to the mean
	std::vector<std::pair<int,int> > splitDists;
	for(int ind=0; ind<nPossibleSplits; ind++)
	{
		int dist = abs(possibleSplits[ind] - meanSplit);
		int origIndex = ind;
		splitDists.push_back(std::pair<int,int>(dist,origIndex));
	}

	// Sort splitDists
	Util::sortVectorOfPairs(splitDists);

	// For every possible split, split the block
	// and determine if the two resulting blocks satisfy the sanity constraints
	// If they do, set the appropriate outgoing variables and leave the function
	bool foundUsablePoint = false;

	int originalBounds[3][2];
	blk->getBounds(originalBounds);

	// std::cout<<"Original bounds in direction "<<examinedDir<<": "<<originalBounds[examinedDir][LOW]<<"\t"<<originalBounds[examinedDir][HIGH]<<std::endl;

	for(int ind=0; ind<nPossibleSplits; ind++)
	{
		// Get the attempted split position from the possibleSplits vector,
		// given the sorting
		int attemptedSplit = possibleSplits[splitDists[ind].second];

		// std::cout<<"Attempted split: "<<attemptedSplit<<std::endl;

		// Get the boundaries for the new blocks
		int newBoundsLeft[3][2], newBoundsRight[3][2];
		blk->getSplitBounds(examinedDir, attemptedSplit, newBoundsLeft, newBoundsRight);


	// std::cout<<"Left block bounds in direction "<<examinedDir<<": "<<newBoundsLeft[examinedDir][LOW]<<"\t"<<newBoundsLeft[examinedDir][HIGH]<<std::endl;
	// std::cout<<"Rght block bounds in direction "<<examinedDir<<": "<<newBoundsRight[examinedDir][LOW]<<"\t"<<newBoundsRight[examinedDir][HIGH]<<std::endl;


		// Create the new blocks
		Block blkLeft(NULL, newBoundsLeft);
		Block blkRight(NULL, newBoundsRight);


		// Shrink the new blocks
		this->shrinkBlock(&blkLeft);
		this->shrinkBlock(&blkRight);

		// Get the efficiency of the original block
		// TODO: Pass this to the function instead of computing it every time
		double efficiencyOld = this->computeBlockEfficiency(blk);

		// Check each block for sanity
		bool isGoodLeft	= this->checkBlockSanity(&blkLeft, this->efficiencyAlpha*efficiencyOld);
		bool isGoodRight	= this->checkBlockSanity(&blkRight, this->efficiencyAlpha*efficiencyOld);

		// std::cout<<"Left sanity: "<<isGoodLeft<<std::endl;
		// std::cout<<"Rght sanity: "<<isGoodRight<<std::endl;

		if(isGoodLeft && isGoodRight)
		{
			// std::cout<<"In here!"<<std::endl;
			foundUsablePoint = true;
			splitInd = attemptedSplit;
			splitDir = examinedDir;
			break;
		}
	}

	// Leave the function
	return foundUsablePoint;

}



bool SignatureBerger::checkBlockSanity(Block* blk, double efficiencyMin)
{
	// Set output variable
	bool isSane = true;

	// Get number of cells in each direction
	int nCells[3];
	blk->getNumberOfCells(nCells);

	// Determine if the number of cells is greater than or equal to the minimum number of cells
	// for every valid direction (1:nDim)
	// Failure to meet this criterion in any direction invalidates the block
	bool isLargerThanMinimum = true;
	if(!this->ignoreMinimumSize)
	{
		for(int dim = 0; dim < this->nDim; dim++)
		{
			isLargerThanMinimum = isLargerThanMinimum && ( nCells[dim] >= this->minCellsPerDir[dim] );
		}
	}

	// Determine if the efficiency of the block is greater than the threshold
	bool isEfficient;
	// isEfficient = this->computeBlockEfficiency(blk) >= this->efficiencyThreshold;
	isEfficient = this->computeBlockEfficiency(blk) >= efficiencyMin;

	isSane = isLargerThanMinimum && isEfficient;

	return isSane;
}

void SignatureBerger::print(std::ostream& out, int indent)
{

	out<<std::string(indent,'\t')<<"Berger clustering information"<<std::endl;
	out<<std::string(indent+1,'\t')<<"Number of blocks (total): "<<this->numBlocks<<std::endl;
	out<<std::string(indent+1,'\t')<<"Block information:"<<std::endl;

	for(std::list<Block*>::iterator b = this->allBlocks.begin(); b != this->allBlocks.end(); ++b)
	{
		(*b)->print(out,indent+2);
	}

}



/*=================================================

	Print block information to a file

	Blk#|iMin|iMax|jMin|jMax|kMin|kMax|Par#|Child1#|Child2#

=================================================*/
void SignatureBerger::printBlocks(std::ostream& out)
{
	std::string delim = "\t";
	std::stringstream header;
	header<<"BLK"<<delim;
	header<<"iMin"<<delim;
	header<<"iMax"<<delim;
	header<<"jMin"<<delim;
	header<<"jMax"<<delim;
	header<<"kMin"<<delim;
	header<<"kMax"<<delim;
	header<<"PAR"<<delim;
	header<<"CHILD1"<<delim;
	header<<"CHILD2"<<delim;

	// Write the header to the file
	out<<header.str()<<std::endl;

	// Loop over all blocks and write the required information
	std::list<Block*>::iterator blkIt;
	for(blkIt = this->allBlocks.begin(); blkIt != this->allBlocks.end(); ++blkIt)
	{
		int bnds[3][2];
		(*blkIt)->getBounds(bnds);

		std::stringstream line;

		// Write block ID
		line<<(*blkIt)->getID()<<delim;

		// Write bounds
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<2; j++)
			{
				line<<bnds[i][j]<<delim;
			}
		}

		// Get parent/child info
		int parID = -1;
		int childIDs[2];

		parID = (*blkIt)->getParentID();
		(*blkIt)->getChildIDs(childIDs);

		// Write parent/child info
		line<<parID<<delim;
		line<<childIDs[0]<<delim;
		line<<childIDs[1];

		// Write line to file
		out<<line.str()<<std::endl;

	}


}


/*=================================================

	Print block connectivity to file

	BlkID|# of neighbors|Neighbor IDs

=================================================*/
void SignatureBerger::printConnectivity(std::ostream& out)
{



	std::string delim = "\t";
	std::stringstream header;
	header<<"BLK"<<delim;
	header<<"NUMNEIGH"<<delim;
	header<<"NEIGHBORS"<<delim;

	// Write the header to the file
	out<<header.str()<<std::endl;

	// Loop over all Nodes
	std::list<Graph::Node<Block*,double>*>::iterator nodeIt;
	for(	nodeIt = this->blockConnectivity.nodeBegin(); 
		nodeIt != this->blockConnectivity.nodeEnd(); ++nodeIt)
	{
		std::stringstream line;
			

		int blkID = (*nodeIt)->getValue()->getID();
		int numNeighbors = (*nodeIt)->getNumberOfEdges();

		line<<blkID<<delim;
		line<<numNeighbors<<delim;

		// Loop over all connected Nodes
		std::list<Graph::Edge<Block*,double>*>::iterator edgeIt;

		for(	edgeIt = (*nodeIt)->begin(); 
			edgeIt != (*nodeIt)->end(); ++edgeIt)
		{
			Graph::Node<Block*,double>* neighbor;
			neighbor = (*edgeIt)->getOtherNode(*nodeIt);
			int neighborID = neighbor->getValue()->getID();

			line<<neighborID<<delim;
		}

		// write line to output
		out<<line.str()<<std::endl;

	}


}

}; // end namespace
#endif
