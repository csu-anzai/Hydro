#include "Block.h"
#include <ostream>
#include <iostream>
#include <string>

namespace Clustering
{


void Block::setBounds(int bounds[3][2])
{
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<2; j++)
		{
			this->bounds[i][j] = bounds[i][j];
		}
	}
}


void Block::setBounds(int dir, int dirBounds[2])
{
	this->bounds[dir][0] = dirBounds[0];
	this->bounds[dir][1] = dirBounds[1];
}



void Block::setBounds(int dir, int dirMin, int dirMax)
{
	this->bounds[dir][0] = dirMin;
	this->bounds[dir][1] = dirMax;
}

void Block::setBounds(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax)
{
	this->setBounds(IAXIS,iMin,iMax);
	this->setBounds(JAXIS,jMin,jMax);
	this->setBounds(KAXIS,kMin,kMax);
}


void Block::getBounds(int dir, int& dirMin, int& dirMax)
{
	dirMin = this->bounds[dir][LOW];
	dirMax = this->bounds[dir][HIGH];
}


void Block::getBounds(int bounds[3][2])
{
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<2; j++)
		{
			bounds[i][j] = this->bounds[i][j];
		}
	}
}


int Block::getNumberOfCells(int dir)
{
	return bounds[dir][HIGH] - bounds[dir][LOW] + 1;
}
void Block::getNumberOfCells(int numCells[3])
{
	numCells[IAXIS] = this->getNumberOfCells(IAXIS);
	numCells[JAXIS] = this->getNumberOfCells(JAXIS);
	numCells[KAXIS] = this->getNumberOfCells(KAXIS);
}

void Block::getSplitBounds(int splitDir, int splitIndex,
							int boundsLeft[3][2], int boundsRight[3][2])
{
	// Initialize new bounds
	this->getBounds(boundsLeft);
	this->getBounds(boundsRight);

	// Set the proper bound information for the new blocks in the direction we're examining
	boundsLeft[splitDir][HIGH] = boundsLeft[splitDir][LOW] + splitIndex;
	boundsRight[splitDir][LOW] = boundsLeft[splitDir][HIGH] + 1;
}



bool Block::addChild(Block* child)
{
	this->children.push_back(child);
	this->m_isLeaf = false;
	return true;
}


void Block::print(std::ostream& out, int indent)
{
	bool isRoot = this->isRoot();
	bool isLeaf = this->isLeaf();
	int parentID = -1;


	if(!isRoot)
	{
		parentID = this->parent->getID();
	}

	out<<std::string(indent,'\t')<<"Block "<<ID<<std::endl;
	out<<std::string(indent+1,'\t')<<std::boolalpha<<"is Root: "<<isRoot<<std::endl;
	out<<std::string(indent+1,'\t')<<std::boolalpha<<"is Leaf: "<<isLeaf<<std::endl;
	

	out<<std::string(indent+1,'\t')<<"Parent ID: "<<parentID<<std::endl;

	out<<std::string(indent+1,'\t')<<"Child IDs: ";
	if(!isLeaf)
	{
		for(std::list<Block*>::iterator c = this->children.begin(); c != this->children.end(); ++c)
		{
			out<<"\t"<<(*c)->getID();
		}
	}
	else
	{
		out<<"\t"<<-1;
	}
	out<<std::endl;


}


}; // end namespace