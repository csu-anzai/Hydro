#ifndef LAVA_CLUSTERING_BLOCK_H
#define LAVA_CLUSTERING_BLOCK_H


#include <list>
#include <ostream>

#include "../includes/LAVAconstants.h"


namespace Clustering
{

class Block
{

private:
	
	// Pointers for the parent block and possible child blocks
	Block* parent;
	std::list<Block*> children;
	bool m_isLeaf;
	bool m_isRoot;
	bool m_isShrunk;

	// Index bounds in the image
	// Rows = dimension
	// Cols = low/high
	int bounds[3][2];

	// Identifying number for writing/reading
	int ID;

public:

	Block()
	{
		m_isRoot = true;
		m_isLeaf = true;
		m_isShrunk = false;
		this->setParent(NULL);
		this->setBounds(IAXIS,0,0);
		this->setBounds(JAXIS,0,0);
		this->setBounds(KAXIS,0,0);
	}

	Block(Block* par)
	{
		m_isRoot = false;
		m_isLeaf = true;
		m_isShrunk = false;
		this->setParent(par);
		this->setBounds(IAXIS,0,0);
		this->setBounds(JAXIS,0,0);
		this->setBounds(KAXIS,0,0);
	}

	Block(Block* par, int bounds[3][2])
	{
		m_isRoot = false;
		m_isLeaf = true;
		m_isShrunk = false;
		this->setParent(par);
		this->setBounds(IAXIS, bounds[IAXIS][LOW], bounds[IAXIS][HIGH]);
		this->setBounds(JAXIS, bounds[JAXIS][LOW], bounds[JAXIS][HIGH]);
		this->setBounds(KAXIS, bounds[KAXIS][LOW], bounds[KAXIS][HIGH]);
	}


	void setParent(Block* par) 
	{ 
		this->parent = par; 
		if(par!=NULL)
		{
			m_isRoot = false;
		}
		else
		{
			m_isRoot = true;
		}
	}

	Block* getParent() { return this->parent; }


	bool isLeaf(){ return m_isLeaf; }
	bool isRoot(){ return m_isRoot; }
	bool isShrunk(){ return m_isShrunk; }
	void setShrunk(bool s = true) { this->m_isShrunk = s; }

	void setBounds(int bounds[3][2]);
	void setBounds(int dir, int dirBounds[2]);
	void setBounds(int dir, int dirMin, int dirMax);
	void setBounds(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax);

	void getBounds(int bounds[3][2]);
	void getBounds(int dir, int* dirBounds);
	void getBounds(int dir, int& dirMin, int& dirMax);

	void setID(int id){ this->ID = id; }
	void getID(int& id){ id = this->ID; }
	int  getID(){ return this->ID; }

	int getParentID()
	{
		if(m_isRoot)
		{
			return -1;
		}
		else
		{
			return this->getParent()->getID();
		}
	}


	int getChildIDs(int childIDs[2])
	{
		if(m_isLeaf)
		{
			childIDs[0] = -1;
			childIDs[1] = -1;
		}
		else
		{
			int count = 0;	
			for(std::list<Block*>::iterator c = this->children.begin(); c != this->children.end(); ++c)
			{
				childIDs[count] = (*c)->getID();
				count++;
			}
		}
	}


	int getNumberOfCells(int dir);
	void getNumberOfCells(int numCells[3]);
	void getSplitBounds(int splitDir, int splitIndex, 
						int boundsLeft[3][2], int boundsRight[3][2]);

	bool addChild(Block* child);



	void print(std::ostream& out, int indent=0);

};

}; // end namespace

#endif