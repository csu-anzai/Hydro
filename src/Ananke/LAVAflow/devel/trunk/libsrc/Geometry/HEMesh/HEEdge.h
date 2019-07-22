#ifndef LAVA_HEEDGE_H
#define LAVA_HEEDGE_H

// #include "HEVertex.h"
// #include "HEFace.h"
#include <cstdlib>
#include <ostream>

class HEVertex;
class HEFace;

class HEEdge
{
public:
	HEEdge()
	{
		this->vertex = NULL;
		this->next = NULL;
		this->prev = NULL;
		this->left = NULL;
		this->pair = NULL;
		isOrphan = true;
	}
	

	bool isOrphan;
	HEVertex *vertex; // Vertex the edge STARTS at (outgoing)
	HEEdge *next, *prev; // next and previous edges (ordered CCW)
	HEEdge *pair; // The corresponding edge lying on the OTHER face the total edge bounds
	HEFace *left; // Face this Edge belongs to

	void print(std::ostream& out, int indent=0)
	{
		if(this==NULL)
			return;
		
		out<<std::string(indent,'\t')<<"HEEdge information:"<<std::endl;
		out<<std::string(indent+1,'\t')<<"Edge Pointer: "<<this<<std::endl;
		out<<std::string(indent+1,'\t')<<"Vertex Pointer: "<<this->vertex<<std::endl;
		out<<std::endl;
	}

	void reset()
	{
		this->left = NULL;
		this->prev = NULL;
		this->next = NULL;
		if(this->pair->left==NULL)
		{
			isOrphan = true;
		}
	}

	bool deleteEdge();
};

#endif