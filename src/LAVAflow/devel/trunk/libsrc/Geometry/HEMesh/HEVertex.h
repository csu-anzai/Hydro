#ifndef LAVA_HEVERTEX_H
#define LAVA_HEVERTEX_H


// #include "HEEdge.h"
#include "../Primitives/Vec3.h"
#include <vector>
#include <ostream>
#include <string>

class HEEdge;

class HEVertex
{
public:
	HEVertex()
	{
		isOrphan = false;

	}
	~HEVertex(){}


	/* data */
	bool isOrphan;
	std::vector<HEEdge*> edgesLeaving; // Edge that STARTS at this vertex (outgoing)
	Vec3 pos;

	void getUnassociatedEdges(std::vector<HEEdge*>& unassocEdges);
	
	void addLeavingEdge(HEEdge* e)
	{
		isOrphan = false;

		std::vector<HEEdge* >::iterator eIt;
		for(eIt=this->edgesLeaving.begin(); eIt < this->edgesLeaving.end(); ++eIt)
		{
			if(e == *eIt)
				return;
		}

		this->edgesLeaving.push_back(e);
	}

	bool removeLeavingEdge(HEEdge* e)
	{

		std::vector<HEEdge* >::iterator eIt;
		for(eIt=this->edgesLeaving.begin(); eIt < this->edgesLeaving.end(); ++eIt)
		{
			if(e == *eIt)
			{
				this->edgesLeaving.erase(eIt);
				break;
			}
		}

		if(this->edgesLeaving.empty())
		{
			this->isOrphan = true;
		}

		return this->isOrphan;

	}

	void print(std::ostream& out, int indent=0)
	{
		out<<std::string(indent,'\t')<<"HEVertex information:"<<std::endl;
		out<<std::string(indent+1,'\t')<<"Vertex Pointer: "<<this<<std::endl;
		out<<std::string(indent+1,'\t')<<"Number of edges leaving: "<<this->edgesLeaving.size()<<std::endl;
		
		out<<std::string(indent+1,'\t')<<"Vertex Position: "<<std::endl;
		out<<std::string(indent+2,'\t');
		this->pos.print(out);

		out<<std::string(indent+1,'\t')<<"Edge Pointers: "<<std::endl;
		for (int i = 0; i < this->edgesLeaving.size(); ++i)
		{
			out<<std::string(indent+2,'\t')<<this->edgesLeaving[i]<<std::endl;
		}
		out<<std::endl;

	}

	bool deleteVertex()
	{
		// No cleanup needs to happen
		return isOrphan;
	}

};

#endif