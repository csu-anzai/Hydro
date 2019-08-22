#include "HEVertex.h"
#include "HEEdge.h"
#include "HEFace.h"

void HEVertex::getUnassociatedEdges(std::vector<HEEdge*>& unassocEdges)
{
	std::vector<HEEdge* >::iterator eIt;
	for(eIt=this->edgesLeaving.begin(); eIt < this->edgesLeaving.end(); ++eIt)
	{
		if((*eIt)->left == NULL)
		{
			unassocEdges.push_back(*eIt);
		}
	}
}