#include "HEEdge.h"
#include "HEVertex.h"

bool HEEdge::deleteEdge()
{
	if(this==NULL)
	{
		return true;
	}
	// if(this->pair->left!=NULL || this->left!=NULL)
	// {
	// 	return false;
	// }
	// if(this->vertex != NULL)
	// {
		this->vertex->removeLeavingEdge(this);
	// }

	return true;
}