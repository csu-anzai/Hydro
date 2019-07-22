#include "HEFace.h"
#include "HEEdge.h"
#include "HEVertex.h"
#include "../Primitives/Vec3.h"
#include "../Primitives/Edge3.h"
#include "../../Math/LAVAMath.h"
#include <vector>
#include <ostream>

void HEFace::getAllEdges(std::vector<HEEdge*>& edgeList)
{
	HEEdge* e = this->edge;
	do
	{	
		edgeList.push_back(e);
		e = e->next;
	}while(e!=this->edge);
}


void HEFace::print(std::ostream& out, int indent)
{
	std::vector<HEEdge*> edgeList;
	this->getAllEdges(edgeList);

	for(std::vector<HEEdge*>::iterator i = edgeList.begin(); i != edgeList.end(); ++i)
	{
		(*i)->vertex->print(out);	
	}

}

bool HEFace::deleteFace()
{
	// Loop over all edges and set their face pointer to NULL
	if(this == NULL)
	{
		return true;
	}
	HEEdge* e = this->edge;
	HEEdge* eTmp;
	do
	{	
		eTmp = e;
		e = e->next;
		eTmp->reset();
	}while(e!=this->edge);
	this->edge = NULL;

	return true;
}

void HEFace::computeNormal()
{
	// NOTE:
	// THIS WILL FAIL WHEN THE TWO EDGES ARE PARALLEL
	// WATCH OUT FOR THIS DISASTER LATER
	//	LOOK FOR NEWELL'S ALGORITHM

	// p1->p2
	Edge3 A(this->edge->vertex->pos, this->edge->pair->vertex->pos);
	// p1->p3
	Edge3 B(this->edge->vertex->pos, this->edge->next->pair->vertex->pos);

	this->normal = Math::crossProduct(A.getNormalizedRay(),B.getNormalizedRay());



}

Vec3 HEFace::getNormal()
{
	return this->normal;
}