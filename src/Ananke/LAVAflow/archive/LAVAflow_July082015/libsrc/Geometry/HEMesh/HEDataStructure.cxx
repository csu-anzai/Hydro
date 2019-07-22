#include "HEDataStructure.h"
#include <iostream>
#include <vector>
#include <set>

#include "HEEdge.h"
#include "HEVertex.h"
#include "../Primitives/Vec3.h"
#include "../Geometry.h"
#include "../../includes/LAVAconstants.h"

#define DEBUG

HEDataStructure* HEDataStructure::generateBrick(double xbnds[2], double ybnds[2],
										   		 double zbnds[2])
{
	double bndBox[3][2];

	bndBox[0][0] = xbnds[0];
	bndBox[0][1] = xbnds[1];
	bndBox[1][0] = ybnds[0];
	bndBox[1][1] = ybnds[1];
	bndBox[2][0] = zbnds[0];
	bndBox[2][1] = zbnds[1];

	return HEDataStructure::generateBrick(bndBox);
}

HEDataStructure* HEDataStructure::generateBrick(double bnds[6])
{
	double bndBox[3][2];

	bndBox[0][0] = bnds[0];
	bndBox[0][1] = bnds[1];
	bndBox[1][0] = bnds[2];
	bndBox[1][1] = bnds[3];
	bndBox[2][0] = bnds[4];
	bndBox[2][1] = bnds[5];

	return HEDataStructure::generateBrick(bndBox);
}

HEDataStructure* HEDataStructure::generateBrick(double xmin, double xmax,
										   		 double ymin, double ymax,
										   		 double zmin, double zmax)
{
	double bndBox[3][2];

	bndBox[0][0] = xmin;
	bndBox[0][1] = xmax;
	bndBox[1][0] = ymin;
	bndBox[1][1] = ymax;
	bndBox[2][0] = zmin;
	bndBox[2][1] = zmax;


	return HEDataStructure::generateBrick(bndBox);
}

HEDataStructure* HEDataStructure::generateBrick(double bndBox[3][2])
{
	int nVerts = 8;
	int nFaces = 6;
	int bbToVerts[][3] =	{{0,0,0},
						 {1,0,0},
				 		 {1,1,0},
				 		 {0,1,0},
				 		 {0,0,1},
				 		 {1,0,1},
				 		 {1,1,1},
				 		 {0,1,1}};

	int faceToVerts[][4] = {{3,2,1,0},
     					  {4,5,6,7},
     					  {0,1,5,4},
     					  {2,3,7,6},
     					  {0,4,7,3},
     					  {1,2,6,5}};

     HEDataStructure* he = new HEDataStructure;

	// create vertices
	for (int i = 0; i < nVerts; ++i)
	{
		Vec3 vertex;
		for (int j = 0; j < 3; ++j)
		{
			// cout<<bbToVerts[i][j]<<endl;
			vertex[j] = bndBox[j][bbToVerts[i][j]];
		}

		he->addVertex(&vertex);
		// he.vertices.back()->print(cout);
	}

	for (int i = 0; i < nFaces; ++i)
	{
		std::vector<int> vertInds(faceToVerts[i],faceToVerts[i]+4);
		he->addFace(vertInds);
	}

	// A brick should be embedded in 3d
	he->embeddingDimension = 3;

	return he;
}

HEDataStructure* HEDataStructure::generateRectangle(double bnds[6])
{
	double bndBox[3][2];

	bndBox[0][0] = bnds[0];
	bndBox[0][1] = bnds[1];
	bndBox[1][0] = bnds[2];
	bndBox[1][1] = bnds[3];
	bndBox[2][0] = bnds[4];
	bndBox[2][1] = bnds[5];

	return HEDataStructure::generateRectangle(bndBox);
}

HEDataStructure* HEDataStructure::generateRectangle(double bndBox[3][2])
{
	int nVerts = 4;
	int nFaces = 1;
	int bbToVerts[][3] =	{{0,0,0},
						 {1,0,0},
				 		 {1,1,0},
				 		 {0,1,0}};

	int faceToVerts[][4] = {{0,1,2,3}};

     HEDataStructure* he = new HEDataStructure;

	// create vertices
	for (int i = 0; i < nVerts; ++i)
	{
		Vec3 vertex;
		for (int j = 0; j < 3; ++j)
		{
			// cout<<bbToVerts[i][j]<<endl;
			vertex[j] = bndBox[j][bbToVerts[i][j]];
		}

		he->addVertex(&vertex);
		// he.vertices.back()->print(cout);
	}

	for (int i = 0; i < nFaces; ++i)
	{
		std::vector<int> vertInds(faceToVerts[i],faceToVerts[i]+4);
		he->addFace(vertInds);
	}

	// A rectangle should be embedded in 2d
	he->embeddingDimension = 2;

	return he;
}






int HEDataStructure::eulerCharacteristic()
{
	int nVertices = this->vertices.size();
	int nEdges = this->edges.size()/2;
	int nFaces = this->faces.size();

	return nVertices - nEdges + nFaces;
}

bool HEDataStructure::isClosed()
{
	// a surface will be closed if every edge is associated with a face
	bool isClosed = true;

	for (std::vector<HEEdge*>::iterator i = this->edges.begin(); i != this->edges.end(); ++i)
	{
		if((*i)->left == NULL)
		{
			isClosed = false;
			break;
		}
	}

	return isClosed;
}

void HEDataStructure::computeNormals()
{
	for (std::vector<HEFace*>::iterator i = this->faces.begin(); i != this->faces.end(); ++i)
	{
		(*i)->computeNormal();
	}
}

bool HEDataStructure::clipWithPlane(Vec3 planePoint, Vec3 planeNormal)
{
	// Get vertices infront of plane
	std::vector<HEVertex*> vertsInFront;
	for(int i=0; i<this->vertices.size(); ++i)
	{
		if(Geometry::isPointInfrontOfPlane(this->vertices[i]->pos,planePoint, planeNormal) )
		{
			vertsInFront.push_back(this->vertices[i]);
		}
	}

	if(vertsInFront.empty())
	{
		std::cerr<<"[HEDataStructure::clipWithPlane] No vertices infront of plane!"<<std::endl;
	}

	// Split all the faces with a plane
	std::vector<HEVertex*> splitVerts;
	this->splitFacesWithPlane(&planePoint, &planeNormal, splitVerts);

#ifdef DEBUG
	for(int i=0; i<splitVerts.size(); i++)
	{
		std::cout<<"-------------------------"<<std::endl;
		std::cout<<"     Split vertex "<<i<<std::endl;
		std::cout<<"-------------------------"<<std::endl;
		splitVerts[i]->print(std::cout,1);

	}
#endif

	this->writeMesh(std::cout);

	// Check to see if we have any split vertices
	// We could have none in the case that the split plane
	// is coplanar with the boundary of the polytope
	// If we don't have any, go ahead and return
	// skipping the deletion and patching part
	if(splitVerts.empty())
	{
		return true;
	}
	else
	{

		for(int i=0; i<vertsInFront.size(); ++i)
		{
			vertsInFront[i]->print(std::cout,1);
			this->deleteConnectionsToVertex(vertsInFront[i]);
		}

		this->writeMesh(std::cout);


		for(int i=0; i<this->edges.size(); ++i)
		{
			// !! HEVertex *vertex; // Vertex the edge STARTS at (outgoing)
			// !! HEEdge *next, *prev; // next and previous edges (ordered CCW)
			// !! HEEdge *pair; // The corresponding edge lying on the OTHER face the total edge bounds
			// !! HEFace *left; // Face this Edge belongs to
			HEEdge* e = this->edges[i];
			if(e->left == NULL)
			{
				std::cout<<"Edge "<<i<<" has no face"<<std::endl;
				e->deleteEdge();
			}
		}

		// Patch the hole
		if(this->embeddingDimension>2)
		{
			this->patchHole(splitVerts[0]);
		}
	}

	return true;

}


bool HEDataStructure::patchHole(HEVertex* vStart)
{

	std::vector<HEVertex*> vertList;
	HEVertex* vCurrent, *vNext;
	HEEdge* E;

	// vertList.push_back(vStart);
	vCurrent = vStart;

	do
	{
		std::vector<HEEdge*> unassocEdges;
		vCurrent->getUnassociatedEdges(unassocEdges);

		// If we have more or less than 1 unassociated edge, this is a
		// poorly defined problem, so error
		if(unassocEdges.size() != 1)
		{
			std::cerr<<"[HEDataStructure::patchHole] For the current vertex there are "
					<<unassocEdges.size()<<" unassociated edges."<<std::endl;
			return false;
		}

		E = unassocEdges[0];

		vertList.push_back(vCurrent);
		vNext = E->pair->vertex;
		vCurrent = vNext;

	}while(vNext != vStart);


	// for(int i=0; i<vertList.size(); i++)
	// {
	// 	std::cout<<"-------------------------"<<std::endl;
	// 	std::cout<<"     Patch vertex "<<i<<std::endl;
	// 	std::cout<<"-------------------------"<<std::endl;
	// 	vertList[i]->print(std::cout,1);

	// }

	return this->addFace(vertList);
}





void HEDataStructure::splitFacesWithPlane(Vec3* planePoint, Vec3* planeNormal, std::vector<HEVertex*>& usedSplitVertices)
{

	std::multimap<HEFace*, HEVertex*> faceToVertex;

	// For every pair of edges, determine any cut vertices
	for (std::vector<HEEdge*>::iterator eIt = uniqueEdges.begin(); eIt != uniqueEdges.end(); ++eIt)
	{
		HEVertex *hevStart, *hevEnd;
		Vec3 ptStart, ptEnd;

		hevStart = (*eIt)->vertex;
		hevEnd = (*eIt)->pair->vertex;
		ptStart = hevStart->pos;
		ptEnd = hevEnd->pos;

		Edge3 line(ptStart,ptEnd);
		Vec3* ptIntersect = new Vec3;
		bool isIntersecting = false;

		// std::cout<<"-------------------------"<<std::endl;
		// std::cout<<"Edge information"<<std::endl;
		// std::cout<<"-------------------------"<<std::endl;
		// line.print(std::cout,1);
		// std::cout<<"\tfaceLeft = "<<(*eIt)->left;;
		// std::cout<<"\tfaceRght = "<<(*eIt)->pair->left<<std::endl;

		isIntersecting = Geometry::intersectionLinePlane(&line,
												    planePoint,
												    planeNormal,
												    &ptIntersect);


		//If it intersects, we should split the edge
		if(isIntersecting)
		{

			HEVertex* vertIntersect;
			HEFace *faceRight, *faceLeft;
			faceLeft = (*eIt)->left;
			faceRight = (*eIt)->pair->left;

			// Just because we're intersecting, doesn't mean we're good to go
			// We may actually be intersecting at one of the tips of the edge
			// If that's the case, we do not need to actually split the edge
			//	and should not split the edge (& create a new vertex)
			// 	We only need to record the proper vertex to use for
			//	face splitting.

			// Get the maximum norm between the starting vertex, ending vertex,
			// and intersection point in order to scale these reasonably
			double normMax = 0.0;
			double normInter	= ptIntersect->norm();
			double normStart	= hevStart->pos.norm();
			double normEnd	= hevEnd->pos.norm();
			normMax = normInter;
			if(normStart>normMax)
				normMax = normStart;
			if(normEnd>normMax)
				normMax = normEnd;


			double distFromStart	= (*ptIntersect - hevStart->pos).norm()/normMax;
			double distFromEnd	= (*ptIntersect - hevEnd->pos).norm()/normMax;

			// std::cout<<"Distance from Start/End: "<<distFromStart<<"\t"<<distFromEnd<<std::endl;

			if(distFromStart < ERRMAX)
			{
				vertIntersect = hevStart;
				faceToVertex.insert(std::pair<HEFace*, HEVertex*>(faceLeft,vertIntersect));
			}
			else if(distFromEnd < ERRMAX)
			{
				vertIntersect = hevEnd;
				faceToVertex.insert(std::pair<HEFace*, HEVertex*>(faceRight,vertIntersect));
			}
			else
			{
				vertIntersect = this->splitEdge(*eIt,ptIntersect);
				faceToVertex.insert(std::pair<HEFace*, HEVertex*>(faceRight,vertIntersect));
				faceToVertex.insert(std::pair<HEFace*, HEVertex*>(faceLeft,vertIntersect));
			}

		}
		else
		{
			continue;
		}
	}

	// For every face, split it if required
	std::vector<HEFace*> facesOld(this->faces);
	for (std::vector<HEFace*>::iterator fIt = facesOld.begin(); fIt != facesOld.end(); ++fIt)
	{
		std::multimap<HEFace*,HEVertex*>::iterator vIt;
  		std::pair<std::multimap<HEFace*,HEVertex*>::iterator,
  				 std::multimap<HEFace*,HEVertex*>::iterator> range;
		std::vector<HEVertex*> splitVerts;

#ifdef DEBUG
		std::cout<<"-------------------------"<<std::endl;
		std::cout<<"Splitting face"<<std::endl;
		std::cout<<"-------------------------"<<std::endl;
#endif
		// Loop over all vertices related to the current face
		range = faceToVertex.equal_range(*fIt);
		for (vIt = range.first; vIt != range.second; ++vIt)
		{
			splitVerts.push_back((*vIt).second);

		}

		// Only perform the edge split if there are 2 vertices in the list
		// There may be 1 vertex. This is expected due to the way we construct
		// the Face->SplitVertex map
		// THERE SHOULD NEVER BE MORE THAN 2
		// If there is, something's wrong
		std::set<HEVertex*> uniqueSplitVertices;
		int nSplitVerts = splitVerts.size();
		if(nSplitVerts==2)
		{
			this->splitFace(*fIt, splitVerts[0], splitVerts[1]);

			// Add the vertices actually used to the unique set of split vertices
			uniqueSplitVertices.insert(splitVerts[0]);
			uniqueSplitVertices.insert(splitVerts[1]);
		}
		else if(nSplitVerts > 2)
		{
			std::cerr<<"[HEDataStructure::splitFacesWithPlane] (WARNING) Found "<<nSplitVerts<<" split vertices. Should only have 2"<<std::endl;
		}


		// Copy over the unique set into the returned vector
		std::set<HEVertex*>::iterator sIt;
		for(sIt=uniqueSplitVertices.begin(); sIt!=uniqueSplitVertices.end(); ++sIt)
		{
			usedSplitVertices.push_back(*sIt);
		}

	}

}


bool HEDataStructure::deleteConnectionsToVertex(HEVertex* V)
{
	std::vector<HEEdge*> edgesLeavingOld(V->edgesLeaving);
	// First, delete all faces touching V
	int counter = 0;
	for (std::vector<HEEdge*>::iterator eIt = edgesLeavingOld.begin(); eIt != edgesLeavingOld.end(); ++eIt)
	{

		// std::cout<<"In here  "<<++counter<<std::endl;
		bool status1, status2;
		status1 = this->deleteFace((*eIt)->left);
		status2 = this->deleteFace((*eIt)->pair->left);
		if(!status1||!status2)
		{
			return false;
		}
	}

	// Second, delete all edges touching V
	for (std::vector<HEEdge*>::iterator eIt = edgesLeavingOld.begin(); eIt != edgesLeavingOld.end(); ++eIt)
	{
		bool status1, status2;
		// std::cout<<"*eIt = "<<*eIt<<std::endl;
		HEEdge* E = *eIt;
		HEEdge* pair = (*eIt)->pair;
		status2 = this->deleteEdge(pair);
		status1 = this->deleteEdge(E);
		if(!status1||!status2)
		{
			return false;
		}
	}

	// Finally, delete V
	bool statusV = this->deleteVertex(V);

	if(!statusV)
	{
		return false;
	}

	return true;
}

bool HEDataStructure::deleteFace(HEFace* F)
{
	// std::cout<<"[HEDataStructure::deleteFace] entering function"<<std::endl;
	// Correct the linkages pointing to the face
	if(!F->deleteFace())
	{
		return false;
	}

	// std::cout<<"about to delete from faces vector"<<std::endl;
	// Remove the face from our face list
	for (std::vector<HEFace*>::iterator i = this->faces.begin(); i != this->faces.end(); ++i)
	{
		// std::cout<<"*i | F"<<*i<<"\t"<<F<<std::endl;
		if( (*i) == F)
		{
			// std::cout<<"Deleting from faces vector"<<std::endl;
			this->faces.erase(i);
			break;
		}
	}
	// std::cout<<"before delete | F = "<<F<<std::endl;
	delete F;
	F = NULL;
	// std::cout<<"after delete | F = "<<F<<std::endl;

	// std::cout<<"[HEDataStructure::deleteFace] leaving function"<<std::endl;
	return true;
}

bool HEDataStructure::deleteEdge(HEEdge* E)
{

	// std::cout<<"[HEDataStructure::deleteEdge] entering function"<<std::endl;
	// if(E == NULL)
	// {
	// 	return true;
	// }


	// Correct the linkages pointing to the edge
	// if(E->pair != NULL)
	// {
	// 	if(!E->pair->deleteEdge())
	// 	{
	// 		return false;
	// 	}
	// }
	if(!E->deleteEdge())
	{
		return false;
	}

	// std::cout<<"about to delete from edges vector"<<std::endl;
	// Remove the edge from our edge lists
	std::vector<HEEdge*>::iterator i;
	for (i = this->edges.begin(); i != this->edges.end(); ++i)
	{
		// std::cout<<"*i | E"<<"\t"<<*i<<"\t"<<E<<std::endl;
		if( (*i) == E)
		{
			this->edges.erase(i);
			break;
		}
	}

	for (i = this->uniqueEdges.begin(); i != this->uniqueEdges.end(); ++i)
	{
		// std::cout<<"*i | E"<<"\t"<<*i<<"\t"<<E<<std::endl;
		if( (*i) == E)
		{
			this->uniqueEdges.erase(i);
			break;
		}
	}

	// for (int i=0; i<this->uniqueEdges.size(); i++)
	// {
	// 	if( this->uniqueEdges[i] == E )
	// 	{
	// 		// Move the last element to our current position, overwritting
	// 		// the pointer to E
	// 		// This will allow us to just pop_back after leaving the loop
	// 		this->uniqueEdges[i] = this->uniqueEdges.back();
	// 		break;
	// 	}
	// }
	// this->uniqueEdges.pop_back();


	// delete E->pair;
	// std::cout<<"before delete | E = "<<E<<std::endl;
	delete E;
	E = NULL;
	// std::cout<<"after delete | E = "<<E<<std::endl;
	// *E = NULL;
	// std::cout<<"[HEDataStructure::deleteEdge] leaving function"<<std::endl;
	return true;
}


bool HEDataStructure::deleteVertex(HEVertex* V)
{
	if(!V->deleteVertex())
	{
		return false;
	}

	for (std::vector<HEVertex*>::iterator i = this->vertices.begin(); i != this->vertices.end(); ++i)
	{
		if( (*i) == V)
		{
			this->vertices.erase(i);
			break;
		}
	}

	delete V;

	return true;
}

HEVertex* HEDataStructure::splitEdge(HEEdge* oldRight, const Vec3* V)
{
	// Create new vertex
	int vertInd = this->addVertex(V);
	HEVertex* newVertex = this->vertices[vertInd];

	// Get old pair
	HEEdge* oldLeft = oldRight->pair;


	// Create two new empty edges
	HEEdge *newRight, *newLeft;
	newRight = new HEEdge();
	newLeft  = new HEEdge();

	// Add them to the edge vector
	this->edges.push_back(newRight);
	this->edges.push_back(newLeft);

	// Update the right edges first
	// (Doesn't matter if we do the left first though)
	newRight->next = oldRight->next;
	newRight->prev = oldRight;
	if(this->embeddingDimension>2 || oldRight->next!=NULL)
	{	// we WILL NOT have a closed surface in 2d, so there is no guarantee that ->next is not null
		// The act of splitting will give edges a valid ->next, so in 2d it MAY happen that we should
		// perform this linkage
		// This had better hold true for 3d though
		oldRight->next->prev = newRight;
	}
	oldRight->next = newRight;
	newRight->pair = oldRight->pair;
	newRight->vertex = newVertex;
	newRight->left = oldRight->left;
	oldRight->pair = newLeft;

	// Update the left edges
	newLeft->next = oldLeft->next;
	newLeft->prev = oldLeft;
	if(this->embeddingDimension>2 || oldLeft->next!=NULL)
	{
		oldLeft->next->prev = newLeft;
	}
	oldLeft->next = newLeft;
	newLeft->pair = oldLeft->pair;
	newLeft->vertex = newVertex;
	newLeft->left = oldLeft->left;
	oldLeft->pair = newRight;

	// Let the new vertex know it has edges leaving it
	newVertex->addLeavingEdge(newRight);
	newVertex->addLeavingEdge(newLeft);

	// Done!
	return newVertex;
}

void HEDataStructure::splitFace(HEFace* oldFace, HEVertex* vStart, HEVertex* vEnd)
{
	// Find the TopLeft (TL), TopRight (TR), BottomLeft (BL), and BottomRight(BR)
	// edges
	// TL enters vStart, TR leaves vStart
	// BR enters vEnd, BL leaves vEnd
	// All right components are now part of the new face
	// while left components are part of the old face
	HEEdge *TL, *TR, *BL, *BR;

	std::vector<HEEdge*> edgeList;

	// Get all the edges of the face
	oldFace->getAllEdges(edgeList);
	// Loop over them until we set TL/TR/BL/BR
	for (int i = 0; i < edgeList.size(); ++i)
	{
		HEEdge* current = edgeList[i];

		// Check if we're TR
		if(current->vertex == vStart)
		{
			TR = current;
			TL = current->prev;
		}

		// Check if we're BL
		if(current->vertex == vEnd)
		{
			BL = current;
			BR = current->prev;
		}
	}

	// Create the new half edges that split the face
	// NOTE: addHalfEdge updates the vertex edge lists
	// for vStart and vEnd, so we don't need to do that
	// explicitly
	int leftInd = this->addHalfEdge(vStart,vEnd);
	HEEdge* newLeft = this->edges[leftInd];
	HEEdge* newRight = newLeft->pair;

	// update new edge linkages
	newLeft->next = BL;
	newLeft->prev = TL;
	newRight->next = TR;
	newRight->prev = BR;
	// Update face information
	// There's a chance that the edge information
	// stored in oldFace will correspond to an edge
	// in the new face we create
	// So, we forcibly reset the edge pointer in the old face
	// This lesson was a painful one to learn
	newLeft->left = oldFace;
	oldFace->edge = newLeft;


	// update old edge linkages
	// NOTE: pairs do not need to change here, since we're inside of a face
	// and not crossing borders
	TL->next = newLeft;
	BL->prev = newLeft;
	TR->prev = newRight;
	BR->next = newRight;

	// Create the new face and add it to the global list
	HEFace* newFace = new HEFace();
	this->faces.push_back(newFace);

	// Set the edge of the face to be equal to newRight
	newFace->edge = newRight;
	newRight->left = newFace;

	// Loop over all of the edges in newFace and update their face pointer
	HEEdge* tmp = newFace->edge;
	do
	{
		tmp->left = newFace;
		tmp = tmp->next;
	}while(tmp!=newFace->edge);

	do
	{
		tmp->left = newFace;
		tmp = tmp->next;
	}while(tmp!=newFace->edge);


	// Done!
}


void HEDataStructure::writeMesh(std::ostream& out)
{

	// Construct face/edge/vertex integer maps for writing an understantable file
	std::map<HEFace*,int>   faceMap;
	std::map<HEEdge*,int>   edgeMap;
	std::map<HEVertex*,int> vertMap;
	for(int i=0; i<this->faces.size(); ++i)
		faceMap[this->faces[i]] = i;
	for(int i=0; i<this->edges.size(); ++i)
		edgeMap[this->edges[i]] = i;
	for(int i=0; i<this->vertices.size(); ++i)
		vertMap[this->vertices[i]] = i;

	// Write header information
	out<<"nFaces\tnHalfEdges\tnVertices"<<std::endl;
	out<<this->faces.size()<<"\t"<<this->edges.size()<<"\t"
	   <<this->vertices.size()<<std::endl;


	// Write face-edge linkages
	std::vector<HEEdge*> edgeList;
	for(int i=0; i<this->faces.size(); ++i)
	{
		edgeList.clear();
		this->faces[i]->getAllEdges(edgeList);

		out<<faceMap[this->faces[i]]+1;
		for (int j = 0; j < edgeList.size(); ++j)
		{
			out<<"\t"<<edgeMap[edgeList[j]]+1;
		}
		out<<std::endl;
	}


	// Write edge-vertex linkages
	// Only write them if the edge is part of a face
	for(int i=0; i<this->edges.size(); ++i)
	{
		HEEdge *curr, *next;
		curr = this->edges[i];
		if(curr->next == NULL)
		// if(curr->left == NULL)
		{
			out<<edgeMap[curr]+1;
			out<<"\t"<<-1<<"\t"<<-1;
			out<<std::endl;
		}
		else
		{
			next = curr->next;

			out<<edgeMap[curr]+1;
			out<<"\t"<<vertMap[curr->vertex]+1<<"\t"<<vertMap[next->vertex]+1;
			out<<std::endl;
		}
	}

	// Write vertex information
	for(int i=0; i<this->vertices.size(); ++i)
	{
		HEVertex* curr = this->vertices[i];
		out<<vertMap[curr]+1;
		out<<"\t"<<curr->pos[0];
		out<<"\t"<<curr->pos[1];
		out<<"\t"<<curr->pos[2];
		out<<std::endl;
	}
}

void HEDataStructure::printVertices(std::ostream& out, int indent=0)
{
	for(int i=0; i<this->vertices.size(); ++i)
	{
		this->vertices[i]->print(out,indent);
	}
}

void HEDataStructure::printFaces(std::ostream& out, int indent=0)
{
	for(int i=0; i<this->faces.size(); ++i)
	{

		out<<std::string(indent,'\t')<<"HEFace information:"<<std::endl;

		HEEdge* e = this->faces[i]->edge;
		do
		{
			out<<std::string(indent+1,'\t')<<"Vertex Position: "<<std::endl;
			out<<std::string(indent+2,'\t');
			e->vertex->pos.print(out);
			out<<std::endl;
			e = e->next;

		}while(e!=this->faces[i]->edge);
	}
}

int HEDataStructure::addVertex(const Vec3* vertGeom)
{
	this->vertices.push_back(new HEVertex);
	this->vertices.back()->pos = vertGeom;
	return this->vertices.size()-1;
}



int HEDataStructure::addHalfEdge(int vStart, int vEnd)
{

	// Get vertex pointer
	HEVertex* vertStart = this->vertices[vStart];
	HEVertex* vertEnd = this->vertices[vEnd];
	return this->addHalfEdge(vertStart, vertEnd);
}

int HEDataStructure::addHalfEdge(HEVertex* vertStart, HEVertex* vertEnd)
{

	// Create new edge
	HEEdge* e1 = new HEEdge;
	HEEdge* e2 = new HEEdge;
	e1->vertex = vertStart;
	e2->vertex = vertEnd;

	// Link to edges leaving the chosen vertex
	vertStart->addLeavingEdge(e1);
	vertEnd->addLeavingEdge(e2);

	// Set them as pairs
	e1->pair = e2;
	e2->pair = e1;

	// Add new edge to list
	this->edges.push_back(e1);
	this->edges.push_back(e2);

	// Add one of the new half-edges to uniqueEdges
	// This stores the information for a geometric query on
	// lines comprising the edges. Eliminates having to
	// get single HEEdges algorithmically later
	uniqueEdges.push_back(e1);

	// return index to vertex
	return this->edges.size()-2; // return the first edge
}

bool HEDataStructure::addFace(std::vector<int> vertexList)
{
	std::vector<HEVertex*> vertexPtrList;
	for(int i=0; i<vertexList.size(); i++)
	{
		vertexPtrList.push_back(this->vertices[vertexList[i]]);
	}

	return this->addFace(vertexPtrList);
}

bool HEDataStructure::addFace(std::vector<HEVertex*> vertexList)
{
	// std::cout<<"----------------------------------------"<<std::endl;
	// std::cout<<"			Adding Face 				   "<<std::endl;
	// std::cout<<"----------------------------------------"<<std::endl;
	// Create face
	HEFace* F = new HEFace;

	// Loop over all vertices and find/create half-edges
	std::vector<HEEdge*> edgePtrs;
	std::vector<HEEdge*>::iterator eIt;
	for(int i=0; i<vertexList.size(); i++)
	{
		HEVertex *vStart, *vEnd;
		vStart 	= vertexList[i];
		vEnd    	= vertexList[(i+1)%vertexList.size()];
		// vStart 	= this->vertices[vStartInd];
		// vEnd 	= this->vertices[vEndInd];

		bool foundNeededEdge = false;
		std::vector<HEEdge*>::iterator eIt;
		for(eIt = vStart->edgesLeaving.begin(); eIt < vStart->edgesLeaving.end(); ++eIt)
		{
			if((*eIt)==NULL)
			{
				// std::cout<<"*eIt == NULL"<<std::endl;
				continue;
			}
			if((*eIt)->prev == NULL)
			{
				// std::cout<<"*eIt->prev == NULL"<<std::endl;
				continue;
			}
			if((*eIt)->prev->vertex == vEnd)
			{
				// std::cout<<"Found the edge I needed!"<<std::endl;
				foundNeededEdge = true;
				edgePtrs.push_back((*eIt)->prev->pair);
				break;
			}
		}

		// If we didn't find the pair of an already created edge
		// then we create a new edge pair
		if(!foundNeededEdge)
		{
			// Create half-edge starting at this vertex
			int edgeInd = this->addHalfEdge(vStart, vEnd);
			edgePtrs.push_back(this->edges[edgeInd]);
		}

	}

	// connect the edges together
	for(int i=0; i<edgePtrs.size(); ++i)
	{
		HEEdge* eCurr = edgePtrs[i];
		HEEdge* eNext = edgePtrs[(i+1)%edgePtrs.size()];
		HEEdge* ePrev = edgePtrs[(i-1+edgePtrs.size())%edgePtrs.size()];

		// Connect the current edge to the previous edge
		eCurr->prev = ePrev;
		// Connect the current edge to the next edge
		eCurr->next = eNext;

		// Connect to face
		eCurr->left = F;
	}

	// Connect the face to an edge
	F->edge = edgePtrs[0];

	// Compute the normal
	F->computeNormal();

	// push onto face vector
	this->faces.push_back(F);



	return true;

}
