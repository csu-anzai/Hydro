#ifndef LAVA_HEDATASTRUCTURE_H
#define LAVA_HEDATASTRUCTURE_H

#include "HEVertex.h"
#include "HEEdge.h"
#include "HEFace.h"
#include <vector>
#include <ostream>
#include <string>
#include <map>
#include <algorithm>
#include "../Primitives/Vec3.h"
#include "../Primitives/Edge3.h"
#include "../Geometry.h"

// template <class T, int N>
// class HEVertex;

// template <class T, int N>
// class HEEdge;

// template <class T, int N>
// class HEFace;

class HEDataStructure
{
public:
	HEDataStructure(){}
	~HEDataStructure(){}

	/* data */
	std::vector<HEVertex*> vertices;
	std::vector<HEEdge*>   edges;
	std::vector<HEFace*>   faces;
	std::vector<HEEdge*>   uniqueEdges;
	int embeddingDimension;
	
	void printVertices(std::ostream& out, int indent);
	void printFaces(std::ostream& out, int indent);
	
	int addVertex(const Vec3* vertGeom);
	int addHalfEdge(int vStart, int vEnd);
	int addHalfEdge(HEVertex* vertStart, HEVertex* vertEnd);
	bool addFace(std::vector<int> vertexList);
	bool addFace(std::vector<HEVertex*> vertexList);
	void writeMesh(std::ostream& out);
	
	HEVertex* splitEdge(HEEdge* oldRight, const Vec3* V);
	void splitFace(HEFace* oldFace, HEVertex* vStart, HEVertex* vEnd);
	void splitFacesWithPlane(Vec3* planePoint, Vec3* planeNormal, std::vector<HEVertex*>& newVertices);

	bool deleteFace(HEFace* F);
	bool deleteEdge(HEEdge* E);
	bool deleteVertex(HEVertex* V);
	bool deleteConnectionsToVertex(HEVertex* V);
	bool patchHole(HEVertex* vStart);
	bool clipWithPlane(Vec3 planePoint, Vec3 planeNormal);
	void computeNormals();


	int eulerCharacteristic();
	bool isClosed();

	// Factory Methods
	static HEDataStructure* generateBrick(double bnds[6]);
	static HEDataStructure* generateBrick(double xbnds[2], double ybnds[2], double zbnds[2]);
	static HEDataStructure* generateBrick(double bndBox[3][2]);
	static HEDataStructure* generateBrick(double xmin, double xmax, double ymin, double ymax,
									  double zmin, double zmax);

	static HEDataStructure* generateRectangle(double bnds[6]);
	static HEDataStructure* generateRectangle(double bndBox[3][2]);
};



typedef HEDataStructure HalfEdgeMesh;
typedef HEDataStructure Polytope;


#endif