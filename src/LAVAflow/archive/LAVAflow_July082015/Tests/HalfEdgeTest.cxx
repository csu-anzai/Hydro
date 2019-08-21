#include <iostream>
#include <fstream>
#include <vector>
#include "../libsrc/Geometry/Geometry.h"
#include "../libsrc/Geometry/Primitives/Vec3.h"
#include "../libsrc/Geometry/HEMesh/HEDataStructure.h"
using namespace std;

void generateBrick(HalfEdgeMesh& he, double bndBox[3][2])
{
	//bndBox[3][2]
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

	// create vertices
	for (int i = 0; i < nVerts; ++i)
	{
		Vec3 vertex;
		for (int j = 0; j < 3; ++j)
		{
			// cout<<bbToVerts[i][j]<<endl;
			vertex[j] = bndBox[j][bbToVerts[i][j]];
		}

		he.addVertex(&vertex);
		// he.vertices.back()->print(cout);
	}

	for (int i = 0; i < nFaces; ++i)
	{
		vector<int> vertInds(faceToVerts[i],faceToVerts[i]+4);
		he.addFace(vertInds);
	}

	// Vec3 splitVec(0.5,0.0,0.0);
	// he.splitEdge(he.faces[2]->edge, &splitVec);

	// he.splitFace(he.faces[0],he.vertices[0],he.vertices[2]);


	Vec3 planePoint(0.8,0.8,0.8);
	Vec3 planeNormal(1,1,1);


	// // Get vertices infront of plane
	// std::vector<HEVertex*> vertsInFront;
	// for(int i=0; i<he.vertices.size(); ++i)
	// {
	// 	if(Geometry::isPointInfrontOfPlane(he.vertices[i]->pos,planePoint, planeNormal) )
	// 	{
	// 		vertsInFront.push_back(he.vertices[i]);
	// 	}
	// }

	// // Split all the faces with a plane
	// std::vector<HEVertex*> splitVerts;
	// he.splitFacesWithPlane(&planePoint, &planeNormal, splitVerts);

	// for(int i=0; i<vertsInFront.size(); ++i)
	// {
	// 	// vertsInFront[i]->print(std::cout,1);
	// 	he.deleteConnectionsToVertex(vertsInFront[i]);
	// }

	// // Patch the hole
	// he.patchHole(splitVerts[0]);

	he.clipWithPlane(planePoint, planeNormal);

	// std::cout<<"Number of faces = "<<he.faces.size()<<std::endl;

	std::cout<<"Euler Characteristic = "<<he.eulerCharacteristic()<<std::endl;
	std::cout<<"Closed surface? "<<std::boolalpha<<he.isClosed()<<std::endl;


}

void simpleTest(HalfEdgeMesh& he)
{
	Vec3 a(0,0,0);
	Vec3 b(1,0,0);
	Vec3 c(1,1,0);
	Vec3 d(0,1,0);
	Vec3 e(2,0,0);
	Vec3 f(2,1,0);
	Vec3 g(3,0.5,0.5);
	int f1[] = {0,1,2,3};
	int f2[] = {1,4,5,2};
	int f3[] = {4,6,5};
	vector<int> vertInds1(f1, f1+4);
	vector<int> vertInds2(f2, f2+4);
	vector<int> vertInds3(f3, f3+3);

	he.addVertex(&a);
	he.addVertex(&b);
	he.addVertex(&c);
	he.addVertex(&d);
	he.addVertex(&e);
	he.addVertex(&f);
	he.addVertex(&g);

	he.addFace(vertInds1);
	he.addFace(vertInds2);
	he.addFace(vertInds3);
}

int main(void)
{
	double bndBox[][2] = {{0,1},
						{0,1},
						{0,1}};
	ofstream meshFile;

	// Polytope* volume = Polytope::generateBrick(bndBox);

	double xmin = 0, xmax = 1.0;
	double ymin = 0, ymax = 1.0;
	double zmin = 0, zmax = 1.0;

	// Polytope* volume = Polytope::generateBrick(xmin,xmax,ymin,ymax,zmin,zmax);
	Polytope* volume = Polytope::generateRectangle(bndBox);
	meshFile.open("half_edge_mesh_0.txt");
	volume->writeMesh(meshFile);
	meshFile.close();

	volume->printVertices(cout,2);
	volume->writeMesh(cout);

	Vec3 planePoint(0.5,0.5,0.5);
	Vec3 planeNormal(1,0,0);

	volume->clipWithPlane(planePoint, planeNormal);

	volume->printVertices(cout,2);
	meshFile.open("half_edge_mesh_1.txt");
	volume->writeMesh(meshFile);
	meshFile.close();

	Vec3 planePoint2(0.3,0.4,0.0);
	Vec3 planeNormal2(1,1,0.0);

	// volume->clipWithPlane(planePoint2, planeNormal2);

	meshFile.open("half_edge_mesh_2.txt");
	volume->writeMesh(meshFile);
	meshFile.close();




	std::cout<<"Euler Characteristic = "<<volume->eulerCharacteristic()<<std::endl;
	std::cout<<"Closed surface? "<<std::boolalpha<<volume->isClosed()<<std::endl;


	meshFile.open("half_edge_mesh.txt");
	volume->writeMesh(meshFile);
	meshFile.close();

	volume->writeMesh(cout);

	delete volume;

}
