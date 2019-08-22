#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/Geometry/Geometry.h"
#include "../libsrc/Geometry/Primitives/Vec3.h"
#include "../libsrc/Geometry/Primitives/Edge3.h"

#include <iostream>
#include <string>


using namespace std;

int main(int argc, char** argv)
{
	string fileName = "/home/tah09e/code/workspace/LAVAflow/Tests/files/oneCell_2d.vtk";
	
	Mesh m(1,2,fileName);

	//m.getDataSet()->PrintSelf(std::cout,vtkIndent(0));


	// Vec3* p0 = new Vec3(0.,0.,0.);
	// Vec3* p1 = new Vec3(0.,0.,2.);
	// Vec3* pPlane = new Vec3(-2.,4.5,3.0);
	// Vec3* planeNormal = new Vec3(0.,0.,1.0);
	// Vec3* pIntersection = new Vec3;
	// Edge3* line = new Edge3(*p0,*p1);
	// bool doesIntersect = Geometry::intersectionLinePlane(line, pPlane, planeNormal,  &pIntersection);


	Vec3 p0(0.,0.,0.);
	Vec3 p1(0.,0.,2.);
	Vec3 pPlane(-2.,4.5,3.0);
	Vec3 planeNormal(0.,0.,1.0);
	Vec3* pIntersection = new Vec3;
	Edge3 line(p0,p1);
	bool doesIntersect = Geometry::intersectionLinePlane(&line, &pPlane, &planeNormal,  &pIntersection);


	cout<<"Intersecting? "<< doesIntersect<<endl;
	cout<<"Intersection pt: \n\t";
	pIntersection->print(cout);
	cout<<endl;

	// Geometry::clipPolygonUsingPlane();
	
}
