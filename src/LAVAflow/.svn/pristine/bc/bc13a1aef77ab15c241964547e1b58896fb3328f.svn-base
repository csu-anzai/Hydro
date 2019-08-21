#ifndef SURFACEAVGPLANE_OPERATOR_H
#define SURFACEAVGPLANE_OPERATOR_H

#include <vector>
#include <string>
#include <vtkSmartPointer.h>
#include "../SurfaceAverageOperator.h"

class vtkDataSet;

class SurfaceAveragePlaneOperator : public SurfaceAverageOperator
{
	double direction[3];
	double startPt[3], endPt[3];
	int fixedAxis;
	int coordSys;
	bool isAxisAligned;
	       
	void generic_init(std::vector<std::string> varNames, bool initSurfacePoints);
	       
public:
	SurfaceAveragePlaneOperator(int nSurfaces,  double* direction, double* startPt, double* endPt, int CS, std::vector<std::string> varNames);
	SurfaceAveragePlaneOperator(int nSurfaces,  int fixedDir, double startPtFixedAxis, double endPtFixedAxis, int CS, std::vector<std::string> varNames);
	SurfaceAveragePlaneOperator(vtkSmartPointer<vtkDataArray> coordinates, int fixedDir, int coordSys, std::vector<std::string> varNames);
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds);


	void addVariable(const char* varName, double initialValue = 0.0, bool ignoreVariableList = false);

	void printSurfaces(std::ostream& out);

	static int cellIntersectsBox(double planeNormal[3], double planePoint[3], double cellEdges[6]);


};

#endif
