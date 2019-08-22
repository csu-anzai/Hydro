#ifndef SHELLAVGPLANE_OPERATOR_H
#define SHELLAVGPLANE_OPERATOR_H

#include <vector>
#include <string>
#include <vtkSmartPointer.h>
#include "../ShellAverageOperator.h"

class vtkDataSet;

class ShellAveragePlaneOperator : public ShellAverageOperator
{
	double direction[3];
	double startPt[3], endPt[3];
	int fixedAxis;
	int coordSys;
	bool isAxisAligned;


	void generic_init(std::vector<std::string> varNames, bool initBinCoordinates);
	       
public:
	ShellAveragePlaneOperator(int nBins,  double* direction, double* startPt, double* endPt, int CS, std::vector<std::string> varNames);
	ShellAveragePlaneOperator(int nBins,  int fixedDir, double startPtFixedAxis, double endPtFixedAxis, int CS, std::vector<std::string> varNames);
	ShellAveragePlaneOperator(vtkSmartPointer<vtkDataArray> coordinates, int fixedDir, int coordSys, std::vector<std::string> varNames);
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds);

	void addVariable(const char* varName);

	void printShells(std::ostream& out);

	static int cellIntersectsBox(double planeNormal[3], double planePoint[3], double cellEdges[6]);
};

#endif
