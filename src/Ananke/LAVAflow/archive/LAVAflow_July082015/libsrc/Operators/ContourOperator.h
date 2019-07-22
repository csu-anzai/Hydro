/*
 * author: Andrew Young
 * brief: The following is a wrapper for the vtkContouring class. The
 * objective is to read in any dataset and output a polydataset.
 */
#ifndef CONTOUR_OPERATOR_H
#define CONTOUR_OPERATOR_H


#include "BaseOperator.h"

//class vtkContourFilter;
class vtkMarchingContourFilter;

class ContourOperator : public BaseOperator
{
	//vtkSmartPointer<vtkContourFilter> contourFilter;
	vtkSmartPointer<vtkMarchingContourFilter> contourFilter;
public:
	ContourOperator();
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds);
	void setIsoLevels(unsigned int numLevels, const double* levels);
        //assumes levels is the appropriate size..
	void getIsoLevels(double* levels);
	void setArrayName(int vtkFieldAssociation, const char* name);
};

#endif
