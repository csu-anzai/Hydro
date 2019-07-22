#ifndef CURVATURE_OPERATOR_H
#define CURVATURE_OPERATOR_H

#include "BaseOperator.h"

class vtkDataSet;
class vtkCurvatures;

class CurvatureOperator : public BaseOperator
{
	vtkSmartPointer<vtkCurvatures> curvatures;
public:
	CurvatureOperator();
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds);
	void useMeanCurvature();
	void useMaximumCurvature();
	void useMinimumCurvature();
};

#endif
