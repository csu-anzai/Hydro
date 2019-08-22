#ifndef SURFACEAREA_OPERATOR_H
#define SURFACEAREA_OPERATOR_H

#include "BaseOperator.h"

class vtkDataSet;

class SurfaceAreaOperator : public BaseOperator
{
	double surfaceArea;
public:
	SurfaceAreaOperator();
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds);
	double getSurfaceArea() const;
};

#endif
