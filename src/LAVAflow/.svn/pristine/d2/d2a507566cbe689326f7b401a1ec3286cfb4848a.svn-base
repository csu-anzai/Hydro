#ifndef MULTIPLY_OPERATOR_H
#define MULTIPLY_OPERATOR_H

#include "BaseOperator.h"

class vtkDataSet;

class MultiplyOperator : public BaseOperator
{
	double scale;
public:
	MultiplyOperator();
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds);
	void setScale(double s);
	double getScale();
};

#endif
