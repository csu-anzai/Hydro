#ifndef OPERATOR_H
#define OPERATOR_H


#include <vtkSmartPointer.h>
#include <vtkDataSet.h>


class BaseOperator
{
	public:
	BaseOperator(){}
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds)=0;
};

#endif
