#ifndef OPERATOR_H
#define OPERATOR_H

class Operator
{
	virtual vtkDataSet* process(vtkDataSet* ds);
};

#endif
