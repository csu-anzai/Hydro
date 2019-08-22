#ifndef NULL_SELECTOR_H
#define NULL_SELECTOR_H

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkDataObject.h>
#include <vtkDataObjectTypes.h>

class NullSelector
{
public:
	virtual vtkSmartPointer<vtkDataSet> select(vtkDataSet* ds)
	{
		vtkDataObject* obj = vtkDataObjectTypes::NewDataObject(ds->GetClassName());
		vtkDataSet* newDS = vtkDataSet::SafeDownCast(obj);
		vtkSmartPointer<vtkDataSet> retval = newDS;
		retval->DeepCopy(ds);
		return retval;
	}
};

#endif
