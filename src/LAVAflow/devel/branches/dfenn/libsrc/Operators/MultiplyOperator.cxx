#include <vtkType.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDataObject.h>
#include <vtkDataObjectTypes.h>

#include "BaseOperator.h"
#include "MultiplyOperator.h"

MultiplyOperator::MultiplyOperator()
{
	scale=1.0f;
}

vtkSmartPointer<vtkDataSet> MultiplyOperator::process(const vtkSmartPointer<vtkDataSet> ds)
{
	vtkDataObject* obj = vtkDataObjectTypes::NewDataObject(ds->GetClassName());
	vtkDataSet* newDS = vtkDataSet::SafeDownCast(obj);
	vtkSmartPointer<vtkDataSet> retval = newDS;
	retval->DeepCopy(ds);
	//ds->GetCellData()->PrintSelf(cout,vtkIndent(0));
	//vtkCellData* cellData=retval->GetCellData();
	//cellData->PrintSelf(cout,vtkIndent(0));
	vtkPointData* pointData = retval->GetPointData();
	vtkCellData* cellData = retval->GetCellData();
	//for(vtkIdType i = 0; i<cellData->GetCellCount(); i++)
	//{
		//vtkDataArray* da = cellData->GetScalars();	
	if(cellData)
	{
		vtkDataArray* da = cellData->GetScalars();	
		if(da)
		{
			for(vtkIdType j=0;j<da->GetNumberOfTuples(); j++)
			{
				for(vtkIdType k = 0; k<da->GetNumberOfComponents(); k++)
				{
					da->SetComponent(j,k,da->GetComponent(j,k)*scale);
				}
			}
		}
	}
	if(pointData)
	{
		vtkDataArray* da = pointData->GetScalars();	
		if(da)
		{
			for(vtkIdType j=0;j<da->GetNumberOfTuples(); j++)
			{
				for(vtkIdType k = 0; k<da->GetNumberOfComponents(); k++)
				{
					da->SetComponent(j,k,da->GetComponent(j,k)*scale);
				}
			}
		}
	}
	return retval;
}
void MultiplyOperator::setScale(double scale)
{
	this->scale = scale;
}

double MultiplyOperator::getScale()
{
	return scale;
}
