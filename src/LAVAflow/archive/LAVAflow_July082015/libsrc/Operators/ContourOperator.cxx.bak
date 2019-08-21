//#include <vtkContourFilter.h>
#include <vtkMarchingContourFilter.h>

#include "BaseOperator.h"
#include "ContourOperator.h"

ContourOperator::ContourOperator()
{
	//contourFilter=vtkSmartPointer<vtkContourFilter>::New();
	contourFilter=vtkSmartPointer<vtkMarchingContourFilter>::New();
}

vtkSmartPointer<vtkDataSet> ContourOperator::process(const vtkSmartPointer<vtkDataSet> ds)
{
	//When we begin dealing with AMR we may need to do the following check
	//if(vtkOverlappingAMR::SafeDownCast(ds))
	//{
	// do DualContouring...
	//}
	contourFilter->SetInputData(ds);
	//contourFilter->PrintSelf(cout,vtkIndent(5));
	contourFilter->Update();
	return contourFilter->GetOutput();
}
void ContourOperator::getIsoLevels(double* levels)
{
	contourFilter->GetValues(levels);
}
void ContourOperator::setIsoLevels(unsigned int numLevels, const double* levels)
{
	contourFilter->SetNumberOfContours(numLevels);
	for(unsigned int i = 0; i<numLevels; i++)
		contourFilter->SetValue(i,levels[i]);
}

void ContourOperator::setArrayName(int vtkFieldAssociation, const char* name)
{
	contourFilter->SetInputArrayToProcess(0,0,0,vtkFieldAssociation,name);
}
