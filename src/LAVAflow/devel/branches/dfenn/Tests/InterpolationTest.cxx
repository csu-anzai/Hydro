#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>



#include <iostream>
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/includes/LAVAconstants.h"



int main(void)
{
	// Create a grid		
	vtkSmartPointer<vtkRectilinearGrid> grid =
	vtkSmartPointer<vtkRectilinearGrid>::New();

	grid->SetDimensions(2,3,1);

	vtkSmartPointer<vtkDoubleArray> xArray =
	vtkSmartPointer<vtkDoubleArray>::New();
	xArray->InsertNextValue(0.0);
	xArray->InsertNextValue(2.0);

	vtkSmartPointer<vtkDoubleArray> yArray =
	vtkSmartPointer<vtkDoubleArray>::New();
	yArray->InsertNextValue(0.0);
	yArray->InsertNextValue(1.0);
	yArray->InsertNextValue(2.0);

	vtkSmartPointer<vtkDoubleArray> zArray =
	vtkSmartPointer<vtkDoubleArray>::New();
	zArray->InsertNextValue(0.0);

	grid->SetXCoordinates(xArray);
	grid->SetYCoordinates(yArray);
	grid->SetZCoordinates(zArray);


	int nDim = 2;
	Mesh data(nDim,CS_CART);
	data.setDataSet(grid);

	data.addVariable("test",1.32);

	double pt[3];

	pt[0] = 1.0;
	pt[1] = 0.5;
	pt[2] = 0.0;
	std::cout<<data.interpolate(pt,"test")<<std::endl;


	pt[0] = 1.0;
	pt[1] = 1.5;
	pt[2] = 0.0;

	std::cout<<data.interpolate(pt,"test")<<std::endl;


}