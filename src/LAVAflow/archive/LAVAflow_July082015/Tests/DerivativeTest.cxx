#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Operators/SurfaceAreaOperator.h"
#include "../libsrc/Operators/shellAveraging/ShellAverageOperator.h"
#include "../libsrc/Operators/shellAveraging/planeShells/ShellAveragePlaneOperator.h"
#include "../libsrc/Operators/surfaceAveraging/planarSurface/SurfaceAveragePlaneOperator.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/includes/LAVAconstants.h"

#include <string>
#include <vtkSphereSource.h>
#include <vtkRegularPolygonSource.h>
#include <vtkSphericalTransform.h>
#include <vtkTransformFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkRectilinearGridGeometryFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>


// void operator/(const vtkSmartPointer<vtkDataArray>& a, const vtkSmartPointer<vtkDataArray>& b)
// {
// 	std::cout<<"In division operator"<<std::endl;
// }

int main(int argc, char** argv)
{

	std::string filename = "../Tests/files/averaging/structpts_constant_2d.vtk";
	// Read in data
	Mesh* newMesh = new Mesh(2, CS_CART, filename);

	newMesh->firstDerivative(XAXIS, "vari2","vari2Dx");
	newMesh->firstDerivative(XAXIS, "vari2Dx","vari2Dxx");

	for(int i=0; i<newMesh->getDataSet()->GetNumberOfCells(); i++)
	{
		std::cout<<"d2/dx2 = "<<newMesh->getDataArray("vari2Dxx")->GetTuple1(i)<<std::endl;
	}

	newMesh->firstDerivative(YAXIS, "vari2","vari2Dy");
	newMesh->firstDerivative(XAXIS, "vari2Dy","vari2Dyy");

	for(int i=0; i<newMesh->getDataSet()->GetNumberOfCells(); i++)
	{
		std::cout<<"d2/dy2 = "<<newMesh->getDataArray("vari2Dyy")->GetTuple1(i)<<std::endl;
	}


	return 0;
}
