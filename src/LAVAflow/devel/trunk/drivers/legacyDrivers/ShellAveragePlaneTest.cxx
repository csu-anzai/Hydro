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

	int nBins = 11;
	double dir[3], startPt[3], endPt[3];
	std::string filename = "../Tests/files/averaging/structpts_constant_2d_spherical.vtk";
    // Plane selector geometry
	startPt[0] =  0;
	startPt[1] =  0;
	startPt[2] =  0;

	endPt[0] = 1;
	endPt[1] = 1;
	endPt[2] = 0;

	dir[0] = 1;
	dir[1] = 0;
	dir[2] = 0;

	int currentCS = CS_SPHERE;
	int currAxis   = XAXIS;



	// Read in data
	Mesh* newMesh = new Mesh(2, currentCS, filename);


	//Average data
	std::vector<std::string> varNames;
	varNames.push_back("vari1");
	ShellAveragePlaneOperator shellAvg(nBins,currAxis,startPt[currAxis],endPt[currAxis], currentCS,varNames);


	shellAvg.process(newMesh->getDataSet());


	std::cout<<"------------------------------------------------------------"<<std::endl;
	std::cout<<"Bin Data"<<std::endl;
	shellAvg.printShells(std::cout);
	std::cout<<"------------------------------------------------------------"<<std::endl;



	for(int b=0; b<shellAvg.bins->GetNumberOfCells();b++)
	{
		std::cout<<shellAvg.bins->GetCellData()->GetArray("vari1")->GetTuple1(b)<<std::endl;
	}

	//Print cellToBin mapping

	std::cout<<"------------------------------------------------------------"<<std::endl;
	shellAvg.printCellToBin(std::cout);
	std::cout<<"------------------------------------------------------------"<<std::endl;
	std::cout<<"    Subtracted Values     "<<std::endl;
	// Initialize data for the new, average subtracted field
	// newMesh->addVariable("vari1Prime");


	// vtkSmartPointer<vtkDataArray> newData =  shellAvg.removeAverageFromData(newMesh->getDataSet()->GetCellData()->GetArray("vari1"), "test");
	// newMesh->storeData(newData,"test");
	// newMesh->storeData(2*newMesh->getDataArray("vari1"),"2vari1");
	newMesh->storeData(shellAvg.removeAverageFromData(newMesh->getDataArray("vari1"), "vari1"),"test");

	for(int i=0; i<newMesh->getDataSet()->GetNumberOfCells(); i++)
	{
		// std::cout<<i<<"\t"<<newData->GetTuple1(i)<<std::endl;
		std::cout<<i<<"\t"<<newMesh->getDataArray("test")->GetTuple1(i)<<std::endl;
	}




	int nSurfaces = 3;
	double surfStart	= 1.0/6.0;
	double surfEnd	= 1.0/6.0 + 2.0/3.0;
	std::vector<std::string> surfaceVars;
	surfaceVars.push_back("vari1");
	SurfaceAveragePlaneOperator surfaceAvg(nSurfaces,currAxis,surfStart, surfEnd, currentCS,surfaceVars);	

	surfaceAvg.process(newMesh->getDataSet());

	surfaceAvg.printSurfaces(std::cout);





	// std::cout<<"------------------------------------------------------------"<<std::endl;
	// std::cout<<"    Checking dataSet operators Values     "<<std::endl;

	// vtkSmartPointer<vtkDataArray> divResult = 1*newData*1 - 1/origData + newData*origData + 0;


	// for(int i=0; i<newMesh->getDataSet()->GetNumberOfCells(); ++i)
	// {
	// 	std::cout<<newData->GetTuple1(i)<<"\t"<<origData->GetTuple1(i)<<"\t"<<divResult->GetTuple1(i)<<std::endl;
	// }

	// newMesh->getDataArray("vari1")/newMesh	->getDataArray("vari2");

	// std::cout<<"----------------------------------------------------------"<<std::endl;

	// ds->PrintSelf(std::cout,vtkIndent(0));

	// std::cout<<"----------------------------------------------------------"<<std::endl;





	// Render raw data
#if 0
	vtkSmartPointer<vtkCellDataToPointData> c2p = vtkSmartPointer<vtkCellDataToPointData>::New();
	c2p->SetInputData(ds);
	c2p->Update();
	vtkSmartPointer<vtkDataSet> pointDS=c2p->GetOutput();

	vtkSmartPointer<vtkRectilinearGridGeometryFilter> geometryFilter = vtkSmartPointer<vtkRectilinearGridGeometryFilter>::New();
	geometryFilter->SetInputData(pointDS);
	geometryFilter->Update();
	vtkSmartPointer<vtkPolyData> polydata = geometryFilter->GetOutput();

	Renderer renderer(800,600);
	renderer.renderPseudoColor(polydata,0,2);
	renderer.interact();
#endif

	//Cleanup

	return 0;
}
