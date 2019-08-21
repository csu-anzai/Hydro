#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Operators/SurfaceAreaOperator.h"
#include "../libsrc/Operators/shellAveraging/ShellAverageOperator.h"
#include "../libsrc/Operators/shellAveraging/planeShells/ShellAveragePlaneOperator.h"
#include "../libsrc/Operators/surfaceAveraging/planarSurface/SurfaceAveragePlaneOperator.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/includes/LAVAconstants.h"
#include "../libsrc/Utilities/LAVAUtil.h"

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
#include <vtkXMLRectilinearGridWriter.h>

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>

// void operator/(const vtkSmartPointer<vtkDataArray>& a, const vtkSmartPointer<vtkDataArray>& b)
// {
// 	std::cout<<"In division operator"<<std::endl;
// }

int main(int argc, char** argv)
{

	std::string filenameIn, fileNameOutBase, scrFileName;
	int currentCS = CS_SPHERE;
	int currAxis   = XAXIS;
	int nDim = 2;
	double MSOLAR = 1.9891e33; // solar mass in grams
	double GRAVCONST = 6.67259e-8; // cm^3/g/s
	double unitFOE = 1.e51;
	double explosionCriterion = 1.e48;
	bool writeOutput = true;



	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		// filenameIn = "../Tests/files/averaging/b155d2a3s123SG1MHwr2LB.plt.0010.vtk";

		filenameIn = "/data1/sne/HOTB/2d/b163d2a3s231SG1MHWr2LB/b163d2a3s231SG1MHWr2LB.plt.0092.vtk";
		scrFileName = "/data1/sne/HOTB/2d/b163d2a3s231SG1MHWr2LB/scr_b163d2a3s231SG1MHWr2LB";

		// nDim = 3;
		// filenameIn = "/data1/sne/HOTB/3d/b157d3a3s503SG1MHWr3LB/b157d3a3s503SG1MHWr3LB.plt.0092.vtk";
		// scrFileName = "/data1/sne/HOTB/3d/b157d3a3s503SG1MHWr3LB/scr_b157d3a3s503SG1MHWr3LB";

		fileNameOutBase = filenameIn;
	}
	else if(argc==4)
	{	
		filenameIn = argv[1];
		fileNameOutBase  = argv[2];
		scrFileName = argv[3];
		std::cout<<"Processing file \""<<filenameIn<<"\""<<std::endl;
	}
	else
	{
		std::cout<<"Incorrect number of arguments passed! Exiting."<<std::endl;
		return -1;
	}

	// Read in data
	Mesh data(nDim, currentCS, filenameIn);

	// data.getDataSet()->PrintSelf(std::cout,vtkIndent(0));

	// Get simulation information
	double simTime = data.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);

	std::cout<<"File information:"<<std::endl;
	std::cout<<"\t             Time: "<<simTime<<std::endl;


	// Create data arrays of coordinates (full 2d/3d arrays of cell centers)
	data.createCoordinateDataArrays("r","phi","theta");



	/* ==================================================

		Compute the strain rate tensor
	
	================================================== */
	data.addVariable("sinPhi");
	data.addVariable("cosPhi");
	data.addVariable("cotPhi");

	if(nDim==2)
	{
		data.addVariable("velz",0.0);
	}

	vtkSmartPointer<vtkDataArray> uR = data["velx"];
	vtkSmartPointer<vtkDataArray> uP = data["vely"];
	vtkSmartPointer<vtkDataArray> uT = data["velz"];
	vtkSmartPointer<vtkDataArray> r = data["r"];
	vtkSmartPointer<vtkDataArray> phi = data["phi"];


	for(int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{
		data["sinPhi"]->SetTuple1(i, sin(phi->GetTuple1(i)));
		data["cosPhi"]->SetTuple1(i, cos(phi->GetTuple1(i)));
		data["cotPhi"]->SetTuple1(i, tan(vtkMath::Pi()/2.0 - phi->GetTuple1(i)));
	}


	// Compute Srr
	data.firstDerivative(XAXIS, "velx", "Srr");

	// Compute Spp
	data.firstDerivative(YAXIS, "vely", "DUpDp");
	data.storeData((data["DUpDp"]+data["velx"])/data["r"],"Spp");


	// Compute Stt
	if(nDim==3)
	{
		data.firstDerivative(ZAXIS, "velz", "DUtDt");
	}
	else
	{
		data.addVariable("DUtDt",0.0);
	}

	data.storeData(	(data["DUtDt"] + data["sinPhi"]*data["velx"] + data["cosPhi"]*data["vely"])/(data["r"]*data["sinPhi"]), "Stt");


	// Compute Srp
	data.firstDerivative(YAXIS,"velx","DUrDp");
	data.firstDerivative(XAXIS,"vely","DUpDr");
	data.storeData(	0.5*(data["DUrDp"]/data["r"] + data["DUpDr"] - data["vely"]/data["r"]), "Srp");


	// Compute Spt
	if(nDim==3)
	{
		data.firstDerivative(ZAXIS, "vely", "DUpDt");
	}
	else
	{
		data.addVariable("DUpDt",0.0);
	}
	data.firstDerivative(YAXIS,"velz","DUtDp");
	data.storeData(	0.5*(( (data["DUpDt"]/data["sinPhi"]) + data["DUtDp"] - data["velz"]*data["cotPhi"])/data["r"]), "Spt");

	std::cout<<"Finished computing Spt"<<std::endl;


	// Compute Srt
	if(nDim==3)
	{
		data.firstDerivative(ZAXIS, "velx", "DUrDt");
	}
	else
	{
		data.addVariable("DUrDt",0.0);
	}
	data.firstDerivative(XAXIS,"velz","DUtDr");
	data.storeData(	0.5*(data["DUrDt"]/(data["r"]*data["sinPhi"]) + data["DUtDr"] - data["velz"]/data["r"]), "Srt");

	std::cout<<"Finished computing Srt"<<std::endl;


	/* ==================================================

		Compute the energy dissipation rate
	
	================================================== */
	// data.storeData(data["Srr"] + data["Spp"] + data["Stt"] + (2.0*(data["Srp"] + data["Srt"] + data["Spt"])), "energyDissipation");
	data.addVariable("energyDissipation",0.0);

	for(int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{

		double edis = 0.0;

		edis += data["Srr"]->GetTuple1(i)*data["Srr"]->GetTuple1(i);
		edis += data["Spp"]->GetTuple1(i)*data["Spp"]->GetTuple1(i);
		edis += data["Stt"]->GetTuple1(i)*data["Stt"]->GetTuple1(i);


		edis += 4.0*data["Srp"]->GetTuple1(i)*data["Srp"]->GetTuple1(i);
		edis += 4.0*data["Srt"]->GetTuple1(i)*data["Srt"]->GetTuple1(i);
		edis += 4.0*data["Spt"]->GetTuple1(i)*data["Spt"]->GetTuple1(i);

		data["energyDissipation"]->SetTuple1(i,edis);

	}

	vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
	writer->SetFileName("energydissipation.vtr");
#if VTK_MAJOR_VERSION <= 5
  	writer->SetInputConnection(data.getDataSet()->GetProducerPort());
#else
	writer->SetInputData(data.getDataSet());
#endif
	writer->Write();



	std::cout<<"Finished"<<std::endl;
}



