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
#include <vtkXMLRectilinearGridWriter.h>

#include <fstream>
#include <algorithm>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int main(int argc, char** argv)
{

	std::string filenameIn, fileNameOutBase;
	int currentCS = CS_SPHERE;
	int currAxis   = XAXIS;
	int nDim = 2;
	double MSOLAR = 1.9891e33; // solar mass in grams
	double GRAVCONST = 6.67259e-8; // cm^3/g/s
	double unitFOE = 1.e51;
	double explosionCriterion = 1.e48;


	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		// filenameIn = "../Tests/files/averaging/b155d2a3s123SG1MHwr2LB.plt.0010.vtk";
		filenameIn = "/home/tah09e/data/sne/data/HOTB/2d/b155d2a3s123SG1MHwr1LB/b155d2a3s123SG1MHwr1LB.plt.0020.vtk";
		fileNameOutBase = filenameIn;



		filenameIn = "/data1/sne/HOTB/3d/b157d3a3s421SG1MHWr3LB.plt.0100.vtk";
		nDim = 3;

		fileNameOutBase = filenameIn;

	}
	else if(argc==3)
	{
		filenameIn = argv[1];
		fileNameOutBase  = argv[2];
	}

	
	std::cout<<"Processing file \""<<filenameIn<<"\""<<std::endl;

	// Read in data
	Mesh data(nDim, currentCS, filenameIn);

	// Get simulation information
	double simTime = data.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);
	std::cout<<"Simulation time: "<<simTime<<std::endl;

	// Get domain and mesh information
	double 	physicalBounds[6]; // Physical space bounds
	int 		meshDimensions[3]; // Number of cells in each dimension
	data.getDataSet()->GetBounds(physicalBounds);
	data.getDataDimension(meshDimensions);

	// std::cout<<"Mesh dimensions:"<<std::endl;
	// std::cout<<meshDimensions[0]<<"\t"<<meshDimensions[1]<<"\t"<<meshDimensions[2]<<std::endl;


	double 	trueVolume = 4.0/3.0*vtkMath::Pi()*(	physicalBounds[1]*physicalBounds[1]*physicalBounds[1] - 
												physicalBounds[0]*physicalBounds[0]*physicalBounds[0] );

	double 	simVolume = 1.0;
	simVolume *= (1.0/3.0)*physicalBounds[1]*physicalBounds[1]*physicalBounds[1] - physicalBounds[0]*physicalBounds[0]*physicalBounds[0];
	simVolume *= -cos(physicalBounds[3]) + cos(physicalBounds[2]);
	if(nDim>2)
	{
		simVolume *= physicalBounds[5] - physicalBounds[4];
	}
	else
	{
		simVolume *= 2.0*vtkMath::Pi();
	}

	// double volumeCorrection = trueVolume/simVolume;
	double volumeCorrection = 2.0/(-cos(physicalBounds[3]) + cos(physicalBounds[2]));

	data.computeCellVolumes("cellVolumeTmp");
	data.storeData(volumeCorrection*data["cellVolumeTmp"],"cellVolume");

	// Create data arrays of coordinates (full 2d/3d arrays of cell centers)
	data.createCoordinateDataArrays("r","phi","theta");

	// Get radial coordinates of points (bounds of cells)
	vtkSmartPointer<vtkDataArray> radCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();




	/*====================================================================
		
		Compute the Reynolds stress tensor components such that
		R_ij = v_i*v'_j

		c.f. Murphy, Dolence, & Burrows (2012)
	=======================================================================*/

	//Average data
	std::vector<std::string> varNames;
	varNames.push_back("velx");
	varNames.push_back("vely");
	varNames.push_back("dens");
	if(nDim>2)
	{
		varNames.push_back("velz");
	}

	ShellAveragePlaneOperator velocityAverager(	radCoords,
												currAxis,
												currentCS,
												varNames		);

	// Process the shell averages
	velocityAverager.process(data.getDataSet());


	data.storeData(velocityAverager.removeAverageFromData(data["velx"],"velx"),"velxPrime");
	data.storeData(velocityAverager.removeAverageFromData(data["vely"],"vely"),"velyPrime");


	data.storeData(data["velxPrime"]*data["velxPrime"], "RxxTmpBothPrimed");

	data.storeData(data["velx"]*data["velxPrime"]*data["dens"], "RxxTmp");
	data.storeData(data["velx"]*data["velyPrime"]*data["dens"], "RxyTmp");
	data.storeData(data["vely"]*data["velxPrime"]*data["dens"], "RyxTmp");
	data.storeData(data["vely"]*data["velyPrime"]*data["dens"], "RyyTmp");

	if(nDim>2)
	{

		data.storeData(velocityAverager.removeAverageFromData(data["velz"],"velz"),"velzPrime");
		
		data.storeData(data["velz"]*data["velzPrime"]*data["dens"], "RzzTmp");

		data.storeData(data["velx"]*data["velzPrime"]*data["dens"], "RxzTmp");
		data.storeData(data["vely"]*data["velzPrime"]*data["dens"], "RyzTmp");
		data.storeData(data["velz"]*data["velyPrime"]*data["dens"], "RzyTmp");
		data.storeData(data["velz"]*data["velxPrime"]*data["dens"], "RzxTmp");
	}

	if(nDim==2)
	{	
		data.storeData(data["dens"]*(data["velx"]*data["velx"] + data["vely"]*data["vely"]),"rampres");
	}
	else if(nDim==3)
	{
		data.storeData(data["dens"]*(data["velx"]*data["velx"] + data["vely"]*data["vely"] + data["velz"]*data["velz"]),"rampres");	
	}
	data.storeData(data["dens"]*data["RxxTmp"],"turbrampresField");
	data.storeData(data["pres"] + data["rampres"] + data["turbrampresField"] ,"totalpres");


	//Average data
	std::vector<std::string> varNames2;
	varNames2.push_back("RxxTmp");
	varNames2.push_back("RxyTmp");
	varNames2.push_back("RyxTmp");
	varNames2.push_back("RyyTmp");
	varNames2.push_back("dens");
	varNames2.push_back("pres");
	varNames2.push_back("rampres");
	varNames2.push_back("totalpres");

	if(nDim>2)
	{
		varNames2.push_back("RzzTmp");
		varNames2.push_back("RxzTmp");
		varNames2.push_back("RzxTmp");
		varNames2.push_back("RyzTmp");
		varNames2.push_back("RzyTmp");
	}


	ShellAveragePlaneOperator shellAvg(	radCoords,
										currAxis,
										currentCS,
										varNames2	);


	shellAvg.process(data.getDataSet());


	// Divide Rij by the background density
	shellAvg.storeData(shellAvg["RxxTmp"]/shellAvg["dens"],"Rxx");
	shellAvg.storeData(shellAvg["RxyTmp"]/shellAvg["dens"],"Rxy");
	shellAvg.storeData(shellAvg["RyxTmp"]/shellAvg["dens"],"Ryx");
	shellAvg.storeData(shellAvg["RyyTmp"]/shellAvg["dens"],"Ryy");

	shellAvg.storeData(shellAvg["Rxx"]*shellAvg["dens"],"turbrampres");

	if(nDim>2)
	{
		shellAvg.storeData(shellAvg["RzzTmp"]/shellAvg["dens"],"Rzz");
		shellAvg.storeData(shellAvg["RxzTmp"]/shellAvg["dens"],"Rxz");
		shellAvg.storeData(shellAvg["RzxTmp"]/shellAvg["dens"],"Rzx");
		shellAvg.storeData(shellAvg["RyzTmp"]/shellAvg["dens"],"Ryz");
		shellAvg.storeData(shellAvg["RzyTmp"]/shellAvg["dens"],"Rzy");
	}


	// ====================================================================
		
	// 	Compute the shell averages of 
	// 	(1)	rho*v^2
	// 	(2)	pres
	// 	(3)	rho*R_xx
	// 	(4)	Rxx, Ryy, Rxy

	// =======================================================================

	// // Create data
	// data.storeData(data["dens"]*(data["velx"]*data["velx"] + data["vely"]*data["vely"]),"rampres");
	// data.storeData(data["dens"]*data["Rxx"],"turbrampres");




	// write to file
	std::string shellOutputFileName = fileNameOutBase + ".reynoldsstresses";

	std::cout<<"Writing file: "<<shellOutputFileName<<std::endl;

	std::ofstream shellOutput;
	shellOutput.open(shellOutputFileName.c_str());
	shellAvg.printShells(shellOutput);
	shellOutput.close();



//	vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
//	writer->SetFileName("reynoldsstresses.vtr");
//#if VTK_MAJOR_VERSION <= 5
//  	writer->SetInputConnection(data.getDataSet()->GetProducerPort());
//#else
//	writer->SetInputData(data.getDataSet());
//#endif
//	writer->Write();





	std::cout<<"Finished"<<std::endl;

}
