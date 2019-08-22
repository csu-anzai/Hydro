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
#include <cmath>
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

//#define WRITEVTK

// void operator/(const vtkSmartPointer<vtkDataArray>& a, const vtkSmartPointer<vtkDataArray>& b)
// {
// 	std::cout<<"In division operator"<<std::endl;
// }

int main(int argc, char** argv)
{

	std::string filenameIn, fileNameOutBase;
	int currentCS = CS_SPHERE;
	int currAxis   = XAXIS;
	int nDim = 2;
	double GRAVCONST = 6.67259e-8; // cm^3/g/s
	
	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		filenameIn = "../Tests/files/averaging/b155d2a3s123SG1MHwr2LB.plt.0010.vtk";

		filenameIn = "/data1/sne/HOTB/3d/b157d3a3s421SG1MHWr3LB.plt.0100.vtk";
		nDim = 3;


		filenameIn = "b163d2a3s530SG1MHWr2LB.plt.0085.vtk";
		nDim = 2;

		fileNameOutBase = filenameIn;
	}
	else if(argc==3)
	{
		filenameIn = argv[1];
		fileNameOutBase  = argv[2];
		std::cout<<"Processing file \""<<filenameIn<<"\""<<std::endl;
	}

	// Read in data
	Mesh data(nDim, currentCS, filenameIn);

	data.getDataSet()->PrintSelf(std::cout,vtkIndent(0));

	// Get simulation information
	double simTime = data.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);
	double pointMass = data.getDataSet()->GetFieldData()->GetArray("PMASS")->GetTuple1(0);
	std::cout<<"Simulation time = "<<simTime<<std::endl;

	// Get domain and mesh information
	double 	physicalBounds[6]; // Physical space bounds
	int 	meshDimensions[3]; // Number of points in each dimension

	data.getDataSet()->GetBounds(physicalBounds);
	vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetDimensions(meshDimensions); // These are the dimensions of points, not cells. dimCells=dsDim-1

	std::cout<<"Physical dimensions: "<<std::endl;
	std::cout<<"\t"<<physicalBounds[0]<<"\t"<<physicalBounds[1]<<std::endl;
	std::cout<<"\t"<<physicalBounds[2]<<"\t"<<physicalBounds[3]<<std::endl;
	std::cout<<"\t"<<physicalBounds[4]<<"\t"<<physicalBounds[5]<<std::endl;

	std::cout<<"Mesh size: "<<std::endl;
	std::cout<<"\t"<<meshDimensions[0]<<"\t"<<meshDimensions[1]<<"\t"<<meshDimensions[2]<<std::endl;

	double trueVolume = 1.0;
	trueVolume *=	1.0/3.0;
	trueVolume *=	physicalBounds[1]*physicalBounds[1]*physicalBounds[1] -
					physicalBounds[0]*physicalBounds[0]*physicalBounds[0];
	trueVolume *=	-cos(physicalBounds[3]) + cos(physicalBounds[2]);
	if(nDim>2)
	{
		trueVolume *= physicalBounds[5]-physicalBounds[4];
	}
	else
	{
		trueVolume *= 2.0*vtkMath::Pi();
	}
	std::cout<<"True volume: "<<trueVolume<<std::endl;

	data.computeCellVolumes();

	double totalVolume = data.sum("volume");

	double volCorrection = trueVolume/totalVolume;


	// Create data arrays of coordinates (full 2d/3d arrays of cell centers)
	std::cout<<"Creating coordinate arrays..."; std::cout.flush();
	data.createCoordinateDataArrays("r","phi","theta");
	std::cout<<"done!"<<std::endl;

	// Get radial coordinates of points (bounds of cells)
	vtkSmartPointer<vtkDataArray> radCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();

	std::cout<<"Creating prereq variables"<<std::endl;
	// Create new variables
	vtkSmartPointer<vtkDataArray> kinetic;
	if(nDim==2)
	{
		kinetic = 0.5*(data["velx"]*data["velx"] + data["vely"]+data["vely"]);
	}
	else if(nDim==3)
	{
		kinetic = 0.5*(data["velx"]*data["velx"] + data["vely"]+data["vely"] + data["velz"]*data["velz"]);
	}

	vtkSmartPointer<vtkDataArray> eint = (1.0/(data["gammae"]-1.0))*data["pres"]/data["dens"];
	vtkSmartPointer<vtkDataArray> enthalpy = eint + data["pres"]/data["dens"];
	vtkSmartPointer<vtkDataArray> radMom = data["velx"]*data["dens"];
	vtkSmartPointer<vtkDataArray> volWork = data["pres"]/data["dens"] + kinetic;
	data.storeData(kinetic,"kinetic");
	data.storeData(eint,"eint");
	data.storeData(enthalpy,"enthalpy");
	data.storeData(radMom,"radMom");
	data.storeData(volWork,"volWork");

	// Compute the divergence of the velocity field
	std::cout<<"Computing the divergence of the velocity field"<<std::endl;
	if(nDim==2)
	{
		vtkSmartPointer<vtkDataArray> rUr = data["r"]*data["velx"];
		data.storeData(rUr,"xUx");
		data.firstDerivative(XAXIS,"xUx","xUxDx");
		data.firstDerivative(YAXIS,"vely","UyDy");
		vtkSmartPointer<vtkDataArray> divU = data["xUxDx"]/data["r"] + data["UyDy"]/data["r"];
		data.storeData(divU,"divU");
	}
	else if(nDim==3)
	{

		data.addVariable("sinPhi");
		vtkSmartPointer<vtkDataArray> sinPhi = data["sinPhi"];
		vtkSmartPointer<vtkDataArray> phi = data["phi"];
		for(int i=0; i<data.getDataSet()->GetCellData()->GetNumberOfTuples(); i++)
		{
			double sphi = sin(phi->GetTuple1(i));
			sinPhi->SetTuple1(i,sphi);
		}

		vtkSmartPointer<vtkDataArray> r2Ur = data["r"]*data["r"]*data["velx"];
		data.storeData(r2Ur,"x2Ux");
		data.firstDerivative(XAXIS,"x2Ux","x2UxDx");


		vtkSmartPointer<vtkDataArray> sinUy = sinPhi*data["vely"];
		data.storeData(sinUy,"sinUy");

		data.firstDerivative(YAXIS,"sinUy","sinUyDy");


		data.firstDerivative(ZAXIS,"velz","UzDz");

		vtkSmartPointer<vtkDataArray> divU = data["x2UxDx"]/(data["r"]*data["r"])+ (data["sinUyDy"]+data["UzDz"])/(sinPhi*data["r"]);
		data.storeData(divU,"divU");
	}


	// Gravity terms required
	data.storeData( GRAVCONST*pointMass/data["r"], "gpotPM" );
	data.storeData(0.5*(data["gpot"]-data["gpotPM"]),"gravity");
	data.firstDerivative(XAXIS,"gravity","gpotDr");

	//Average data
	std::cout<<"Averaging data"<<std::endl;
	std::vector<std::string> varNames;
	varNames.push_back("dens");
	varNames.push_back("velx");
	varNames.push_back("vely");
	varNames.push_back("pres");
	varNames.push_back("temp");
	varNames.push_back("enthalpy");
	varNames.push_back("kinetic");
	// varNames.push_back("shock");
	varNames.push_back("qenr");

	ShellAveragePlaneOperator shellAvg(	radCoords,
										currAxis,
										currentCS,
										varNames		);

	// Process the shell averages
	shellAvg.process(data.getDataSet());


	// Create kernels for flux equations
	std::cout<<"Creating flux kernels"<<std::endl;

	data.storeData(shellAvg.removeAverageFromData(data["enthalpy"],"enthalpy"),"enthalpyPrime");
	data.storeData(shellAvg.removeAverageFromData(data["kinetic"],"kinetic"),"kineticPrime");
	data.storeData(shellAvg.removeAverageFromData(data["pres"],"pres"),"presPrime");

	vtkSmartPointer<vtkDataArray> FcKernel = data["velx"]*data["dens"]*data["enthalpyPrime"];
	vtkSmartPointer<vtkDataArray> FkKernel = data["velx"]*data["dens"]*data["kineticPrime"];
	vtkSmartPointer<vtkDataArray> FpKernel = -1*data["velx"]*data["presPrime"];


	
	vtkSmartPointer<vtkDataArray> PaKernel = -1*data["velx"]*data["gpotDr"]*
								 shellAvg.removeAverageFromData(data["dens"],"dens");


	vtkSmartPointer<vtkDataArray> PpKernel = data["divU"]*
								 shellAvg.removeAverageFromData(data["pres"],"pres");



	data.storeData(FcKernel,"FcKernel");
	data.storeData(FkKernel,"FkKernel");
	data.storeData(FpKernel,"FpKernel");
	data.storeData(PaKernel,"PaKernel");
	data.storeData(PpKernel,"PpKernel");

	// Create surface points integrator
	vtkSmartPointer<vtkDoubleArray> surfacePoints = vtkSmartPointer<vtkDoubleArray>::New();
	for(int i=0; i<radCoords->GetNumberOfTuples()-1; i++)
	{
		double binLeft, binRight, surfPt;
		binLeft		= radCoords->GetTuple1(i);
		binRight 	= radCoords->GetTuple1(i+1);

		surfPt = 0.5*(binRight*binRight*binRight + binLeft*binLeft*binLeft);
		surfPt = pow(surfPt,1.0/3.0);

		surfacePoints->InsertNextValue(surfPt);
	}

	// Integrate the averaged qenr from the outside radius (with initial condition = 0)
	// to obtain the kernel for Fr
	shellAvg.addVariable("FrKernelIntegration");
	vtkSmartPointer<vtkDataArray> qenrAvg = shellAvg.getDataArray("qenr");
	vtkSmartPointer<vtkDataArray> FrKernelIntegration = shellAvg.getDataArray("FrKernelIntegration");
	int nBins = shellAvg.getNumberOfBins();

	// Set initial condition
	FrKernelIntegration->SetTuple1(nBins-1,0.0);
	for(int i=nBins-2; i>=0; i--)
	{
		// i+1 is actually i-1, since we're going backward
		double f = 0.5*(qenrAvg->GetTuple1(i+1) + qenrAvg->GetTuple1(i));
		double h = surfacePoints->GetTuple1(i) - surfacePoints->GetTuple1(i+1);
		double value = FrKernelIntegration->GetTuple1(i+1) + h*f;
		FrKernelIntegration->SetTuple1(i,value);
		// std::cout<<"FrKernelValue = "<<value<<std::endl;
	}

	// Create full field populated with -FrKernelIntegration
	data.addVariable("FrKernelTmp");
	vtkSmartPointer<vtkDataArray> FrKernel = shellAvg.removeAverageFromData(data["FrKernelTmp"],"FrKernelIntegration");
	data.storeData(FrKernel,"FrKernel");

	std::vector<std::string> surfaceVars;
	surfaceVars.push_back("FcKernel");
	surfaceVars.push_back("FkKernel");
	surfaceVars.push_back("FpKernel");
	surfaceVars.push_back("FrKernel");
	surfaceVars.push_back("PaKernel");
	surfaceVars.push_back("PpKernel");

	surfaceVars.push_back("radMom");
	surfaceVars.push_back("enthalpy");
	surfaceVars.push_back("kinetic");
	surfaceVars.push_back("gravity");
	surfaceVars.push_back("volWork");

	surfaceVars.push_back("pres");
	surfaceVars.push_back("divU");
	surfaceVars.push_back("velx");
	surfaceVars.push_back("dens");
	surfaceVars.push_back("gpotDr");


	std::cout<<"Surface averaging..."<<std::endl;
	SurfaceAveragePlaneOperator surfaceAvg(	surfacePoints,
											currAxis,
											currentCS,
											surfaceVars	);	

	surfaceAvg.process(data.getDataSet());
	

	// Create F_E, F_{E,K}, and P_{E,K} from surfaced averaged quantities
	vtkSmartPointer<vtkDataArray> Fe = surfaceAvg.getDataArray("radMom")*(surfaceAvg.getDataArray("enthalpy") 
																		+ surfaceAvg.getDataArray("kinetic") 
																		+ surfaceAvg.getDataArray("gravity"));
	surfaceAvg.storeData(Fe,"Fe");

	vtkSmartPointer<vtkDataArray> Fek = surfaceAvg.getDataArray("radMom")*surfaceAvg.getDataArray("volWork");
	surfaceAvg.storeData(Fek,"Fek");
	
	vtkSmartPointer<vtkDataArray> Pek = surfaceAvg.getDataArray("radMom")*(surfaceAvg.getDataArray("pres")*surfaceAvg.getDataArray("divU") 
																		+ surfaceAvg.getDataArray("velx")*surfaceAvg.getDataArray("dens")*surfaceAvg.getDataArray("gpotDr"));
	surfaceAvg.storeData(Pek,"Pek");


	// Determine the total internal energy and the total kinetic energy in the system
	data.storeData(data["dens"]*data["volume"]*volCorrection,"mass");
	data.storeData(data["eint"]*data["mass"],"eintMassWeighted");
	data.storeData(data["kinetic"]*data["mass"],"ekinMassWeighted");
	double E = data.sum("eintMassWeighted");
	double Ek = data.sum("ekinMassWeighted");

	// Store the total internal and kinetic energies in the shell file (constant as a function of radius)
	surfaceAvg.addVariable("E",E);
	surfaceAvg.addVariable("Ek",Ek);

	// print to file
	std::string shellOutputFileName = fileNameOutBase + ".shell";
	std::ofstream shellOutput;
	shellOutput.open(shellOutputFileName.c_str());
	shellOutput<<"simTime\t"<<simTime<<std::endl;
	shellAvg.printShells(shellOutput);
	shellOutput.close();


	std::string surfOutputFileName = fileNameOutBase + ".surf";
	std::ofstream surfOutput;
	surfOutput.open(surfOutputFileName.c_str());
	surfOutput<<"simTime\t"<<simTime<<std::endl;
	surfaceAvg.printSurfaces(surfOutput);
	surfOutput.close();

#ifdef WRITEVTK
	
		vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
		writer->SetFileName("flux.vtr");
	#if VTK_MAJOR_VERSION <= 5
	  	writer->SetInputConnection(data.getDataSet()->GetProducerPort());
	#else
		writer->SetInputData(data.getDataSet());
	#endif
		writer->Write();

#endif

	std::cout<<"Finished"<<std::endl;

	return 0;
}
