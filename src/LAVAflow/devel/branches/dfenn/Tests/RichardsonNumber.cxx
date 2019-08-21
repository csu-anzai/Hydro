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
#include <cmath>
#include <complex>

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

		// filenameIn = "/data1/sne/HOTB/2d/b163d2a3s231SG1MHWr2LB/b163d2a3s231SG1MHWr2LB.plt.0092.vtk";
		filenameIn = "/data1/sne/HOTB/2d/b163d2a3s231SG1MHWr2LB/b163d2a3s231SG1MHWr2LB.plt.0120.vtk";
		scrFileName = "/data1/sne/HOTB/2d/b163d2a3s231SG1MHWr2LB/scr_b163d2a3s231SG1MHWr2LB";

		// nDim = 3;
		// filenameIn = "/data1/sne/HOTB/3d/b157d3a3s503SG1MHWr3LB/b157d3a3s503SG1MHWr3LB.plt.0092.vtk";
		// scrFileName = "/data1/sne/HOTB/3d/b157d3a3s503SG1MHWr3LB/scr_b157d3a3s503SG1MHWr3LB";

		fileNameOutBase = filenameIn;
	}
	else if(argc==3)
	{	
		filenameIn = argv[1];
		fileNameOutBase  = argv[2];
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
	double pointMass = data.getDataSet()->GetFieldData()->GetArray("PMASS")->GetTuple1(0);

	std::cout<<"File information:"<<std::endl;
	std::cout<<"\t             Time: "<<simTime<<std::endl;


	// Create data arrays of coordinates (full 2d/3d arrays of cell centers)
	data.createCoordinateDataArrays("r","phi","theta");



	/* ==================================================

		Compute the gradients of pressure and density
	
	================================================== */
	data.addVariable("sinPhi");

	if(nDim==2)
	{
		data.addVariable("velz",0.0);
	}

	vtkSmartPointer<vtkDataArray> phi = data["phi"];
	for(int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{
		data["sinPhi"]->SetTuple1(i, sin(phi->GetTuple1(i)));
	}


	// Compute radial components of gradients
	data.firstDerivative(XAXIS,"dens","gradDens_r");
	data.firstDerivative(XAXIS,"pres","gradPres_r");
	data.firstDerivative(XAXIS,"velx","gradVelx_r");
	data.firstDerivative(XAXIS,"vely","gradVely_r");

	if(nDim==3)
	{
		data.firstDerivative(XAXIS,"velz","gradVelz_r");
	}
	else
	{
		data.addVariable("gradVelz_r",0.0);
	}

	// Compute phi components of gradients
	data.firstDerivative(YAXIS, "dens", "DdensDphi");
	data.firstDerivative(YAXIS, "pres", "DpresDphi");
	
	data.storeData(data["DdensDphi"]/data["r"],"gradDens_phi");
	data.storeData(data["DpresDphi"]/data["r"],"gradPres_phi");

	// Compute theta components of gradients
	if(nDim==3)
	{
		data.firstDerivative(ZAXIS, "dens", "DdensDtheta");
		data.firstDerivative(ZAXIS, "pres", "DpresDtheta");
		data.storeData(data["DdensDtheta"]/(data["sinPhi"]*data["r"]),"gradDens_theta");
		data.storeData(data["DpresDtheta"]/(data["sinPhi"]*data["r"]),"gradPres_theta");
	}
	else
	{
		data.addVariable("gradDens_theta",0.0);
		data.addVariable("gradPres_theta",0.0);
	}


	/* ==================================================

		Compute the square of the Brunt-Vaisala frequency

		N^2 = (gradP/rho)dot(gradRho/rho-gradP/(gammac*P))

		See Tassoul (78), Eqn. 52
	
	================================================== */
	vtkSmartPointer<vtkDataArray> Nsqr_r;
	vtkSmartPointer<vtkDataArray> Nsqr_p;
	vtkSmartPointer<vtkDataArray> Nsqr_t;

	Nsqr_r = (data["gradPres_r"]/data["dens"])*((data["gradDens_r"]/data["dens"]) - (data["gradPres_r"]/(data["gammac"]*data["pres"])));

	Nsqr_p = (data["gradPres_phi"]/data["dens"])*(data["gradDens_phi"]/data["dens"] - data["gradPres_phi"]/(data["gammac"]*data["pres"]));

	Nsqr_t = (data["gradPres_theta"]/data["dens"])*(data["gradDens_theta"]/data["dens"] - data["gradPres_theta"]/(data["gammac"]*data["pres"]));

	data.storeData(Nsqr_r, "Nsqr_r");
	data.storeData(Nsqr_p, "Nsqr_p");
	data.storeData(Nsqr_r + Nsqr_p + Nsqr_t, "Nsqr");



	/* ==================================================

		Compute the square of the Brunt-Vaisala frequency
		for the alternate description given by Foglizzo, Scheck, & Janka (2006)

		N^2 = |dPhi/dr|*[(dP/dr)/(gammac*P) - (dRho/dr)/rho]

	================================================== */
	data.storeData( GRAVCONST*pointMass/data["r"], "gpotPM" );
	data.storeData(0.5*(data["gpot"]-data["gpotPM"]),"gravity");
	data.firstDerivative(XAXIS,"gravity","gpotDr");

	vtkSmartPointer<vtkDataArray> omega_buoy;
	omega_buoy = (-1.0*data["gravity"])*( data["gradPres_r"]/(data["gammac"]*data["pres"]) - data["gradDens_r"]/data["dens"] );

	data.storeData(omega_buoy,"omegaSqr");

	/* ==================================================

		Decompose the Brunt-Vaisala frequency into its 
		real and imaginary components (take the square root)
	
	================================================== */
	data.addVariable("BVreal");
	data.addVariable("BVimag");
	data.addVariable("Omegareal");
	data.addVariable("Omegaimag");
	for(int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{
		double Nsqr = data["Nsqr"]->GetTuple1(i);
		double Osqr = data["omegaSqr"]->GetTuple1(i);
		std::complex<double> N = sqrt(std::complex<double>(Nsqr,0.0));
		std::complex<double> O = sqrt(std::complex<double>(Osqr,0.0));
		data["BVreal"]->SetTuple1(i,N.real());
		data["BVimag"]->SetTuple1(i,N.imag());
		data["Omegareal"]->SetTuple1(i,O.real());
		data["Omegaimag"]->SetTuple1(i,O.imag());
	}

	/* ==================================================

		Compute the Richardson number as 

		Ri = N^2/((dVelx/dr)^2)
	
	================================================== */

	data.storeData(data["Nsqr"]/((data["gradVelx_r"]*data["gradVelx_r"]) + 0.0001),"Ri");
	// data.storeData(data["Nsqr"]/((data["gradVely_r"]*data["gradVely_r"]) + (data["gradVelz_r"]*data["gradVelz_r"])),"Ri");


	/* ==================================================

		Average the imaginary, real, and squared BV freq,
		and radial velocity

	================================================== */
	// Get radial coordinates of points (bounds of cells)
	vtkSmartPointer<vtkDataArray> radCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();		
	std::vector<std::string> varNames;
	varNames.push_back("Nsqr");
	varNames.push_back("BVreal");
	varNames.push_back("BVimag");
	varNames.push_back("omegaSqr");
	varNames.push_back("Omegareal");
	varNames.push_back("Omegaimag");
	varNames.push_back("Ri");
	varNames.push_back("velx");

	ShellAveragePlaneOperator shellAvg(	radCoords,
										currAxis,
										currentCS,
										varNames		);

	// Process the shell averages
	shellAvg.process(data.getDataSet());


	// Write average quantities
	std::string shellOutputFileName = fileNameOutBase + ".bruntvaisala";
	std::ofstream shellOutput;
	shellOutput.open(shellOutputFileName.c_str());
	shellOutput<<"simTime\t"<<simTime<<std::endl;
	shellAvg.printShells(shellOutput);
	shellOutput.close();


// 	vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
// 	writer->SetFileName("richardson.vtr");
// #if VTK_MAJOR_VERSION <= 5
//   	writer->SetInputConnection(data.getDataSet()->GetProducerPort());
// #else
// 	writer->SetInputData(data.getDataSet());
// #endif
// 	writer->Write();



	std::cout<<"Finished"<<std::endl;
}



