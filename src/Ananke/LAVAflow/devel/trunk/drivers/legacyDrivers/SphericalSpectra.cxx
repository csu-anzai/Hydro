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

#include "../libsrc/Math/LAVAMath.h"

#include <string>
#include <vtkRectilinearGrid.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLRectilinearGridWriter.h>

#include <fstream>
#include <algorithm>
#include <complex>

int main(int argc, char** argv)
{

	std::string filenameIn, filenameExp, fileNameOutBase;
	int currentCS = CS_SPHERE;
	int currAxis   = XAXIS;
	int nDim = 2;
	double MSOLAR = 1.9891e33; // solar mass in grams
	double GRAVCONST = 6.67259e-8; // cm^3/g/s
	double unitFOE = 1.e51;
	double km2cm = 100000;
	int fileNum = 400;

	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		char fnameChar[256], fnameExp[256];
		
		std::sprintf(fnameChar,"%s%04d%s","/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/b163d2a3s530SG1MHWr2LBp.plt.",fileNum,".vtk");
		std::sprintf(fnameExp, "%s%04d%s","/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/output/b163d2a3s530SG1MHWr2LBp.plt.",fileNum,".expstat");
		nDim = 2;


		std::sprintf(fnameChar,"%s%04d%s","/data1/sne/HOTB/3d/b157d3a3s401SG1MHWr3LB/b157d3a3s401SG1MHWr3LB.plt.",100,".vtk");
		std::sprintf(fnameExp, "%s%04d%s","/data1/sne/HOTB/3d/b157d3a3s401SG1MHWr3LB/b157d3a3s401SG1MHWr3LB.plt.",100,".expstat");
		nDim = 3;


		filenameIn = fnameChar;
		filenameExp = fnameExp;
		fileNameOutBase = filenameIn;
	}
	else if(argc==4)
	{
		filenameIn = argv[1];
		fileNameOutBase  = argv[2];
		filenameExp  = argv[3];
	}
	else
	{
		std::cerr<<"Wrong number of arguments!"<<std::endl;
	}

	std::cout<<"Processing file \""<<filenameIn<<"\""<<std::endl;
	
	// Read in data
	Mesh data(nDim, currentCS, filenameIn);

	// Get simulation information
	double simTime = data.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);

	std::cout<<"File information:"<<std::endl;
	std::cout<<"\t             Time: "<<simTime<<std::endl;

	// Get domain and mesh information
	double 	physicalBounds[6]; // Physical space bounds
	int 		meshDimensions[3]; // Number of cells in each dimension
	data.getDataSet()->GetBounds(physicalBounds);
	data.getDataDimension(meshDimensions);

	int totalCells = 1.0; for(int i=0;i<nDim;i++){totalCells *= meshDimensions[i];}

	std::cout<<"Mesh dimensions:"<<std::endl;
	std::cout<<"\t"<<meshDimensions[0]<<"\t"<<meshDimensions[1]<<"\t"<<meshDimensions[2]<<std::endl;
	std::cout<<"\tTotal cells: "<<totalCells<<std::endl;



	// double 	trueVolume = 4.0/3.0*vtkMath::Pi()*(	physicalBounds[1]*physicalBounds[1]*physicalBounds[1] - 
	// 											physicalBounds[0]*physicalBounds[0]*physicalBounds[0] );

	// double 	simVolume = 1.0;
	// simVolume *= (1.0/3.0)*physicalBounds[1]*physicalBounds[1]*physicalBounds[1] - physicalBounds[0]*physicalBounds[0]*physicalBounds[0];
	// simVolume *= -cos(physicalBounds[3]) + cos(physicalBounds[2]);
	// if(nDim>2)
	// {
	// 	simVolume *= physicalBounds[5] - physicalBounds[4];
	// }
	// else
	// {
	// 	simVolume *= 2.0*vtkMath::Pi();
	// }

	// double volumeCorrection = trueVolume/simVolume;
	double volumeCorrection = 2.0/(-cos(physicalBounds[3]) + cos(physicalBounds[2]));
	std::cout<<"Ratio of true volume to simulated volume: "<<volumeCorrection<<std::endl;

	// data.computeCellVolumes("cellVolumeTmp");
	// data.storeData(volumeCorrection*data["cellVolumeTmp"],"cellVolume");

	// Create data arrays of coordinates (full 2d/3d arrays of cell centers)
	data.createCoordinateDataArrays("r","phi","theta");

	// Get the radial coordinates
	// This produces the left(i) and right(i+1) coordinates for cell i
	vtkSmartPointer<vtkDataArray> radCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();


	/*====================================================================
		
		Create the azimuthal kinetic energy
		integral kernel for each cell
		(c.f. Hanke, Marek, Muller, Janka 2012, Eqn 8)

		This will later be weighted by the spherical harmonic to
		compute the energy

	=======================================================================*/
	data.addVariable("aziEKernel");

	for(int i=0; i<data["aziEKernel"]->GetNumberOfTuples(); i++)
	{
		double value = 1.0;
			

		value *= sqrt(data["dens"]->GetTuple1(i));
		value *= data["vely"]->GetTuple1(i);
		// value *= data["cellVolume"]->GetTuple1(i);

		double cellEdges[6];
		data.getDataSet()->GetCellBounds(i,cellEdges);

		value *= sin(0.5*(cellEdges[3]+cellEdges[2]))*(cellEdges[3]-cellEdges[2]);
		if(nDim>2)
		{
			value *= cellEdges[5] - cellEdges[4];
		}
		else
		{
			value *= 2.0*vtkMath::Pi();
		}


		data["aziEKernel"]->SetTuple1(i, value);
	}


	/*====================================================================
		
		Read in the maximum gain radius and the minimum shock radius. 
		This information will be used to determine the proper range 
		of radial coordinates to average over.

	=======================================================================*/
	Util::ScalarList expStat;
	expStat.read(filenameExp);
	double gainRadiusMax  = expStat.getScalar("grMax");
	double shockRadiusMin = expStat.getScalar("shockMin");
	double averagingFraction = 0.75;
	double usableWidth = averagingFraction*(shockRadiusMin-gainRadiusMax);
	double averagingWidth = 30*km2cm;

	if(averagingWidth>usableWidth)
	{
		averagingWidth = usableWidth;
	}

	double averagingMidPt = gainRadiusMax + 0.5*(shockRadiusMin-gainRadiusMax);
	double averagingLower = averagingMidPt - 0.5*averagingWidth;
	double averagingUpper = averagingMidPt + 0.5*averagingWidth;

	std::cout<<"Minimum shock radius [km]: "<<shockRadiusMin/km2cm<<std::endl;
	std::cout<<"Maximum gain  radius [km]: "<<gainRadiusMax/km2cm<<std::endl;
	std::cout<<"Usable Width         [km]: "<<usableWidth/km2cm<<std::endl;
	std::cout<<"Averaging width      [km]: "<<averagingWidth/km2cm<<std::endl;
	std::cout<<"Averaging midpoint   [km]: "<<averagingMidPt/km2cm<<std::endl;
	std::cout<<"Averaging lower      [km]: "<<averagingLower/km2cm<<std::endl;
	std::cout<<"Averaging upper      [km]: "<<averagingUpper/km2cm<<std::endl;

	/*====================================================================
		
		Determine the start and end radial indices that lie in the averaging region

	=======================================================================*/
	bool foundFirst = false, foundLast = false;;
	int startInd = -1, finalInd = -1;

	std::cout<<"Number of tuples in radCoords: "<<radCoords->GetNumberOfTuples()<<std::endl;

	for(int i=0; i<=meshDimensions[0]; i++)
	{
		double rL = radCoords->GetTuple1(i);
		double rR = radCoords->GetTuple1(i+1); 

		// Find starting index
		if( (rR >= averagingLower) && !foundFirst )
		{
			foundFirst = true;
			startInd = i;
		}
		else if( (rR >= averagingUpper) && foundFirst)
		{
			finalInd = i;
			break;
		}
	}

	int radialIndexStart = startInd;
	int nRadialSamples = finalInd-startInd;

	double** binEdges;
	binEdges = new double*[nRadialSamples];
	for(int i=0; i<nRadialSamples; i++)
	{
		binEdges[i] = new double[2];
		binEdges[i][LOW] = radCoords->GetTuple1(i+startInd);
		binEdges[i][HIGH] = radCoords->GetTuple1(i+1+startInd);
	}
	binEdges[0][LOW] = averagingLower;
	binEdges[nRadialSamples-1][HIGH] = averagingUpper;


	std::cout<<"Start index: "<<startInd<<std::endl;
	std::cout<<"Final index: "<<finalInd<<std::endl;


	/*====================================================================
		
		Compute the spectral energy Espec(l) of the azimuthal kinetic energy as a 
		function of multipole order l

		(1)	For a fixed radial index, loop over theta indices and 
				obtain the rays along phi
		(2)	Choose a range of multipole orders, 0<=l<=L
		(3)	For each order, loop over the degree (m=-l:l)
		(4)	Sum the contributions from aziEKernel weighted with the 
				spherical harmonic at that point
		(5)	Square the sum and store it in Espec(l)

	=======================================================================*/

	int LMin = 0;
	int LMax = 100;
	double* Espec = new double[LMax+1];
	double* EspecSingleBin = new double[LMax+1];

	for(int i=0; i<=LMax; i++){ Espec[i] = 0.0; EspecSingleBin[i] = 0.0; }


	double* radialPositions = new double[nRadialSamples];
	

	int startingIndex[3];
	
	for(int s=0; s<nRadialSamples; s++)
	{
		int radialIndex = radialIndexStart + s;
		double rLeft = radCoords->GetTuple1(radialIndex);
		double rRight = radCoords->GetTuple1(radialIndex+1);
		radialPositions[s] = pow(0.5*(rLeft*rLeft*rLeft + rRight*rRight*rRight),1.0/3.0);
		std::cout<<"Evaluation radius [km]: "<<radialPositions[s]/km2cm<<std::endl;

		

		// Move over multipole order
		for(int l=LMin; l<=LMax; l++)
		{
			std::cout<<"l = "<<l<<std::endl;
			EspecSingleBin[l] = 0.0;

			// In 2d, do not consider degrees larger than 0. This is because the only symmetric mode is the m=0 mode.
			// Inclusion of other modes leads to incorrect answers 
			int mMin 	= 0;
			int mMax 	= 0;

			if(nDim>2)
			{
				mMin = -l;
				mMax = l;
			}

			// Loop over degree
			for(int m=mMin; m<=mMax; m++)
			{

				// Perform the integration of the spherical harmonic-weighted function (the actual integral in mathy terms)
				// aziEKernel already contains the sin(phi)*dPhi*dTheta term
				std::complex<double> sum(0.,0.);
				
				for(int thetaInd = 0; thetaInd < (nDim>2 ? meshDimensions[2] : 1); ++thetaInd)
				{


					startingIndex[XAXIS]	= radialIndex;
					startingIndex[YAXIS]	= 0;
					startingIndex[ZAXIS]	= thetaInd;
					std::vector<int> cellIndices = data.getCellsAlongRay(YAXIS, startingIndex);

					// Loop over cells to compute the sum
					for (std::vector<int>::iterator cellIt = cellIndices.begin(); cellIt!= cellIndices.end(); ++cellIt)
					{
						// Get the spherical harmonic
						std::complex<double> sh = Math::sphericalHarmonic(	l,
																			m,
																			data["theta"]->GetTuple1(*cellIt),
																			data["phi"]->GetTuple1(*cellIt));
						
						// Take the conjugate
						sh = std::conj(sh);

						// Add to multipole sum
						std::complex<double> value(data["aziEKernel"]->GetTuple1(*cellIt), 0.);
						sum += sh*value;

					} // cells
				} // theta


				// Perform the norm and square operations
				double energyValue = std::real(sum)*std::real(sum) + std::imag(sum)*std::imag(sum);

				// Store in temporary Espec
				EspecSingleBin[l] += energyValue;

			} // m
		} // l

		// Accumulate weighted sum for average
		for(int l=LMin; l<=LMax; l++)
		{
			Espec[l] += EspecSingleBin[l]*(binEdges[s][HIGH]-binEdges[s][LOW]);
		}
		std::cout<<"binEdges["<<s<<"][HIGH]-binEdges["<<s<<"][LOW]: "<<(binEdges[s][HIGH]-binEdges[s][LOW])/km2cm<<std::endl;
	}

	// Finish computing the average
	for(int l=LMin; l<=LMax; l++)
	{
		Espec[l] /= binEdges[nRadialSamples-1][HIGH]-binEdges[0][LOW];
	}

	// Print Espec

	char foutChar[256];
	std::sprintf(foutChar,"%s%s",fileNameOutBase.c_str(),".spectra");

	std::ofstream output;
	output.open(foutChar);

	output<<"simTime\t"<<simTime<<std::endl;
	output<<"nRadialSamples\t"<<nRadialSamples<<std::endl;
	for(int s=0; s<nRadialSamples; s++)
	{
		output<<radialPositions[s]<<"\t";
	}
	output<<std::endl;

	output<<"L\tE"<<std::endl;
	std::cout<<"L\tE"<<std::endl;
	for(int l=LMin; l<=LMax; l++)
	{
		std::cout<<l<<"\t"<<Espec[l]<<std::endl;
		output<<l<<"\t"<<Espec[l]<<std::endl;
	}
	output.close();



	// Cleanup
	delete [] Espec;
	delete [] EspecSingleBin;
	delete [] radialPositions;
	delete [] binEdges;


// 	vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
// 	writer->SetFileName("sh_debugging.vtr");
// #if VTK_MAJOR_VERSION <= 5
//   	writer->SetInputConnection(data.getDataSet()->GetProducerPort());
// #else
// 	writer->SetInputData(data.getDataSet());
// #endif
// 	writer->Write();


	std::cout<<"Finished"<<std::endl;
}



