#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Operators/SurfaceAreaOperator.h"
#include "../libsrc/Operators/shellAveraging/ShellAverageOperator.h"
#include "../libsrc/Operators/shellAveraging/planeShells/ShellAveragePlaneOperator.h"
#include "../libsrc/Operators/surfaceAveraging/planarSurface/SurfaceAveragePlaneOperator.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/includes/LAVAconstants.h"

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

	std::string filenameIn, fileNameOutBase;
	int currentCS = CS_SPHERE;
	int currAxis   = XAXIS;
	int nDim = 2;
	double MSOLAR = 1.9891e33; // solar mass in grams
	double GRAVCONST = 6.67259e-8; // cm^3/g/s
	double unitFOE = 1.e51;
	double km2cm = 100000;
	int fileNum = 200;
	bool writeOutputFile = false;

	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		// filenameIn = "../Tests/files/averaging/b155d2a3s123SG1MHwr2LB.plt.0010.vtk";
		char fnameChar[256];
		// std::sprintf(fnameChar,"%s%04d%s","/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/b163d2a3s530SG1MHWr2LBp.plt.",fileNum,".vtk");
		std::sprintf(fnameChar,"%s%04d%s","/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/b163d2a3s530SG1MHWr2LBp.plt.",fileNum,".vtk");
		// filenameIn = "/data1/sne/HOTB/2d/b163d2a3s231SG1MHWr2LB/b163d2a3s231SG1MHWr2LB.plt.0040.vtk";
		filenameIn = fnameChar;
		fileNameOutBase = filenameIn;
		writeOutputFile = false;
	}
	else if(argc==3)
	{
		filenameIn = argv[1];
		fileNameOutBase  = argv[2];
		std::cout<<"Processing file \""<<filenameIn<<"\""<<std::endl;
		writeOutputFile = true;
	}

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

	std::cout<<"Physical dimensions:"<<std::endl;
	std::cout<<"\t"<<physicalBounds[0]<<"\t"<<physicalBounds[1]<<std::endl;
	std::cout<<"\t"<<physicalBounds[2]<<"\t"<<physicalBounds[3]<<std::endl;
	std::cout<<"\t"<<physicalBounds[4]<<"\t"<<physicalBounds[5]<<std::endl;
	std::cout<<"Mesh dimensions:"<<std::endl;
	std::cout<<"\t"<<meshDimensions[0]<<"\t"<<meshDimensions[1]<<"\t"<<meshDimensions[2]<<std::endl;
	std::cout<<"\tTotal cells: "<<totalCells<<std::endl;

	// Create data arrays of coordinates (full 2d/3d arrays of cell centers)
	data.createCoordinateDataArrays("r","phi","theta");

	// Get the radial coordinates
	vtkSmartPointer<vtkDataArray> radEdges   = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();
	vtkSmartPointer<vtkDataArray> phiEdges   = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetYCoordinates();
	vtkSmartPointer<vtkDataArray> thetaEdges = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetZCoordinates();


	/*====================================================================

		Determine the shock radius at a given (phi,theta) coordinate
		We'll store this radius in an nPhi x nTheta array
		Also compute the phi and theta coordinates

	=======================================================================*/

	double dPhiConst = phiEdges->GetTuple1(2)-phiEdges->GetTuple1(1);
	int nStart = floor(std::max(0.0,(physicalBounds[2]-0.0)/dPhiConst));//+1;
	int nEnd   = floor(std::max(0.0,(vtkMath::Pi()-physicalBounds[3])/dPhiConst));//+1;

	if(nStart>0) nStart++;
	if(nEnd>0) nEnd++;

	std::cout<<"Number of cells to duplicate at start: "<<nStart<<std::endl;
	std::cout<<"Number of cells to duplicate at end:   "<<nEnd<<std::endl;


	int nPhi = meshDimensions[1];
	int nPhiTotal = nPhi + nStart + nEnd;
	int nTheta = (nDim>2) ? meshDimensions[2] : 1;
	double* phiCoord     = new double[nPhiTotal];
	double* thetaCoord   = new double[nTheta];
	double* deltaPhi     = new double[nPhiTotal];
	double* deltaTheta   = new double[nTheta];
	double** shockRadius = new double*[nPhiTotal];
	for(int i=0; i<nPhiTotal; i++)
	{
		shockRadius[i] = new double[nTheta];
	}


	double phiStart = 0.0;
	double phiEnd   = vtkMath::Pi();
	double dPhi = (phiEnd-phiStart)/double(nPhi);
	for(int i=0; i<nPhi; i++)
	{
		phiCoord[i+nStart] = 0.5*(phiEdges->GetTuple1(i+1)+phiEdges->GetTuple1(i));
		deltaPhi[i+nStart] = phiEdges->GetTuple1(i+1)-phiEdges->GetTuple1(i);
		// std::cout<<"phiCoord["<<i<<"]: "<<phiCoord[i]<<std::endl;
	}

	for(int i=0; i<nStart; i++)
	{
		phiCoord[i] = phiCoord[nStart] - (nStart-i)*dPhiConst;
		deltaPhi[i] = dPhiConst;
	}

	for(int i=0; i<nEnd; i++)
	{
		int ind = i + nPhi + nStart;
		phiCoord[ind] = phiCoord[nPhi+nStart-1] + i*dPhiConst;
		deltaPhi[ind] = dPhiConst;
	}

	if(nDim>2)
	{
		for(int i=0; i<nTheta; i++)
		{
			thetaCoord[i] = 0.5*(thetaEdges->GetTuple1(i+1)+thetaEdges->GetTuple1(i));
			deltaTheta[i] = thetaEdges->GetTuple1(i+1)-thetaEdges->GetTuple1(i);
		}
	}
	else
	{
		thetaCoord[0] = 0.0;
		deltaTheta[0] = 2.0*vtkMath::Pi();
	}

	vtkSmartPointer<vtkDataArray> shock = data["shock"];
	vtkSmartPointer<vtkDataArray> radius = data["r"];
	std::vector<double> shockPositions;
	std::vector<int> cellIndicesOutsideShock;
	int startingIndex[3];


	for(int p = 0; p<nPhi; ++p)
	{
		for(int t = 0; t<nTheta; ++t)
		{

			// std::cout<<"p = "<<p<<"\tt= "<<t<<std::endl;
			
			startingIndex[XAXIS]	= 0;
			startingIndex[YAXIS]	= p;
			startingIndex[ZAXIS]	= t;
			std::vector<int> cellIndices = data.getCellsAlongRay(XAXIS, startingIndex);

			bool foundFirst = false;
			double shockPosAvg = 0.0;
			int nShocked = 0;

			for (std::vector<int>::reverse_iterator i = cellIndices.rbegin(); i != cellIndices.rend(); ++i)
			{
				bool isShocked = shock->GetTuple1(*i) > 0.5; // The data is either 0 or 1, so you need to pick something reasonable

				// If it's shocked, add the radial coordinate to the sum
				if(isShocked)
				{
					if(!foundFirst)
					{
						foundFirst = true;
					}

					// Increase the sum
					shockPosAvg += radius->GetTuple1(*i);

					// Increase the number of cells to use in the average
					nShocked++;
				}
				else
				{
					// If we've already found the first element, then not being shocked means we left the shock
					// Therefore, we should break the loop
					if(foundFirst)
					{
						break;
					}

				}

			}

			// Finish computing the average by dividing by the number of cells at the intersection of the shock
			// and the ray
			if(nShocked > 0)
			{
				shockPosAvg /= double(nShocked);
			}

			// Store the shock position
			shockRadius[p+nStart][t] = shockPosAvg;

		}
	}


	// Reflect interior shock positions to fill in axis gap
	for(int t=0; t<nTheta; t++)
	{
		for(int p=0; p<nStart; p++)
		{
			shockRadius[nStart-p-1][t] = shockRadius[nStart+p][t];
		}	


		for(int p=0; p<nEnd; p++)
		{
			shockRadius[nStart+nPhi+p][t] = shockRadius[nStart+nPhi-p-1][t];
		}	
	}

	// for(int i=0; i<nPhiTotal; i++)
	// {
	// 	std::cout<<i<<"\t"<<phiCoord[i]<<"\t"<<deltaPhi[i]<<"\t"<<shockRadius[i][0]<<std::endl;
	// }


	// data.addVariable("aziEKernel");

	// for(int i=0; i<data["aziEKernel"]->GetNumberOfTuples(); i++)
	// {
	// 	double value = 1.0;
			

	// 	value *= sqrt(data["dens"]->GetTuple1(i));
	// 	value *= data["vely"]->GetTuple1(i);
	// 	// value *= data["cellVolume"]->GetTuple1(i);

	// 	double cellEdges[6];
	// 	data.getDataSet()->GetCellBounds(i,cellEdges);

	// 	value *= -cos(physicalBounds[3]) + cos(physicalBounds[2]);
	// 	if(nDim>2)
	// 	{
	// 		value *= physicalBounds[5] - physicalBounds[4];
	// 	}
	// 	else
	// 	{
	// 		value *= 2.0*vtkMath::Pi();
	// 	}


	// 	data["aziEKernel"]->SetTuple1(i, value);

	// }

	/*====================================================================
		
		Compute the spherical harmonic decomposition

	=======================================================================*/

	int LMin = 0;
	int LMax = 3;
	double* A = new double[LMax+1];

	for(int i=0; i<=LMax; i++){ A[i] = 0.0; }


	// Move over multipole order
	for(int l=LMin; l<=LMax; l++)
	{	
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
			std::complex<double> sum(0.,0.);	
			
			for(int p = 0; p<nPhiTotal; ++p)
			{
				for(int t = 0; t<nTheta; ++t)
				{

					/*-----------------------

						Compute the integral for this mode

						Use the trapezoidal rule for the phi-direction
					
					---------------------- */	
					// double phiLeft  = phiCoord[p]-0.5*deltaPhi[p];
					// double phiRight = phiCoord[p]+0.5*deltaPhi[p];

					// // Get the spherical harmonics at the edges
					// std::complex<double> shLeft = Math::sphericalHarmonic(	l,
					// 							 							m,
					// 														thetaCoord[t],
					// 														phiLeft);
					// std::complex<double> shRight = Math::sphericalHarmonic(	l,
					// 							 							m,
					// 														thetaCoord[t],
					// 														phiRight);
					
					// // Take the conjugate
					// shLeft  = std::conj(shLeft);
					// shRight = std::conj(shRight);

					// // Add to multipole sum
					// sum += 0.5*shockRadius[p][t]*(std::sin(phiLeft)*shLeft + std::sin(phiRight)*shRight)*deltaPhi[p]*deltaTheta[t];
	


					/*-----------------------

						Compute the integral for this mode

						Use the midpt rule. No appreciable gain by using the trapezoidal rule
					
					---------------------- */	
					std::complex<double> sh = Math::sphericalHarmonic(	l,
												 						m,
																		thetaCoord[t],
																		phiCoord[p]);
					
					// Take the conjugate
					sh  = std::conj(sh);

					// Add to multipole sum
					sum += shockRadius[p][t]*sh*sin(phiCoord[p])*deltaPhi[p]*deltaTheta[t];

				}
			}
			
			// Square the a_lm sum
			double energy = std::real(sum)*std::real(sum) + std::imag(sum)*std::imag(sum);
			
			// Update order coefficient
			A[l] += energy;
		
		}
	}

	// Take the square of the A[l] summation (over a_lm)
	for(int l=0; l<=LMax; l++)
	{
		A[l] = sqrt(A[l]);
		std::cout<<"A["<<l<<"]: "<<A[l]<<std::endl;
	}


	if(writeOutputFile)
	{
		char foutChar[256];
		std::sprintf(foutChar,"%s%s",fileNameOutBase.c_str(),".shockdecomp");
		std::ofstream output;
		output.open(foutChar);
		output.setf(std::ios::scientific,std::ios::floatfield);
		output.precision(15);
		output<<"simTime\t"<<simTime<<std::endl;
		// output<<"L\treal(A_L)\timag(A_L)"<<std::endl;
		// for(int l=LMin; l<=LMax; l++)
		// {
		// 	output<<l<<"\t"<<std::real(A[l])<<"\t"<<std::imag(A[l])<<std::endl;
		// }
		output<<"L\tA_L"<<std::endl;
		for(int l=LMin; l<=LMax; l++)
		{
			output<<l<<"\t"<<A[l]<<std::endl;
		}
		output.close();
	}

	// Cleanup
	delete [] A;
	delete [] phiCoord;
	delete [] thetaCoord;
	delete [] deltaPhi;
	delete [] deltaTheta;
	delete [] shockRadius;

	// vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
	// writer->SetFileName("sh_debugging.vtr");
	// writer->SetInputData(data.getDataSet());
	// writer->Write();


	std::cout<<"Finished"<<std::endl;
}



