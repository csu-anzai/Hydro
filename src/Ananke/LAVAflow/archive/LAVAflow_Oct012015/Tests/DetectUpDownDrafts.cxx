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
		filenameIn = "/home/tah09e/data/sne/data/HOTB/2d/b155d2a3s123SG1MHwr1LB/b155d2a3s123SG1MHwr1LB.plt.0006.vtk";
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
	double pointMass = data.getDataSet()->GetFieldData()->GetArray("PMASS")->GetTuple1(0);
	double pointGrav = data.getDataSet()->GetFieldData()->GetArray("PMGRV")->GetTuple1(0);

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




	/*====================================================================
		
		Compute the angular distribution of mass weighted velocity
		integrated along radial rays in the gain region. 
		(1)	Get all the cells along a radial ray
		(1.1)	Moving outward from the star to the shock, compute the 
				sum of radial momentum and normalize this by the total 
				volume in the ray.
	=======================================================================*/

	double radMin = 0;

	int nPhi = meshDimensions[1];
	int nTheta = nDim>2 ? meshDimensions[2] : 1;

	double** radialMomentumIntegrals;
	radialMomentumIntegrals = new double*[nPhi];
	for(int i=0; i<nPhi; i++)
	{
		radialMomentumIntegrals[i] = new double[nTheta]; 
		memset(radialMomentumIntegrals[i], 0.0, nTheta);
	}

	data.storeData(data["velx"]*data["dens"]*data["cellVolume"],"radMom");
	std::string variableForAnalysis = "radMom";

	int startingIndex[3];

	for(int p = 0; p<nPhi; ++p)
	{
		for(int t = 0; t<nTheta; ++t)
		{
			startingIndex[XAXIS]	= 0;
			startingIndex[YAXIS]	= p;
			startingIndex[ZAXIS]	= t;
			std::vector<int> cellIndices = data.getCellsAlongRay(XAXIS, startingIndex);

			double raySum = 0.0, rayVolume = 0.0;
			bool foundFirst = false;
			bool foundLast  = false;
			bool useCell     = false;
			int nCellsUsed = 0;

			for (std::vector<int>::reverse_iterator i = cellIndices.rbegin(); i != cellIndices.rend(); ++i)
			{
				
				bool isShocked = data["shock"]->GetTuple1(*i) > 0.5; // The data is either 0 or 1, so you need to pick something reasonable

				// If the cell is shocked
				if(isShocked)
				{
					if(!foundFirst)
					{
						foundFirst = true;
					}

				}
				else
				{
					// If we've already found the first element, then not being shocked means we left the shock
					// We should start using all cells from here on out
					if(foundFirst && !foundLast)
					{
						foundLast = true;
						useCell  = true;
					}

				}

				if(data["qenr"]->GetTuple1(*i) < 0.0)
				{
					break;
				}

				if(data["r"]->GetTuple1(*i)<radMin)
				{
					break;
				}


				if(useCell)
				{
					raySum += data[variableForAnalysis.c_str()]->GetTuple1(*i);
					rayVolume += data["dens"]->GetTuple1(*i)*data["cellVolume"]->GetTuple1(*i); 
					nCellsUsed++;
				}

			}

			// Store results
			radialMomentumIntegrals[p][t] = raySum/rayVolume;
			// radialMomentumIntegrals[p][t] = raySum;

			if(nCellsUsed == 0)
			{
				std::cout<<"No cells used for a ray!"<<std::endl;
			}

		}
	}




	/*====================================================================
		
		Determine the number and fractional size (fractional solid angle) of up and down drafts.
		
		In 2D, this consists of finding the zeros, and distance between zeros,
		of radialMomentumIntegrals as a function of phi

		In 3D, the problem is expanded to finding the area (in phi-theta space)
		enclosed by zero level sets

		2D algorithm outline
		(1)	Starting at p=0, determine if the first ray is a down draft or up draft
			(isDownDraft). Set phiInterOld = phi[p]
		(2)	Increase p until a switch in sign occurs. Flag this by switching the state
			of isDownDraft. Interpolate to find the zero-intersect. Store this in 
			phiInterNew
		(3)	Store (phiInterOld-phiInterNew) in the appropriate solidAngle vector
		(4)	Set phiInterOld=phiInterNew
		(5)	Continue until all rays have been exhausted. When we reach the last ray
			use phiInterNew = phi[end]
	

	=======================================================================*/

	// Get coordinate arrays

	vtkSmartPointer<vtkDataArray> phiCoords	= vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetYCoordinates();
	vtkSmartPointer<vtkDataArray> thetaCoords	= vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetZCoordinates();

	std::vector<double>	solidAngleDownDraft;
	std::vector<double>	solidAngleUpDraft;
	std::vector<double>	radmomDownDraft;
	std::vector<double>	radmomUpDraft;
	if(nDim==2)
	{


		int phiInd = 0;
		int thetaInd = 0;
		bool isDownDraft = false;
		double valuePrev = 0.0, valueCurr = 0.0;
		double intersectionOld = 0.0, intersectionNew = 0.0;
		double coordPrev = 0.0, coordCurr = 0.0;
		
		double saRay = 0.0;
		double radMomSum = 0.0;
		int radMomRaysUsed = 0;


		// Determine whether the first point is a down draft
		if(radialMomentumIntegrals[phiInd][thetaInd] <= 0.0)
		{
			isDownDraft = true;
		}
		else
		{
			isDownDraft = false;
		}

		// std::cout<<"Starting with a down draft? "<<isDownDraft<<std::endl;

		valuePrev = radialMomentumIntegrals[phiInd][thetaInd];
		intersectionOld = phiCoords->GetTuple1(phiInd);
		coordPrev = 0.5*(phiCoords->GetTuple1(0)+phiCoords->GetTuple1(1));

		for(phiInd=1; phiInd<nPhi; phiInd++)
		{	
			// Get the current value
			valueCurr = radialMomentumIntegrals[phiInd][thetaInd];

			// Get current coordinate
			coordCurr = 0.5*(phiCoords->GetTuple1(phiInd)+phiCoords->GetTuple1(phiInd+1));

			// Determine if there is a sign change
			bool signChanged = sgn(valuePrev) != sgn(valueCurr);

			// If there is a sign change, perform the procedure outlined above
			// Otherwise, keep going
			if(signChanged)
			{


				// Determine the intersection point via linear interpolation
				double slope = (valueCurr-valuePrev)/(coordCurr-coordPrev);

				intersectionNew = -1.0*valuePrev/slope + coordPrev;

				// Compute the solid angle spanned by the structure
				double sa = -cos(intersectionNew)+cos(intersectionOld);



				// Handle the average momentum for this draft
				double prevSum = radMomSum;
				double currSum = 0.0;
				double prevL = phiCoords->GetTuple1(phiInd-1);
				double prevR = phiCoords->GetTuple1(phiInd);
				double currL = prevR;
				double currR = phiCoords->GetTuple1(phiInd+1);

				if(intersectionNew < prevR)
				{
					// handle previous sum
					double phiBndPrev = std::max(intersectionOld, prevL);
					saRay = -cos(intersectionNew) + cos(phiBndPrev);
					prevSum += saRay*valuePrev;


					// handle current sum
					saRay = -cos(currL) + cos(intersectionNew);
					currSum += saRay*valuePrev;

				}
				else
				{
					// handle previous sum
					// nothing to handle with current sum. It will be taken care of either in the 
					// else statement or in the next intersection test
					double phiBndPrev = std::max(intersectionOld, prevL);
					saRay = -cos(prevR) + cos(phiBndPrev);
					prevSum += saRay*valuePrev;

					saRay = -cos(intersectionNew) + cos(prevR);
					prevSum += saRay*valueCurr;

				}

				// Set radMomSum equal to currSum
				radMomSum = currSum;




				// Store the solid angle and radial momentum in the proper vector
				if(isDownDraft)
				{
					solidAngleDownDraft.push_back(sa);
					radmomDownDraft.push_back(prevSum/sa);
				}
				else
				{
					solidAngleUpDraft.push_back(sa);
					radmomUpDraft.push_back(prevSum/sa);
				}


				// Set the proper old value for the intersection
				intersectionOld = intersectionNew;

				// Flip the isDownDraft flag
				isDownDraft = !isDownDraft;
			}
			else
			{
				// Add the previous ray's contribution to the solid angle weighted average of momentum
				double phiBndPrev = std::max(intersectionOld, phiCoords->GetTuple1(phiInd-1));
				saRay = -cos(phiCoords->GetTuple1(phiInd)) + cos(phiBndPrev);
				radMomSum += saRay*valuePrev;
			}

			// Update the "old" values
			valuePrev 	= valueCurr;
			coordPrev 	= coordCurr;
		}


		// Handle the last structure, since we're going to hit the end 
		// of the data and no more sign changes will occur
		double saFinal = -cos(phiCoords->GetTuple1(nPhi)) + cos(intersectionOld);
		radMomSum += (-cos(phiCoords->GetTuple1(nPhi)) + cos(std::max(intersectionOld, phiCoords->GetTuple1(nPhi-1))))*valueCurr;
		if(isDownDraft)
		{
			solidAngleDownDraft.push_back(saFinal);
			radmomDownDraft.push_back(radMomSum/saFinal);
		}
		else
		{
			solidAngleUpDraft.push_back(saFinal);
			radmomUpDraft.push_back(radMomSum/saFinal);
		}


		// Normalize the solid angles by the solid angle spanned in 2d
		double saDomain = -cos(phiCoords->GetTuple1(nPhi)) + cos(phiCoords->GetTuple1(0));
		for(int i=0; i<solidAngleDownDraft.size(); i++)
		{
			solidAngleDownDraft[i] /= saDomain;
		}
		for(int i=0; i<solidAngleUpDraft.size(); i++)
		{
			solidAngleUpDraft[i] /= saDomain;
		}


		double saSum = 0.0;

		// std::cout<<"Down drafts"<<std::endl;
		// for(int i=0; i<solidAngleDownDraft.size(); i++)
		// {
		// 	std::cout<<"\t"<<solidAngleDownDraft[i]<<"\t"<<radmomDownDraft[i]<<std::endl;
		// 	saSum += solidAngleDownDraft[i];
		// }

		// std::cout<<"Up drafts"<<std::endl;
		// for(int i=0; i<solidAngleUpDraft.size(); i++)
		// {
		// 	std::cout<<"\t"<<solidAngleUpDraft[i]<<"\t"<<radmomUpDraft[i]<<std::endl;
		// 	saSum += solidAngleUpDraft[i];
		// }

		// std::cout<<"Total SA: "<<saSum<<std::endl;

		std::cout<<" Number of down drafts: "<<solidAngleDownDraft.size()<<std::endl;
		std::cout<<" Number of up drafts:   "<<solidAngleUpDraft.size()<<std::endl;


	}
	else if(nDim>2)
	{
		std::cout<<"[ERROR] Dimension "<<nDim<<" does not support calculating solid angle fractions! Please add this functionality"<<std::endl;
	}
	else
	{
		std::cout<<"[ERROR] Dimension "<<nDim<<" does not support calculating solid angle fractions!"<<std::endl;
	}




	/*====================================================================
		
		Write phi-theta ray data to file
		(1)	Header consists of
			a) Simulation time
			b) ...
		(2)	Line after header consists of theta coordinates
		(3) Subsequent lines consist of [phi Coordinate|ray0|ray1|...|rayN]
	=======================================================================*/
	std::ofstream output;
	output.open((fileNameOutBase + ".rayRadialMomenta").c_str());


	// Set flags
	output.precision(15);
	output.setf(std::ios::scientific,std::ios::floatfield);


	// Write header
	output<<"simtime\t"<<simTime<<std::endl;



	// Write theta coordinates
	// if(nTheta==1)
	// {
	// 	output<<thetaCoords->GetTuple1(0)<<std::endl;
	// }
	// else
	// {
	// 	for(int t = 0; t<nTheta; ++t)
	// 	{
	// 		double thetaCoord = 0.5*(thetaCoords->GetTuple1(t) + thetaCoords->GetTuple1(t+1));
	// 		output<<thetaCoord<<"\t";
	// 	}
	// 	output<<std::endl;
	// }


	for(int p = 0; p<nPhi; ++p)
	{
		double phiCoord = 0.5*(phiCoords->GetTuple1(p) + phiCoords->GetTuple1(p+1));

		// Write phiCoord at start of row
		output<<phiCoord;

		for(int t = 0; t<nTheta; ++t)
		{
			output<<"\t"<<radialMomentumIntegrals[p][t];
		}

		output<<std::endl;
	}


	output.close();


	/*====================================================================
		
		Write down/up draft information to file
		(1)	Header consists of 
			a) Simulation time
			b) Number of down drafts (Nd)
			c) Number of up drafts (Nu)
			with each value on a separate line
		(2)	Next Nd lines consist of [solid angle, radMom] for down drafts
		(2)	Next Nu lines consist of [solid angle, radMom] for up drafts
	=======================================================================*/
	std::ofstream draftOutput;
	draftOutput.open((fileNameOutBase + ".draftSolidAngles").c_str());


	// Set flags
	draftOutput.precision(15);
	draftOutput.setf(std::ios::scientific,std::ios::floatfield);


	// Write header
	draftOutput<<"simtime\t"<<simTime<<std::endl;
	draftOutput<<"Ndowndrafts\t"<<solidAngleDownDraft.size()<<std::endl;
	draftOutput<<"Nupdrafts\t"<<solidAngleUpDraft.size()<<std::endl;

	// write down draft data
	for(int i=0; i<solidAngleDownDraft.size(); i++)
	{
		draftOutput<<solidAngleDownDraft[i]<<"\t"<<radmomDownDraft[i]<<std::endl;
	}


	// write up draft data
	for(int i=0; i<solidAngleUpDraft.size(); i++)
	{
		draftOutput<<solidAngleUpDraft[i]<<"\t"<<radmomUpDraft[i]<<std::endl;
	}

	draftOutput.close();









	// vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
	// writer->SetFileName("draft_debug.vtr");
	// writer->SetInputData(data.getDataSet());
	// writer->Write();

	std::cout<<"Finished"<<std::endl;

}