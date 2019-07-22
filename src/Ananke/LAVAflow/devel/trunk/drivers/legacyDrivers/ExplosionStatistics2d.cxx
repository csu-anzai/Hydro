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
	int nDim = 3;
	double MSOLAR = 1.9891e33; // solar mass in grams
	double GRAVCONST = 6.67259e-8; // cm^3/g/s
	double unitFOE = 1.e51;
	double explosionCriterion = 1.e48;
	bool writeOutput = false;

	// Create scalar list utility to store data we want to output
	Util::ScalarList dataList;




	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		// filenameIn = "../Tests/files/averaging/b155d2a3s123SG1MHwr2LB.plt.0010.vtk";

		filenameIn = "/data1/sne/HOTB/2d/b163d2a3s231SG1MHWr2LB/b163d2a3s231SG1MHWr2LB.plt.0092.vtk";
		scrFileName = "/data1/sne/HOTB/2d/b163d2a3s231SG1MHWr2LB/scr_b163d2a3s231SG1MHWr2LB";

		nDim = 3;
		filenameIn = "/data1/sne/HOTB/3d/b157d3a3s401SG1MHWr3LB/b157d3a3s401SG1MHWr3LB.plt.0100.vtk";
		scrFileName = "/data1/sne/HOTB/3d/b157d3a3s401SG1MHWr3LB/scr_b157d3a3s401SG1MHWr3LB";

		// nDim = 3;
		// filenameIn = "/data1/sne/HOTB/3d/b157d3a3s503SG1MHWr3LB/b157d3a3s503SG1MHWr3LB.plt.0092.vtk";
		// scrFileName = "/data1/sne/HOTB/3d/b157d3a3s503SG1MHWr3LB/scr_b157d3a3s503SG1MHWr3LB";

		// fileNameOutBase = filenameIn;
		fileNameOutBase = "expstats2dtest";

		writeOutput = true;
	}
	else if(argc==4)
	{	
		filenameIn = argv[1];
		fileNameOutBase  = argv[2];
		scrFileName = argv[3];
		std::cout<<"Processing file \""<<filenameIn<<"\""<<std::endl;
		writeOutput = true;
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
	int simCycle = data.getDataSet()->GetFieldData()->GetArray("CYCLE")->GetTuple1(0);
	double pointMass = data.getDataSet()->GetFieldData()->GetArray("PMASS")->GetTuple1(0);
	double pointGrav = data.getDataSet()->GetFieldData()->GetArray("PMGRV")->GetTuple1(0);

	std::cout<<"File information:"<<std::endl;
	std::cout<<"\t             Time: "<<simTime<<std::endl;
	std::cout<<"\t   Point Mass [g]: "<<pointMass<<std::endl;
	std::cout<<"\tPoint Mass [Msun]: "<<pointMass/MSOLAR<<std::endl;
	std::cout<<"\t Point Mass current [Msun]: "<<pointGrav/MSOLAR<<std::endl;

	// data.getDataSet()->GetFieldData()->PrintSelf(std::cout,vtkIndent(0));

	// Store
	dataList.addScalar(simTime,"simtime");
	dataList.addScalar(simCycle,"cycle");
	dataList.addScalar(pointMass,"pointmass");


	// Get domain and mesh information
	double 	physicalBounds[6]; // Physical space bounds
	int 		meshDimensions[3]; // Number of cells in each dimension
	data.getDataSet()->GetBounds(physicalBounds);
	data.getDataDimension(meshDimensions);

	std::cout<<"Mesh dimensions:"<<std::endl;
	std::cout<<meshDimensions[0]<<"\t"<<meshDimensions[1]<<"\t"<<meshDimensions[2]<<std::endl;

	int nPhi = (nDim>1 ? meshDimensions[1] : 1);
	int nTheta = (nDim>2 ? meshDimensions[2] : 1);

	double 	trueVolume = 4.0/3.0*vtkMath::Pi()*(	physicalBounds[1]*physicalBounds[1]*physicalBounds[1] - 
												physicalBounds[0]*physicalBounds[0]*physicalBounds[0] );

	double 	simVolume = 1.0;
	simVolume *= (1.0/3.0)*physicalBounds[1]*physicalBounds[1]*physicalBounds[1] - physicalBounds[0]*physicalBounds[0]*physicalBounds[0];
	
	if(nDim>1)
	{
		simVolume *= -cos(physicalBounds[3]) + cos(physicalBounds[2]);	
	}
	else
	{
		simVolume *= 2.0;
	}
	if(nDim>2)
	{
		simVolume *= physicalBounds[5] - physicalBounds[4];
	}
	else
	{
		simVolume *= 2.0*vtkMath::Pi();
	}

	// double volumeCorrection = trueVolume/simVolume;
	double volumeCorrection = 1.0;
	if(nDim>1)
	{
		2.0/(-cos(physicalBounds[3]) + cos(physicalBounds[2]));
	}
	std::cout<<"Ratio of true volume to simulated volume: "<<volumeCorrection<<std::endl;

	data.computeCellVolumes("cellVolumeTmp");
	data.storeData(volumeCorrection*data["cellVolumeTmp"],"cellVolume");

	std::cout<<"True volume/Corrected volume: "<<trueVolume/data.sum("cellVolume")<<std::endl;


	// Create mass variable
	std::clog<<"Creating mass variable..."; std::clog.flush();
	data.storeData(data["dens"]*data["cellVolume"],"mass");
	std::clog<<"done!"<<std::endl;

	// Create data arrays of coordinates (full 2d/3d arrays of cell centers)
	std::cout<<"Creating coordinate data arrays...";std::cout.flush();
	data.createCoordinateDataArrays("r","phi","theta");
	std::cout<<"done"<<std::endl;
	vtkSmartPointer<vtkDataArray> radCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();

	/*====================================================================
		
		Compute shock statistics for the model. 
		(1)	Get all the cells along a radial ray
		(1.1)	Moving from the outside in, find the first band of 
				"shocked" cells
		(1.2)	Determine the average radius for these cells
		(1.3)	Store this information as the shock position for
				this ray
		(2)	Aggregate these shock positions for every ray (given a phi-theta pair)
		(2.1)	Compute the min, max, mean

	=======================================================================*/

	std::clog<<"Computing shock statistics..."; std::clog.flush();

	vtkSmartPointer<vtkDataArray> shock = data["shock"];
	vtkSmartPointer<vtkDataArray> radius = data["r"];
	std::vector<double> shockPositions;
	std::vector<int> cellIndicesOutsideShock;
	int startingIndex[3];

	std::clog<<"shock positions..."; std::clog.flush();

	for(int p = 0; p<nPhi; ++p)
	{
		for(int t = 0; t< nTheta; ++t)
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
						cellIndicesOutsideShock.push_back(*(i-1));
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

			// Add the ray's shock position to the vector
			shockPositions.push_back(shockPosAvg);
		}
	}


	std::clog<<"min/max/mean..."; std::clog.flush();

	// Compute the min, max, and mean shock position in the domain
	double shockPosMin = physicalBounds[1];
	double shockPosMax = physicalBounds[0];
	double shockPosMean = 0.0;
	for(int i=0; i<shockPositions.size(); i++)
	{
		double sp = shockPositions[i];
		// Max
		if(sp > shockPosMax)
		{
			shockPosMax = sp;
		}

		// Min
		if(sp < shockPosMin)
		{
			shockPosMin = sp;
		}

		// Average sum
		shockPosMean += sp;
	}

	// Finish computing the average
	if(shockPositions.size() > 0)
	{
		shockPosMean /= double(shockPositions.size());
	}

	// Compute the sample standard variance (s^2 = sqrt[1/(N-1) Sum[(x-xavg)^2]])
	double shockPosVar = 0.0;
	for(int i=0; i<shockPositions.size(); i++)
	{
		double sp = shockPositions[i];
		shockPosVar += (sp-shockPosMean)*(sp-shockPosMean);
	}

	// Finish computing the variance
	if(shockPositions.size()-1 > 0)
	{
		shockPosVar = sqrt(shockPosVar/double(shockPositions.size()-1));
	}


	std::clog<<"done!"<<std::endl;

	std::cout<<"Shock statistics"<<std::endl;
	std::cout<<"\t Mean: "<<shockPosMean<<std::endl;
	std::cout<<"\t  Min: "<<shockPosMin<<std::endl;
	std::cout<<"\t  Max: "<<shockPosMax<<std::endl;
	std::cout<<"\t  Var: "<<shockPosVar<<std::endl;


	// Store
	dataList.addScalar(shockPosMean,"shockMean");
	dataList.addScalar(shockPosMin,"shockMin");
	dataList.addScalar(shockPosMax,"shockMax");
	dataList.addScalar(shockPosVar,"shockVar");


	/*====================================================================
		
		Create the gain region mask
		(1) Moving inward along radial rays, determine where we cross the shock
		(2) Begin flagging cells as in the GR
		(3) When qenr<0, we've left the gain region. Store the innermost cell
			of the GR in a list so we can compute the MFR later


	=======================================================================*/

	std::clog<<"Computing gain region mask..."; std::clog.flush();

	data.addVariable("mask");
	vtkSmartPointer<vtkDataArray> qenr = data["qenr"];
	vtkSmartPointer<vtkDataArray> mask = data["mask"];
	std::vector<int> cellIndicesOutsideGR;

	for(int p = 0; p<meshDimensions[1]; ++p)
	{
		for(int t = 0; t< (nDim>2 ? meshDimensions[2] : 1); ++t)
		{

			// std::cout<<"p = "<<p<<"\tt= "<<t<<std::endl;
			
			startingIndex[XAXIS]	= 0;
			startingIndex[YAXIS]	= p;
			startingIndex[ZAXIS]	= t;
			std::vector<int> cellIndices = data.getCellsAlongRay(XAXIS, startingIndex);

			bool foundFirst = false;
			double shockPosAvg = 0.0;
			int nShocked = 0;
			bool belowShock = false;

			for (std::vector<int>::reverse_iterator i = cellIndices.rbegin(); i != cellIndices.rend(); ++i)
			{
				
				bool isShocked = shock->GetTuple1(*i) > 0.5; // The data is either 0 or 1, so you need to pick something reasonable
				bool isQenrPositive = qenr->GetTuple1(*i) > 0.0;

				// This section of if statements determines if we've made it inside the shock
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
					// Therefore, we should start considering that we're in the GR
					if(foundFirst)
					{
						belowShock = true;
					}

				}

				// This section of if statements determines if we're above the GR (along this ray)
				// If we are, set the mask flag to 1.0
				// If this evaluates to false, that means we've 
				if(belowShock) // If we're below the shock
				{
					if(isQenrPositive) // And qenr is positive
					{
						// Then we should flag this cell as part of the mask
						mask->SetTuple1(*i,1.0);
					}
					else // and qenr is negative
					{
						// Then we've left the GR
						// We should save this cell and then break the loop
						cellIndicesOutsideGR.push_back(*i);
						break;
					}
				}

			}
		}
	}
	std::clog<<"done!"<<std::endl;


	double Mgain = data.sum("mass","mask");
	std::cout<<"Total species mass [Msun]: "<<data.sum("mass")/MSOLAR<<std::endl;
	std::cout<<"Gain region mass:"<<std::endl;
	std::cout<<"\t    Mass [g]: "<<Mgain<<std::endl;
	std::cout<<"\t Mass [Msun]: "<<Mgain/MSOLAR<<std::endl;


	// Store
	dataList.addScalar(Mgain,"massgainregion");



	/*====================================================================
		
		Compute statistics about the gain radius
		(1) Min, Mean, Max, Var

	=======================================================================*/
	// Compute the min, max, and mean shock position in the domain
	double grPosMin = physicalBounds[1];
	double grPosMax = physicalBounds[0];
	double grPosMean = 0.0;

	for (std::vector<int>::iterator i = cellIndicesOutsideGR.begin(); i != cellIndicesOutsideGR.end(); ++i)
	{
		double sp = radius->GetTuple1(*i);
		// Max
		if(sp > grPosMax)
		{
			grPosMax = sp;
		}

		// Min
		if(sp < grPosMin)
		{
			grPosMin = sp;
		}

		// Average sum
		grPosMean += sp;
	}

	// Finish computing the average
	if(cellIndicesOutsideGR.size() > 0)
	{
		grPosMean /= double(cellIndicesOutsideGR.size());
	}

	// Compute the sample standard variance (s^2 = sqrt[1/(N-1) Sum[(x-xavg)^2]])
	double grPosVar = 0.0;
	for (std::vector<int>::iterator i = cellIndicesOutsideGR.begin(); i != cellIndicesOutsideGR.end(); ++i)
	{
		double sp = radius->GetTuple1(*i);
		grPosVar += (sp-grPosMean)*(sp-grPosMean);
	}

	// Finish computing the variance
	if(cellIndicesOutsideGR.size()-1 > 0)
	{
		grPosVar = sqrt(grPosVar/double(cellIndicesOutsideGR.size()-1));
	}


	std::clog<<"done!"<<std::endl;

	std::cout<<"Gain region statistics"<<std::endl;
	std::cout<<"\t Mean: "<<grPosMean<<std::endl;
	std::cout<<"\t  Min: "<<grPosMin<<std::endl;
	std::cout<<"\t  Max: "<<grPosMax<<std::endl;
	std::cout<<"\t  Var: "<<grPosVar<<std::endl;


	// Store
	dataList.addScalar(grPosMean,"grMean");
	dataList.addScalar(grPosMin,"grMin");
	dataList.addScalar(grPosMax,"grMax");
	dataList.addScalar(grPosVar,"grVar");


	/*====================================================================
		
		Compute the mass flow rate at the shock front. 

		(1)	Compute individual mass flow rates for each ray,
				using the cell directly after the shock (moving outward)

	=======================================================================*/
	double MFR = 0.0;
	for (std::vector<int>::iterator i = cellIndicesOutsideShock.begin(); i != cellIndicesOutsideShock.end(); ++i)
	{
		// Determine the area at the volume centroid of the cell that the radial flow sees
		double cellEdges[6];
		data.getDataSet()->GetCellBounds(*i,cellEdges);
		double r = data["r"]->GetTuple1(*i);
		double cellArea = 1.0;
		cellArea *= r*r; // Fixed radius
		if(nDim>1)
		{
			cellArea *= -cos(cellEdges[3]) + cos(cellEdges[2]); // Phi contribution
		}
		else
		{
			cellArea *= 2.0;
		}
		if(nDim>2)
		{
			cellArea *= cellEdges[5]-cellEdges[4]; // theta contribution
		}
		else
		{
			cellArea *= 2.0*vtkMath::Pi();
		}

		// Correct for the full 4pir^2 area

		cellArea *= volumeCorrection;

		MFR += -1*data["dens"]->GetTuple1(*i)*data["velx"]->GetTuple1(*i)*cellArea;
	}



	std::cout<<"Mass flow rate information at the shock"<<std::endl;
	std::cout<<"\t    MFR [g/s]: "<<MFR<<std::endl;
	std::cout<<"\t MFR [Msun/s]: "<<MFR/MSOLAR<<std::endl;

	// Store
	dataList.addScalar(MFR,"mfrShock");

	/*====================================================================
		
		Compute the mass flow rate from the lower GR boundary 

		(1)	Compute individual mass flow rates for each ray,
				using the cell directly below the GR 

	=======================================================================*/
	double mfrGR = 0.0;
	for (std::vector<int>::iterator i = cellIndicesOutsideGR.begin(); i != cellIndicesOutsideGR.end(); ++i)
	{
		// Determine the area at the volume centroid of the cell that the radial flow sees
		double cellEdges[6];
		data.getDataSet()->GetCellBounds(*i,cellEdges);
		double r = data["r"]->GetTuple1(*i);
		double cellArea = 1.0;
		cellArea *= r*r; // Fixed radius
		if(nDim>1)
		{
			cellArea *= -cos(cellEdges[3]) + cos(cellEdges[2]); // Phi contribution
		}
		else
		{
			cellArea *= 2.0;
		}
		if(nDim>2)
		{
			cellArea *= cellEdges[5]-cellEdges[4]; // theta contribution
		}
		else
		{
			cellArea *= 2.0*vtkMath::Pi();
		}

		// Correct for the full 4pir^2 area

		cellArea *= volumeCorrection;

		mfrGR += data["dens"]->GetTuple1(*i)*data["velx"]->GetTuple1(*i)*cellArea;
	}



	std::cout<<"Mass flow rate information at the gain radius"<<std::endl;
	std::cout<<"\t    MFR [g/s]: "<<mfrGR<<std::endl;
	std::cout<<"\t MFR [Msun/s]: "<<mfrGR/MSOLAR<<std::endl;

	// Store
	dataList.addScalar(mfrGR,"mfrGR");

	/*====================================================================
		
		Compute the specific binding energy

		(1)	First, create a new variable to hold "gpot" - point gravity
				"gpot" only holds the gravity of the material outside
				the PNS
		(2)	Compute the specific kinetic energy
		(3)	Compute the specific internal energy
		(4)	Compute the specific binding energy as 
				internal + kinetic + grav (note: grav is negative)

	=======================================================================*/

				
	data.addVariable("gpotPMCorners");
	for(int i=0; i<data["gpotPMCorners"]->GetNumberOfTuples(); i++)
	{

		double cellEdges[6];
		data.getDataSet()->GetCellBounds(i,cellEdges);
		double rMin = cellEdges[0];
		double rMax = cellEdges[1];

		double gravPot = 0.0;

		// 	Determine the average gravitational potential from the corners 
		// 	In spherical geometry, the radius is the only component that matters
		// 	Therefore, 
		//		(2d) => (2*GM/r_min + 2*GM/r_min)/(4 corners)
		//		(3d) => (4*GM/r_min + 4*GM/r_max)/(8 corners)
		//	Both cases reduce to
		//		(GM/r_min + GM/r_max)/2

		gravPot = 0.5*GRAVCONST*pointMass*((rMin+rMax)/(rMin*rMax));

		// Store the value
		data["gpotPMCorners"]->SetTuple1(i,gravPot);

	}


	// Compute the total gravitational energy
	data.storeData( GRAVCONST*pointMass/data["r"], "gpotPM" );
	// data.storeData( 0.5*(data["gpot"] - data["gpotPM"]), "egrav" );
	data.storeData( 0.5*(data["gpot"] - data["gpotPMCorners"]), "egrav" );
	// data.storeData( data["gpot"], "egrav" );

	// Compute the kinetic energy
	if(nDim==2)
	{
		data.storeData( 0.5*(data["velx"]*data["velx"] + data["vely"]*data["vely"]), "ekin" );
	}
	else if(nDim==3)
	{
	data.storeData( 0.5*(data["velx"]*data["velx"] + data["vely"]*data["vely"] + data["velz"]*data["velz"]), "ekin" );
	}

	// Compute the internal energy
	data.storeData( (1.0/(data["gammae"]-1.0))*data["pres"]/data["dens"], "eint" );

	// Compute the thermal + kinetic energy
	// data.storeData( data["ekin"] + data["eint"], "eInt+Kin" );
	// data.storeData( data["eInt+Kin"]*data["dens"]*data["cellVolume"], "eInt+KinMassWeighted" );


	// Compute the binding energy
	data.storeData( data["ekin"] + data["eint"] + data["egrav"], "ebind" );


	// Compute the explosion energy for the entire star
	data.storeData( data["ebind"]*data["dens"]*data["cellVolume"], "ebindMassWeighted");

	// data.storeData( data["ekin"]*data["dens"]*data["cellVolume"], "ekinMassWeighted");
	// data.storeData( data["eint"]*data["dens"]*data["cellVolume"], "eintMassWeighted");
	// data.storeData( data["egrav"]*data["dens"]*data["cellVolume"], "egravMassWeighted");



	// Create masks to filter positive and negative binding energy summations
	// data.storeData( data["ebind"]>0, "ebindPositiveMask");
	data.storeData( data["ebind"]<0, "ebindNegativeMask");
	// data.storeData( data["ebindPositiveMask"]*data["mask"], "ebindPositiveMaskGR");
	data.storeData( data["ebindNegativeMask"]*data["mask"], "ebindNegativeMaskGR");



	/*
	double ekinStar = data.sum("ekinMassWeighted");
	double eintStar = data.sum("eintMassWeighted");
	double egravStar = data.sum("egravMassWeighted");
	double ekinGR = data.sum("ekinMassWeighted","mask");
	double eintGR = data.sum("eintMassWeighted","mask");
	double egravGR = data.sum("egravMassWeighted","mask");

	double bindingNegStar = data.sum("ebindMassWeighted","ebindNegativeMask");
	double bindingPosStar = data.sum("ebindMassWeighted","ebindPositiveMask");
	double bindingNegGR = data.sum("ebindMassWeighted","ebindNegativeMaskGR");
	double bindingPosGR = data.sum("ebindMassWeighted","ebindPositiveMaskGR");
	
	double expEnergyTest = 0.0;
	for(int i=0; i<data["ebindMassWeighted"]->GetNumberOfTuples(); i++)
	{
		double value = data["ebindMassWeighted"]->GetTuple1(i);
		if(value > 0.0)
		{
			expEnergyTest += value;
		}
	}
	


	std::cout.setf(std::ios::scientific,std::ios::floatfield);
	std::cout.setf(std::ios::showpos);
	std::cout.precision(15);
	std::cout<<"\tEntire star:"<<std::endl;
	std::cout<<"\t\t   Total kinetic energy: "<<ekinStar/unitFOE<<std::endl;
	std::cout<<"\t\t  Total internal energy: "<<eintStar/unitFOE<<std::endl;
	std::cout<<"\t\t   Total int+kin energy: "<<data.sum("eInt+KinMassWeighted")/unitFOE<<std::endl;
	std::cout<<"\t\t      Total grav energy: "<<egravStar/unitFOE<<std::endl;
	std::cout<<"\t\t           Total energy: "<<data.sum("ebindMassWeighted")/unitFOE<<std::endl;
		std::cout<<"\t\t       Explosion energy: "<<data.sum("ebindMassWeighted","ebindPositiveMask")/unitFOE<<std::endl;
	std::cout<<"\t\tNegative binding energy: "<<data.sum("ebindMassWeighted","ebindNegativeMask")/unitFOE<<std::endl;


	// std::cout<<"\tGain region:"<<std::endl;
	// std::cout<<"\t\t   Total kinetic energy: "<<ekinGR/unitFOE<<std::endl;
	// std::cout<<"\t\t  Total internal energy: "<<eintGR/unitFOE<<std::endl;
	// std::cout<<"\t\t      Total grav energy: "<<egravGR/unitFOE<<std::endl;
	// std::cout<<"\t\t       Explosion energy: "<<data.sum("ebindMassWeighted","ebindPositiveMaskGR")/unitFOE<<std::endl;
	// std::cout<<"\t\tNegative binding energy: "<<data.sum("ebindMassWeighted","ebindNegativeMaskGR")/unitFOE<<std::endl;
	

	double eKinHOTB = data.getDataSet()->GetFieldData()->GetArray("EKIN")->GetTuple1(0);
	double eGravHOTB = data.getDataSet()->GetFieldData()->GetArray("EGRAV")->GetTuple1(0);
	double eIntHOTB = data.getDataSet()->GetFieldData()->GetArray("EINT")->GetTuple1(0);
	double eBindPHOTB = data.getDataSet()->GetFieldData()->GetArray("EBINDP")->GetTuple1(0);
	double eBindNHOTB = data.getDataSet()->GetFieldData()->GetArray("EBINDM")->GetTuple1(0);

	std::cout<<" eKinHOTB: "<<eKinHOTB<<std::endl;
	std::cout<<" eIntHOTB: "<<eIntHOTB<<std::endl;
	std::cout<<"eGravHOTB: "<<eGravHOTB<<std::endl;
	std::cout<<" Error in kinetic: "<<(10.*ekinStar/unitFOE-eKinHOTB)/eKinHOTB*100.0<<std::endl;
	std::cout<<"Error in internal: "<<(10.*eintStar/unitFOE-eIntHOTB)/eIntHOTB*100.0<<std::endl;
	std::cout<<"    Error in grav: "<<(10.*egravStar/unitFOE-eGravHOTB)/eGravHOTB*100.0<<std::endl;
	// std::cout<<"expEnergyTest: "<<expEnergyTest<<std::endl;


	std::cout.unsetf(std::ios::floatfield);
	std::cout.unsetf(std::ios::showpos);
	std::cout.precision(5);

	// Print if the explosion is successful based on the current explosion energy in the entire star
	double explosionEnergy = data.sum("ebindMassWeighted","ebindPositiveMask");
	std::cout<<"Explosion energy: "<<explosionEnergy/unitFOE<<std::endl;
	std::cout<<"Explosion occured? "<< (explosionEnergy>explosionCriterion) <<std::endl;

	
	// Store
	dataList.addScalar(explosionEnergy,"exposionenergy");
	dataList.addScalar(explosionCriterion,"explosioncriterion");
	dataList.addScalar(bindingNegStar,"bindingnegstar");
	dataList.addScalar(bindingPosStar,"bindingposstar");
	dataList.addScalar(bindingNegGR,"bindingneggr");
	dataList.addScalar(bindingPosGR,"bindingposgr");

	dataList.addScalar(ekinStar,"ekinstar");
	dataList.addScalar(eintStar,"eintstar");
	dataList.addScalar(egravStar,"egravstar");
	dataList.addScalar(ekinGR,"ekingr");
	dataList.addScalar(eintGR,"eintgr");
	dataList.addScalar(egravGR,"egravgr");


	dataList.addScalar(double(explosionEnergy>explosionCriterion),"explosionstarted");

	

	// Store
	dataList.addScalar(ekinGR,"ekingr");
	dataList.addScalar(eintGR,"eintgr");
	dataList.addScalar(egravGR,"egravgr");
	*/


	// Get HOTB calculated energies from the VTK file
	double eKinHOTB = data.getDataSet()->GetFieldData()->GetArray("EKIN")->GetTuple1(0);
	double eGravHOTB = data.getDataSet()->GetFieldData()->GetArray("EGRAV")->GetTuple1(0);
	double eIntHOTB = data.getDataSet()->GetFieldData()->GetArray("EINT")->GetTuple1(0);
	double eBindPHOTB = data.getDataSet()->GetFieldData()->GetArray("EBINDP")->GetTuple1(0);
	double eBindNHOTB = data.getDataSet()->GetFieldData()->GetArray("EBINDM")->GetTuple1(0);


	dataList.addScalar(eKinHOTB,"ekin");
	dataList.addScalar(eIntHOTB,"eint");
	dataList.addScalar(eGravHOTB,"egrav");
	dataList.addScalar(eBindPHOTB,"ebindp");
	dataList.addScalar(eBindNHOTB,"ebindn");

	// Parse explosion time and theshold from scr file

	double explosionTime = -1.0;

	std::ifstream scrFID;
	scrFID.open(scrFileName.c_str());

	if(scrFID.is_open())
	{
		std::string line;

		// Read file line-by-line until "Neutrino luminosities" is found.
		while(scrFID.good())
		{
			std::getline(scrFID, line);

			if(line.find("[WRITE_INTEGRALS] Energy threshold for explosion")!=std::string::npos)
			{
				std::stringstream ss;

				// Find the position of the equal sign
				// The explosion time is pos+1:end
				size_t eqPos = line.find("=");
				ss<<line.substr(eqPos+1);
				ss>>explosionTime;

				// We're done here
				break;
			}

		}

		scrFID.close();
	}
	else
	{
		std::cerr<<"[ERROR] Problem opening \""<<scrFileName<<"\". Neutrino luminosities will not be available."<<std::endl;
	}

	// Store
	dataList.addScalar(explosionTime,"exptime");

	std::cout<<"Explosion Time: "<<explosionTime<<std::endl;





	/*====================================================================
		
		Compute the partition of mass between down flows (velx<0)
		and up drafts (velx>0)


	=======================================================================*/

	data.storeData( (data["velx"]>0.0)*data["mask"], "updraftMask");
	data.storeData( (data["velx"]<0.0)*data["mask"], "downdraftMask");

	double massUpDraft = data.sum("mass","updraftMask");
	double massDownDraft = data.sum("mass","downdraftMask");
	double totalDraftMass = massUpDraft + massDownDraft;

	std::cout<<"Draft information: "<<std::endl;
	std::cout<<"\t  Up Draft Mass [Msun]: "<<massUpDraft/MSOLAR<<std::endl;
	std::cout<<"\tDown Draft Mass [Msun]: "<<massDownDraft/MSOLAR<<std::endl;
	std::cout<<"\t     Total Mass [Msun]: "<<totalDraftMass/MSOLAR<<std::endl;
	std::cout<<"\t    Up Draft Mass Frac: "<<massUpDraft/totalDraftMass<<std::endl;
	std::cout<<"\t  Down Draft Mass Frac: "<<massDownDraft/totalDraftMass<<std::endl;


	// Store
	dataList.addScalar(massUpDraft,"massupdraft");
	dataList.addScalar(massDownDraft,"massdowndraft");



	/*====================================================================
		
		Compute the radial momentum in the gain region


	=======================================================================*/
	data.storeData(data["velx"]*data["dens"],"radmomDens");

	double radmomDens = data.sum("radmomDens","mask");
	double radmomDensUp = data.sum("radmomDens","updraftMask");
	double radmomDensDown = data.sum("radmomDens","downdraftMask");

	// Store
	dataList.addScalar(radmomDens,"radmomDens");
	dataList.addScalar(radmomDensUp,"radmomDensUp");
	dataList.addScalar(radmomDensDown,"radmomDensDown");



	/*====================================================================
		
		Compute the advective and heating timescales

		(1)	Compute the advective timescale using Tadv = Mgain/MFR@shock
		(2)	Heating timescale is the ration between the Sum[mass weighted ebind]
				and the Sum[mass weighted qenr]

	=======================================================================*/

	// Compute the advective timescale
	double timescaleAdv = Mgain/MFR;


	// Compute the de-volume weighted net heating rate
	data.storeData( data["qenr"]*data["cellVolume"], "qenrMassWeighted");

	// Compute the total net heating rate
	double netHeatingRate = data.sum("qenrMassWeighted","mask");

	// Compute the heating timescale
	double timescaleHeat = fabs(data.sum("ebindMassWeighted","ebindNegativeMaskGR"))/netHeatingRate;

	// Compute the efficiency of heating
	double heatingEfficiency = timescaleAdv/timescaleHeat;


	std::cout<<"Timescale information"<<std::endl;
	std::cout<<"\tNet Heating Rate [1e51 erg/s]: "<<netHeatingRate/unitFOE<<std::endl;
	std::cout<<"\tAdvection [s]: "<<timescaleAdv<<std::endl;
	std::cout<<"\t  Heating [s]: "<<timescaleHeat<<std::endl;
	std::cout<<"\t   Efficiency: "<<heatingEfficiency<<std::endl;

	// Store
	dataList.addScalar(netHeatingRate,"netheatingrate");
	dataList.addScalar(timescaleAdv,"timescaleadv");
	dataList.addScalar(timescaleHeat,"timescaleheat");
	dataList.addScalar(heatingEfficiency,"heatingefficiency");

	/*====================================================================
		
		Compute the antesonic condition specified by Pejcha (2012).

		max(sound_speed^2/local_escape_velocity^2)

		We take the approach of Mueller (2012, Rel Exp Models of ccSNe)
		and calculate both the point-wise maximum in the post-shock region,
		as well as averaging the local sound speed over spherical shells and then
		taking the maximum.

		The escape velocity is given by sqrt(2GM/r)
		The sound speed is given by sqrt(gammac*pres/dens)

	=======================================================================*/

	data.storeData(data["gammac"]*data["pres"]/data["dens"],"sound_speedSqr");
	data.storeData((2.0*pointMass*GRAVCONST)/data["r"],"escape_velocitySqr");
	data.storeData(data["sound_speedSqr"]/data["escape_velocitySqr"],"antesonicPointwise");


	// Average the sound speed
	std::vector<std::string> antesonicVars;
	antesonicVars.push_back("sound_speedSqr");

	ShellAveragePlaneOperator antesonicAvg(	radCoords,
											XAXIS,
											CS_SPHERE,
											antesonicVars );
	antesonicAvg.process(data.getDataSet());

	vtkSmartPointer<vtkDataArray> soundSpeedSqrAvg = antesonicAvg.getDataArray("sound_speedSqr");

	// Compute the pointwise maximum
	double antesonicPWMax = 0.0;
	double antesonicPWMaxRadius = 0.0;

	for(int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{
		// Gain region check
		if(data.getDataArray("mask")->GetTuple1(i)<0.5)
		{
			continue;
		}

		double tmp = data.getDataArray("antesonicPointwise")->GetTuple1(i);
		if(tmp>antesonicPWMax)
		{
			antesonicPWMax = tmp;
			antesonicPWMaxRadius = data.getDataArray("r")->GetTuple1(i);
		}

	}

	// Compute the average maximum
	double antesonicAvgMax = 0.0;
	double antesonicAvgMaxRadius = 0.0;

	for(int i=0; i<antesonicAvg.getBins()->GetNumberOfCells(); i++)
	{
		double cellBounds[6];
		antesonicAvg.getBins()->GetCellBounds(i,cellBounds);
		double cellRadius = std::pow(0.5*(cellBounds[0]*cellBounds[0]*cellBounds[0] + 
										  cellBounds[1]*cellBounds[1]*cellBounds[1]),1.0/3.0);

		// Ensure we're inside the GR
		if(cellRadius < grPosMax || cellRadius > shockPosMin)
		{
			continue;
		}

		double tmp = antesonicAvg.getDataArray("sound_speedSqr")->GetTuple1(i)/(2.0*pointMass*GRAVCONST/cellRadius);

		if(tmp > antesonicAvgMax)
		{
			antesonicAvgMax = tmp;
			antesonicAvgMaxRadius = cellRadius;
		}
	}

	// Store data
	dataList.addScalar(antesonicPWMax,"antesonic_pointwise_max");
	dataList.addScalar(antesonicPWMaxRadius,"antesonic_pointwise_rad");
	dataList.addScalar(antesonicAvgMax,"antesonic_average_max");
	dataList.addScalar(antesonicAvgMaxRadius,"antesonic_average_rad");

	std::cout<<"Antesonic condition:"<<std::endl;
	std::cout<<"\tPointwise: "<<antesonicPWMax<<" ("<<antesonicPWMaxRadius<<")"<<std::endl;
	std::cout<<"\t Averaged: "<<antesonicAvgMax<<" ("<<antesonicAvgMaxRadius<<")"<<std::endl;


	/*====================================================================
		
		Parse neutrino luminosities from scr file

		In the scr file, the neutrino luminosities are found in a section
		that looks like:
		
		BEGIN VERBATIM
		
		---------------------
		Neutrino luminosities
		---------------------
		L0_tot     (total neutrino luminosity)
		L0_nu      (electron neu) (electron antineu) (muon neu) (muon antineu) (tau neu) (tau antineu)

		END VERBATIM

		We really only need electron neutrinos, but might as well get them all

	=======================================================================*/

	double luminosityTotal = 0.0;
	double luminosities[6]; memset(luminosities,0.0,6);

	// std::ifstream scrFID;
	scrFID.open(scrFileName.c_str());

	if(scrFID.is_open())
	{
		std::string line;

		// Read file line-by-line until "Neutrino luminosities" is found.
		while(scrFID.good())
		{
			std::getline(scrFID, line);
			// std::cout<<line<<std::endl;

			if(line.compare("Neutrino luminosities")==0)
			{
				std::string identifier;
				std::stringstream ss;
				// Discard next "-------" line
				std::getline(scrFID, line);

				// Read total luminosity
				ss.str(""); ss.clear();
				std::getline(scrFID, line);
				ss<<line;
				ss>>identifier>>luminosityTotal;

				// Read individual luminosities
				ss.str(""); ss.clear();
				std::getline(scrFID, line);
				ss<<line;
				ss>>identifier;
				for(int i=0; i<6; i++)
				{
					double tmp = 0.0;
					ss>>tmp;
					luminosities[i] = tmp;
				}

				// We're done here
				break;
			}

		}

		scrFID.close();
	}
	else
	{
		std::cerr<<"[ERROR] Problem opening \""<<scrFileName<<"\". Neutrino luminosities will not be available."<<std::endl;
	}


	// Store
	dataList.addScalar(luminosityTotal,"lumtotal");
	dataList.addScalar(luminosities[0],"lumelec");
	dataList.addScalar(luminosities[1],"lumelecanti");
	dataList.addScalar(luminosities[2],"lummuon");
	dataList.addScalar(luminosities[3],"lummuonanti");
	dataList.addScalar(luminosities[4],"lumtau");
	dataList.addScalar(luminosities[5],"lumtauanti");


	// Print
	std::cout<<"Neutrino luminosities: "<<std::endl;
	std::cout<<"\t"<<luminosities[0]<<std::endl;
	std::cout<<"\t"<<luminosities[1]<<std::endl;
	std::cout<<"\t"<<luminosities[2]<<std::endl;
	std::cout<<"\t"<<luminosities[3]<<std::endl;
	std::cout<<"\t"<<luminosities[4]<<std::endl;
	std::cout<<"\t"<<luminosities[5]<<std::endl;
	std::cout<<"\t"<<luminosityTotal<<std::endl;




	/*====================================================================
		
		Write scalar data stored in the ScalarList to file

	=======================================================================*/
	if(writeOutput)
	{
		std::cout<<"Writing scalars to expstat file..."<<std::endl;
		std::cout<<(fileNameOutBase + ".expstat").c_str()<<std::endl;
		
		std::ofstream scalarOutput;
		scalarOutput.open((fileNameOutBase + ".expstat").c_str());
		dataList.print(scalarOutput);
		scalarOutput.close();
	}


	// vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
	// writer->SetFileName("rectilineargrid.vtr");
	// writer->SetInputData(data.getDataSet());
	// writer->Write();



	std::cout<<"Finished"<<std::endl;
}



