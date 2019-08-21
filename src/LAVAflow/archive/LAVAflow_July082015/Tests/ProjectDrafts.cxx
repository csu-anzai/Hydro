#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/includes/LAVAconstants.h"
#include "../libsrc/Utilities/LAVAUtil.h"
#include "../libsrc/Clustering/FloodFill/FloodFill.h"

#include <string>
#include <vtkRectilinearGrid.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
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
		filenameIn = "/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/b163d2a3s530SG1MHWr2LBp.plt.0200.vtk";
		fileNameOutBase = "/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/output/b163d2a3s530SG1MHWr2LBp.plt.0200";
		nDim = 2;
		
		// filenameIn = "/data1/sne/HOTB/3d/b157d3a3s401SG1MHWr3LB/b157d3a3s401SG1MHWr3LB.plt.0100.vtk";
		// nDim = 3;

		// fileNameOutBase = filenameIn;
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
	std::cout<<"Creating coordinate data arrays..."; std::cout.flush();
	data.createCoordinateDataArrays("r","phi","theta");
	std::cout<<"done"<<std::endl;

	vtkSmartPointer<vtkDataArray> radEdges   = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();
	vtkSmartPointer<vtkDataArray> phiEdges   = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetYCoordinates();
	vtkSmartPointer<vtkDataArray> thetaEdges = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetZCoordinates();


	int nPhi = meshDimensions[1];
	int nTheta = nDim>2 ? meshDimensions[2] : 1;

	std::cout<<"nTheta: "<<nTheta<<std::endl;
	std::cout<<"nPhi: "<<nPhi<<std::endl;

/*====================================================================
		
		Compute the minimum and maximum extents of the gain region for every radial ray
		(1) Moving inward along radial rays, determine where we cross the shock
		(2) Begin flagging cells as in the GR
		(3) When qenr<0, we've left the gain region. Store the innermost cell
			of the GR in a list so we can compute the MFR later


	=======================================================================*/

	double** gainRadius;
	double** shockRadius;
	gainRadius  = new double*[nPhi];
	shockRadius = new double*[nPhi];
	for(int i=0; i<nPhi; i++)
	{
		gainRadius[i]  = new double[nTheta]; 
		shockRadius[i] = new double[nTheta];
		for(int j=0; j<nTheta; j++)
		{
			gainRadius[i][j] = 0.0;
			shockRadius[i][j] = 0.0;
		} 
	}

	vtkSmartPointer<vtkDataArray> shock = data["shock"];
	vtkSmartPointer<vtkDataArray> qenr  = data["qenr"];

	for(int p = 0; p<nPhi; ++p)
	{
		for(int t = 0; t<nTheta; ++t)
		{

			// std::cout<<"p = "<<p<<"\tt= "<<t<<std::endl;
			
			int startingIndex[3];
			startingIndex[XAXIS]	= 0;
			startingIndex[YAXIS]	= p;
			startingIndex[ZAXIS]	= t;
			std::vector<int> cellIndices = data.getCellsAlongRay(XAXIS, startingIndex);

			bool foundFirst = false;
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
					if(foundFirst && !belowShock)
					{
						shockRadius[p][t] = data["r"]->GetTuple1(*i);
						belowShock = true;
					}

				}

				// This section of if statements determines if we're above the GR (along this ray)
				// If we are, set the mask flag to 1.0
				// If this evaluates to false, that means we've 
				if(belowShock) // If we're below the shock
				{
					if(!isQenrPositive) // And qenr is positive
					{
						gainRadius[p][t] = data["r"]->GetTuple1(*(i));
						break;
					}
				}

			}
		}
	}


	/*====================================================================
		
		Smooth the radius arrays. This will remove "jagged" behavior 
		near the two boundaries

	=======================================================================*/



	/*====================================================================
		
		Compute the angular distribution of mass weighted velocity
		integrated along radial rays in the gain region. 
		(1)	Get all the cells along a radial ray
		(1.1)	Moving outward from the star to the shock, compute the 
				sum of radial momentum and normalize this by the total 
				volume in the ray.
	=======================================================================*/

	// Create VTK dataset
	vtkSmartPointer<vtkRectilinearGrid> projDataSet = vtkSmartPointer<vtkRectilinearGrid>::New();
	
	projDataSet->SetDimensions(nPhi+1,nTheta+1,1); // 2D rectilinear grid
	projDataSet->SetXCoordinates(vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetYCoordinates());
	projDataSet->SetYCoordinates(vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetZCoordinates());

	std::vector<double> grWidthFractionalOffset;
	grWidthFractionalOffset.push_back(0.0);
	grWidthFractionalOffset.push_back(0.1);
	grWidthFractionalOffset.push_back(0.2);
	grWidthFractionalOffset.push_back(0.3);
	grWidthFractionalOffset.push_back(0.4);

	double radMin = 0;
	double** radialMomentumIntegrals;
	int** image;

	radialMomentumIntegrals = new double*[nPhi];
	image = new int*[nPhi];
	for(int i=0; i<nPhi; i++)
	{
		radialMomentumIntegrals[i] = new double[nTheta]; 
		image[i] = new int[nTheta]; 
		for(int j=0; j<nTheta; j++)
		{
			radialMomentumIntegrals[i][j] = 0.0;
			image[i][j] = Clustering::FloodFill::UNUSED;
		}
	}


	// Create variables 
	data.addVariable("sgnVelx",0.0);
	for(int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{
		data["sgnVelx"]->SetTuple1(i,double(sgn(data["velx"]->GetTuple1(i))));
	}
	if(nDim<3)
	{
		data.storeData(data["velx"]*data["velx"] + data["vely"]*data["vely"], "kinetic");
	}
	else
	{
		data.storeData(data["velx"]*data["velx"] + data["vely"]*data["vely"]+ data["velz"]*data["velz"], "kinetic");
	}

	data.storeData(data["velx"]*data["dens"]*data["cellVolume"],"radMom");
	data.storeData(data["sgnVelx"]*data["dens"]*data["kinetic"]*data["cellVolume"],"signedKinetic");

	std::vector<std::string> analysisVariables;
	analysisVariables.push_back("radMom");
	analysisVariables.push_back("signedKinetic");


	std::string variableForAnalysis = "signedKinetic";


	char clustFileName[512];
	std::sprintf(clustFileName, "%s.cluster", fileNameOutBase.c_str());
	std::cout<<"Opening output file "<<clustFileName<<std::endl;
	std::ofstream clusterOutput;
	clusterOutput.open(clustFileName);
	clusterOutput<<"simTime"<<"\t"<<simTime<<std::endl;
	int startingIndex[3];


	for (std::vector<double>::iterator  widthFracIt= grWidthFractionalOffset.begin(); widthFracIt != grWidthFractionalOffset.end(); ++widthFracIt)
	{
	


		for(std::vector<std::string>::iterator varIt = analysisVariables.begin(); varIt != analysisVariables.end(); ++varIt)
		{

			std::string variableForAnalysis = *varIt;
			std::cout<<"Processing variable: "<<variableForAnalysis<<std::endl;

			for(int p = 0; p<nPhi; ++p)
			{
				for(int t = 0; t<nTheta; ++t)
				{

					// Reset image
					image[p][t] = Clustering::FloodFill::UNUSED;

					startingIndex[XAXIS]	= 0;
					startingIndex[YAXIS]	= p;
					startingIndex[ZAXIS]	= t;
					std::vector<int> cellIndices = data.getCellsAlongRay(XAXIS, startingIndex);

					double raySum = 0.0, rayVolume = 0.0;
					bool foundFirst = false;
					bool foundLast  = false;
					bool useCell     = false;
					int nCellsUsed = 0;

					double rGR = gainRadius[p][t];
					double rSH = shockRadius[p][t];
					double width = rSH-rGR;
					double offset = (*widthFracIt)*width;
					double rMin = rGR+offset;
					double rMax = rSH-offset;

						// std::cout<<"rGR: "<<rGR<<std::endl;
						// std::cout<<"rSH: "<<rSH<<std::endl;
						// std::cout<<"width: "<<width<<std::endl;
						// std::cout<<"offset: "<<offset<<std::endl;
						// std::cout<<"rMin: "<<rMin<<std::endl;
						// std::cout<<"rMax: "<<rMax<<std::endl;

					for (std::vector<int>::reverse_iterator i = cellIndices.rbegin(); i != cellIndices.rend(); ++i)
					{
						useCell = false;
						double cellRad = data["r"]->GetTuple1(*i);
						if(cellRad > rMin && cellRad < rMax)
						{
							useCell = true;
						}


						if(useCell)
						{
							raySum += data[variableForAnalysis.c_str()]->GetTuple1(*i);
							// rayVolume += data["dens"]->GetTuple1(*i)*data["cellVolume"]->GetTuple1(*i); 
							nCellsUsed++;
						}

					}

					// Store results
					// radialMomentumIntegrals[p][t] = raySum/rayVolume;
					radialMomentumIntegrals[p][t] = raySum;

					// Update image
					if(radialMomentumIntegrals[p][t]>0.0)
					{
						image[p][t] = Clustering::FloodFill::TARGET;					
					}

					if(nCellsUsed == 0)
					{
						std::cout<<"No cells used for a ray!"<<std::endl;
					}

				}
			}

			// Create flood fill clustering object
			std::cout<<"Creating FloodFill object..."; std::cout.flush();
			Clustering::FloodFill floodClusterObj(image, nPhi, nTheta);
			floodClusterObj.setBoundaryType(JAXIS,LOW,Clustering::FloodFill::FIXED);
			floodClusterObj.setBoundaryType(JAXIS,HIGH,Clustering::FloodFill::FIXED);
			floodClusterObj.setBoundaryType(KAXIS,LOW,Clustering::FloodFill::PERIODIC);
			floodClusterObj.setBoundaryType(KAXIS,HIGH,Clustering::FloodFill::PERIODIC);
			std::cout<<"done"<<std::endl;
			std::cout<<"Processing flood fill..."; std::cout.flush();
			floodClusterObj.process();
			std::cout<<"done"<<std::endl;


			// Compute the total solid angle occupied by each cluster
			// This involves looping over each cell in every cluster, computing the solid angle, and adding 
			// to the cluster sum
			int nClusters = floodClusterObj.getNumberOfClusters();

			if(nClusters==0)
			{
				break;
			}


			double* clusterSolidAngle = new double[nClusters];

			for(int i=0; i<nClusters; i++)
			{
				clusterSolidAngle[i] = 0.0;
				std::list<Clustering::FloodFill::cellIndices>* clustIndices = floodClusterObj.getClusterCellIndices(i);
				std::list<Clustering::FloodFill::cellIndices>::iterator cellIt;
				for (cellIt = clustIndices->begin(); cellIt != clustIndices->end(); ++cellIt)
				{
					double localSolidAngle = 0.0;
					double phiMin, phiMax, thetaMin, thetaMax;
					int phiInd   = (*cellIt).i;
					int thetaInd = (*cellIt).j;

					// Get the proper phi and theta bounds for the dimension
					phiMin = phiEdges->GetTuple1(phiInd);
					phiMax = phiEdges->GetTuple1(phiInd+1);

					if(nDim>2)
					{
						thetaMin = thetaEdges->GetTuple1(thetaInd);
						thetaMax = thetaEdges->GetTuple1(thetaInd+1);
					}
					else
					{
						thetaMin = 0.0;
						thetaMax = 2.0*vtkMath::Pi();
					}

					// Compute the local solid angle contribution
					localSolidAngle = (thetaMax-thetaMin)*(-std::cos(phiMax)+std::cos(phiMin));

					// Add to the cluster solid angle
					clusterSolidAngle[i] += localSolidAngle;
				}
			}
	 
			// Write the cluster solid angles to a file
			char variableID[256];
			sprintf(variableID, "%s_%0.2f",variableForAnalysis.c_str(), *widthFracIt);
			clusterOutput<<variableID<<"\t"<<nClusters<<std::endl;
			for(int i=0; i<nClusters; i++)
			{
				clusterOutput<<i<<"\t"<<clusterSolidAngle[i]<<std::endl;
			}

			double saSum = 0.0;
			for(int i=0; i<nClusters; i++)
			{
				saSum += clusterSolidAngle[i];
			}
			std::cout<<"Sum of cluster solid angles: "<<saSum<<std::endl;


			delete[] clusterSolidAngle;
		}

#if 0
		char varName[256];
		std::sprintf(varName,"projection_wf%0.2f",*widthFracIt);
		vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
		tmpArray->SetName(varName);
		std::sprintf(varName,"projection_wf%0.2f_clusters",*widthFracIt);
		vtkSmartPointer<vtkDoubleArray> clustArray  = vtkSmartPointer<vtkDoubleArray>::New();
		clustArray->SetName(varName);


			
		// Store data in vtk file
		for(int p=0; p<nPhi; p++)
		{
			for(int t=0; t<nTheta; t++)
			{
				tmpArray->InsertNextValue(0.0);
				clustArray->InsertNextValue(0.0);
			}
		}


		for(int p = 0; p<nPhi; ++p)
		{
			for(int t = 0; t<nTheta; ++t)
			{

				int ind = Util::sub2ind(p, t, 0, nPhi, nTheta);
				tmpArray->SetTuple1(ind,radialMomentumIntegrals[p][t]);
				clustArray->SetTuple1(ind,floodClusterObj.getPixelValue(p,t));

			}
		}

		projDataSet->GetCellData()->AddArray(tmpArray);
		projDataSet->GetCellData()->AddArray(clustArray);
#endif
	}

	// Close output file
	clusterOutput.close();

	


#if 0
	vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
	writer->SetFileName("draft_debug.vtr");
#if VTK_MAJOR_VERSION <= 5
  	writer->SetInputConnection(projDataSet->GetProducerPort());
#else
	writer->SetInputData(projDataSet);
#endif
	writer->Write();
#endif

	delete[] radialMomentumIntegrals;
	delete[] gainRadius;
	delete[] shockRadius;

	std::cout<<"Finished"<<std::endl;

}