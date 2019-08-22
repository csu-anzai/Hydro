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
#include "../libsrc/Clustering/SignatureBerger/SignatureBerger.h"

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

	// Create scalar list utility to store data we want to output
	Util::ScalarList dataList;




	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		// filenameIn = "../Tests/files/averaging/b155d2a3s123SG1MHwr2LB.plt.0010.vtk";
		filenameIn = "/home/tah09e/data/sne/data/HOTB/2d/b163d2a3s231SG1MHWr2LB/b163d2a3s231SG1MHWr2LB.plt.0092.vtk";
		// filenameIn = "/home/tah09e/data/sne/data/HOTB/2d/b163d2a3s231SG1MHWr2LB/b163d2a3s231SG1MHWr2LB.plt.0011.vtk";
		scrFileName = "/home/tah09e/data/sne/data/HOTB/2d/b163d2a3s231SG1MHWr2LB/scr_b163d2a3s231SG1MHWr2LB";
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

	// Get simulation information
	double simTime = data.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);

	std::cout<<"File information:"<<std::endl;
	std::cout<<"\t             Time: "<<simTime<<std::endl;

	// data.getDataSet()->GetFieldData()->PrintSelf(std::cout,vtkIndent(0));

	// Get domain and mesh information
	double 	physicalBounds[6]; // Physical space bounds
	int 		meshDimensions[3]; // Number of cells in each dimension
	data.getDataSet()->GetBounds(physicalBounds);
	data.getDataDimension(meshDimensions);

	std::cout<<"Mesh dimensions:"<<std::endl;
	std::cout<<meshDimensions[0]<<"\t"<<meshDimensions[1]<<"\t"<<meshDimensions[2]<<std::endl;


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
	std::cout<<"Ratio of true volume to simulated volume: "<<volumeCorrection<<std::endl;

	data.computeCellVolumes("cellVolumeTmp");
	data.storeData(volumeCorrection*data["cellVolumeTmp"],"cellVolume");

	std::cout<<"True volume/Corrected volume: "<<trueVolume/data.sum("cellVolume")<<std::endl;

	// Create data arrays of coordinates (full 2d/3d arrays of cell centers)
	data.createCoordinateDataArrays("r","phi","theta");




	/*====================================================================
		
		Create Gain Region Mask 
		This flags whether a cell is inside the gain region or 
		outside the gain region.

		The gain region is defined as the region bounded by:
		(start)	The point at which the lateral average of qenr 
				transitions from negative to positive
		(end)	The maximum of the shock radius and the point at 
				which the lateral average of qenr goes negative

		(1)	Average qenr
		(2)	Compute the radial positions where qenr goes from 
			(a)	negative to positive (qenrNP)
			(a)	positive to negative (qenrPN)
		(3)	For every ray
			(1)	For every cell, determine if the cell lies inside the
				gain region

				If it does, set the mask to 1.0
				Otherwise, set the mask to 0.0

	=======================================================================*/

	double qenrNP = 0.0;
	double qenrPN = 0.0;
	bool foundNP = false, foundPN = false;

	double gainRadius = 0.0;

	// Average qenr
	std::vector<std::string> shellVars;
	shellVars.push_back("qenr");

	vtkSmartPointer<vtkDataArray> radCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();
	ShellAveragePlaneOperator shellAvg(	radCoords,
										XAXIS,
										CS_SPHERE,
										shellVars	);
	shellAvg.process(data.getDataSet());

	vtkSmartPointer<vtkDataArray> qenrAvg = shellAvg.getDataArray("qenr");

	// Find qenrNP
	// For every bin in the interior (not on the edge) of the domain,
	// look at the previous and next bin values.
	// If the next is positive and the previous negative, interpolate
	// to find the coordinate where qenr becomes zero
	double yPrev = 0.0;
	double yNext = 1.0;
	for(int i=0; i<qenrAvg->GetNumberOfTuples()-1; i++)
	{
		yPrev = qenrAvg->GetTuple1(i);
		yNext = qenrAvg->GetTuple1(i+1);


		if( ((yPrev < 0.0) && (yNext > 0.0)) || ((yPrev > 0.0) && (yNext < 0.0)) )
		{
			// Get the volumetric centers of the previous and next bins
			double xPrevBndL = radCoords->GetTuple1(i);
			double xPrevBndR = radCoords->GetTuple1(i+1);
			double xNextBndL = radCoords->GetTuple1(i+1);
			double xNextBndR = radCoords->GetTuple1(i+2);
			double xPrev = pow(0.5*(xPrevBndL*xPrevBndL*xPrevBndL + xPrevBndR*xPrevBndR*xPrevBndR), 1.0/3.0);
			double xNext = pow(0.5*(xNextBndL*xNextBndL*xNextBndL + xNextBndR*xNextBndR*xNextBndR), 1.0/3.0);
			double slope = (yNext-yPrev)/(xNext-xPrev);
			double xCrossing = -yPrev/slope + xPrev;

			// qenr is going from negative to positive => qenrNP = xCrossing
			if( (yPrev < 0.0) && (yNext > 0.0) && !foundNP)
			{
				qenrNP = xCrossing;
				foundNP = true;
			}
			// qenr is going from positive to negative => qenrPN = xCrossing
			else if( (yPrev > 0.0) && (yNext < 0.0) && !foundPN )
			{
				qenrPN = xCrossing;
				foundPN = true;
			}
		}
	}

	// Set the gain radius as the transition from negative to positive
	gainRadius = qenrNP; 

	// Store
	dataList.addScalar(gainRadius,"gainradius");

	// Construct mask
	data.addVariable("maskGR");
	vtkSmartPointer<vtkDataArray> mask = data["maskGR"];
	int startingIndex[3];

	for(int p = 0; p<meshDimensions[1]; ++p)
	{
		for(int t = 0; t< (nDim>2 ? meshDimensions[2] : 1); ++t)
		{
			startingIndex[XAXIS]	= 0;
			startingIndex[YAXIS]	= p;
			startingIndex[ZAXIS]	= t;
			std::vector<int> cellIndices = data.getCellsAlongRay(XAXIS, startingIndex);

			bool foundFirst = false, foundLast = false;
			double shockPos = 0.0;

			for (std::vector<int>::reverse_iterator i = cellIndices.rbegin(); i != cellIndices.rend(); ++i)
			{
				
				bool isShocked = data["shock"]->GetTuple1(*i) > 0.5; // The data is either 0 or 1, so you need to pick something reasonable

				// If it's shocked, add the radial coordinate to the sum
				if(isShocked)
				{
					if(!foundFirst)
						foundFirst = true;
				}
				else
				{
					if(foundFirst && !foundLast)
					{
						foundLast = true;
						shockPos = data["r"]->GetTuple1(*i);
						// std::cout<<"shock pos = "<<shockPos<<std::endl;
					}
				}

				// Get cell coordinates
				double cellEdges[6];
				data.getDataSet()->GetCellBounds(*i,cellEdges);

				double xMin = cellEdges[0];
				double xMax = cellEdges[1];
				double maskValue = 0.0;

				if(xMax > qenrNP)
				{
					if(xMin < std::max(qenrPN, shockPos))
					{
						maskValue = 1.0;
					}
					else
					{
						maskValue = 0.0;
					}
				}
				else
				{
					maskValue = 0.0;
				}

				mask->SetTuple1(*i,maskValue);
			}
		}
	}



	/*====================================================================
		
		Create velocity masks as the intersection of the 
		gain region mask and velocity clips

	=======================================================================*/
	double velxMin = -5.e7;
	double velxMax = +1.e8;
	data.storeData( data["velx"]<velxMin, "maskTmpVelxNeg");
	data.storeData( data["velx"]>velxMax, "maskTmpVelxPos");
	data.storeData( data["maskTmpVelxNeg"]*data["maskGR"], "maskVelxNeg");
	data.storeData( data["maskTmpVelxPos"]*data["maskGR"], "maskVelxPos");



	/*====================================================================
		
		Apply the image erosion operator to kill small links in 
		the mask

	=======================================================================*/

	// Create rectangular structure element
	int nIStruct = 10;
	int nJStruct = 1;
	int seWidth[] = {(nIStruct-1)/2, (nJStruct-1)/2};
	bool** structElem;
	structElem = new bool*[nIStruct];
	for(int i=0; i<nIStruct; i++)
	{
		structElem[i] = new bool[nJStruct];
		memset(structElem[i],false,nJStruct);
	}

	for(int i=0; i<nIStruct; i++)
	{
		for(int j=0; j<nJStruct; j++)
		{
			structElem[i][j] = true;
		}
	}
	
	std::cout<<"Structure element size: "<<nIStruct<<"\t"<<nJStruct<<std::endl;
	std::cout<<"Structure element width: "<<seWidth[0]<<"\t"<<seWidth[1]<<std::endl;


	data.addVariable("maskVelxNegEroded");
	
	// Set the eroded mask equal to the original mask
	for(int i=0; i<data["maskVelxNeg"]->GetNumberOfTuples(); i++)
	{
		data["maskVelxNegEroded"]->SetTuple1(i,data["maskVelxNeg"]->GetTuple1(i));
	}

	// Loop over every "pixel" in the mesh
	// for(int i=0; i<meshDimensions[0]; i++)
	// {
	// 	for(int j=0; j<meshDimensions[1]; j++)
	// 	{
	int indexMax = 0;
	for(int i=seWidth[0]; i<meshDimensions[0]-seWidth[0]; i++)
	{
		for(int j=seWidth[1]; j<meshDimensions[1]-seWidth[1]; j++)
		{

			bool allPixelsTrue = true;

			// Loop over Subset of image inside the structure element
			for(int ii=0; ii<nIStruct; ii++)
			{
				for(int jj=0; jj<nJStruct; jj++)
				{

					// Handle conditions near boundaries
					// TODO: Implement periodic, reflecting, wall
					// 		and handle the cases therein
					// Right now, the outer for loops wont allow this

					int iIndex = i-seWidth[0]+ii;
					int jIndex = j-seWidth[1]+jj;
					int kIndex = 0;

					int index = Util::sub2ind(	iIndex,
												jIndex,
												kIndex,
												meshDimensions[0],
												meshDimensions[1]);

					// Determine if the mask is off AND if the structure element
					// requires this cell
					double maskValue = data["maskVelxNeg"]->GetTuple1(index);
					// std::cout<<"maskValue "<<maskValue<<std::endl;
					bool isMaskOn = maskValue > 0.5;
					// std::cout<<isMaskOn<<std::endl;
					if( !isMaskOn )//&& structElem[ii][jj])
					{
						// std::cout<<"In here!"<<std::endl;
						allPixelsTrue = false;
						break;
					}

				}
			}

			// if(allPixelsTrue)
			// {
			// 			std::cout<<"untouchable pixel i/j: "
			// 						<<i<<"\t"
			// 						<<j<<"\t"<<std::endl;
			// }

			// If not all pixels in the SE were true, set all of them to false
			// in the erosion mask
			if(!allPixelsTrue)
			{
				// Loop over Subset of image inside the structure element
				for(int ii=0; ii<nIStruct; ii++)
				{
					for(int jj=0; jj<nJStruct; jj++)
					{

						// std::cout<<"i/j/ii/jj: "
						// 			<<i<<"\t"
						// 			<<j<<"\t"
						// 			<<ii<<"\t"
						// 			<<jj<<"\t"<<std::endl;

						// Handle conditions near boundaries
						// TODO: Implement periodic, reflecting, wall
						// 		and handle the cases therein
						// Right now, the outer for loops wont allow this

						int iIndex = i-seWidth[0]+ii;
						int jIndex = j-seWidth[1]+jj;
						int kIndex = 0;

						int index = Util::sub2ind(	iIndex,
													jIndex,
													kIndex,
													meshDimensions[0],
													meshDimensions[1]);

						// std::cout<<"Index: "<<index<<std::endl;

						// Set mask values to false
						data["maskVelxNegEroded"]->SetTuple1(index,0.0);

					}
				}
			}


			
		}
	}




	/*	----------------------------------
			
			Create SignatureBerger 
			clusting object

		----------------------------------- */

	std::cout<<"Mesh dimensions: "<<meshDimensions[0]<<"\t"<<meshDimensions[1]<<"\t"<<meshDimensions[2]<<std::endl;

	Clustering::SignatureBerger cluster(nDim, meshDimensions);

	// Attach the data array
	cluster.setData(data["maskVelxNeg"]);

	// Set the minimum number of cells
	cluster.setMinCellsPerDim(5,5,3);

	cluster.setEfficiency(0.9);

	cluster.process();

	cluster.computeClusters();

	// Print information
	int nClusters = cluster.getNumberOfClusters();
	int nBlocksFromConstruction = cluster.getNumberOfLeafBlocks();
	int nBlocksFromGraph = cluster.getConnectivityGraph()->getNumberOfNodes();

	std::cout<<"Number of clusters: "<<nClusters<<std::endl;
	std::cout<<"Number of leaf blocks from construction: "<<nBlocksFromConstruction<<std::endl;
	std::cout<<"Number of leaf blocks from graph: "<<nBlocksFromGraph<<std::endl;

	// Write block information to file
	std::ofstream osBlockInfo;
	osBlockInfo.open("sigberger.blkinfo");
	cluster.printBlocks(osBlockInfo);
	osBlockInfo.close();

	// Write block connectivity to file
	std::ofstream osBlockConn;
	osBlockConn.open("sigberger.blkconn");
	cluster.printConnectivity(osBlockConn);
	osBlockConn.close();


	/*	----------------------------------
			
			Create grid variable which
			plots the block ID for every leaf

		----------------------------------- */

	data.addVariable("leafID");

	// Set the values in leafID
	std::list<Clustering::Block*>* blkList = cluster.getBlockList();
	std::list<Clustering::Block*>::iterator blk;

	for(blk = blkList->begin(); blk != blkList->end(); ++blk)
	{
		if(!(*blk)->isLeaf())
			continue;

		int blkID = (*blk)->getID();
		// std::cout<<"leaf block ID: "<<blkID<<std::endl;
		int bnds[3][2];
		(*blk)->getBounds(bnds);
		for(int k=bnds[KAXIS][LOW]; k<=bnds[KAXIS][HIGH]; k++)
		{
			for(int j=bnds[JAXIS][LOW]; j<=bnds[JAXIS][HIGH]; j++)
			{
				for(int i=bnds[IAXIS][LOW]; i<=bnds[IAXIS][HIGH]; i++)
				{
					int index = Util::sub2ind(i,j,k,meshDimensions[0],meshDimensions[1]);
					data["leafID"]->SetTuple1(index,blkID);
				}
			}
		}
	}




	vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
	writer->SetFileName("clustering.vtr");
#if VTK_MAJOR_VERSION <= 5
  	writer->SetInputConnection(data.getDataSet()->GetProducerPort());
#else
	writer->SetInputData(data.getDataSet());
#endif
	writer->Write();
}