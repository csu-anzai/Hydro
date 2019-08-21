#include <iostream>
#include <fstream>
#include <cmath>

#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkUniformGrid.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLImageDataWriter.h>

#include "../libsrc/Clustering/SignatureBerger/SignatureBerger.h"
#include "../libsrc/Utilities/LAVAUtil.h"

int main(void)
{

	int nDim = 3;
	int nI = 100;
	int nJ = 100;
	int nK = 100;

	double xMin = -1.0, xMax = 1.0;
	double yMin = -1.0, yMax = 1.0;
	double zMin = -1.0, zMax = 1.0;

	/*	----------------------------------
			Create example grid

			Create a 3d rectilinear grid from
			on [-1,1]x[-1,1]x[0,1] with
			nIxnJxnK cells
		----------------------------------- */
#if 0
	vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();

	// Set grid point dimension
	grid->SetDimensions(nI+1,nJ+1,nK+1);

	vtkSmartPointer<vtkDoubleArray> xArray = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> yArray = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> zArray = vtkSmartPointer<vtkDoubleArray>::New();
	
	// Create coordinates
	for(int i=0; i<nI+1; i++)
	{
		double x = xMin + (xMax-xMin)*double(i)/double(nI);
		xArray->InsertNextValue(x);
	}

	for(int i=0; i<nJ+1; i++)
	{
		double y = yMin + (yMax-yMin)*double(i)/double(nI);
		yArray->InsertNextValue(y);
	}

	for(int i=0; i<nK+1; i++)
	{
		double z = zMin + (zMax-zMin)*double(i)/double(nI);
		zArray->InsertNextValue(z);
	}

	// Set coordinates
	grid->SetXCoordinates(xArray);
	grid->SetYCoordinates(yArray);
	grid->SetZCoordinates(zArray);

	// Create data array
	vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
	tmpArray->SetName("mask");
	for(int i=0; i<nI*nJ*nK; i++)
	{
		tmpArray->InsertNextValue(0.0);
	}
	grid->GetCellData()->AddArray(tmpArray);
#endif	

	/*	----------------------------------
			Create a uniform grid
		----------------------------------- */

	vtkSmartPointer<vtkUniformGrid> grid = vtkSmartPointer<vtkUniformGrid>::New();

	// Set grid point dimension
	grid->SetDimensions(nI+1,nJ+1,nK+1);

	// Set spacing to 1
	grid->SetSpacing(double(1),double(1),double(1));// Create data array

	vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
	tmpArray->SetName("mask");
	for(int i=0; i<nI*nJ*nK; i++)
	{
		tmpArray->InsertNextValue(0.0);
	}
	grid->GetCellData()->AddArray(tmpArray);






	/*	----------------------------------
			Create a hollow sphere and 
			clip it with y==0
		----------------------------------- */
#if 0
	double radiusInner = 0.3;
	double radiusOuter = 0.5;
	double clipY = 0.0;

	for(int i=0; i<nI*nJ*nK; i++)
	{
		double cellEdges[6];
		grid->GetCellBounds(i,cellEdges);
		double xCell = 0.5*(cellEdges[0]+cellEdges[1]);
		double yCell = 0.5*(cellEdges[2]+cellEdges[3]);
		double zCell = 0.5*(cellEdges[4]+cellEdges[5]);
		double radCell = std::sqrt(xCell*xCell + yCell*yCell + zCell*zCell);

		if(radCell <= radiusOuter && radCell >= radiusInner && yCell<=clipY)
		{
			grid->GetCellData()->GetArray("mask")->SetTuple1(i,1.0);
		}
	}
#endif

	/*	----------------------------------
			Create a brick 
			5 	cells in the x
			2 	cells in the y
			1 	cells in the z
		----------------------------------- */
#if 0
	for(int i=0; i<nI; i++)
	{
		for(int j=0; j<nJ; j++)
		{
			for(int k=0; k<nK; k++)
			{
				int index = Util::sub2ind(i,j,k,nI,nJ);

				if(i < 5 && j<2 && k<1)
					grid->GetCellData()->GetArray("mask")->SetTuple1(index,1.0);

			}
		}
		
	}
#endif


	/*	----------------------------------
			Create 2 disjoint bricks
		----------------------------------- */
#if 0
	for(int i=0; i<nI; i++)
	{
		for(int j=0; j<nJ; j++)
		{
			for(int k=0; k<nK; k++)
			{
				int index = Util::sub2ind(i,j,k,nI,nJ);

				if(i < 4 && j<2 && k<1)
					grid->GetCellData()->GetArray("mask")->SetTuple1(index,1.0);

				if(i > 5 && j>3 && k>2)
					grid->GetCellData()->GetArray("mask")->SetTuple1(index,1.0);
			}
		}
		
	}
#endif

		/*	----------------------------------
			Create a bowl and sphere
		----------------------------------- */
#if 1
	double diameter = sqrt(nI*nI + nJ*nJ + nK*nK);
	double rInner = diameter/5.0;
	double rOuter = diameter/4.0;
	double rSphere = 15.0;
	double yClip = nJ/2.0;
	double origin[] = {80.,80.,80.};
	for(int i=0; i<nI; i++)
	{
		for(int j=0; j<nJ; j++)
		{
			for(int k=0; k<nK; k++)
			{
				int index = Util::sub2ind(i,j,k,nI,nJ);

				double rBowl = sqrt( (i-nI/2.0)*(i-nI/2.0) + (j-nJ/2.0)*(j-nJ/2.0) + (k-nK/2.0)*(k-nK/2.0));

				if(rBowl >= rInner && rBowl<=rOuter && j<yClip)
					grid->GetCellData()->GetArray("mask")->SetTuple1(index,1.0);



				double r = sqrt( (i-origin[0])*(i-origin[0]) + (j-origin[1])*(j-origin[1]) + (k-origin[2])*(k-origin[2]));

				if(r <= rSphere)
					grid->GetCellData()->GetArray("mask")->SetTuple1(index,1.0);


			}
		}
		
	}
#endif




	/*	----------------------------------
			
			Create SignatureBerger 
			clusting object

		----------------------------------- */

	int nCells[] = {nI, nJ, nK};
	Clustering::SignatureBerger cluster(nDim, nCells);

	// Attach the data array
	cluster.setData(grid->GetCellData()->GetArray("mask"));

	// Set the minimum number of cells
	cluster.setMinCellsPerDim(5,5,5);

	cluster.setEfficiency(0.9);

	cluster.process();

	cluster.computeClusters();

	// Print information
	// cluster.print(std::cout);
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
	vtkSmartPointer<vtkDoubleArray> blkIDArray  = vtkSmartPointer<vtkDoubleArray>::New();
	blkIDArray->SetName("leafID");
	for(int i=0; i<nI*nJ*nK; i++)
	{
		blkIDArray->InsertNextValue(0.0);
	}
	grid->GetCellData()->AddArray(blkIDArray);

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
					int index = Util::sub2ind(i,j,k,nI,nJ);
					grid->GetCellData()->GetArray("leafID")->SetTuple1(index,blkID);
				}
			}
		}
	}



	/*	----------------------------------
			Write grid to vtr file
		----------------------------------- */
#if 0
	vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
	writer->SetFileName("signatureberger.vtr");
	writer->SetInputData(grid);
	writer->Write();
#endif


	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName("signatureberger.vti");
#if VTK_MAJOR_VERSION <= 5
  	writer->SetInputConnection(grid->GetProducerPort());
#else
	writer->SetInputData(grid);
#endif
	writer->Write();

}

