#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>

#include <map>
#include <vector>
#include <sstream>
#include "ShellAveragePlaneOperator.h"

#include "Geometry/Geometry.h"
#include "Geometry/HEMesh/HEDataStructure.h"
#include "Geometry/Primitives/Vec3.h"
#include "includes/LAVAconstants.h"

#define INSIDE 1
#define OUTSIDE 2
#define INTERSECT 3


void ShellAveragePlaneOperator::printShells(std::ostream& out)
{	
	char delim = '\t';
	std::streamsize prec = 15;

	// Create header
	std::stringstream header;
	header<<"binMin"<<delim<<"binMax"<<delim<<"volume";
	for (std::vector<std::string>::iterator i = this->variableNames.begin(); i != this->variableNames.end(); ++i)
	{
		header<<delim<<*i;
	}
	out<<header.str()<<std::endl;


	// Write data lines
	for(int binId = 0; binId < this->nBins; binId++)
	{

		std::stringstream dataLine;

		// Set formatting
		dataLine.precision(prec);
		dataLine.setf(std::ios::scientific,std::ios::floatfield);

		// Get bin extents
		double binEdges[6];
		this->bins->GetCellBounds(binId,binEdges);

		// Write bin extents
		double startDist = sqrt(startPt[0]*startPt[0] + startPt[1]*startPt[1] + startPt[2]*startPt[2]);		
		dataLine<<startDist + binEdges[0]<<delim<<startDist + binEdges[1];

		// Write bin volume
		dataLine<<delim<<this->bins->GetCellData()->GetArray("volume")->GetTuple1(binId);

		// Loop over all variables for this bin
		for (std::vector<std::string>::iterator varIt = this->variableNames.begin(); varIt != this->variableNames.end(); ++varIt)
		{
			dataLine<<delim<<this->bins->GetCellData()->GetArray((*varIt).c_str())->GetTuple1(binId);
			
		}


		// write dataLine to ostream
		out<<dataLine.str()<<std::endl;

	}

}

void ShellAveragePlaneOperator::addVariable(const char* varName)
{
	vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
	
	tmpArray->SetName(varName);

	for(int i=0; i<nBins; i++)
	{
		tmpArray->InsertNextValue(0.0);
	}
		
	this->bins->GetCellData()->AddArray(tmpArray);
}


void ShellAveragePlaneOperator::generic_init(std::vector<std::string> varNames, bool initBinCoordinates = true)
{


	if(initBinCoordinates)
	{
		// Construct the bin cells
		this->bins = vtkSmartPointer<vtkRectilinearGrid>::New();
		vtkSmartPointer<vtkDoubleArray> binEdges  = vtkSmartPointer<vtkDoubleArray>::New();
		
		for(int i=0; i<(nBins+1); i++)
		{	
			double t = float(i)/float(nBins);
			double x = (endPt[0] - startPt[0])*t;
			double y = (endPt[1] - startPt[1])*t;
			double z = (endPt[2] - startPt[2])*t;
			double length = sqrt(x*x + y*y + z*z);
			binEdges->InsertNextValue(length);
		}
		
		
		this->bins->SetDimensions(nBins+1,1,1); // 1D rectilinear grid
		this->bins->SetXCoordinates(binEdges);
	}
	
	// Initialize the data arrays
	// vtkSmartPointer<vtkDoubleArray> avgArray  = vtkSmartPointer<vtkDoubleArray>::New();
	// vtkSmartPointer<vtkDoubleArray> volArray  = vtkSmartPointer<vtkDoubleArray>::New();
	// avgArray->SetName("average");
	// volArray->SetName("volume");
	// for(int i=0; i<nBins; i++)
	// {
	// 	avgArray->InsertNextValue(0.0);
	// 	volArray->InsertNextValue(0.0);
	// }
	// this->bins->GetCellData()->AddArray(avgArray);
	// this->bins->GetCellData()->AddArray(volArray);

	this->variableNames = varNames;
	std::vector<std::string> varNamesTmp(this->variableNames);
	varNamesTmp.push_back("volume");
	for(int v=0; v<varNamesTmp.size(); ++v)
	{
		// std::cout<<"Creating variable "<<varNamesTmp[v]<<std::endl;
		vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
		
		tmpArray->SetName(varNamesTmp[v].c_str());

		for(int i=0; i<nBins; i++)
		{
			tmpArray->InsertNextValue(0.0);
		}
			
		this->bins->GetCellData()->AddArray(tmpArray);
	}



	// Print diagnostics
	// std::cout << "There are " << this->bins->GetNumberOfPoints() 
	// 		<< " points." << std::endl;
	// std::cout << "There are " << this->bins->GetNumberOfCells() 
	// 		<< " cells." << std::endl;
	  
	// std::cout<<"The bin can hold a Cell->Bin map of "<<this->cellToBinMap.max_size()<<" elements."<<std::endl;  
	  
	// std::cout<<"[Leaving] ShellAveragePlaneOperator::ShellAveragePlaneOperator"<<std::endl;
}

ShellAveragePlaneOperator::ShellAveragePlaneOperator(int nBins, double* direction, double* startPt, double* endPt, int coordSys, std::vector<std::string> varNames)
{
	
	// std::cout<<"[Entering] ShellAveragePlaneOperator::ShellAveragePlaneOperator"<<std::endl;
	this->nBins = nBins;

	this->isAxisAligned = false;
	this->fixedAxis = -1;
	this->coordSys = coordSys;

	for(int i=0; i<3; i++)
	{
		this->direction[i]	= direction[i];
		this->startPt[i]		= startPt[i];
		this->endPt[i]		= endPt[i];
	}
	
	this->generic_init(varNames);

}

ShellAveragePlaneOperator::ShellAveragePlaneOperator(int nBins, int fixedDir, double startPtFixedAxis, double endPtFixedAxis, int coordSys, std::vector<std::string> varNames)
{
	
	// std::cout<<"[Entering] ShellAveragePlaneOperator::ShellAveragePlaneOperator"<<std::endl;
	this->nBins = nBins;

	this->isAxisAligned = true;
	this->fixedAxis = fixedDir;
	this->coordSys = coordSys;

	for(int i=0; i<3; i++)
	{
		this->direction[i]	= 0.0;
		this->startPt[i]		= 0.0;
		this->endPt[i]		= 0.0;
	}

	if(fixedDir==XAXIS) //"x"-axis
	{
		this->direction[XAXIS]	= 1.0;
		this->startPt[XAXIS]		= startPtFixedAxis;
		this->endPt[XAXIS]		= endPtFixedAxis;
	}
	else if(fixedDir==YAXIS) //"y"-axis
	{
		this->direction[YAXIS]	= 1.0;
		this->startPt[YAXIS]		= startPtFixedAxis;
		this->endPt[YAXIS]		= endPtFixedAxis;
	}
	else if(fixedDir==ZAXIS) //"z"-axis
	{
		this->direction[ZAXIS]	= 1.0;
		this->startPt[ZAXIS]		= startPtFixedAxis;
		this->endPt[ZAXIS]		= endPtFixedAxis;
	}
	else
	{
		std::cerr<<"ERROR [ShellAveragePlaneOperator] Fixed coordinate axis of "
				 <<fixedDir<<" is not "<<XAXIS<<", "<<YAXIS<<", or "<<ZAXIS<<"."<<std::endl;
		return;
	}
	
	this->generic_init(varNames);
}


ShellAveragePlaneOperator::ShellAveragePlaneOperator(vtkSmartPointer<vtkDataArray> coordinates, int fixedDir, int coordSys, std::vector<std::string> varNames)
{
	
	// std::cout<<"[Entering] ShellAveragePlaneOperator::ShellAveragePlaneOperator"<<std::endl;
	this->nBins = coordinates->GetNumberOfTuples()-1;

	this->isAxisAligned = true;
	this->fixedAxis = fixedDir;
	this->coordSys = coordSys;

	for(int i=0; i<3; i++)
	{
		this->direction[i]	= 0.0;
		this->startPt[i]		= 0.0;
		this->endPt[i]		= 0.0;
	}

	this->direction[fixedDir] 	= 1.0;
	this->startPt[fixedDir]		= coordinates->GetTuple1(0);
	this->endPt[fixedDir]			= coordinates->GetTuple1(nBins);

	// Assign the bin cells
	// First, shift all cells to the "left" so that the first cell starts at 0.0 
	// (zero distance from the starting point)
	vtkSmartPointer<vtkDoubleArray> newCoordinates = vtkSmartPointer<vtkDoubleArray>::New();
	for(int i=0; i<coordinates->GetNumberOfTuples(); i++)
	{
		double coord = coordinates->GetTuple1(i);
		double newCoord = coord - this->startPt[fixedDir];

		newCoordinates->InsertNextValue(newCoord);
	}

	this->bins = vtkSmartPointer<vtkRectilinearGrid>::New();
	this->bins->SetDimensions(nBins+1,1,1); // 1D rectilinear grid
	this->bins->SetXCoordinates(newCoordinates);

	int dim[3];
	this->bins->GetDimensions(dim);
	// std::cout<<"dims = "<<dim[0]<<"\t"<<dim[1]<<"\t"<<dim[2]<<std::endl;
	// std::cout<<"Start and end points: "<<this->startPt[fixedDir]<<"\t"<<this->endPt[fixedDir]<<std::endl;

	this->generic_init(varNames,false);
}


vtkSmartPointer<vtkDataSet> ShellAveragePlaneOperator::process(const vtkSmartPointer<vtkDataSet> ds)
{

	// std::cout<<"[Entering] ShellAveragePlaneOperator::process"<<std::endl;
	

	
	//Determine the dimensionality of the dataset
	int dsDim[3];
	int nDim=1;

	if(ds->GetDataObjectType() == VTK_RECTILINEAR_GRID)
	{
		vtkRectilinearGrid::SafeDownCast(ds)->GetDimensions(dsDim); // These are the dimensions of points, not cells. dimCells=dsDim-1
	}
	if(ds->GetDataObjectType() == VTK_UNIFORM_GRID)
	{
		vtkImageData::SafeDownCast(ds)->GetDimensions(dsDim); // These are the dimensions of points, not cells. dimCells=dsDim-1
	}
	
	// std::cout<<"dsDim = "<<dsDim[0]<<" "<<dsDim[1]<<" "<<dsDim[2]<<std::endl;
	
	if((dsDim[0]-1) < 2)
		nDim = 0;
	else if((dsDim[1]-1) < 2)
		nDim = 1;
	else if((dsDim[2]-1) < 2)
		nDim = 2;
	else
		nDim = 3;
	
	// std::cout<<"nDim = "<<nDim<<std::endl;
	
	if(nDim<2)
	{
		//error
		std::cerr<<"[ERROR] (ShellAveragePlaneOperator::process) Dimension of grid is less than 2! Averaging makes no sense."<<std::endl;
		//TODO: abort
	}



	vtkSmartPointer<vtkDataArray> binVolumeArray = this->bins->GetCellData()->GetArray("volume");
	// Find which bin the cell is in
	// Assumes bins are axis aligned for now
	// This loops over all bins and determines if the grid cell is outside, fully inside, or partially inside
	for(int b=0; b<this->bins->GetNumberOfCells();b++)
	// for(int b=0; b<1;b++)
	{
		double binEdges[6];
		this->bins->GetCellBounds(b,binEdges);
		
		
		double planeMinNormal[3], planeMaxNormal[3], planeMinPoint[3], planeMaxPoint[3], directionLength;
		directionLength = sqrt(this->direction[0]*this->direction[0] + this->direction[1]*this->direction[1] + this->direction[2]*this->direction[2]);
		planeMinNormal[0]	=	this->direction[0]/directionLength;
		planeMinNormal[1]	=	this->direction[1]/directionLength;
		planeMinNormal[2]	=	this->direction[2]/directionLength;
		
		planeMaxNormal[0] 	= 	-planeMinNormal[0];
		planeMaxNormal[1]	= 	-planeMinNormal[1];
		planeMaxNormal[2]	=	-planeMinNormal[2];

		planeMinPoint[0]		= 	binEdges[0]*this->direction[0]/directionLength + startPt[0];
		planeMinPoint[1]		= 	binEdges[0]*this->direction[1]/directionLength + startPt[1];
		planeMinPoint[2]		= 	binEdges[0]*this->direction[2]/directionLength + startPt[2];
		
		planeMaxPoint[0]		= 	binEdges[1]*this->direction[0]/directionLength + startPt[0];
		planeMaxPoint[1]		= 	binEdges[1]*this->direction[1]/directionLength + startPt[1];
		planeMaxPoint[2]		= 	binEdges[1]*this->direction[2]/directionLength + startPt[2];

		//std::cout<<"Plane normal: "<<planeNormal[0]<<"  "<<planeNormal[1]<<"  "<<planeNormal[2]<<std::endl;
		// std::cout<<"Plane points:  "<<std::endl;
		// std::cout<<"\t"<<planeMinPoint[0]<<"  "<<planeMaxPoint[0]<<std::endl;
		// std::cout<<"\t"<<planeMinPoint[1]<<"  "<<planeMaxPoint[1]<<std::endl;
		// std::cout<<"\t"<<planeMinPoint[2]<<"  "<<planeMaxPoint[2]<<std::endl;

		
		//Loop over all cells in the dataset
		for(int i=0; i<ds->GetNumberOfCells(); i++)
		// for(int i=0; i<3; i++)
		{
			// Diagnostics
			double cellEdges[6];
			ds->GetCellBounds(i,cellEdges);

			// Shrink the cell slightly to get around floating point geometry problems.
			double cellEdgesTmp[6];
			cellEdgesTmp[0]	= (1.0 + ERRMAX) * cellEdges[0];
			cellEdgesTmp[2]	= (1.0 + ERRMAX) * cellEdges[2];
			cellEdgesTmp[4]	= (1.0 + ERRMAX) * cellEdges[4];
			cellEdgesTmp[1]	= (1.0 - ERRMAX) * cellEdges[1];
			cellEdgesTmp[3]	= (1.0 - ERRMAX) * cellEdges[3];
			cellEdgesTmp[5]	= (1.0 - ERRMAX) * cellEdges[5];
		

			int InOutIntersectMin = ShellAveragePlaneOperator::cellIntersectsBox(planeMinNormal, planeMinPoint, cellEdgesTmp);
			int InOutIntersectMax = ShellAveragePlaneOperator::cellIntersectsBox(planeMaxNormal, planeMaxPoint, cellEdgesTmp);
			
			
			// There are 4 cases
			// 1) cell registers as In for both planes. The cell is completely enclosed in the bin
			// 2) cell registers as In-Out or Out-In. The cell is complete outside the bin
			// 3) cell registers as Intersected for one plane and In for another. Need to weight by the volume clipped by the intersecting plane.
			// 4) cell registers as Intersected for both planes. Need to weight by the volume between the planes
//http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
//http://www.cs.berkeley.edu/~jfc/mirtich/massProps.html
			// Case 1 - Cell completely inside

			bool includeCell = false;
			double cellVolume = 0.0;				
			// If we have axis-aligned bins, we do not have to resort to complicated polytope integrals
			if(this->isAxisAligned)
			{
				double cellEdgesClipped[6];
				cellEdgesClipped[0]	= cellEdges[0];
				cellEdgesClipped[1]	= cellEdges[1];
				cellEdgesClipped[2]	= cellEdges[2];
				cellEdgesClipped[3]	= cellEdges[3];
				cellEdgesClipped[4]	= cellEdges[4];
				cellEdgesClipped[5]	= cellEdges[5];

				// Case 1 Cell is fully inside bin
				// we do not need to do anything except flag that we should compute properties
				// for this cell
				if( (InOutIntersectMin==INSIDE) && (InOutIntersectMax==INSIDE))
				{
					includeCell = true;
				}
				// Case 3&4 Cell is intersected
				// Need to trim the edges accordingly. Since this is axis aligned,
				// we need to only change the values in cellEdgesClipped
				else if( (InOutIntersectMin==INTERSECT) || (InOutIntersectMax==INTERSECT) )
				{
					includeCell = true;

					// If the minimum plane of the bin clips the cell
					// we need to set the minimum value for the fixed axis of cellEdgesClipped
					if(InOutIntersectMin == INTERSECT)
					{
						cellEdgesClipped[this->fixedAxis] = planeMinPoint[this->fixedAxis];
					}

					// If the maximum plane of the bin clips the cell
					// we need to set the maximum value for the fixed axis of cellEdgesClipped
					if(InOutIntersectMax == INTERSECT)
					{
						cellEdgesClipped[this->fixedAxis+1] = planeMaxPoint[this->fixedAxis];
					}
				}

				// If the cell is inside or clipped, 
				// we need to actually compute its contribution to the bin
				if(includeCell)
				{

					// Compute volume (or area in 2d) of cell
					if(this->coordSys == CS_CART)
					{
						cellVolume = 1.0;
						for(int j=0; j<nDim; j++)
						{
							cellVolume *= cellEdgesClipped[2*j+1]-cellEdgesClipped[2*j];
						}
					}
					else if(this->coordSys == CS_CYL)
					{
						cellVolume = 1.0;
						cellVolume *= cellEdgesClipped[1]*cellEdgesClipped[1] - cellEdgesClipped[0]*cellEdgesClipped[0]; // Radial component
						cellVolume *= 0.5*(cellEdgesClipped[3]-cellEdgesClipped[2]); // Theta component
						if(nDim>2)
						{
							cellVolume *= cellEdgesClipped[5]-cellEdgesClipped[4]; // height component
						}
					}
					else if(this->coordSys == CS_SPHERE)
					{
						cellVolume = 1.0;
						cellVolume *= (1.0/3.0)*(	cellEdgesClipped[1]*cellEdgesClipped[1]*cellEdgesClipped[1] - 
												cellEdgesClipped[0]*cellEdgesClipped[0]*cellEdgesClipped[0]); // Radial component
						
						if(nDim>1)
						{
							cellVolume *= -cos(cellEdgesClipped[3]) + cos(cellEdgesClipped[2]); // Phi component
						}
						else
						{
							cellVolume *= 2.0;
						}
						if(nDim>2)
						{
							cellVolume *= cellEdgesClipped[5] - cellEdgesClipped[4]; // Theta component
						}
						else
						{
							cellVolume *= 2.0*vtkMath::Pi();
						}
					}
					
				}

			}
			// Otherwise, we need to be a bit more generic in our approach
			else
			{
#if 0
				if( (InOutIntersectMin==INSIDE) && (InOutIntersectMax==INSIDE))
				{
					// Compute volume of cell (or area in 2d)
					double volume;				
					volume = 1.0;
					for(int j=0; j<nDim; j++)
					{
						volume *= cellEdges[2*j+1]-cellEdges[2*j];
						//std::cout<<j<<"volume : "<< volume<<std::endl;
					}
					
					// std::cout<<"volume : "<< volume<<std::endl;
					
					double cellValue, binValue, localContribution;

					cellValue = cellDataArray->GetTuple1(i);
					binValue  = binDataArray->GetTuple1(b);
					// std::cout<<"Bin value: "<<binValue<<std::endl;
					// localContribution = cellValue*volume;
					localContribution = volume;
					
					binDataArray->SetTuple1(b,binValue + localContribution);

					// Add to map
					this->cellToBinMap.insert(std::pair<int,int>(i,b));
					this->binToCellMap.insert(std::pair<int,int>(b,i));

				}
				// Case 3&4 Cell is intersected
				// Need to make a polytope, and cut it with one or more cells
				else if( (InOutIntersectMin==INTERSECT) || (InOutIntersectMax==INTERSECT) )
				{
					std::cout<<"Should create a new Polytope"<<std::endl;

					Polytope* cellPoly;
					if(nDim==2)
					{
						std::cout<<"\tCreating rectangle"<<std::endl;
						cellPoly = Polytope::generateRectangle(cellEdges);
					}
					else if(nDim==3)
					{
						std::cout<<"\tCreating brick"<<std::endl;
						cellPoly = Polytope::generateBrick(cellEdges);
					}


					if(InOutIntersectMin == INTERSECT)
					{
						// Remember to flip normal
						std::cout<<"\tClipping polytope with Min plane"<<std::endl;

						Vec3 planePoint(planeMinPoint[0],planeMinPoint[1],planeMinPoint[2]);
						Vec3 planeNormal(-planeMinNormal[0],-planeMinNormal[1],-planeMinNormal[2]);
						cellPoly->clipWithPlane(planePoint, planeNormal);


					}

					if(InOutIntersectMax == INTERSECT)
					{
						// Remember to flip normal
						std::cout<<"\tClipping polytope with Max plane"<<std::endl;
						Vec3 planePoint(planeMaxPoint[0],planeMaxPoint[1],planeMaxPoint[2]);
						Vec3 planeNormal(-planeMaxNormal[0],-planeMaxNormal[1],-planeMaxNormal[2]);
						cellPoly->clipWithPlane(planePoint, planeNormal);
					}

					// Integrating polytope
					std::cout<<"\tIntegrating polytope"<<std::endl;

					// Add to map
					this->cellToBinMap.insert(std::pair<int,int>(i,b));
					this->binToCellMap.insert(std::pair<int,int>(b,i));

				}
				// Case 2
				else
				{
					//std::cout<<"Outside of bin"<<std::endl;
				}
#endif
			}

			if(includeCell)
			{
				double cellValue, binValue, binVolume, localContribution;

				// Update bin variable data
				for(int v=0; v < this->variableNames.size(); v++)
				{
					std::string var = this->variableNames[v];

					// std::cout<<"var: "<<var<<std::endl;

					vtkSmartPointer<vtkDataArray> cellDataArray = ds->GetCellData()->GetArray(var.c_str());
					vtkSmartPointer<vtkDataArray> binDataArray = this->bins->GetCellData()->GetArray(var.c_str());

					// std::cout<<"dataset pointer: "<<ds<<std::endl;
					// std::cout<<"celldata pointer: "<<cellDataArray<<std::endl;
					// std::cout<<i<<std::endl;

					cellValue	= cellDataArray->GetTuple1(i);
					binValue	= binDataArray->GetTuple1(b);
					localContribution = cellValue*cellVolume;
				
					// Update bin
					binDataArray->SetTuple1(b,binValue + localContribution);
				}


				// Update bin volume data
				binVolume 	= binVolumeArray->GetTuple1(b);
				binVolumeArray->SetTuple1(b, binVolume + cellVolume);

#if 1
				// update cell volume fraction
				// Right now, this adds the raw volume
				// we'll normalize later to obtain the true volume fraction
				this->setVolumeFraction(i, b, cellVolume);

				// Add to map
				this->cellToBinMap.insert(std::pair<int,int>(i,b));
				this->binToCellMap.insert(std::pair<int,int>(b,i));
#endif
			}
			
		}
	
	}
	
	// Divide the bin values by the bin volume to obtain the actual average
	double totalVolume = 0.0;
	for(int v=0; v < this->variableNames.size(); v++)
	{
		std::string var = this->variableNames[v];
		vtkSmartPointer<vtkDataArray> binDataArray = this->bins->GetCellData()->GetArray(var.c_str());

		for(int b=0; b<this->bins->GetNumberOfCells();b++)
		{
			double binValue, binVolume;

			binValue		= binDataArray->GetTuple1(b);
			binVolume	= binVolumeArray->GetTuple1(b);
			// std::cout<<"Bin volume: "<<binVolume<<std::endl;
			if(binVolume > ERRMAX)
			{
				binDataArray->SetTuple1(b,binValue/binVolume);
			}

			if(v==0)
				totalVolume += binVolume;
		}

	}
	// std::cout<<"Total volume: "<<totalVolume<<std::endl;
#if 1
	// Update the raw cell volumes stored in the volumeFraction variable to be true volume fractions
	for(int cellInd=0; cellInd<ds->GetNumberOfCells(); cellInd++)
	{
		double totalCellVolume = 0.0;
		std::map<int,double>::iterator it;

		for(it = this->cellVolumeFraction[cellInd].begin(); it != this->cellVolumeFraction[cellInd].end(); ++it)
		{
			totalCellVolume += (*it).second;
		}


		if(totalCellVolume < 1.e-10)
		{
			std::cout<<"[ShellAveragePlaneOperator::process] WARNING: totalCellVolume is tiny!\n";
			std::cout.setf(std::ios_base::scientific,std::ios_base::floatfield);
			double cellEdges[6];
			ds->GetCellBounds(cellInd,cellEdges);

			std::cout<<"\tCell Bounding Box:\n";
			std::cout<<"\t\t"<<cellEdges[0]<<"\t"<<cellEdges[1]<<std::endl;
			std::cout<<"\t\t"<<cellEdges[2]<<"\t"<<cellEdges[3]<<std::endl;
			std::cout<<"\t\t"<<cellEdges[4]<<"\t"<<cellEdges[5]<<std::endl;

			std::cout<<"\ttotalCellVolume = "<<totalCellVolume<<std::endl;
			std::cout<<"\t\tRadial component = "<<	cellEdges[1]*cellEdges[1]*cellEdges[1] - 
													cellEdges[0]*cellEdges[0]*cellEdges[0]<<std::endl;
			std::cout<<"\t\tPhi component = "<<-cos(cellEdges[3]) + cos(cellEdges[2])<<std::endl;

			std::cout.unsetf(std::ios_base::scientific);
		}

		for(it = this->cellVolumeFraction[cellInd].begin(); it != this->cellVolumeFraction[cellInd].end(); ++it)
		{
			this->cellVolumeFraction[cellInd][(*it).first] = (*it).second/totalCellVolume;
			// std::cout<<"volume fraction: "<<this->cellVolumeFraction[i][(*it).first]<<std::endl;
		}

	}
#endif

	// std::cout<<"[Leaving] ShellAveragePlaneOperator::process"<<std::endl;
	return this->bins;

}

int ShellAveragePlaneOperator::cellIntersectsBox(double planeNormal[3], double planePoint[3], double cellEdges[6])
{
	// This assumes the grid cell is axis aligned (an AAB - Axis Aligned Box)
	// Ref: http://zach.in.tu-clausthal.de/teaching/cg_literatur/lighthouse3d_view_frustum_culling/index.html
	double pVert[3], nVert[3];
	
	// Determine the positive vertex (pVert)
	pVert[0] = cellEdges[0]; //xmin
	pVert[1] = cellEdges[2]; //ymin
	pVert[2] = cellEdges[4]; //zmin
	
	if(planeNormal[0] >= 0)
		pVert[0] = cellEdges[1];
	if(planeNormal[1] >= 0)
		pVert[1] = cellEdges[3];
	if(planeNormal[2] >= 0)
		pVert[2] = cellEdges[5];
	
	// Determine the negative vertex (pVert)
	nVert[0] = cellEdges[1]; //xmax
	nVert[1] = cellEdges[3]; //ymax
	nVert[2] = cellEdges[5]; //zmax
	
	if(planeNormal[0] >= 0)
		nVert[0] = cellEdges[0];
	if(planeNormal[1] >= 0)
		nVert[1] = cellEdges[2];
	if(planeNormal[2] >= 0)
		nVert[2] = cellEdges[4];
				
	
	//std::cout<<"pVert = "<<pVert[0]<<"  "<<pVert[1]<<"  "<<pVert[2]<<std::endl;
	//std::cout<<"nVert = "<<nVert[0]<<"  "<<nVert[1]<<"  "<<nVert[2]<<std::endl;
	
	int InOutIntersect = INSIDE; //inside = 1 | outside = 2 | intersect = 3
	
	double distanceFromPlane;
	distanceFromPlane = planeNormal[0]*(pVert[0]-planePoint[0]) + planeNormal[1]*(pVert[1]-planePoint[1]) + planeNormal[2]*(pVert[2]-planePoint[2]);
	// std::cout<<"distanceFromPlane positive pt: "<<distanceFromPlane<<std::endl;
	if(distanceFromPlane<0)// && abs(distanceFromPlane)>ERRMAX)
	{
		InOutIntersect = OUTSIDE; // outside
	}
	else
	{
		distanceFromPlane = planeNormal[0]*(nVert[0]-planePoint[0]) + planeNormal[1]*(nVert[1]-planePoint[1]) + planeNormal[2]*(nVert[2]-planePoint[2]);
		// std::cout<<"distanceFromPlane negative pt: "<<distanceFromPlane<<std::endl;

		if(distanceFromPlane<0)// || abs(distanceFromPlane)>ERRMAX)
		{
			InOutIntersect = INTERSECT; // intersect
		}
	}
	
	// if(InOutIntersect==INSIDE)
	// {
	// 	std::cout<<"\tCell is inside"<<std::endl;
	// }
	// if(InOutIntersect==OUTSIDE)
	// {
	// 	std::cout<<"\tCell is outside"<<std::endl;
	// }
	// if(InOutIntersect==INTERSECT)
	// {
	// 	std::cout<<"\tCell is intersecting"<<std::endl;
	// }
	
	return InOutIntersect;
	
}
















