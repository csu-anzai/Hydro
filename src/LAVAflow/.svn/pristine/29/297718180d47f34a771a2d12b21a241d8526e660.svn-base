#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkMath.h>
#include <map>
#include <sstream>
#include <string>
#include "SurfaceAveragePlaneOperator.h"

#include "Geometry/Geometry.h"
#include "Geometry/HEMesh/HEDataStructure.h"
#include "Geometry/Primitives/Vec3.h"
#include "includes/LAVAconstants.h"

#define INSIDE 1
#define OUTSIDE 2
#define INTERSECT 3



void SurfaceAveragePlaneOperator::printSurfaces(std::ostream& out)
{	
	char delim = '\t';
	std::streamsize prec = 15;
	
	// Create header
	std::stringstream header;
	header<<"surfacePt"<<delim<<"area";
	// header<<"surfacePt";
	for (std::vector<std::string>::iterator i = this->variableNames.begin(); i != this->variableNames.end(); ++i)
	{
		header<<delim<<*i;
		// std::cout<<"*i="<<*i<<std::endl;
	}
	out<<header.str()<<std::endl;

	// Write data lines
	for(int surfaceId = 0; surfaceId < this->nSurfaces; surfaceId++)
	{

		std::stringstream dataLine;

		// Set formatting
		dataLine.precision(prec);
		dataLine.setf(std::ios::scientific,std::ios::floatfield);

		// Get surface extents		
		double surfacePt[3];
		this->surfaces->GetPoint(surfaceId,surfacePt);

		// Write surface extents
		double startDist = sqrt(startPt[0]*startPt[0] + startPt[1]*startPt[1] + startPt[2]*startPt[2]);
		dataLine<<startDist + surfacePt[0];

		// Write surface area
		dataLine<<delim<<this->surfaces->GetPointData()->GetArray("area")->GetTuple1(surfaceId);

		// Loop over all variables for this surface
		for (std::vector<std::string>::iterator varIt = this->variableNames.begin(); varIt != this->variableNames.end(); ++varIt)
		{
			// std::cout<<"Writing variable "<<*varIt<<" for surface#"<<surfaceId<<std::endl;
			dataLine<<delim<<this->surfaces->GetPointData()->GetArray((*varIt).c_str())->GetTuple1(surfaceId);
			
		}


		// write dataLine to ostream
		out<<dataLine.str()<<std::endl;

	}
}

void SurfaceAveragePlaneOperator::generic_init(std::vector<std::string> varNames, bool initShellCoordinates = false)
{
	if(initShellCoordinates)
	{
		// Construct the surface cells
		this->surfaces = vtkSmartPointer<vtkRectilinearGrid>::New();
		vtkSmartPointer<vtkDoubleArray> surfaceEdges  = vtkSmartPointer<vtkDoubleArray>::New();
		
		for(int i=0; i<nSurfaces; i++)
		{	
			double t = float(i)/float(nSurfaces-1);
			double x = (endPt[0] - startPt[0])*t;
			double y = (endPt[1] - startPt[1])*t;
			double z = (endPt[2] - startPt[2])*t;
			double length = sqrt(x*x + y*y + z*z);
			surfaceEdges->InsertNextValue(length);
		}
		
		
		this->surfaces->SetDimensions(nSurfaces,1,1); // 1D rectilinear grid
		this->surfaces->SetXCoordinates(surfaceEdges);
	}

	// Initialize the data
	// this->variableNames = varNames;
	// std::vector<std::string> varNamesTmp(this->variableNames);
	// varNamesTmp.push_back("area");
	// for(int v=0; v<varNamesTmp.size(); ++v)
	this->addVariable("area",0.0,true);
	for(int v=0; v<varNames.size(); ++v)
	{
		// // std::cout<<"Creating variable "<<varNamesTmp[v]<<std::endl;
		// vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
		
		// tmpArray->SetName(varNamesTmp[v].c_str());

		// for(int i=0; i<nSurfaces; i++)
		// {
		// 	tmpArray->InsertNextValue(0.0);
		// }
			
		// this->surfaces->GetPointData()->AddArray(tmpArray);
		this->addVariable(varNames[v].c_str(),-1.0,false);
	}

	
	// Print diagnostics
	// std::cout << "There are " << this->surfaces->GetNumberOfPoints() 
	// 		<< " points." << std::endl;
	// std::cout << "There are " << this->surfaces->GetNumberOfCells() 
	// 		<< " cells." << std::endl;

	// vtkSmartPointer<vtkDataArray> data = this->surfaces->GetPointData()->GetArray("average");
	// for(int i=0; i<this->surfaces->GetNumberOfPoints();i++)
	// {
	// 	double point[3];
	// 	this->surfaces->GetPoint(i,point);
	// 	std::cout<<"Point "<<i<<" is located at by "<<std::endl;
	// 	std::cout<<"\t"<<point[0]<<" , "<<point[1]<<" , "<<point[2]<<std::endl;
	//     std::cout<<"Point "<<i<<" has value: "<<data->GetTuple1(i)<<std::endl;
	// }
  
	  
	// std::cout<<"The surface can hold a Cell->Surface map of "<<this->cellToSurfaceMap.max_size()<<" elements."<<std::endl;  
	  
	// std::cout<<"[Leaving] SurfaceAveragePlaneOperator::SurfaceAveragePlaneOperator"<<std::endl;
}

SurfaceAveragePlaneOperator::SurfaceAveragePlaneOperator(int nSurfaces, double* direction, double* startPt, double* endPt, int coordSys, std::vector<std::string> varNames)
{
	
	// std::cout<<"[Entering] SurfaceAveragePlaneOperator::SurfaceAveragePlaneOperator"<<std::endl;
	this->nSurfaces = nSurfaces;

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

SurfaceAveragePlaneOperator::SurfaceAveragePlaneOperator(int nSurfaces, int fixedDir, double startPtFixedAxis, double endPtFixedAxis, int coordSys, std::vector<std::string> varNames)
{
	
	// std::cout<<"[Entering] SurfaceAveragePlaneOperator::SurfaceAveragePlaneOperator"<<std::endl;
	this->nSurfaces = nSurfaces;

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
		std::cerr<<"ERROR [SurfaceAveragePlaneOperator] Fixed coordinate axis of "
				 <<fixedDir<<" is not "<<XAXIS<<", "<<YAXIS<<", or "<<ZAXIS<<"."<<std::endl;
		return;
	}
	
	this->generic_init(varNames);
}

SurfaceAveragePlaneOperator::SurfaceAveragePlaneOperator(vtkSmartPointer<vtkDataArray> coordinates, int fixedDir, int coordSys, std::vector<std::string> varNames)
{
	
	// std::cout<<"[Entering] SurfaceAveragePlaneOperator::SurfaceAveragePlaneOperator"<<std::endl;
	this->nSurfaces = coordinates->GetNumberOfTuples();

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
	this->endPt[fixedDir]			= coordinates->GetTuple1(nSurfaces-1);
	
	// Assign the surface points
	// First, shift all points to the "left" so that the first cell starts at 0.0 
	// (zero distance from the starting point)
	vtkSmartPointer<vtkDoubleArray> newCoordinates = vtkSmartPointer<vtkDoubleArray>::New();
	for(int i=0; i<coordinates->GetNumberOfTuples(); i++)
	{
		double coord = coordinates->GetTuple1(i);
		double newCoord = coord - this->startPt[fixedDir];

		newCoordinates->InsertNextValue(newCoord);
	}

	this->surfaces = vtkSmartPointer<vtkRectilinearGrid>::New();
	this->surfaces->SetDimensions(nSurfaces,1,1); // 1D rectilinear grid
	this->surfaces->SetXCoordinates(newCoordinates);

	this->generic_init(varNames,false);
}


vtkSmartPointer<vtkDataSet> SurfaceAveragePlaneOperator::process(const vtkSmartPointer<vtkDataSet> ds)
{

	// std::cout<<"[Entering] SurfaceAveragePlaneOperator::process"<<std::endl;
	


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
		std::cerr<<"[ERROR] (SurfaceAveragePlaneOperator::process) Dimension of grid is less than 2! Averaging makes no sense."<<std::endl;
		//TODO: abort
	}

	vtkSmartPointer<vtkDataArray> surfaceAreaArray = this->surfaces->GetPointData()->GetArray("area");
	// Find which surface the cell is in
	// Assumes surfaces are axis aligned for now
	// This loops over all surfaces and determines if the grid cell is outside, fully inside, or partially inside
	for(int b=0; b<this->surfaces->GetNumberOfPoints();b++)
	// for(int b=0; b<1;b++)
	{
		double surfacePt[3];
		this->surfaces->GetPoint(b,surfacePt);
		
		
		double planeNormal[3], planePoint[3], directionLength;
		directionLength = sqrt(this->direction[0]*this->direction[0] + this->direction[1]*this->direction[1] + this->direction[2]*this->direction[2]);
		planeNormal[0]	=	this->direction[0]/directionLength;
		planeNormal[1]	=	this->direction[1]/directionLength;
		planeNormal[2]	=	this->direction[2]/directionLength;
				
		planePoint[0]		= 	surfacePt[0]*this->direction[0]/directionLength + startPt[0];
		planePoint[1]		= 	surfacePt[1]*this->direction[1]/directionLength + startPt[1];
		planePoint[2]		= 	surfacePt[2]*this->direction[2]/directionLength + startPt[2];

		//std::cout<<"Plane normal: "<<planeNormal[0]<<"  "<<planeNormal[1]<<"  "<<planeNormal[2]<<std::endl;
		// std::cout<<"Plane points:  "<<std::endl;
		// std::cout<<"\t"<<planePoint[0]<<"  "<<planePoint[1]<<"  "<<planePoint[2]<<std::endl;

		
		//Loop over all cells in the dataset
		for(int i=0; i<ds->GetNumberOfCells(); i++)
		// for(int i=0; i<3; i++)
		{
			// Diagnostics
			double cellEdges[6];
			ds->GetCellBounds(i,cellEdges);

			// Shrink the cell slightly to get around floating point geometry problems.
			double cellEdgesTmp[6];
			cellEdgesTmp[0]	= cellEdges[0]	+ ERRMAX;
			cellEdgesTmp[2]	= cellEdges[2]	+ ERRMAX;
			cellEdgesTmp[4]	= cellEdges[4]	+ ERRMAX;
			cellEdgesTmp[1]	= cellEdges[1]	- ERRMAX;
			cellEdgesTmp[3]	= cellEdges[3]	- ERRMAX;
			cellEdgesTmp[5]	= cellEdges[5]	- ERRMAX;
		

			int InOutIntersect = SurfaceAveragePlaneOperator::cellIntersectsBox(planeNormal, planePoint, cellEdgesTmp);
			
			
			// There are 2 cases
			// 1) cell is intersected by plane
			// 2) cell is not intersected by plane (could be in or out)

			bool includeCell = false;
			double cellArea = 0.0;				
			// If we have axis-aligned surfaces, we do not have to resort to complicated polytope integrals
			if(this->isAxisAligned)
			{
				//Case 1 - Cell is intersected by the plane
				if( InOutIntersect==INTERSECT )
				{
					includeCell = true;
				}

				// If the cell is intersected, we need to compute the area and area-weighted contribution
				if(includeCell)
				{

					// Compute area (or arc length in 1d) of surface
					if(this->coordSys == CS_CART)
					{
						cellArea = 1.0;
						for(int j=0; j<nDim; j++)
						{
							if(j!=this->fixedAxis)
							{
								cellArea *= cellEdges[2*j+1]-cellEdges[2*j];
							}
						}
					}
					else if(this->coordSys == CS_CYL)
					{
						cellArea = 1.0;
						if(this->fixedAxis==XAXIS)
						{
							cellArea *= planePoint[0]; // Fixed radius
							cellArea *= cellEdges[3]-cellEdges[2]; // Theta contribution
							if(nDim>2)
							{
								cellArea *= cellEdges[5]-cellEdges[4]; // height contribution
							}
						}
						else if(this->fixedAxis==YAXIS)
						{
							cellArea *= cellEdges[1]-cellEdges[0]; // R contribution
							if(nDim>2)
							{
								cellArea *= cellEdges[5]-cellEdges[4]; // height contribution
							}
						}
						else if(this->fixedAxis==ZAXIS)
						{
							cellArea *= cellEdges[1]*cellEdges[1] - cellEdges[0]*cellEdges[0]; // Radial component
							cellArea *= 0.5*(cellEdges[3]-cellEdges[2]); // Theta component						
						}
						else
						{
							cellArea = 0.0;
							std::cerr<<"God damnit"<<std::endl;
						}


					}
					else if(this->coordSys == CS_SPHERE)
					{
						cellArea = 1.0;
						if(this->fixedAxis==XAXIS)
						{
							cellArea *= planePoint[0]*planePoint[0]; // Fixed radius
							cellArea *= -cos(cellEdges[3]) + cos(cellEdges[2]); // Phi contribution
							if(nDim>2)
							{
								cellArea *= cellEdges[5]-cellEdges[4]; // theta contribution
							}
							else
							{
								cellArea *= 2.0*vtkMath::Pi();
							}
						}
						else if(this->fixedAxis==YAXIS)
						{
							std::cout<<"WARNING : Using spherical area with fixedAxis = YAXIS (Phi). This may not be correct."<<std::endl;
							cellArea *= abs(cos(planePoint[1])); // Phi contribution
							cellArea *= cellEdges[1]*cellEdges[1]-cellEdges[0]*cellEdges[0]; // R contribution
							if(nDim>2)
							{
								cellArea *= cellEdges[5]-cellEdges[4]; // theta contribution
							}
						}
						else if(this->fixedAxis==ZAXIS)
						{
							std::cout<<"WARNING : Using spherical area with fixedAxis = ZAXIS (Theta). This may not be correct."<<std::endl;
							cellArea *= cellEdges[1]*cellEdges[1] - cellEdges[0]*cellEdges[0]; // Radial component
							cellArea *= 0.5*(cellEdges[3]-cellEdges[2]); // Theta component						
						}
						else
						{
							cellArea = 0.0;
							std::cerr<<"God damnit"<<std::endl;
						}
					}

					// double cellValue, surfaceValue, surfaceArea, localContribution;

					// cellValue	= cellDataArray->GetTuple1(i);
					// surfaceValue	= surfaceDataArray->GetTuple1(b);
					// surfaceArea 	= surfaceAreaArray->GetTuple1(b);
					// localContribution = cellValue*cellArea;
					
					// // Update surface
					// surfaceDataArray->SetTuple1(b,	surfaceValue + localContribution);
					// surfaceAreaArray->SetTuple1(b,	surfaceArea + cellArea);

					// // Add to map
					// this->cellToSurfaceMap.insert(std::pair<int,int>(i,b));
					// this->surfaceToCellMap.insert(std::pair<int,int>(b,i));
				}

			}
			// Otherwise, we need to be a bit more generic in our approach
			else
			{
				/* ... */
			}
			
			if(includeCell)
			{
				double cellValue, surfValue, surfArea, localContribution;

				// Update surf variable data
				for(int v=0; v < this->variableNames.size(); v++)
				{
					std::string var = this->variableNames[v];
// std::cout<<"var = "<<var<<std::endl;
					vtkSmartPointer<vtkDataArray> cellDataArray = ds->GetCellData()->GetArray(var.c_str());
					vtkSmartPointer<vtkDataArray> surfDataArray = this->surfaces->GetPointData()->GetArray(var.c_str());


					cellValue	= cellDataArray->GetTuple1(i);
					surfValue	= surfDataArray->GetTuple1(b);
					localContribution = cellValue*cellArea;
				
					// Update surf
					surfDataArray->SetTuple1(b,surfValue + localContribution);
				}


				// Update surf volume data
				surfArea 	= surfaceAreaArray->GetTuple1(b);
				surfaceAreaArray->SetTuple1(b, surfArea + cellArea);

				// update cell volume fraction
				// Right now, this adds the raw volume
				// we'll normalize later to obtain the true volume fraction
				// this->setVolumeFraction(i, b, cellVolume);

#if 1
				// Add to map
				this->cellToSurfaceMap.insert(std::pair<int,int>(i,b));
				this->surfaceToCellMap.insert(std::pair<int,int>(b,i));
#endif
			}

		}
	
	}
	
	// Divide the surface values by the surface area to obtain the actual average
	double totalArea = 0.0;

	for(int v=0; v < this->variableNames.size(); v++)
	{
		std::string var = this->variableNames[v];
		vtkSmartPointer<vtkDataArray> surfaceDataArray = this->surfaces->GetPointData()->GetArray(var.c_str());

		for(int b=0; b<this->surfaces->GetNumberOfPoints();b++)
		{
			double surfaceValue, surfaceArea;

			surfaceValue	= surfaceDataArray->GetTuple1(b);
			surfaceArea	= surfaceAreaArray->GetTuple1(b);
			totalArea += surfaceArea;
			// std::cout<<"Surface area: "<<surfaceArea<<std::endl;
			if(surfaceArea > ERRMAX)
			{
				surfaceDataArray->SetTuple1(b,surfaceValue/surfaceArea);
			}
		}
	}
	// std::cout<<"Total area: "<<totalArea<<std::endl;

	// std::cout<<"[Leaving] SurfaceAveragePlaneOperator::process"<<std::endl;
	return this->surfaces;

}

int SurfaceAveragePlaneOperator::cellIntersectsBox(double planeNormal[3], double planePoint[3], double cellEdges[6])
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


void SurfaceAveragePlaneOperator::addVariable(const char* variableName, double initialValue, bool ignoreVariableList)
{	

	int numPts = this->surfaces->GetNumberOfPoints();
	vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
	tmpArray->Allocate(numPts);
	tmpArray->SetName(variableName);
	tmpArray->SetNumberOfComponents(1);
	tmpArray->SetNumberOfTuples(numPts);
	// tmpArray->PrintSelf(std::cout,vtkIndent(0));

	double* tmpPtr = (double*)tmpArray->GetVoidPointer(0);

	for(int i=0; i<numPts; i++)
	{
		// tmpArray->SetTuple1(i,initialValue);
		tmpPtr[i] = initialValue;
	}
	
	this->surfaces->GetPointData()->AddArray(tmpArray);
	if(!ignoreVariableList)
	{
		this->variableNames.push_back(variableName);
	}
}














