#include "SurfaceAverageOperator.h"

#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <sstream>

void SurfaceAverageOperator::printSurfaces(std::ostream& out)
{	
	char delim = '\t';
	std::streamsize prec = 15;

	// Create header
	std::stringstream header;
	header<<"surfacePt"<<delim<<"area";
	for (std::vector<std::string>::iterator i = this->variableNames.begin(); i != this->variableNames.end(); ++i)
	{
		header<<delim<<*i;
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
		dataLine<<surfacePt[0];

		// Write surface area
		dataLine<<delim<<this->surfaces->GetPointData()->GetArray("area")->GetTuple1(surfaceId);

		// Loop over all variables for this surface
		for (std::vector<std::string>::iterator varIt = this->variableNames.begin(); varIt != this->variableNames.end(); ++varIt)
		{
			dataLine<<delim<<this->surfaces->GetPointData()->GetArray((*varIt).c_str())->GetTuple1(surfaceId);
			
		}


		// write dataLine to ostream
		out<<dataLine.str()<<std::endl;

	}
}