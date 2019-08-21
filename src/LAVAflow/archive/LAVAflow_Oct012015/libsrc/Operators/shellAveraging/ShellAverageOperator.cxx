#include "ShellAverageOperator.h"


#include <sstream>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkCellData.h>



vtkSmartPointer<vtkDataArray> ShellAverageOperator::removeAverageFromData(vtkSmartPointer<vtkDataArray> origData, 
																		const char* variableName)
{

	// Initialize the new dataarray
	vtkSmartPointer<vtkDoubleArray> newData  = vtkSmartPointer<vtkDoubleArray>::New();
	int numCells = origData->GetNumberOfTuples();
	for(int i=0; i<numCells; i++)
	{
		newData->InsertNextValue(0.0);
	}

	// Loop through all cells and remove the appropriate volume-weighted average values
	std::multimap<int,int>::iterator cellIt, binIt;
	std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> binRangeIt;
	int cellInd, binInd;

	for (cellIt = this->cellToBinMap.begin(); cellIt != this->cellToBinMap.end(); ++cellIt)
	{
		double origValue=0, avgValue=0, newValue = 0;

		// Get the original value
		cellInd = (*cellIt).first;
		origValue = origData->GetTuple1(cellInd);
		newValue = origValue;

		// Get the bin range that this cell spans
		binRangeIt = this->cellToBinMap.equal_range((*cellIt).first);

		// For every bin, add the volume-weighted contribution
		for (binIt=binRangeIt.first; binIt!=binRangeIt.second; ++binIt)
		{
			// Get the average value from the bin
			binInd = (*binIt).second;
			avgValue = this->bins->GetCellData()->GetArray(variableName)->GetTuple1(binInd)*this->getVolumeFraction(cellInd, binInd);

			// Update the new value
			newValue -= avgValue;
		}

		// Store the new value
		newData->SetTuple1(cellInd,newValue);
	}

	return newData;
}



void ShellAverageOperator::printShells(std::ostream& out)
{	
	char delim = '\t';
	
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

		// Get bin extents		
		double binEdges[6];
		this->bins->GetCellBounds(binId,binEdges);

		// Write bin extents
		dataLine<<binEdges[0]<<delim<<binEdges[1];

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
