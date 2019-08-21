#ifndef SURFACEAVG_OPERATOR_H
#define SURFACEAVG_OPERATOR_H

#include <map>
#include <vector>
#include <string>
#include <ostream>

#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>


class vtkDataSet;

class SurfaceAverageOperator
{
	// Rectilinear dataset storing the surfacened data
	
public:
	vtkSmartPointer<vtkRectilinearGrid> surfaces;
	int nSurfaces;
	std::multimap<int,int> cellToSurfaceMap;
	std::multimap<int,int> surfaceToCellMap;
	std::vector<std::string> variableNames;

	SurfaceAverageOperator(){};
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds) = 0;
	vtkSmartPointer<vtkRectilinearGrid> getSurfaces(){return this->surfaces;} 

	vtkSmartPointer<vtkDataArray> getDataArray(const char* variableName){ return this->surfaces->GetPointData()->GetArray(variableName);}
	virtual void addVariable(const char* variableName, double initialValue = 0.0, bool ignoreVariableList = false) = 0;

	virtual void printSurfaces(std::ostream& out);

	void printCellToSurface(std::ostream& out)
	{
		std::multimap<int,int>::iterator cellIt, surfIt;
		std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> surfRangeIt;
		out<<"\nCell\tSurfaces"<<std::endl;
		for (cellIt = this->cellToSurfaceMap.begin(); cellIt != this->cellToSurfaceMap.end(); ++cellIt)
		{
			surfRangeIt = this->cellToSurfaceMap.equal_range((*cellIt).first);

			out<<(*cellIt).first;
			for (surfIt=surfRangeIt.first; surfIt!=surfRangeIt.second; ++surfIt)
			{
				out<<"\t"<<(*surfIt).second;
			}
			out<<endl;
			
		}
	}

	void printSurfaceToCell(std::ostream& out)
	{
		std::multimap<int,int>::iterator surfIt, cellIt;
		std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> cellRangeIt;
		out<<"\nSurface\tCells"<<std::endl;
		for (surfIt = this->surfaceToCellMap.begin(); surfIt != this->surfaceToCellMap.end(); ++surfIt)
		{
			cellRangeIt = this->surfaceToCellMap.equal_range((*surfIt).first);

			out<<(*surfIt).first;
			for (cellIt=cellRangeIt.first; cellIt!=cellRangeIt.second; ++cellIt)
			{
				out<<"\t"<<(*cellIt).second;
			}
			out<<endl;
			
		}
	}


	void storeData(vtkSmartPointer<vtkDataArray> dataIn, const char* variableName)
	{
		this->variableNames.push_back(variableName);
		dataIn->SetName(variableName);
		this->surfaces->GetPointData()->AddArray(dataIn);	
	}

};

#endif
