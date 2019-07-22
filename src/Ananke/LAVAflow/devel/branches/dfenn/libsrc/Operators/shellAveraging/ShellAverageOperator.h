#ifndef SHELLAVG_OPERATOR_H
#define SHELLAVG_OPERATOR_H

#include <map>
#include <vector>
#include <string>
#include <ostream>

#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>

class vtkDataSet;


class ShellAverageOperator
{
	// Rectilinear dataset storing the binned data
	
public:
	vtkSmartPointer<vtkRectilinearGrid> bins;
	int nBins;
	std::multimap<int,int> cellToBinMap;
	std::multimap<int,int> binToCellMap;
	std::vector<std::map<int,double> > cellVolumeFraction;
	std::vector<std::string> variableNames;

	ShellAverageOperator(){};
	virtual vtkSmartPointer<vtkDataSet> process(const vtkSmartPointer<vtkDataSet> ds) = 0;
	vtkSmartPointer<vtkRectilinearGrid> getBins(){return this->bins;} 
	int getNumberOfBins(){return nBins;}

	vtkSmartPointer<vtkDataArray> removeAverageFromData(vtkSmartPointer<vtkDataArray> origData, const char* variableName);


	vtkSmartPointer<vtkDataArray> getDataArray(const char* variableName){ return this->bins->GetCellData()->GetArray(variableName);}
	virtual void addVariable(const char* variableName) = 0;

	void storeData(vtkSmartPointer<vtkDataArray> dataIn, const char* variableName)
	{
		this->variableNames.push_back(variableName);
		dataIn->SetName(variableName);
		this->bins->GetCellData()->AddArray(dataIn);	
	}

	// Operators
	vtkSmartPointer<vtkDataArray> operator[](const char* variableName)
	{
		return this->getDataArray(variableName);
	}

	void setVolumeFraction(int cellInd, int binInd, double volumeFraction)
	{
		// If the cell ind exceeds the size of the list of cells, we need to increase the size
		// In general, cellInd could be greater than 1 past the end, so we need to actually resize
		// instead of just pushing_back
		if(cellInd >= this->cellVolumeFraction.size())
		{			
			this->cellVolumeFraction.resize(cellInd+1);
		}

		// Now, add/update the pair of <binId, volumeFraction> to the map
		// First, check if it's in the map. If so, update the value
		if(this->cellVolumeFraction[cellInd].count(binInd)>0)
		{
			this->cellVolumeFraction[cellInd][binInd] = volumeFraction;
		}
		// Otherwise, we need to insert the pair
		else
		{
			this->cellVolumeFraction[cellInd].insert(std::pair<int,double>(binInd, volumeFraction));
		}
	}

	double getVolumeFraction(int cellInd, int binInd)
	{
		// Make sure cellInd is actually allocated in the vector
		if(cellInd >= this->cellVolumeFraction.size() && cellInd > 0)
		{
			std::cerr<<"[ShellAverageOperator::getVolumeFraction] WARNING: CellInd is not in the vector of volume fractions."<<std::endl;
			return 0.0;
		}

		// Make sure binInd is in the map
		if(this->cellVolumeFraction[cellInd].count(binInd) == 0)
		{
			std::cerr<<"[ShellAverageOperator::getVolumeFraction] WARNING: The requested bin is not associated with the requested cell in the mapping."<<std::endl;
			return 0.0;
		}

		return this->cellVolumeFraction[cellInd][binInd];

	}

	void printCellToBin(std::ostream& out)
	{
		int prevCell = -1;
		std::multimap<int,int>::iterator cellIt, binIt;
		std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> binRangeIt;
		out<<"\nCell\tBins"<<std::endl;
		for (cellIt = this->cellToBinMap.begin(); cellIt != this->cellToBinMap.end(); ++cellIt)
		{

			if((*cellIt).first == prevCell)
				continue;

			binRangeIt = this->cellToBinMap.equal_range((*cellIt).first);

			out<<(*cellIt).first;
			for (binIt=binRangeIt.first; binIt!=binRangeIt.second; ++binIt)
			{
				out<<"\t"<<(*binIt).second;
			}
			out<<endl;

			prevCell = (*cellIt).first;
			
		}
	}

	void printBinToCell(std::ostream& out)
	{
		std::multimap<int,int>::iterator binIt, cellIt;
		std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> cellRangeIt;
		out<<"\nBin\tCells"<<std::endl;
		for (binIt = this->binToCellMap.begin(); binIt != this->binToCellMap.end(); ++binIt)
		{
			cellRangeIt = this->binToCellMap.equal_range((*binIt).first);

			out<<(*binIt).first;
			for (cellIt=cellRangeIt.first; cellIt!=cellRangeIt.second; ++cellIt)
			{
				out<<"\t"<<(*cellIt).second;
			}
			out<<endl;
			
		}
	}

	virtual void printShells(std::ostream& out);


};

#endif
