#include "Particles.h"

#include "../Readers/VTKReader.h"


Particles::Particles(int coordSys, const char* fileName)
{
	this->coordinateSystem = coordSys;
	this->readVTKDataSet(fileName);
}
Particles::Particles(int coordSys, std::string fileName)
{
	this->coordinateSystem = coordSys;
    	this->readVTKDataSet(fileName);
}

Particles::Particles(int coordSys)
{
	this->coordinateSystem = coordSys;
}
Particles::~Particles(){;}



int Particles::readVTKDataSet(std::string fileName)
{
	return readVTKDataSet(fileName.c_str());
}
int Particles::readVTKDataSet(const char* fileName)
{
	// Create VTK reader to read unstructured grid
	VTKUnstructuredReader reader;
	vtkDataSet* dsNew = reader.readFile(fileName);

	// Error check
	if(dsNew==NULL)
	{
		std::cerr<<"[Particles | readVTKDataSet] Error reading file \""<<fileName<<"\""<<std::endl;
		return -1;
	}

	this->dataSet = dsNew;

	// Set the number of particles in the dataset
	this->numberOfParticles = this->dataSet->GetPointData()->GetNumberOfTuples();

	return 1;
	
}


void Particles::addVariable(const char* variableName, double initialValue)
{	
	int numPoints = this->dataSet->GetNumberOfPoints();

	vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
	tmpArray->SetName(variableName);
	tmpArray->SetNumberOfComponents(1);
	tmpArray->SetNumberOfTuples(numPoints);
	for(int i=0; i<numPoints; i++)
	{
		tmpArray->SetTuple1(i,initialValue);
	}
	
	this->dataSet->GetPointData()->AddArray(tmpArray);
}

void Particles::storeData(vtkSmartPointer<vtkDataArray> dataIn, const char* variableName)
{
	dataIn->SetName(variableName);
	this->dataSet->GetPointData()->AddArray(dataIn);	
}

void Particles::createCoordinateDataArrays()
{
	this->createCoordinateDataArrays("x","y","z");
}

void Particles::createCoordinateDataArrays(const char* xName, const char* yName, const char* zName)
{

	// std::clog<<"[Particles::createCoordinateDataArrays] Entering function"<<std::endl;

	int nPts = this->getDataSet()->GetNumberOfPoints();

	// std::cout<<"[Particles::createCoordinateDataArrays] nPts = "<<nPts<<std::endl;

	vtkSmartPointer<vtkDoubleArray> xArray = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> yArray = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> zArray = vtkSmartPointer<vtkDoubleArray>::New();

	// std::cout<<"[Particles::createCoordinateDataArrays] Allocating arrays...";
	xArray->SetNumberOfComponents(1);
	xArray->SetNumberOfTuples(nPts);
	yArray->SetNumberOfComponents(1);
	yArray->SetNumberOfTuples(nPts);
	zArray->SetNumberOfComponents(1);
	zArray->SetNumberOfTuples(nPts);
	// std::cout<<"done!"<<std::endl;

	// std::cout<<"[Particles::createCoordinateDataArrays] Setting coordinates...";
	// std::cout.flush();
	double ptLocation[3];
	for(int i=0; i<nPts; i++)
	{
		this->dataSet->GetPoint(i,ptLocation);

		xArray->SetTuple1(i,ptLocation[0]);
		yArray->SetTuple1(i,ptLocation[1]);
		zArray->SetTuple1(i,ptLocation[2]);
	}
	// std::cout<<"done!"<<std::endl;


	// std::cout<<"[Particles::createCoordinateDataArrays] Storing coordinate arrays...";
	this->storeData(xArray, xName);
	this->storeData(yArray, yName);
	this->storeData(zArray, zName);
	// std::cout<<"done!"<<std::endl;

	// std::clog<<"[Particles::createCoordinateDataArrays] Leaving function"<<std::endl;
}



double Particles::sum(const char* varName, const char* maskName)
{
	double sum = 0.0;

	bool useMask = false;

	if(maskName!=NULL)
	{
		useMask = true;
	}

	// Loop through the mesh and compute the sum

	if(useMask)
	{	
		vtkSmartPointer<vtkDataArray> array = this->getDataArray(varName);
		vtkSmartPointer<vtkDataArray> mask = this->getDataArray(maskName);
		for(int i=0; i<array->GetNumberOfTuples(); i++)
		{
			if(mask->GetTuple1(i)>ERRMAX)
			{
				sum += array->GetTuple1(i);
			}
		}
	}
	else
	{
		vtkSmartPointer<vtkDataArray> array = this->getDataArray(varName);
		for(int i=0; i<array->GetNumberOfTuples(); i++)
		{
			sum += array->GetTuple1(i);
		}
	}

	return sum;
}

double Particles::mean(const char* varName, const char* maskName)
{
	double sum = 0.0;
	double numberOfSamples = 0.0;

	sum = this->sum(varName, maskName);
	if(maskName!=NULL)
	{
		numberOfSamples = this->sum(maskName);
	}
	else
	{
		numberOfSamples = double(this->numberOfParticles);
	}


	return sum/numberOfSamples;
}

double Particles::variance(const char* varName, const char* maskName)
{

	double mean = this->mean(varName,maskName);
	int numberOfSamples = 0;
	bool useMask = false;
	double variance = 0.0;

	if(maskName!=NULL)
	{
		useMask = true;
	}

	// Loop through the mesh and compute the sum

	if(useMask)
	{	
		vtkSmartPointer<vtkDataArray> array = this->getDataArray(varName);
		vtkSmartPointer<vtkDataArray> mask = this->getDataArray(maskName);
		for(int i=0; i<array->GetNumberOfTuples(); i++)
		{
			if(mask->GetTuple1(i)>ERRMAX)
			{
				double value = array->GetTuple1(i);
				variance += (value-mean)*(value-mean);
				numberOfSamples++;
			}
		}
	}
	else
	{
		vtkSmartPointer<vtkDataArray> array = this->getDataArray(varName);
		for(int i=0; i<array->GetNumberOfTuples(); i++)
		{
			double value = array->GetTuple1(i);
			variance += (value-mean)*(value-mean);
			numberOfSamples++;
		}
	}

	return variance/double(numberOfSamples-1);

}