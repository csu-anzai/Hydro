#ifndef LAVA_MESH_H
#define LAVA_MESH_H

/*
	Wrapper for VTK backend
*/

// LAVAflow includes
#include "../includes/LAVAconstants.h"
#include "../Readers/VTKReader.h"
// VTK includes
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>

//STD includes
#include <string>
#include <iostream>
#include <vector>	

class Mesh
{
public:
	Mesh(int dim, int coordSys);
	Mesh(int dim, int coordSys, std::string fileName);
	Mesh(int dim, int coordSys, const char* fileName);
	Mesh(int dim, int coordSys, const char* fileName, std::vector<std::string> scalarNames);
	~Mesh();

	/* data */
private:
	// Dimensionality of the mesh
	int dimension;
	int dataDim[3];
	// Coordinate system of the mesh
	int coordinateSystem;
	bool coordArraysCreated;
	char iCoordName[256], jCoordName[256], kCoordName[256];
	vtkSmartPointer<vtkDoubleArray> iCoordArray;
	vtkSmartPointer<vtkDoubleArray> jCoordArray;
	vtkSmartPointer<vtkDoubleArray> kCoordArray;
	// vtkDataSet comprising the backend mesh
	vtkSmartPointer<vtkDataSet> dataSet;
	double** cellBounds;




public:

	// Accessors
	void setDimension(int dim)	{ this->dimension = dim; }
	int  getDimension()        { return this->dimension; }
	void getDimension(int& dim){ dim = this->dimension; }

	void getDataDimension(int dim[])
	{
		for(int i=0; i<3; i++)
		{
			dim[i] = this->dataDim[i];
		}
	}

	void setCoordinateSystem(int cs) { this->coordinateSystem = cs; }
	int  getCoordinateSystem()       { return this->coordinateSystem; }
	void getCoordinateSystem(int& cs){ cs = this->coordinateSystem; }

	void setDataSet(vtkSmartPointer<vtkDataSet> ds) { this->dataSet = ds; }
	vtkSmartPointer<vtkDataSet>  getDataSet()       { return this->dataSet; }
	void getDataSet(vtkSmartPointer<vtkDataSet> ds) { ds = this->dataSet; }

	vtkSmartPointer<vtkDataArray> getDataArray(const char* variableName){ return this->dataSet->GetCellData()->GetArray(variableName);}


	// Utility functions
	int readVTKDataSet(std::string fileName);
	// int readVTKDataSet(const char* fileName);
	int readVTKDataSet(const char* fileName, std::vector<std::string>* scalarNames = NULL);

	void addVariable(const char* variableName, double initialValue = 0.0);
	void storeData(vtkSmartPointer<vtkDataArray> dataIn, const char* variableName);
	void createCoordinateDataArrays();
	void createCoordinateDataArrays(const char* xName, const char* yName, const char* zName);

	void computeCellVolumes();
	void computeCellVolumes(const char* volName);


	void getCellCenter(int ind, double center[3]);

	double sum(const char* varName, const char* maskName = NULL);


	virtual std::vector<int> getCellsAlongRay(int axis, double startCoord[3]){;}
	virtual std::vector<int> getCellsAlongRay(int axis, int startInd[3]);

	// Operators
	vtkSmartPointer<vtkDataArray> operator[](const char* variableName)
	{
		return this->getDataArray(variableName);
	}

	// Derivatives
	void firstDerivative(int axis, const char* variableName, const char* derivName);

	// Interpolation
	double interpolate(double pt[3], const char* varName, int interpOrder = 0);

	// Populator to store the cell bounds without having to call GetCellBounds
	void populateCellBounds();



};








vtkSmartPointer<vtkDataArray> operator /(const vtkSmartPointer<vtkDataArray>& a, 
										  const vtkSmartPointer<vtkDataArray>& b);

vtkSmartPointer<vtkDataArray> operator *(const vtkSmartPointer<vtkDataArray>& a, 
										  const vtkSmartPointer<vtkDataArray>& b);

vtkSmartPointer<vtkDataArray> operator +(const vtkSmartPointer<vtkDataArray>& a, 
										  const vtkSmartPointer<vtkDataArray>& b);

vtkSmartPointer<vtkDataArray> operator -(const vtkSmartPointer<vtkDataArray>& a, 
										  const vtkSmartPointer<vtkDataArray>& b);

vtkSmartPointer<vtkDataArray> operator*(double scalar, const vtkSmartPointer<vtkDataArray>& A);
vtkSmartPointer<vtkDataArray> operator*(const vtkSmartPointer<vtkDataArray>& A, double scalar);
vtkSmartPointer<vtkDataArray> operator/(double scalar, const vtkSmartPointer<vtkDataArray>& A);
vtkSmartPointer<vtkDataArray> operator/(const vtkSmartPointer<vtkDataArray>& A, double scalar);
vtkSmartPointer<vtkDataArray> operator+(double scalar, const vtkSmartPointer<vtkDataArray>& A);
vtkSmartPointer<vtkDataArray> operator+(const vtkSmartPointer<vtkDataArray>& A, double scalar);
vtkSmartPointer<vtkDataArray> operator-(double scalar, const vtkSmartPointer<vtkDataArray>& A);
vtkSmartPointer<vtkDataArray> operator-(const vtkSmartPointer<vtkDataArray>& A, double scalar);


vtkSmartPointer<vtkDataArray> operator >(const vtkSmartPointer<vtkDataArray>& a, double value);
vtkSmartPointer<vtkDataArray> operator >(double value, const vtkSmartPointer<vtkDataArray>& A);
vtkSmartPointer<vtkDataArray> operator <(const vtkSmartPointer<vtkDataArray>& a, double value);
vtkSmartPointer<vtkDataArray> operator <(double value, const vtkSmartPointer<vtkDataArray>& A);

#endif