#ifndef LAVA_PARTICLES_H
#define LAVA_PARTICLES_H


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
#include <vtkPointData.h>

//STD includes
#include <string>
#include <iostream>

class Particles
{
public:
	Particles(int coordSys, std::string fileName);
	Particles(int coordSys, const char* fileName);
	Particles(int coordSys);
	~Particles();

	/* data */
private:
	int numberOfParticles;
	// Coordinate system of the mesh
	int coordinateSystem;
	// vtkDataSet comprising the backend unstructured mesh
	vtkSmartPointer<vtkDataSet> dataSet;




public:

	// Accessors
	void setCoordinateSystem(int cs) { this->coordinateSystem = cs; }
	int  getCoordinateSystem()       { return this->coordinateSystem; }
	void getCoordinateSystem(int& cs){ cs = this->coordinateSystem; }

	int getNumberOfParticles(){ return this->numberOfParticles; }

	void setDataSet(vtkSmartPointer<vtkDataSet> ds) { this->dataSet = ds; }
	vtkSmartPointer<vtkDataSet>  getDataSet()       { return this->dataSet; }
	void getDataSet(vtkSmartPointer<vtkDataSet> ds) { ds = this->dataSet; }

	vtkSmartPointer<vtkDataArray> getDataArray(const char* variableName){ return this->dataSet->GetPointData()->GetArray(variableName);}


	// Utility functions


	int readVTKDataSet(std::string fileName);
	int readVTKDataSet(const char* fileName);


	void getParticlePosition(int ind, double ptLocation[3]){ this->dataSet->GetPoint(ind,ptLocation); }

	void addVariable(const char* variableName, double initialValue = 0.0);
	void storeData(vtkSmartPointer<vtkDataArray> dataIn, const char* variableName);
	void createCoordinateDataArrays();
	void createCoordinateDataArrays(const char* xName, const char* yName, const char* zName);

	double sum(const char* varName, const char* maskName = NULL);
	double mean(const char* varName, const char* maskName = NULL);
	double variance(const char* varName, const char* maskName = NULL);


	// Operators
	vtkSmartPointer<vtkDataArray> operator[](const char* variableName)
	{
		return this->getDataArray(variableName);
	}



};





#endif