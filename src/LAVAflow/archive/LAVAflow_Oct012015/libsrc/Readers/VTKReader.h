#ifndef VTK_READER_H
#define VTK_READER_H

#include "BaseReader.h"
#include <vtkDataSet.h>
#include <vtkRectilinearGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>

#include <vector>
#include <string>


template<class readerType>
class VTKReader
{

private:
	vtkSmartPointer<readerType> reader;

public:
	VTKReader();
	virtual vtkDataSet* readFile(const char* fname);
	virtual vtkDataSet* readFile(const char* fname, std::vector<std::string>* scalarNames);
};

/* ===================================
	Implementation
   =================================== */

template<class readerType>
VTKReader<readerType>::VTKReader()
{
	this->reader=vtkSmartPointer<readerType>::New();
}

template<class readerType>
vtkDataSet* VTKReader<readerType>::readFile(const char* fname)
{
	if(reader)
	{
		reader->SetFileName(fname);
		reader->ReadAllScalarsOn();
		reader->ReadAllVectorsOn(); 
		reader->Update();
		return reinterpret_cast<vtkDataSet*>(reader->GetOutput());
	}
	return NULL;
}

template<class readerType>
vtkDataSet* VTKReader<readerType>::readFile(const char* fname, std::vector<std::string>* scalarNames)
{
	if(reader)
	{
		reader->SetFileName(fname);
		reader->ReadAllScalarsOff();
		reader->ReadAllVectorsOff();
		reader->Update();

		vtkSmartPointer<vtkDataSet> ds;
		// ds->CopyStructure(reader->GetOutput());
		// ds->CopyStructure(reinterpret_cast<vtkDataSet*>(reader->GetOutput()));


		// Set scalars we want to read
		for (std::vector<std::string>::iterator i = scalarNames->begin(); i != scalarNames->end(); ++i)
		{
			reader->SetScalarsName((*i).c_str());
			reader->Update();

			// vtkSmartPointer<vtkAbstractArray> aa = reader->ReadArray("vari3",9,1);
			// ds->GetCellData()->AddArray(aa);
		}

		// reader->ReadAllScalarsOn();
		// reader->ReadAllVectorsOn(); 

		// reader->Update();
		// return reinterpret_cast<vtkDataSet*>(reader->GetOutput());

		return ds;
	}
	return NULL;
}

// Typedefs for different readers
typedef VTKReader<vtkRectilinearGridReader> VTKRectilinearReader;
typedef VTKReader<vtkUnstructuredGridReader> VTKUnstructuredReader;



#endif
