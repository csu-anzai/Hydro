#include <vtkAMRFlashReader.h>
#include <vtkOverlappingAMR.h>
#include <vtkDataArraySelection.h>
#include <vtkDataSet.h>
#include "FlashReader.h"

FlashReader::FlashReader()
{
	reader=vtkSmartPointer<vtkAMRFlashReader>::New();
}

vtkDataSet* FlashReader::readFile(const char* fname)
{
	if(reader)
	{
		vtkOverlappingAMR *amr = NULL;
		reader->SetFileName( fname );
		vtkDataArraySelection* cdsa = reader->GetCellDataArraySelection();
		//cdsa->EnableArray(varName);
		//for now just enable all arrays for reading and read all levels
		cdsa->EnableAllArrays();
		reader->SetMaxLevel(reader->GetNumberOfLevels()-1);
		reader->Update();
		amr = reader->GetOutput();
		//cast to a generic data set.
		return reinterpret_cast<vtkDataSet*>(amr);
	}
	return NULL;
}
