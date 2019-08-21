#ifndef FLASH_READER_H
#define FLASH_READER_H

#include "BaseReader.h"
#include <vtkDataSet.h>
#include <vtkAMRFlashReader.h>
#include <vtkSmartPointer.h>

class FlashReader: public BaseReader
{
vtkSmartPointer<vtkAMRFlashReader> reader;
public:
FlashReader();
virtual vtkDataSet* readFile(const char* fname);
};

#endif
