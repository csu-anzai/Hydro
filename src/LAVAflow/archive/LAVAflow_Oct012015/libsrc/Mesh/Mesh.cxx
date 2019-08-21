#include "Mesh.h"
#include <vtkRectilinearGrid.h>
#include <vtkImageData.h>
#include <vtkMath.h>
#include <vtkCell.h>
#include "../Utilities/LAVAUtil.h"


Mesh::Mesh(int dim, int coordSys, const char* fileName)
{
	// std::clog<<"In Mesh (const char* fileName) constructor"<<std::endl;
	this->dimension = dim;
	this->coordinateSystem = coordSys;
	this->coordArraysCreated = false;
	this->readVTKDataSet(fileName);
	this->cellBounds = NULL;
	this->populateCellBounds();
}
Mesh::Mesh(int dim, int coordSys, std::string fileName)
{
	// std::clog<<"In Mesh (string fileName) constructor"<<std::endl;
	this->dimension = dim;
	this->coordinateSystem = coordSys;
	this->coordArraysCreated = false;
	this->readVTKDataSet(fileName);
	this->cellBounds = NULL;
	this->populateCellBounds();
}

Mesh::Mesh(int dim, int coordSys)
{
	this->dimension = dim;
	this->coordinateSystem = coordSys;
	this->coordArraysCreated = false;
	this->cellBounds = NULL;
	this->populateCellBounds();
}


Mesh::Mesh(int dim, int coordSys, const char* fileName, std::vector<std::string> scalarNames)
{
	// std::clog<<"In Mesh (const char* fileName) constructor"<<std::endl;
	this->dimension = dim;
	this->coordinateSystem = coordSys;
	this->coordArraysCreated = false;
	this->readVTKDataSet(fileName, &scalarNames);
	this->cellBounds = NULL;
	this->populateCellBounds();
}



Mesh::~Mesh()
{
	if(this->cellBounds != NULL)
	{
		delete[] this->cellBounds;
	}
}

void Mesh::populateCellBounds()
{
// 	std::cout<<"Entering populateCellBounds..."<<std::endl;
// 	std::cout<<"\tthis->cellBounds at start: "<<this->cellBounds<<std::endl;
	if(this->cellBounds != NULL)
	{
		delete[] this->cellBounds;
	}
	// std::cout<<"\tPast delete"<<std::endl;

	// Get the number of cells in the mesh
	int nCells = this->dataSet->GetNumberOfCells();
	// std::cout<<"\tNumber of cells: "<<nCells<<std::endl;
	// Allocate and populate cell bounds
	double tmpCellEdges[6];
	this->cellBounds = new double*[nCells];
	for(int i=0; i<nCells; i++)
	{
		this->cellBounds[i] = new double[6];

		this->dataSet->GetCellBounds(i,tmpCellEdges);

		for(int j=0; j<6; j++)
		{
			this->cellBounds[i][j] = tmpCellEdges[j];
		}
	}

	// std::cout<<"Leaving populateCellBounds..."<<std::endl;
}


int Mesh::readVTKDataSet(std::string fileName)
{
	return readVTKDataSet(fileName.c_str());
}
int Mesh::readVTKDataSet(const char* fileName, std::vector<std::string>* scalarNames)
{
	// Create VTK reader to read rectilinear grid
	VTKRectilinearReader reader;
	vtkDataSet* dsNew;
	if(scalarNames==NULL)
	{	
		dsNew = reader.readFile(fileName);
	}
	else
	{
		dsNew = reader.readFile(fileName, scalarNames);
	}

	// Error check
	if(dsNew==NULL)
	{
		std::cerr<<"[Mesh | readVTKDataSet] Error reading file \""<<fileName<<"\""<<std::endl;
		return -1;
	}

	this->dataSet = dsNew;



	if(this->dataSet->GetDataObjectType() == VTK_RECTILINEAR_GRID)
	{
		vtkRectilinearGrid::SafeDownCast(this->dataSet)->GetDimensions(this->dataDim); // These are the dimensions of points, not cells. dimCells=dsDim-1
	}
	if(this->dataSet->GetDataObjectType() == VTK_UNIFORM_GRID)
	{
		vtkImageData::SafeDownCast(this->dataSet)->GetDimensions(this->dataDim); // These are the dimensions of points, not cells. dimCells=dsDim-1
	}

	this->dataDim[0]	-= 1;
	this->dataDim[1]	-= 1;
	this->dataDim[2]	-= 1;

	return 0;
	
}


void Mesh::addVariable(const char* variableName, double initialValue)
{	

	int numCells = this->dataSet->GetNumberOfCells();
	vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
	tmpArray->Allocate(numCells);
	tmpArray->SetName(variableName);
	tmpArray->SetNumberOfComponents(1);
	tmpArray->SetNumberOfTuples(numCells);
	// tmpArray->PrintSelf(std::cout,vtkIndent(0));

	double* tmpPtr = (double*)tmpArray->GetVoidPointer(0);

	for(int i=0; i<numCells; i++)
	{
		// tmpArray->SetTuple1(i,initialValue);
		tmpPtr[i] = initialValue;
	}
	
	this->dataSet->GetCellData()->AddArray(tmpArray);
}

void Mesh::storeData(vtkSmartPointer<vtkDataArray> dataIn, const char* variableName)
{
	dataIn->SetName(variableName);
	this->dataSet->GetCellData()->AddArray(dataIn);	
}

void Mesh::createCoordinateDataArrays()
{
	this->createCoordinateDataArrays("x","y","z");
}

void Mesh::createCoordinateDataArrays(const char* xName, const char* yName, const char* zName)
{

	// std::clog<<"[Mesh::createCoordinateDataArrays] Entering function"<<std::endl;

	vtkSmartPointer<vtkDoubleArray> xArray = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> yArray = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> zArray = vtkSmartPointer<vtkDoubleArray>::New();

	int nPts = this->getDataSet()->GetNumberOfCells();
	xArray->SetNumberOfComponents(1);
	xArray->SetNumberOfTuples(nPts);
	yArray->SetNumberOfComponents(1);
	yArray->SetNumberOfTuples(nPts);
	zArray->SetNumberOfComponents(1);
	zArray->SetNumberOfTuples(nPts);

	double* xArrayPtr = (double*)xArray->GetVoidPointer(0);
	double* yArrayPtr = (double*)yArray->GetVoidPointer(0);
	double* zArrayPtr = (double*)zArray->GetVoidPointer(0);

	double cellBounds[6];
	double cellCenter[3];

	for(int i=0; i<nPts; i++)
	{
		// this->getDataSet()->GetCellBounds(i,cellBounds);
		for(int j=0; j<6; j++)
		{
			cellBounds[j] = this->cellBounds[i][j];
		}

		// x
		if(this->coordinateSystem == CS_CART)
		{
			cellCenter[0] = 0.5*(cellBounds[0]+cellBounds[1]);
		}
		else if(this->coordinateSystem == CS_SPHERE || this->coordinateSystem == CS_CYL)
		{
			double R = 0.0;
			if(this->coordinateSystem == CS_CYL)
			{
				R = 0.5*(cellBounds[1]*cellBounds[1] + cellBounds[0]*cellBounds[0]);
				R = pow(R,0.5);
			}
			else if(this->coordinateSystem == CS_SPHERE)
			{
				R = 0.5*(cellBounds[1]*cellBounds[1]*cellBounds[1] + cellBounds[0]*cellBounds[0]*cellBounds[0]);
				R = pow(R,1.0/3.0);
			}
			cellCenter[0] = R;
		}
		// y-z
		cellCenter[1] = 0.5*(cellBounds[2]+cellBounds[3]);
		cellCenter[2] = 0.5*(cellBounds[4]+cellBounds[5]);

		// xArray->SetTuple1(i,cellCenter[0]);
		// yArray->SetTuple1(i,cellCenter[1]);
		// zArray->SetTuple1(i,cellCenter[2]);

		xArrayPtr[i] = cellCenter[0];
		yArrayPtr[i] = cellCenter[1];
		zArrayPtr[i] = cellCenter[2];
	}

	this->storeData(xArray, xName);
	this->storeData(yArray, yName);
	this->storeData(zArray, zName);

	// Save coordinate axes names
	strcpy(this->iCoordName,xName);
	strcpy(this->jCoordName,yName);
	strcpy(this->kCoordName,zName);

	// Flag that we've created the coordinate arrays
	this->coordArraysCreated = true;

	// Store a copy of the pointer to the coordinate arrays
	this->iCoordArray = xArray;
	this->jCoordArray = yArray;
	this->kCoordArray = zArray;

	// std::clog<<"[Mesh::createCoordinateDataArrays] Leaving function"<<std::endl;
}



void Mesh::computeCellVolumes()
{
	this->computeCellVolumes("volume");
}
void Mesh::computeCellVolumes(const char* volName)
{

	// Add array to the dataset
	this->addVariable(volName);

	// Get the pointer to the new array
	vtkSmartPointer<vtkDataArray> volume = this->getDataArray(volName);

	double* volumePtr = (double*)volume->GetVoidPointer(0);
	double cellEdges[6];
	
	// Loop over all cells and compute their volumes
	for(int i=0; i<this->getDataSet()->GetNumberOfCells(); i++)
	{

		// Get cell boundaries
		// this->getDataSet()->GetCellBounds(i, cellEdges);
		for(int j=0; j<6; j++)
		{
			cellEdges[j] = this->cellBounds[i][j];
		}

		// Compute volume (or area in 2d) of cell
		double cellVolume = 1.0;
		if(this->getCoordinateSystem() == CS_CART)
		{
			cellVolume = 1.0;
			for(int j=0; j<this->dimension; j++)
			{
				cellVolume *= cellEdges[2*j+1]-cellEdges[2*j];
			}
		}
		else if(this->getCoordinateSystem() == CS_CYL)
		{
			cellVolume = 1.0;
			cellVolume *= cellEdges[1]*cellEdges[1] - cellEdges[0]*cellEdges[0]; // Radial component
			cellVolume *= 0.5*(cellEdges[3]-cellEdges[2]); // Theta component
			if(this->dimension>2)
			{
				cellVolume *= cellEdges[5]-cellEdges[4]; // height component
			}
		}
		else if(this->getCoordinateSystem() == CS_SPHERE)
		{
			cellVolume = 1.0;
			cellVolume *= (1.0/3.0)*(	cellEdges[1]*cellEdges[1]*cellEdges[1] - 
									cellEdges[0]*cellEdges[0]*cellEdges[0]); // Radial component
			
			if(this->dimension>1)
			{
				cellVolume *= -cos(cellEdges[3]) + cos(cellEdges[2]); // Phi component
			}
			else
			{
				cellVolume *= 2.0;
			}
			if(this->dimension>2)
			{
				cellVolume *= cellEdges[5] - cellEdges[4]; // Theta component
			}
			else
			{
				cellVolume *= 2.0*vtkMath::Pi();
			}
		}
		else
		{
			std::cerr<<"[Mesh::computeCellVolumes] ERROR: Coordinate system not supported!"<<std::endl;
			return;
		}

		// Set the cell volume
		volume->SetTuple1(i,cellVolume);

	}
}


void Mesh::firstDerivative(int axis, const char* variableName, const char* derivName)
{

	// std::cout<<"[Mesh::firstDerivative] Entering function."<<std::endl;

	// Get some information about the mesh
	int nDim		= this->dimension;
	int coordSys	= this->getCoordinateSystem();
	int dataDim[3];   this->getDataDimension(dataDim);

	// std::cout<<"nDim = "<<nDim<<std::endl;
	// std::cout<<"coordSys = "<<coordSys<<std::endl;
	// std::cout<<"dataDim = "<<dataDim[0]<<"\t"<<dataDim[1]<<"\t"<<dataDim[2]<<"\t"<<std::endl;


	// Get the dataArray we want to take the derivative of
	vtkSmartPointer<vtkDataArray> data = this->getDataArray(variableName);

	// Add new variable
	this->addVariable(derivName);
	vtkSmartPointer<vtkDataArray> derivative = this->getDataArray(derivName);

	int k, kMin, kMax;
	int j, jMin, jMax;
	int i, iMin, iMax;

	iMin = 0;
	iMax = dataDim[0];
	jMin = 0;
	jMax = nDim>1 ? dataDim[1] : 1;
	kMin = 0;
	kMax = nDim>2 ? dataDim[2] : 1;

	int indexRange[3][2];
	indexRange[IAXIS][LOW]	= iMin;
	indexRange[IAXIS][HIGH]	= iMax;
	indexRange[JAXIS][LOW]	= jMin;
	indexRange[JAXIS][HIGH]	= jMax;
	indexRange[KAXIS][LOW]	= kMin;
	indexRange[KAXIS][HIGH]	= kMax;

	std::cout<<"Index range: "<<std::endl;
	std::cout<<indexRange[IAXIS][LOW]<<"\t"<<indexRange[JAXIS][LOW]<<"\t"<<indexRange[KAXIS][LOW]<<std::endl;
	std::cout<<indexRange[IAXIS][HIGH]<<"\t"<<indexRange[JAXIS][HIGH]<<"\t"<<indexRange[KAXIS][HIGH]<<std::endl;

	double cellBounds[6];
	double cellCenterP[3], cellCenterM[3], cellCenterC[3];
	double fP, fM, fC;

	int cellIndC[3], cellIndP[3], cellIndM[3];

	// std::cout<<"i"<<"\t"<<"j"<<"\t"<<"k"<<"\t"<<"x"<<"\t"<<"y"<<"\t"<<"z"<<std::endl;
	// std::cout<<"i"<<"\t"<<"j"<<"\t"<<"k"<<"\t"<<"fP"<<"\t"<<"fM"<<"\t"<<"deriv"<<std::endl;
	for(k = indexRange[KAXIS][LOW]; k < indexRange[KAXIS][HIGH]; ++k)
	{
		for(j = indexRange[JAXIS][LOW]; j < indexRange[JAXIS][HIGH]; ++j)
		{
			for(i = indexRange[IAXIS][LOW]; i < indexRange[IAXIS][HIGH]; ++i)
			{

				cellIndP[IAXIS] = i;	cellIndP[JAXIS] = j;	cellIndP[KAXIS] = k;
				cellIndM[IAXIS] = i;	cellIndM[JAXIS] = j;	cellIndM[KAXIS] = k;
				cellIndC[IAXIS] = i;	cellIndC[JAXIS] = j;	cellIndC[KAXIS] = k;
				
				cellIndP[axis] += 1;
				cellIndM[axis] -= 1;

				// Shift for backward difference
				if(cellIndP[axis]>=(indexRange[axis][HIGH]-1))
				{
					cellIndP[axis] -= 1;
					cellIndM[axis] -= 1;
					cellIndC[axis] -= 1;
				}
				// Shift for forward difference
				else if(cellIndM[axis]<=indexRange[axis][LOW])
				{
					cellIndP[axis] += 1;
					cellIndM[axis] += 1;
					cellIndC[axis] += 1;
				}
				// Leave things alone for centered difference
				else if(cellIndM[axis]>=indexRange[axis][LOW] && cellIndP[axis]<indexRange[axis][HIGH])
				{
					// std::cout<<"Central difference available."<<std::endl;
				}
				else
				{
					std::cout<<"This shouldn't be happening."<<std::endl;
				}

				int ind	 = Util::sub2ind(i, j, k, dataDim[0], dataDim[1]);
				int indP = Util::sub2ind(cellIndP[0], cellIndP[1], cellIndP[2], dataDim[0], dataDim[1]);
				int indM = Util::sub2ind(cellIndM[0], cellIndM[1], cellIndM[2], dataDim[0], dataDim[1]);
				int indC = Util::sub2ind(cellIndC[0], cellIndC[1], cellIndC[2], dataDim[0], dataDim[1]);

				// std::cout<<"ind\tindM\tindC\tindP"<<std::endl;
				// std::cout<<ind<<"\t"<<indM<<"\t"<<indC<<"\t"<<indP<<std::endl;

				if(ind<0 || indP<0 || indM<0 || indC<0)
				{

					std::cout<<"ind\tindM\tindC\tindP"<<std::endl;
					std::cout<<ind<<"\t"<<indM<<"\t"<<indC<<"\t"<<indP<<std::endl;
					std::cout<<"i\tj\tk"<<std::endl;
					std::cout<<i<<"\t"<<j<<"\t"<<k<<std::endl;
				}


				this->getCellCenter(indP, cellCenterP);
				fP = data->GetTuple1(indP);


				this->getCellCenter(indM, cellCenterM);
				fM = data->GetTuple1(indM);


				this->getCellCenter(indC, cellCenterC);
				fC = data->GetTuple1(indC);




				double dist = cellCenterP[axis] - cellCenterM[axis];

				if(fabs(dist)<ERRMAX)
				{
					std::cout<<"[Mesh::firstDerivative] Distance = "<<dist<<std::endl;
					std::cout<<"CellCenterP["<<axis<<"] = "<<cellCenterP[axis]<<std::endl;
					std::cout<<"CellCenterM["<<axis<<"] = "<<cellCenterM[axis]<<std::endl;
					std::cout<<"CellIndP = "<<cellIndP[IAXIS]<<"\t"<<cellIndP[JAXIS]<<"\t"<<cellIndP[KAXIS]<<std::endl;
					std::cout<<"CellIndM = "<<cellIndM[IAXIS]<<"\t"<<cellIndM[JAXIS]<<"\t"<<cellIndM[KAXIS]<<std::endl;
				}

				double deriv = (fP-fM)/dist;
				// double deriv = (fP - 2*fC + fM)/(distSqr);

				// Store derivative
				derivative->SetTuple1(ind,deriv);
			}
		}
	}

	// std::cout<<"[Mesh::firstDerivative] Leaving function."<<std::endl;
}

void Mesh::getCellCenter(int ind, double center[3])
{


	// If we've already computed the coordinate arrays, just retrieve the data from the arrays
	if(this->coordArraysCreated)
	{
		center[0] = this->iCoordArray->GetTuple1(ind);
		center[1] = this->jCoordArray->GetTuple1(ind);
		center[2] = this->kCoordArray->GetTuple1(ind);
	}
	// Otherwise, we need to compute the center
	else
	{
		// Get the cell bounds
		double cellBounds[6];
		vtkSmartPointer<vtkDataSet> ds = this->getDataSet();
		// this->getDataSet()->GetCellBounds(ind,cellBounds);
		ds->GetCellBounds(ind,cellBounds);
		
		switch(this->coordinateSystem)
		{
			case(CS_CART):
				center[0] = 0.5*(cellBounds[0]+cellBounds[1]);
				center[1] = 0.5*(cellBounds[2]+cellBounds[3]);
				center[2] = 0.5*(cellBounds[4]+cellBounds[5]);
				break;

			case(CS_CYL):
				center[0] = sqrt(0.5*(cellBounds[0]*cellBounds[0] + cellBounds[1]*cellBounds[1]));
				center[1] = 0.5*(cellBounds[2]+cellBounds[3]);
				center[2] = 0.5*(cellBounds[4]+cellBounds[5]);
				break;

			case(CS_SPHERE):
				center[0] = pow(0.5*(cellBounds[0]*cellBounds[0]*cellBounds[0] + cellBounds[1]*cellBounds[1]*cellBounds[1]),1.0/3.0);
				center[1] = 0.5*(cellBounds[2]+cellBounds[3]);
				center[2] = 0.5*(cellBounds[4]+cellBounds[5]);
				break;

			default:
				std::cerr<<"[Mesh::getCellCenter] ERROR: Coordinate system "<<this->coordinateSystem<<" not implemented!"<<std::endl;
				break;
		}
	}

}

double Mesh::sum(const char* varName, const char* maskName)
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

std::vector<int> Mesh::getCellsAlongRay(int axis, int startInd[3])
{
	std::vector<int> cellIndices;

	int k, kMin, kMax;
	int j, jMin, jMax;
	int i, iMin, iMax;

	// Handle dimensionality
	iMin = 0;
	iMax = this->dataDim[0];
	jMin = 0;
	jMax = this->dimension>1 ? this->dataDim[1] : 1;
	kMin = 0;
	kMax = this->dimension>2 ? this->dataDim[2] : 1;

	// Set the range to the entire domain
	int indexRange[3][2];
	indexRange[IAXIS][LOW]	= iMin;
	indexRange[IAXIS][HIGH]	= iMax;
	indexRange[JAXIS][LOW]	= jMin;
	indexRange[JAXIS][HIGH]	= jMax;
	indexRange[KAXIS][LOW]	= kMin;
	indexRange[KAXIS][HIGH]	= kMax;

	// Trim this to only look at cells along the axis specified, starting at the startInd
	// I'd rather be verbose here than do some fancy loop trick to set this data
	if(axis == XAXIS)
	{
		indexRange[JAXIS][LOW] = startInd[JAXIS];
		indexRange[JAXIS][HIGH] = startInd[JAXIS]+1;
		indexRange[KAXIS][LOW] = startInd[KAXIS];
		indexRange[KAXIS][HIGH] = startInd[KAXIS]+1;
	}
	else if(axis == YAXIS)
	{
		indexRange[IAXIS][LOW] = startInd[IAXIS];
		indexRange[IAXIS][HIGH] = startInd[IAXIS]+1;
		indexRange[KAXIS][LOW] = startInd[KAXIS];
		indexRange[KAXIS][HIGH] = startInd[KAXIS]+1;

	}
	else if(axis == ZAXIS)
	{
		indexRange[IAXIS][LOW] = startInd[IAXIS];
		indexRange[IAXIS][HIGH] = startInd[IAXIS]+1;
		indexRange[JAXIS][LOW] = startInd[JAXIS];
		indexRange[JAXIS][HIGH] = startInd[JAXIS]+1;

	}
	else
	{
		std::cerr<<"I shouldn't be here!"<<std::endl;
		return cellIndices;
	}

	for(k = indexRange[KAXIS][LOW]; k < indexRange[KAXIS][HIGH]; ++k)
	{
		for(j = indexRange[JAXIS][LOW]; j < indexRange[JAXIS][HIGH]; ++j)
		{
			for(i = indexRange[IAXIS][LOW]; i < indexRange[IAXIS][HIGH]; ++i)
			{

				int ind	= k*dataDim[0]*dataDim[1] + j*dataDim[0] + i;
				cellIndices.push_back(ind);
			}
		}
	}

	return cellIndices;
}


double Mesh::interpolate(double pt[3], const char* varName, int interpOrder)
{
	double value = 0.0;

	double tol2 = 0.001;
	int subId;
	double pcoords[3], weights[8];
	int cellIndex = -1;

	// Get the index of the cell this point lies in
	cellIndex = this->dataSet->FindCell(pt,NULL,-1,tol2,subId,pcoords,weights);


	if(cellIndex<0)
	{
		std::cerr<<"[Mesh::interpolate] ERROR: Sample point not inside mesh!"<<std::endl;
		std::cerr<<"[Mesh::interpolate]\tSample point: ( "<<pt[0]<<" , "<<pt[0]<<" , "<<pt[0]<<" )"<<std::endl;
		return -1e99;
	}

	// Determine interpolation scheme
	switch(interpOrder)
	{
		case(0):
			value = this->getDataArray(varName)->GetTuple1(cellIndex);
			break;
		
		default:
			std::cerr<<"[Mesh::interpolate] ERROR: No implemented interpolation scheme for order = "<<interpOrder<<std::endl;

	}

	return value;
}










vtkSmartPointer<vtkDataArray> operator /(const vtkSmartPointer<vtkDataArray>& A, const vtkSmartPointer<vtkDataArray>& B)
{
	vtkSmartPointer<vtkDataArray> C = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuplesA, nTuplesB, nTuples;

	nTuplesA = A->GetNumberOfTuples();
	nTuplesB = B->GetNumberOfTuples();

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);

	// std::cout<<"[operator /] Number of tuples for A and B: "<<nTuplesA<<"\t"<<nTuplesB<<std::endl;

	if(nTuplesA == nTuplesB)
	{
		nTuples = nTuplesA;

		// Set the number of tuples for the result
		C->Allocate(nTuples);
		C->SetNumberOfTuples(nTuples);
		C->SetNumberOfComponents(1);
		double* cPtr = (double*)C->GetVoidPointer(0);

		// Do the division
		for(int i=0; i<nTuples; i++)
		{
			double a = A->GetTuple1(i);
			double b = B->GetTuple1(i);

			// double a = aPtr[i];
			// double b = bPtr[i];
			// std::cout<<"a = "<<a<<"\t";
			// std::cout<<"b = "<<b<<std::endl;
			if(fabs(b)!=0.0)
			{
				// C->SetTuple1(i,a/b);
				cPtr[i] = a/b;
			}
			else if(fabs(a) == 0.0)
			{
				// C->SetTuple1(i,0.0);
				cPtr[i] = 0.0;
			}
			else
			{
				// C->SetTuple1(i,1e99);
				cPtr[i] = 1e99;
				std::cerr<<"[operator/] Divison by zero!"<<std::endl;
				std::cerr<<"\tc=a/b: a = "<<a<<"\t"<<"b = "<<b<<std::endl;
			}
		}

	}
	else
	{
		/* ... */
		std::cerr<<"ERROR [Mesh operator/] Differing number of tuples between datasets!"<<std::endl;
	}

	return C;
}



vtkSmartPointer<vtkDataArray> operator*(const vtkSmartPointer<vtkDataArray>& A, const vtkSmartPointer<vtkDataArray>& B)
{
	vtkSmartPointer<vtkDataArray> C = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuplesA, nTuplesB, nTuples;

	nTuplesA = A->GetNumberOfTuples();
	nTuplesB = B->GetNumberOfTuples();

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);
	// std::cout<<"[operator /] Number of tuples for A and B: "<<nTuplesA<<"\t"<<nTuplesB<<std::endl;


	// vtkSmartPointer<vtkDoubleArray> tmpArray  = vtkSmartPointer<vtkDoubleArray>::New();
	// tmpArray->Allocate(numCells);
	// tmpArray->SetName(variableName);
	// tmpArray->SetNumberOfComponents(1);
	// tmpArray->SetNumberOfTuples(numCells);
	// tmpArray->PrintSelf(std::cout,vtkIndent(0));


	if(nTuplesA == nTuplesB)
	{
		nTuples = nTuplesA;

		// Set the number of tuples for the result
		C->Allocate(nTuples);
		C->SetNumberOfComponents(1);
		C->SetNumberOfTuples(nTuples);
		
		// C->PrintSelf(std::cout,vtkIndent(1));

		double* cPtr = (double*)C->GetVoidPointer(0);
		
		// Do the multiplication
		for(int i=0; i<nTuples; i++)
		{
			C->SetTuple1(i,A->GetTuple1(i)*B->GetTuple1(i));
			// cPtr[i] = aPtr[i]*bPtr[i];
			// cPtr[i] = 1.0;
		}

	}
	else
	{
		/* ... */
		std::cerr<<"ERROR [Mesh operator*] Differing number of tuples between datasets!"<<std::endl;
	}

	return C;
}


vtkSmartPointer<vtkDataArray> operator+(const vtkSmartPointer<vtkDataArray>& A, const vtkSmartPointer<vtkDataArray>& B)
{
	vtkSmartPointer<vtkDataArray> C = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuplesA, nTuplesB, nTuples;

	nTuplesA = A->GetNumberOfTuples();
	nTuplesB = B->GetNumberOfTuples();

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);
	// std::cout<<"[operator /] Number of tuples for A and B: "<<nTuplesA<<"\t"<<nTuplesB<<std::endl;

	if(nTuplesA == nTuplesB)
	{
		nTuples = nTuplesA;

		// Set the number of tuples for the result
		C->SetNumberOfTuples(nTuples);
		C->SetNumberOfComponents(1);
		double* cPtr = (double*)C->GetVoidPointer(0);

		// Do the addition
		for(int i=0; i<nTuples; i++)
		{
			C->SetTuple1(i,A->GetTuple1(i)+B->GetTuple1(i));
			// cPtr[i] = aPtr[i]+bPtr[i];
		}

	}
	else
	{
		/* ... */
		std::cerr<<"ERROR [Mesh operator+] Differing number of tuples between datasets!"<<std::endl;
	}

	return C;
}


vtkSmartPointer<vtkDataArray> operator-(const vtkSmartPointer<vtkDataArray>& A, const vtkSmartPointer<vtkDataArray>& B)
{
	vtkSmartPointer<vtkDataArray> C = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuplesA, nTuplesB, nTuples;

	nTuplesA = A->GetNumberOfTuples();
	nTuplesB = B->GetNumberOfTuples();

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);
	// std::cout<<"[operator /] Number of tuples for A and B: "<<nTuplesA<<"\t"<<nTuplesB<<std::endl;

	if(nTuplesA == nTuplesB)
	{
		nTuples = nTuplesA;

		// Set the number of tuples for the result
		C->SetNumberOfTuples(nTuples);
		C->SetNumberOfComponents(1);
		double* cPtr = (double*)C->GetVoidPointer(0);
	
		// Do the subtraction
		for(int i=0; i<nTuples; i++)
		{
			C->SetTuple1(i,A->GetTuple1(i)-B->GetTuple1(i));
			// cPtr[i] = aPtr[i] - bPtr[i];
		}

	}
	else
	{
		/* ... */
		std::cerr<<"ERROR [Mesh operator-] Differing number of tuples between datasets!"<<std::endl;
	}

	return C;
}


vtkSmartPointer<vtkDataArray> operator*(double scalar, const vtkSmartPointer<vtkDataArray>& A)
{
	vtkSmartPointer<vtkDataArray> B = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	B->SetNumberOfTuples(nTuples);
	B->SetNumberOfComponents(1);

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);

	// Do the scalar multiplication
	for(int i=0; i<nTuples; i++)
	{
		B->SetTuple1(i,scalar*A->GetTuple1(i));
		// bPtr[i] = scalar*aPtr[i];
	}

	return B;
}

vtkSmartPointer<vtkDataArray> operator*(const vtkSmartPointer<vtkDataArray>& A, double scalar)
{
	vtkSmartPointer<vtkDataArray> B = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	B->SetNumberOfTuples(nTuples);
	B->SetNumberOfComponents(1);

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);

	// Do the scalar multiplication
	for(int i=0; i<nTuples; i++)
	{
		B->SetTuple1(i,scalar*A->GetTuple1(i));
		// bPtr[i] = scalar*aPtr[i];
	}

	return B;
}


vtkSmartPointer<vtkDataArray> operator/(const vtkSmartPointer<vtkDataArray>& A, double scalar)
{
	vtkSmartPointer<vtkDataArray> B = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	B->SetNumberOfTuples(nTuples);
	B->SetNumberOfComponents(1);

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);
	// Do the scalar multiplication
	for(int i=0; i<nTuples; i++)
	{
		B->SetTuple1(i,A->GetTuple1(i)/scalar);
		// bPtr[i] = aPtr[i]/scalar;
	}

	return B;
}



vtkSmartPointer<vtkDataArray> operator/(double scalar, const vtkSmartPointer<vtkDataArray>& A)
{
	vtkSmartPointer<vtkDataArray> B = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	B->SetNumberOfTuples(nTuples);
	B->SetNumberOfComponents(1);
	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);

	// Do the scalar multiplication
	for(int i=0; i<nTuples; i++)
	{
		B->SetTuple1(i,scalar/A->GetTuple1(i));
		// bPtr[i] = scalar/aPtr[i];
	}

	return B;
}


vtkSmartPointer<vtkDataArray> operator+(double scalar, const vtkSmartPointer<vtkDataArray>& A)
{
	vtkSmartPointer<vtkDataArray> B = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	B->SetNumberOfTuples(nTuples);
	B->SetNumberOfComponents(1);

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);

	// Do the scalar multiplication
	for(int i=0; i<nTuples; i++)
	{
		B->SetTuple1(i,scalar+A->GetTuple1(i));
		// bPtr[i] = aPtr[i] + scalar;
	}

	return B;
}


vtkSmartPointer<vtkDataArray> operator+(const vtkSmartPointer<vtkDataArray>& A, double scalar)
{
	vtkSmartPointer<vtkDataArray> B = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	B->SetNumberOfTuples(nTuples);
	B->SetNumberOfComponents(1);

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);

	// Do the scalar multiplication
	for(int i=0; i<nTuples; i++)
	{
		B->SetTuple1(i,scalar+A->GetTuple1(i));
		// bPtr[i] = scalar + aPtr[i];
	}

	return B;
}


vtkSmartPointer<vtkDataArray> operator-(const vtkSmartPointer<vtkDataArray>& A, double scalar)
{
	vtkSmartPointer<vtkDataArray> B = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	B->SetNumberOfTuples(nTuples);
	B->SetNumberOfComponents(1);

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);

	// Do the scalar multiplication
	for(int i=0; i<nTuples; i++)
	{
		B->SetTuple1(i,A->GetTuple1(i)-scalar);
		// bPtr[i] = aPtr[i]-scalar;
	}

	return B;
}


vtkSmartPointer<vtkDataArray> operator-(double scalar, const vtkSmartPointer<vtkDataArray>& A)
{
	vtkSmartPointer<vtkDataArray> B = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	B->SetNumberOfTuples(nTuples);
	B->SetNumberOfComponents(1);

	double* aPtr = (double*)A->GetVoidPointer(0);
	double* bPtr = (double*)B->GetVoidPointer(0);
	
	// Do the scalar multiplication
	for(int i=0; i<nTuples; i++)
	{
		B->SetTuple1(i,scalar-A->GetTuple1(i));
		// bPtr[i] = scalar-aPtr[i];
	}

	return B;
}



vtkSmartPointer<vtkDataArray> operator >(const vtkSmartPointer<vtkDataArray>& A, double value)
{
	vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	mask->SetNumberOfTuples(nTuples);
	mask->SetNumberOfComponents(1);

	// Do the scalar multiplication
	double a;
	for(int i=0; i<nTuples; i++)
	{
		a = A->GetTuple1(i);
		if(a>value)
		{
			mask->SetTuple1(i,1.0);
		}
		else
		{
			mask->SetTuple1(i,0.0);
		}
	}

	return mask;
}


vtkSmartPointer<vtkDataArray> operator >(double value, const vtkSmartPointer<vtkDataArray>& A)
{
	vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	mask->SetNumberOfTuples(nTuples);
	mask->SetNumberOfComponents(1);

	// Do the scalar multiplication
	double a;
	for(int i=0; i<nTuples; i++)
	{
		a = A->GetTuple1(i);
		if(value>a)
		{
			mask->SetTuple1(i,1.0);
		}
		else
		{
			mask->SetTuple1(i,0.0);
		}
	}

	return mask;
}


vtkSmartPointer<vtkDataArray> operator <(const vtkSmartPointer<vtkDataArray>& A, double value)
{
	vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	mask->SetNumberOfTuples(nTuples);
	mask->SetNumberOfComponents(1);

	// Do the scalar multiplication
	double a;
	for(int i=0; i<nTuples; i++)
	{
		a = A->GetTuple1(i);
		if(a<value)
		{
			mask->SetTuple1(i,1.0);
		}
		else
		{
			mask->SetTuple1(i,0.0);
		}
	}

	return mask;
}

vtkSmartPointer<vtkDataArray> operator <(double value, const vtkSmartPointer<vtkDataArray>& A)
{
	vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkDoubleArray>::New();
	int nTuples = A->GetNumberOfTuples();

	// Set the number of tuples for the result
	mask->SetNumberOfTuples(nTuples);
	mask->SetNumberOfComponents(1);

	// Do the scalar multiplication
	double a;
	for(int i=0; i<nTuples; i++)
	{
		a = A->GetTuple1(i);
		if(value<a)
		{
			mask->SetTuple1(i,1.0);
		}
		else
		{
			mask->SetTuple1(i,0.0);
		}
	}

	return mask;
}
