#include <sstream>
#include "Histogram.h"


namespace Statistics
{
	void computeMinMax(double* data, int numDataPts, double& minVal, double& maxVal)
	{
		minVal = data[0];
		maxVal = data[0];

		for(int i=1; i<numDataPts; i++)
		{
			if(data[i]<minVal)
				minVal = data[i];
			if(data[i]>maxVal)
				maxVal = data[i];
		}
	}
};

namespace Statistics
{

Histogram::Histogram()
{
	//empty constructor -- nothing really needs to be done, since this is essentially a container for processed histogram objects
	// std::cout<<"[Histogram::Histogram()] Empty constructor"<<std::endl;
}
Histogram::~Histogram()
{

	// std::cout<<"[Histogram::~Histogram()] Entering destructor"<<std::endl;
	// Delete all HistObjs
	while (!this->histObjs.empty())
	{
		HistObj* tmp = this->histObjs.back();
		this->histObjs.pop_back();

		// Free histObj memory
		delete[] tmp->binEdges;
		delete[] tmp->binCount;
		delete tmp;
	}

	// std::cout<<"[Histogram::~Histogram()] Leaving destructor"<<std::endl;
}


// Create a static histogram object
Histogram::HistObj* Histogram::createStaticHistogram(double binRange[2], int numBins,              std::string name, Histogram::spacing_t spacingType)
{
	// Create histogram object
	HistObj* hobj = new HistObj;

	// Set the spacing type
	hobj->spacingType = spacingType;

	// Set the bin range
	hobj->binRange[0] = binRange[0];
	hobj->binRange[1] = binRange[1];

	// Set the number of bins
	hobj->nBins = numBins;


	// Do some sanity checking on combinations of parameters
	if( spacingType == Histogram::SPACING_LOG )
	{
		if(binRange[0]<=0 || binRange[1]<=0)
		{
			std::cerr<<"[Histogram::createStaticHistogram] Fixed data extents are <=0! Exiting.\n\tbinRange[0] = "<<binRange[0]<<"\n\tbinRange[1] = "<<binRange[1]<<std::endl;
			exit(-1);
		}
	}


	// Create the bin edges
	double binWidth = Histogram::computeBinWidth(hobj->binRange, hobj->nBins, spacingType);
	hobj->binWidth = binWidth;
	hobj->binEdges = new double[hobj->nBins+1];
	Histogram::generateBinEdges(hobj->binEdges, hobj->binRange, hobj->nBins, hobj->binWidth, hobj->spacingType);

	// Create the count array
	hobj->binCount = new int[hobj->nBins];
	for(int i=0; i<hobj->nBins; i++)
	{
		hobj->binCount[i] = 0;
	}

	// Set the name
	hobj->name = name;

	// Add to the container
	this->histObjs.push_back(hobj);

	// Return the HistObj
	return hobj;
}

// Create a static histogram object
Histogram::HistObj* Histogram::createStaticHistogram(double binRange[2], double binWidthSuggested, std::string name, Histogram::spacing_t spacingType)
{
	// Create histogram object
	HistObj* hobj = new HistObj;

	// Set the spacing type
	hobj->spacingType = spacingType;

	// Set the bin range
	hobj->binRange[0] = binRange[0];
	hobj->binRange[1] = binRange[1];

	// Do some sanity checking on combinations of parameters
	if( spacingType == Histogram::SPACING_LOG )
	{
		if(binRange[0]<=0 || binRange[1]<=0)
		{
			std::cerr<<"[Histogram::createStaticHistogram] Fixed data extents are <=0! Exiting.\n\tbinRange[0] = "<<binRange[0]<<"\n\tbinRange[1] = "<<binRange[1]<<std::endl;
			exit(-1);
		}
	}

	// Determine the number of bins
	int numBins = Histogram::computeNumBins(hobj->binRange, binWidthSuggested, spacingType);
	if(numBins==0)
	{
		std::cout<<"[Histogram::createStaticHistogram] No bins. Populating with fake data."<<std::endl;
		numBins = 1;
	}

	// Set the number of bins
	hobj->nBins = numBins;

	// Create the bin edges
	double binWidth = Histogram::computeBinWidth(hobj->binRange, hobj->nBins, spacingType);
	hobj->binWidth = binWidth;
	hobj->binEdges = new double[hobj->nBins+1];
	Histogram::generateBinEdges(hobj->binEdges, hobj->binRange, hobj->nBins, hobj->binWidth, hobj->spacingType);

	// Create the count array
	hobj->binCount = new int[hobj->nBins];
	for(int i=0; i<hobj->nBins; i++)
	{
		hobj->binCount[i] = 0;
	}

	// Set the name
	hobj->name = name;

	// Add to the container
	this->histObjs.push_back(hobj);

	// Return the HistObj
	return hobj;
}


// Determine the bin width from the number of bins
double Histogram::computeBinWidth(double binRange[2], int numBins, Histogram::spacing_t spacingType)
{

	double binWidth = -1.0;
	if(spacingType == Histogram::SPACING_LIN)
	{
		binWidth = (binRange[1]-binRange[0])/numBins;
		return binWidth;
	}
	else if(spacingType == Histogram::SPACING_LOG)
	{		
		binWidth = (log10(binRange[1])-log10(binRange[0]))/numBins;
		return binWidth;
	}
	else
	{
		std::cerr<<"[Histogram::computeBinWidth] ERROR: Spacing type not supported. Exiting."<<std::endl;
		exit(-1);
	}

	return binWidth;
}

int Histogram::computeNumBins(double binRange[2], double binWidth, Histogram::spacing_t spacingType)
{
	int numBins = -1;

	if(spacingType == Histogram::SPACING_LIN)
	{
		numBins = ceil((binRange[1]-binRange[0])/binWidth);
		return numBins;
	}
	else if(spacingType == Histogram::SPACING_LOG)
	{		
		numBins = ceil((log10(binRange[1])-log10(binRange[0]))/binWidth);
		return numBins;
	}
	else
	{
		std::cerr<<"[Histogram::computeNumBins] ERROR: Spacing type not supported. Exiting."<<std::endl;
		exit(-1);
	}

	return numBins;
}

// Generate bin edges
void Histogram::generateBinEdges(double* binEdges, double binRange[2], int numBins, double binWidth, Histogram::spacing_t spacingType)
{
	if(spacingType == Histogram::SPACING_LIN)
	{
		// Set the boundary values
		binEdges[0] = binRange[0];
		binEdges[numBins] = binRange[1];

		// Get the edges on the interior
		for(int i=1; i<numBins; i++)
		{
			binEdges[i] = double(i)*binWidth + binRange[0];
		}
	}
	else if(spacingType == Histogram::SPACING_LOG)
	{
		// Set the boundary values
		binEdges[0] = binRange[0];
		binEdges[numBins] = binRange[1];
		
		// Get the edges on the interior
		for(int i=1; i<numBins+1; i++)
		{
			// Determine the exponent of the bin edge
			double exponent = double(i)*binWidth + log10(binRange[0]);
			binEdges[i] = pow(10.0, exponent);
		}
	}
	else
	{
		std::cerr<<"[Histogram::generateBinEdges] ERROR: Spacing type not supported. Exiting."<<std::endl;
		exit(-1);
	}
}


// Update a histogram object (online update)
int Histogram::addDatum(Histogram::HistObj* hobj, double datum)
{
	// Get the bin the datum lies in
	int binInd = Histogram::getBinIndex(hobj, datum);

	// Increase the binCount counter
	hobj->binCount[binInd]++;

	// Return the bin index
	return binInd;
}

// Get the bin index that a given data point will lie in
int Histogram::getBinIndex(Histogram::HistObj* hobj, double datum)
{
	int binInd = -1;
	if(hobj->spacingType == Histogram::SPACING_LIN)
	{
		// Get bin index
		binInd = int((datum - hobj->binRange[0])/hobj->binWidth);
	}
	else if(hobj->spacingType == Histogram::SPACING_LOG)
	{
		double log10Datum = log10(datum);
		binInd = int((log10Datum - log10(hobj->binRange[0]))/hobj->binWidth);
	}
	else
	{
		std::cerr<<"[Histogram::getBinIndex] ERROR: Spacing type not supported. Exiting."<<std::endl;
		exit(-1);
	}

	// If we're outside the range, lump it into the first or last bin
	if(binInd<0)
	{
		binInd = 0;
	}
	if(binInd>=hobj->nBins)
	{
		binInd = hobj->nBins-1;
	}

	return binInd;
}


// Fix the bin width to binWidth
void Histogram::processFixedWidth(double* data, int numDataPts, double binWidth, std::string name, Histogram::spacing_t spacingType)
{

	// Compute the range of the data
	double minVal, maxVal;
	computeMinMax(data, numDataPts, minVal, maxVal);
	double binRange[] = {minVal,maxVal};

	HistObj* hobj = this->createStaticHistogram(binRange, binWidth, name, spacingType);
	for(int d=0; d<numDataPts; d++)
	{
		this->addDatum(hobj, data[d]);
	}

}

void Histogram::processStaticHistogram(double* data, int numDataPts, double binRange[2], int numBins, std::string name, Histogram::spacing_t spacingType)
{

	if(numBins==0)
	{
		std::cerr<<"[Histogram::processStaticHistogram] ERROR: Can't have 0 bins! Exiting."<<std::endl;
		exit(-1);
	}

	HistObj* hobj = this->createStaticHistogram(binRange, numBins, name, spacingType);
	for(int d=0; d<numDataPts; d++)
	{
		this->addDatum(hobj, data[d]);
	}

}


// Compute the number of bins as k=sqrt(numDataPts), and use this to determine the fixed bin width
void Histogram::processSquareRootChoice(double* data, int numDataPts, std::string name)
{
	std::cerr<<"[Histogram::processSquareRootChoice] ERROR: Not implemented."<<std::endl;
	exit(-1);
}


void Histogram::printObject(std::ostream& out)
{
	char delim = '\t';
	std::streamsize prec = 15;
	
	for (std::vector<HistObj*>::iterator i = this->histObjs.begin(); i != this->histObjs.end(); ++i)
	{
		int nBins = (*i)->nBins;

		// Create header
		std::stringstream header;
		header<<(*i)->name<<delim<<nBins;

		// Write header to stream
		out<<header.str()<<std::endl;

		// Write data lines
		for(int binId = 0; binId < nBins; binId++)
		{
			std::stringstream dataLine;

			// Set formatting
			dataLine.precision(prec);
			dataLine.setf(std::ios::scientific,std::ios::floatfield);

			// Create data line
			dataLine<<(*i)->binEdges[binId]<<delim<<(*i)->binEdges[binId+1]<<delim<<(*i)->binCount[binId];

			// write dataLine to ostream
			out<<dataLine.str()<<std::endl;

		}		
	}

}


};