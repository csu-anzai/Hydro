#ifndef LAVA_STATISTICS_HISTOGRAM_H
#define LAVA_STATISTICS_HISTOGRAM_H

#include <list>
#include <ostream>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkDataArray.h>

#include "../../includes/LAVAconstants.h"
#include "HistObj.h"

namespace Statistics
{



void computeMinMax(double* data, int numDataPts, double& minVal, double& maxVal);

class Histogram
{

public:
	enum spacing_t 
	{
		SPACING_LIN = 0, 
		SPACING_LOG = 1
	};

	struct HistObj
	{
		// Number of bins
		int nBins;
		// Bin range
		double binRange[2];
		// Bin spacing
		// NOTE: 	For linear spacingType, this is the bin width in physical coordinates
		//			For log spacingType, this is the bin width in log10(physical coordinates) [i.e. the width of the exponent]
		double binWidth;
		// Spacing type
		spacing_t spacingType;
		// Edges of bins (size: nBins+1)
		double* binEdges;
		// Histogram count for bins (size: nBins)
		int*	binCount;
		// Identifying name
		std::string name;
	};

private:
	std::vector<HistObj*> histObjs;

public:




	Histogram();
	~Histogram();

	HistObj* createStaticHistogram(double binRange[2], int numBins, std::string name, Histogram::spacing_t);
	HistObj* createStaticHistogram(double binRange[2], double binWidthSuggested, std::string name, Histogram::spacing_t);

	// Fix the bin width to binWidth
	void processFixedWidth(double* data, int numDataPts, double binWidth, std::string name, Histogram::spacing_t spacingType = Histogram::SPACING_LIN);
	// Fix the histogram range and the number of bins
	void processStaticHistogram(double* data, int numDataPts, double binRange[2], int numBins, std::string name, Histogram::spacing_t spacingType = Histogram::SPACING_LIN);
	// Compute the number of bins as k=sqrt(numDataPts), and use this to determine the fixed bin width
	void processSquareRootChoice(double* data, int numDataPts, std::string name);
	
	int computeNumBins(double binRange[2], double binWidth, Histogram::spacing_t spacingType);
	double computeBinWidth(double binRange[2], int numBins, Histogram::spacing_t spacingType);
	void generateBinEdges(double* binEdges, double binRange[2], int numBins, double binWidth, Histogram::spacing_t spacingType);
	int addDatum(Histogram::HistObj* hobj, double datum);
	int getBinIndex(Histogram::HistObj* hobj, double datum);

	void printObject(std::ostream& out);


};

};

#endif