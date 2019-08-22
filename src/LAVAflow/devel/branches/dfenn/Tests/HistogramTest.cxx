#include <iostream>
#include <cstdlib>
#include "../libsrc/Statistics/Histogram/Histogram.h"

using namespace Statistics;

int main(void)
{

	Histogram hist;
	int nDataPts = 10000;
	double* data = new double[nDataPts];
	double binRange[] = {1,1001};

	Histogram::HistObj* onlineHist = hist.createStaticHistogram(binRange, 30, "onlineHist", Histogram::SPACING_LOG);

	std::cout<<"onlineHist->spacingType: "<<onlineHist->spacingType<<std::endl;

	for(int i=0; i<nDataPts; i++)
	{
		data[i] = double(rand()%1001+1);

		hist.addDatum(onlineHist, data[i]);
	}



	hist.processFixedWidth(data, nDataPts, 25, "test");


	hist.printObject(std::cout);

	delete[] data;

	std::cout<<"Finished"<<std::endl;
}