#ifndef LAVA_TURBULENCE_TIMEAUTOCORR_H
#define LAVA_TURBULENCE_TIMEAUTOCORR_H

#include <string>
#include <vector>
#include <array>
#include "../Mesh/FlashAmrMesh.h"
#include <memory>
#include "ArrayOperations.h"

class TimeAutoCorr
{
public:
	int numTimeSeps;
	int numVars;
	std::vector<std::string> meshVars;
	std::vector<double> timeSeps;
	std::vector<std::vector<double>> autoCorrelations;

	std::string baseFilename;
	int fileIndexBounds[2];
	std::string *filenames;

	int numSamples;
	
	MPI::Intracomm meshComm;

	// TimeAutoCorr();
	TimeAutoCorr(std::string, const char**, int, int[2], MPI::Intracomm = MPI::COMM_WORLD);
	~TimeAutoCorr();

	void GetAutoCorrelations(int=-1, std::string="");
	virtual std::shared_ptr<FlashIO> LoadData(std::string &, std::vector<std::string> &, std::string="", int=-1)=0;
	virtual void GeneratePoints(std::shared_ptr<FlashIO>)=0;
	virtual int GetQuantityAtPoints(int , std::shared_ptr<FlashIO>, double[] , int = -1, std::string = "")=0;
};


#endif
