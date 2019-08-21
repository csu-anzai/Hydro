#include "TimeAutoCorr.h" 
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <time.h>

using namespace std;



TimeAutoCorr::TimeAutoCorr(std::string baseFilename, const char** varNames, int nVars, int fileIndexBounds[2], MPI::Intracomm meshComm)
{
	this->meshComm = meshComm;
	this->baseFilename = baseFilename;
	this->fileIndexBounds[0] = fileIndexBounds[0];
	this->fileIndexBounds[1] = fileIndexBounds[1];
	numSamples = 0;
	numVars = nVars;

	// mesh variables
	meshVars = CharToVecOfString(varNames, numVars);

	// std::cout << "Character array: " << varNames[0][0] << endl;
	// std::cout << "Mesh vars to load: " << meshVars[0] << endl;
	
	//  get the number of time separations and figure out the filenames
	numTimeSeps = fileIndexBounds[1] - fileIndexBounds[0];
	filenames = new std::string[numTimeSeps+1];

	for (int i=fileIndexBounds[0], index=0; i <= fileIndexBounds[1]; i++, index++)
	{
		// create the filenames
		std::ostringstream fileNumStr;
    	fileNumStr << std::setw( 4 ) << std::setfill( '0' ) << i;
    	std::string filename = baseFilename + fileNumStr.str();
    	
    	filenames[index] = filename;
	}
	
	// reserve space for the time values and AC values
	timeSeps.resize(numTimeSeps);

	
	autoCorrelations.resize(numVars);
	for (auto &meshVar : autoCorrelations)
		meshVar.resize(numTimeSeps);

}

TimeAutoCorr::~TimeAutoCorr()
{
	delete[] filenames;
}



void TimeAutoCorr::GetAutoCorrelations(int maskID, std::string maskIDVar)
{
	// cout << "Getting acs for tag " << maskID << endl;

	// load initial mesh

	// if we're masking (for example, based on particle tag), load the masking variable as well
	std::vector<std::string> allLoadVars = meshVars;

	if (maskID > -1)
		allLoadVars.push_back(maskIDVar);
	
	std::cout << "Loading initial file..." << std::endl;
	std::shared_ptr<FlashIO> initData = LoadData(filenames[0], allLoadVars, maskIDVar, maskID);
	double initSimTime = initData->getTime();

	GeneratePoints(initData);
	
	for (int meshVar=0; meshVar<meshVars.size(); meshVar++)
	{
		std::cout << "Processing mesh var " << meshVars[meshVar] << std::endl;
		int meshVarIndex = initData->findVarIndex(meshVars[meshVar].c_str());
	
		// save the values at each point for the original mesh
		double *initPointVals = new double[numSamples];
		int actualNumSamples = GetQuantityAtPoints(meshVarIndex, initData, initPointVals, maskID, maskIDVar);

		// cout << "Number of particles with ID in initial file: " << actualNumSamples << endl;
		
		// for each separation, load the new mesh and compare with initial mesh
		for (int i = 0; i < numTimeSeps; i++)
		{
			// load second mesh
			std::vector<std::string> varName = {meshVars[meshVar]};

			if (maskID > -1)
				varName.push_back(maskIDVar); 

			std::cout << "Loading file at index " << i+1 << std::endl;

			clock_t t = clock();
			std::shared_ptr<FlashIO> currentData = LoadData(filenames[i+1], varName, maskIDVar, maskID);

			std::cout << "Elapsed time: " << (double)(clock() - t) / CLOCKS_PER_SEC << std::endl;

			// cout << "Current sim time: " << currentData->getTime() << endl;

			// get velocities at points
			double *currentPointVals = new double[actualNumSamples];
			// the zero here is correct, because we're only loading one variable at a time
			actualNumSamples = GetQuantityAtPoints(0, currentData, currentPointVals, maskID, maskIDVar); 

			// TEMPORARY!!!
			// actualNumSamples = 5;

			// cout << "Number of particles with ID in current file: " << actualNumSamples << endl;

			// Average of above normalized to vel^2
			double cumProduct = 0;
			double cumVelPt1Squared = 0;
			double cumVelPt2Squared = 0;
			for (int pointIndex=0; pointIndex<actualNumSamples; pointIndex++)
			{
				// cout << "Initial velocity at " << pointIndex << ": " << initPointVals[pointIndex] << '\n';
				// cout << "Current velocity at " << pointIndex << ": " << currentPointVals[pointIndex] << endl;

				cumProduct += initPointVals[pointIndex]*currentPointVals[pointIndex];
				cumVelPt1Squared += initPointVals[pointIndex]*initPointVals[pointIndex];
				cumVelPt2Squared += currentPointVals[pointIndex]*currentPointVals[pointIndex];
			}

			autoCorrelations[meshVar][i] = cumProduct / (std::sqrt(cumVelPt1Squared) * std::sqrt(cumVelPt2Squared));
			
			if (meshVar == 0)
				timeSeps[i] = currentData->getTime() - initSimTime;

			delete[] currentPointVals;
		}

		delete[] initPointVals;
	}
}


