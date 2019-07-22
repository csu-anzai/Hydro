#include "LagTimeAutoCorr.h"

LagTimeAutoCorr::LagTimeAutoCorr(std::string baseFilename, const char** varNames, int nVars, int fileIndexBounds[2], MPI::Intracomm meshComm)
    : TimeAutoCorr(baseFilename, varNames, nVars, fileIndexBounds, meshComm)
{    
    
    // figure out how many points to generate based on the MPI rank
    // load the first particle file
    // we are assuming that all the particles files contain the same number of particles
    const char** noVars = NULL;
    FlashParticles firstPartFile(filenames[0].c_str(), noVars, 0, PART_SORT_TAG);
    numSamples = firstPartFile.nParticlesLocal;
}

LagTimeAutoCorr::~LagTimeAutoCorr()
{
}

std::shared_ptr<FlashIO> LagTimeAutoCorr::LoadData(std::string &filename, std::vector<std::string> &meshVars, std::string maskVar, int maskID)
{
    // std::cout << "Loading particles" << std::endl;

    std::vector<const char*> varsChar = VecStringToVecChar(meshVars);
    std::shared_ptr<FlashIO> partPtr(new FlashParticles (filename.c_str(), &varsChar[0], meshVars.size(), PART_SORT_TAG, maskVar, maskID));
    double simTime = partPtr->getTime();

    // std::cout << "Sim time in LoadData: " << simTime << std::endl;
    return partPtr;
}


int LagTimeAutoCorr::GetQuantityAtPoints(int meshVar, std::shared_ptr<FlashIO> partFlashIO, double pointVals[], int maskID, std::string maskVar)
{
    std::shared_ptr<FlashParticles> particles = std::dynamic_pointer_cast<FlashParticles>(partFlashIO);
    double** partData = particles->getDataPtr();

    int maskVarID = -1;
    if (maskID > -1)
        maskVarID = particles->findVarIndex(maskVar.c_str());

    int numPts = 0;

    for (int i=0; i<particles->nParticlesLocal; i++)
    {
        if (maskID > -1)
        {
            if (partData[maskVarID][i] == maskID)
            {
                pointVals[numPts] = partData[meshVar][i];
                numPts++;
            }
        }
        else
        {
            numPts++;
            pointVals[i] = partData[meshVar][i];
        }
    }

    return numPts;

    
}
