#ifndef LAVA_TURBULENCE_LAG_TIMEAUTOCORR_H
#define LAVA_TURBULENCE_LAG_TIMEAUTOCORR_H

#include "TimeAutoCorr.h"
#include "FlashParticles.h"
#include "ArrayOperations.h"

class LagTimeAutoCorr : public TimeAutoCorr
{
    
public:
    LagTimeAutoCorr(std::string, const char**, int, int[2], MPI::Intracomm = MPI::COMM_WORLD);
    ~LagTimeAutoCorr();

    virtual int GetQuantityAtPoints(int , std::shared_ptr<FlashIO>, double[] , int = -1, std::string = "");
    // virtual void GetQuantityAtPoints(int , std::shared_ptr<FlashIO>, double[]);
    virtual void GeneratePoints(std::shared_ptr<FlashIO>) {};
    virtual std::shared_ptr<FlashIO> LoadData(std::string &, std::vector<std::string> &, std::string="", int=-1);

};

#endif