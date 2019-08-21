#ifndef LAVA_TURBULENCE_EUL_TIMEAUTOCORR_H
#define LAVA_TURBULENCE_EUL_TIMEAUTOCORR_H

#include "TimeAutoCorr.h"

class EulTimeAutoCorr : public TimeAutoCorr
{
    std::vector<std::array<double, MESH_MDIM>> points;
    int *pointBlockIDs;
    std::vector<std::array<int, MESH_MDIM>> pointZoneIDs;

public:
    EulTimeAutoCorr(std::string, const char**, int, int[2], int, MPI::Intracomm = MPI::COMM_WORLD);
    ~EulTimeAutoCorr();

    virtual int GetQuantityAtPoints(int , std::shared_ptr<FlashIO>, double[] , int = -1, std::string = "");
    virtual void GeneratePoints(std::shared_ptr<FlashIO>);
    virtual std::shared_ptr<FlashIO> LoadData(std::string &, std::vector<std::string> &, std::string="", int=-1);

};

#endif