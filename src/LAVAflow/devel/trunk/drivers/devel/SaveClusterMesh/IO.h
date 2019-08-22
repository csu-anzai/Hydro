#include <string>
#include "libconfig.h++"

class ConfigData
{
  public:

    double thresholdVarCutoff = 1.;
    int minClusterSize = 10000;
    std::string filePath = "/home/wkshop24/LAVAflow/data/wdStir_sol_data/";
    std::string inFilename = "wdStir_hdf5_chk_0017";
    std::string clusterFilename = "markedClusters";
};


int ReadConfigFile(ConfigData &, std::string );
void OutputOnRoot(std::string );
