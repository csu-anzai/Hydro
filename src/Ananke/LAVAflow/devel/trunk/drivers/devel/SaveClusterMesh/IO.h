#include <string>
#include "libconfig.h++"

class ConfigData
{
  public:

    double thresholdVarCutoff = 1.;
    int minClusterSize = 10000;
    std::string filePath = ".";
    std::string inFilename;
    std::string clusterFilename = "markedClusters";
};


int ReadConfigFile(ConfigData &, std::string );
void OutputOnRoot(std::string );
