#include <string>
#include "libconfig.h++"

class ConfigData
{
  public:

    double radiusStep = 1e5;
    std::string filePath = ".";
    std::string inFilename;
    std::string outFilename = "clusterMassDistribution.dat";
};


int ReadConfigFile(ConfigData &, std::string );
void OutputOnRoot(std::string );
