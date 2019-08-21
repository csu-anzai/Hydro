#include <string>
#include <array>
#include "libconfig.h++"

class ConfigData
{
  public:

    std::string inFilePath = "."; //
    std::string inFilename = ""; //
    
    std::string outFilePath = "."; //
    std::string outFilename = ""; //
    

    std::array<double,2> subdomainBoundsX;
    std::array<double,2> subdomainBoundsY;
    std::array<double,2> subdomainBoundsZ;
};


int ReadConfigFile(ConfigData &, std::string );
void OutputOnRoot(std::string );
