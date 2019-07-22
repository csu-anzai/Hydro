#include <string>
#include "libconfig.h++"
#include <vector>


class ConfigData
{
  public:

    std::string filePath         = ".";
    std::vector<std::string> varNames; 

    bool calcEulStats       = false;
    bool calcLagStats       = false;
    bool calcACFuncs        = false;
    bool calcStructFuncs    = false;

    std::string eulBasename      = "flash_hdf5_plt_cnt_";
    int eulInitFileIndex    = 0;
    int eulFinalFileIndex   = 0;
    int nSamplePoints;
    
    std::string lagBasename      = "flash_hdf5_part_";
    int lagInitFileIndex    = 0;
    int lagFinalFileIndex   = 0;
};


int ReadConfigFile(ConfigData &, std::string );
void OutputOnRoot(std::string );
