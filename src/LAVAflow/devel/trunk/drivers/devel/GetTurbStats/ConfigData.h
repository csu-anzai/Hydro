#include <string>


class ConfigData
{
  public:

    string filePath         = ".";

    bool calcEulStats       = false;
    bool calcLagStats       = false;
    bool calcACFuncs        = false;
    bool calcStructFuncs    = false;

    string eulBasename      = "flash_hdf5_plt_cnt_";
    int eulInitFileIndex    = 0;
    int eulFinalFileIndex   = 0;
    int nSamplePoints;

    string lagBasename      = "flash_hdf5_part_";
    int lagInitFileIndex    = 0;
    int lagFinalFileIndex   = 0;
};