#include <string>
#include <array>
#include "libconfig.h++"

class ConfigData
{
  public:

    double unmixedClusterTol = 1e-4;
    double clusterIDTol = 0.2;
    double thresholdVarCutoff = 1.;
    int minClusterSize = 10000;
    int backgroundClusterID = 1;

    std::string filePath = "."; //
    std::string meshBasename = "tburn_hdf5_plt_cnt_"; //
    std::string partBasename = "tburn_hdf5_part_"; //
    std::string originalClusterFilename = "markedClusters"; //
    std::string massTracerDataFilename = "mtClusterAnalysis.dat";
    std::string particleDataFilename = "particleClusterAnalysis.dat";

    std::array<int,2> meshFileIndexBounds; //
    std::array<int,2> partFileIndexBounds; //

    bool trackClustersMassTracers = true;
    bool trackClustersParticlesMT = false;
    bool trackClustersParticles = true;
    bool saveNewClusterMeshes = false;
};


int ReadConfigFile(ConfigData &, std::string );
void OutputOnRoot(std::string );
