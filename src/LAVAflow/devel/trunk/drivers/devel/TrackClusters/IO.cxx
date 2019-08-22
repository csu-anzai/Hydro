#include "IO.h"
#include <iostream>
#include "LavaMPI.h"

using namespace libconfig;

int ReadConfigFile(ConfigData &data, std::string filename)
{
    Config cfg;
    
    // Read the file. If there is an error, report it and exit.
    try
    {
      cfg.readFile(filename.c_str());
    }
    catch(const FileIOException &fioex)
    {
      std::cerr << "I/O error while reading file." << std::endl;
      return(EXIT_FAILURE);
    }
    catch(const ParseException &pex)
    {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                << " - " << pex.getError() << std::endl;
      return(EXIT_FAILURE);
    }

    // Look up values.
    if (!(cfg.lookupValue("unmixedClusterTol", data.unmixedClusterTol) &&
          cfg.lookupValue("clusterIDTol", data.clusterIDTol) &&
          cfg.lookupValue("thresholdVarCutoff", data.thresholdVarCutoff) &&
          cfg.lookupValue("minClusterSize", data.minClusterSize) &&
          cfg.lookupValue("backgroundClusterID", data.backgroundClusterID) &&

          cfg.lookupValue("filePath", data.filePath) &&
          cfg.lookupValue("meshBasename", data.meshBasename) &&
          cfg.lookupValue("partBasename", data.partBasename) &&
          cfg.lookupValue("originalClusterFilename", data.originalClusterFilename) &&
          cfg.lookupValue("massTracerDataFilename", data.massTracerDataFilename) &&
          cfg.lookupValue("particleDataFilename", data.particleDataFilename) &&

          cfg.lookupValue("trackClustersMassTracers", data.trackClustersMassTracers) &&
          cfg.lookupValue("trackClustersParticlesMT", data.trackClustersParticlesMT) &&
          cfg.lookupValue("trackClustersParticles", data.trackClustersParticles) &&
          cfg.lookupValue("saveNewClusterMeshes", data.saveNewClusterMeshes)
          ))
    {
        std::cout << "Could not read in all values" << std::endl;
        return EXIT_FAILURE;
    }

    // make sure that the path ends with a slash
    if (data.filePath.back() != '/')
      data.filePath.push_back('/');

    // std::cout << "File path is: " << data.filePath << std::endl;

    const Setting& meshFileIndexBoundsSetting = cfg.getRoot()["meshFileIndexBounds"];
    
    // std::cout << "Number of items in array: " << meshFileIndexBoundsSettingLength << std::endl;
    for (int i = 0; i < 2; i++)
      data.meshFileIndexBounds[i] = (int) meshFileIndexBoundsSetting[i];

    const Setting& partFileIndexBoundsSetting = cfg.getRoot()["partFileIndexBounds"];
    for (int i = 0; i < 2; i++)
      data.partFileIndexBounds[i] = (int) partFileIndexBoundsSetting[i];

    return EXIT_SUCCESS;
}

void OutputOnRoot(std::string outString)
{
  if (MPI::COMM_WORLD.Get_rank()==0)
    std::cout << outString << std::endl;
}
