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
    if (!(cfg.lookupValue("filePath", data.filePath) &&
          cfg.lookupValue("calcEulStats", data.calcEulStats) &&
          cfg.lookupValue("calcLagStats", data.calcLagStats) &&
          cfg.lookupValue("calcACFuncs", data.calcACFuncs) &&
          cfg.lookupValue("calcStructFuncs", data.calcStructFuncs) &&
          cfg.lookupValue("eulBasename", data.eulBasename) &&
          cfg.lookupValue("eulInitFileIndex", data.eulInitFileIndex) &&
          cfg.lookupValue("eulFinalFileIndex", data.eulFinalFileIndex) &&
          cfg.lookupValue("nSamplePoints", data.nSamplePoints) &&
          cfg.lookupValue("lagBasename", data.lagBasename) &&
          cfg.lookupValue("lagInitFileIndex", data.lagInitFileIndex) &&
          cfg.lookupValue("lagFinalFileIndex", data.lagFinalFileIndex)
          ))
    {
        std::cout << "Could not read in all values" << std::endl;
        return EXIT_FAILURE;
    }

    // make sure that the path ends with a slash
    if (data.filePath.back() != '/')
      data.filePath.push_back('/');

    // std::cout << "File path is: " << data.filePath << std::endl;

    const Setting& varNamesSetting = cfg.getRoot()["varNames"];
    int varNamesSettingLength = varNamesSetting.getLength();

    // std::cout << "Number of items in array: " << varNamesSettingLength << std::endl;
    for (int i = 0; i < varNamesSettingLength; i++)
      data.varNames.push_back((const char*) varNamesSetting[i]);

    return EXIT_SUCCESS;
}

void OutputOnRoot(std::string outString)
{
  if (MPI::COMM_WORLD.Get_rank()==0)
    std::cout << outString << std::endl;
}
