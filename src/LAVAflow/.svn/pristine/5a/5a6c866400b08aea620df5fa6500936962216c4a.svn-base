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
    if (!(cfg.lookupValue("inFilePath", data.inFilePath) &&
          cfg.lookupValue("outFilePath", data.outFilePath) &&
          cfg.lookupValue("inFilename", data.inFilename) &&
          cfg.lookupValue("outFilename", data.outFilename)
          
          ))
    {
        std::cout << "Could not read in all values" << std::endl;
        return EXIT_FAILURE;
    }

    // make sure that the path ends with a slash
    if (data.inFilePath.back() != '/')
      data.inFilePath.push_back('/');

    if (data.outFilePath.back() != '/')
      data.outFilePath.push_back('/');

    // std::cout << "File path is: " << data.filePath << std::endl;

    const Setting& subdomainBoundsXSetting = cfg.getRoot()["subdomainBoundsX"];    
    // std::cout << "Number of items in array: " << meshFileIndexBoundsSettingLength << std::endl;
    for (int i = 0; i < 2; i++)
      data.subdomainBoundsX[i] = (double) subdomainBoundsXSetting[i];

    const Setting& subdomainBoundsYSetting = cfg.getRoot()["subdomainBoundsY"];    
    for (int i = 0; i < 2; i++)
      data.subdomainBoundsY[i] = (double) subdomainBoundsYSetting[i];

    const Setting& subdomainBoundsZSetting = cfg.getRoot()["subdomainBoundsZ"];    
    for (int i = 0; i < 2; i++)
      data.subdomainBoundsZ[i] = (double) subdomainBoundsZSetting[i];

    return EXIT_SUCCESS;
}

void OutputOnRoot(std::string outString)
{
  if (MPI::COMM_WORLD.Get_rank()==0)
    std::cout << outString << std::endl;
}
