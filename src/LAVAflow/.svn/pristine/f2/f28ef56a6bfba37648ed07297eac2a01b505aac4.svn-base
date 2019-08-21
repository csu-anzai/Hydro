#include "EulTimeAutoCorr.h"
#include "LagTimeAutoCorr.h"
#include "LavaMPI.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include "Driver_main.h"
#include "ArrayOperations.h"
#include "GetTurbStats.h"

using namespace std;

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{

    LavaMPI mpiObj(argc, argv);

    
    ConfigData cfgData;
    string cfgFile = "config";
    ReadConfigFile(cfgData, cfgFile);

    OutputOnRoot("\nLoaded parameters from config file.");


    //------------------Do analysis------------------------------------|

    std::vector<std::string> meshVars = cfgData.varNames;

    const int numVars = meshVars.size();
    std::vector<const char*> meshVarConverted = VecStringToVecChar(meshVars);

    if (cfgData.calcEulStats)
    {
      OutputOnRoot("\nPerforming Eulerian autocorrelation analysis.");
      string fullFilename = cfgData.filePath + cfgData.eulBasename;
      int fileIndexBounds[2] = {cfgData.eulInitFileIndex, cfgData.eulFinalFileIndex};
      EulTimeAutoCorr autoCorr(fullFilename, &meshVarConverted[0], numVars, fileIndexBounds, cfgData.nSamplePoints);
      autoCorr.GetAutoCorrelations();
      CollectAndOutput(&autoCorr, cfgData.filePath + "eulerianAutoCorr");
    }
    
    if (cfgData.calcLagStats)
    {
      OutputOnRoot("\nPerforming Lagrangian autocorrelation analysis.");
      string fullFilename = cfgData.filePath + cfgData.lagBasename;
      int fileIndexBounds[2] = {cfgData.lagInitFileIndex, cfgData.lagFinalFileIndex};
      LagTimeAutoCorr autoCorr(fullFilename, &meshVarConverted[0], numVars, fileIndexBounds);
      autoCorr.GetAutoCorrelations();
      CollectAndOutput(&autoCorr, cfgData.filePath + "lagrangianAutoCorr");
    }
        
}
