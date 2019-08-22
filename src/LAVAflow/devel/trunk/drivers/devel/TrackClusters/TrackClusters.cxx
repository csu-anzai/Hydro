#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <array>
#include <numeric>
#include <iomanip>
#include <sstream>
#include "Driver_main.h"
#include "ArrayOperations.h"
#include "FlashAmrMesh.h"
#include "FloodFill.h"
#include "LavaMPI.h"
#include "FlashParticles.h"
#include "IndexOperations.h"
#include "TrackClusters.h"

using namespace std;

std::vector<std::string> GetFilenames(std::string , std::array<int,2>&);



int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{
    // only do on root
    if (MPI::COMM_WORLD.Get_size() > 1)
    {
        cerr << "This is a serial program. Rerun using only one processor." << endl;
        exit(1);
    }

    cout << "Beginning cluster tracking analysis..." << endl << endl;

    ConfigData cfgData;
    string cfgFile = "config";
    ReadConfigFile(cfgData, cfgFile);

    // figure out output file names
    // get the number of time separations and figure out the filenames

    std::string filePath = cfgData.filePath;
    std::string meshBasename = cfgData.meshBasename;
    std::string partBasename = cfgData.partBasename;

    std::string originalClusterFilename = filePath + cfgData.originalClusterFilename;
    std::string baseMeshFilename = filePath + meshBasename;
    std::string basePartFilename = filePath + partBasename;  
    
    std::array<int,2> meshFileIndexBounds = cfgData.meshFileIndexBounds;
    std::array<int,2> partFileIndexBounds = cfgData.partFileIndexBounds;
    
    auto meshFilenames = GetFilenames(baseMeshFilename, meshFileIndexBounds);
    auto partFilenames = GetFilenames(basePartFilename, partFileIndexBounds);

    
    bool trackClustersMassTracers = cfgData.trackClustersMassTracers;
    bool trackClustersParticlesMT = cfgData.trackClustersParticlesMT;
    bool trackClustersParticles = cfgData.trackClustersParticles;

    if (trackClustersMassTracers)
    {
        cout << "Performing mass tracer analysis..." << endl;
        TrackClustersMassTracers(meshFilenames, originalClusterFilename, cfgData);
        cout << "Completed mass tracer analysis." << endl << endl;
    }
    
    if (trackClustersParticlesMT)
    {
        cout << "Performing hydbrid particle / mass tracer analysis..." << endl;
        TrackClustersParticlesMT(partFilenames, meshFilenames, originalClusterFilename, cfgData);
        cout << "Completed hydbrid particle / mass tracer analysis." << endl << endl;
    }

    if (trackClustersParticles)
    {
        cout << "Performing particle analysis..." << endl;
        TrackClustersParticles(originalClusterFilename, cfgData);
        cout << "Completed particle analysis." << endl << endl;
    }

    cout << "Execution complete" << endl;

}

std::vector<std::string> GetFilenames(std::string baseFilename, std::array<int,2> &fileIndexBounds)
{
    int numTimeSeps = fileIndexBounds[1] - fileIndexBounds[0];
    std::vector<std::string> filenames;
    
    for (int i=fileIndexBounds[0]; i <= fileIndexBounds[1]; i++)
    {
        // create the filenames
        std::ostringstream fileNumStr;
        fileNumStr << std::setw( 4 ) << std::setfill( '0' ) << i;
        
        std::string meshFilename = baseFilename + fileNumStr.str();        
        filenames.push_back(meshFilename);
    }

    return filenames;
}
