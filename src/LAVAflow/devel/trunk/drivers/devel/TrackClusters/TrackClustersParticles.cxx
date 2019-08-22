#include "TrackClusters.h"
#include "LagTimeAutoCorr.h"

using namespace std;

const int colWidth = 22;


void TrackClustersParticles(std::string originalClusterFilename, ConfigData &cfgData)
{
    // get cluster unique IDs
    // load the mesh that we're comparing against (the one marking the clusters)
    // use it to get the initial simulation time, as well as the number of clusters
    int nOldClusters;
    double origSimTime;
    std::vector<double> origClusterMasses;
    GetOrigMeshInfo(originalClusterFilename, nOldClusters, origSimTime, origClusterMasses);

    // TEMPORARY!!!
    // nOldClusters = 1;

    std::string fullFilename = cfgData.filePath + cfgData.partBasename;

    // std::vector<std::string> partVars = {"velx"};
    std::vector<std::string> partVars = {"velx", "vely", "velz"};
    std::vector<const char*> meshVarsChar = VecStringToVecChar(partVars);
    string partMaskName = "clst";

    string fullOutFilename = cfgData.filePath+cfgData.particleDataFilename;
    ofstream outFile(fullOutFilename, ios::out);

    // column headers
    const char *vinit[] = {"clusterID", "deltaT", "autoCorr_velx", "autoCorr_vely", "autoCorr_velz"};
    std::vector<std::string> colHeaders(vinit, end(vinit)); // definition

    outFile << std::left;
    for (auto thisColHeader : colHeaders)
        outFile << std::setw(colWidth) << thisColHeader.c_str();

    outFile << '\n';

    // for each cluster ID
    int backGroundClusterID = cfgData.backgroundClusterID;
    int beginLoopIndex = backGroundClusterID+1;
    for (int i=beginLoopIndex; i<nOldClusters+beginLoopIndex; i++)
    {
        cout << "Processing cluster " << i << "..." << endl; 

        // get autocorrelation
        LagTimeAutoCorr autoCorr(fullFilename, &meshVarsChar[0], partVars.size(), cfgData.partFileIndexBounds.data());
        autoCorr.GetAutoCorrelations(i, partMaskName);
        

        int nTimeSeps = autoCorr.timeSeps.size();

        
        for (int t=0; t<nTimeSeps; t++)             
        {
            outFile << std::left << std::setw(colWidth) << i << std::setw(colWidth) << autoCorr.timeSeps[t];
            for (int dim=0; dim<partVars.size(); dim++)
                outFile << std::left << std::setw(colWidth) << autoCorr.autoCorrelations[dim][t];
            
            outFile << endl;
        }

        
    }

    outFile.close();

    cout << "Wrote particle analysis to " << fullOutFilename << '\n';

}