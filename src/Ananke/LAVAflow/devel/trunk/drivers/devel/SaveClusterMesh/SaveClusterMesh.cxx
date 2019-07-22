// This only works for uniformly-refined meshes right now

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include "Driver_main.h"
#include "ArrayOperations.h"
#include "FlashAmrMesh.h"
#include "FloodFill.h"
#include "LavaMPI.h"
#include "IO.h"

using namespace std;

double GetIgnitionTime(double , double , double );
bool MaskCritIgTime(double**** , double, int[], int, int, int);


int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{

    LavaMPI mpiObj(argc, argv);

    ConfigData cfgData;
    string cfgFile = "config";
    ReadConfigFile(cfgData, cfgFile);

    const double thresholdVarCutoff = cfgData.thresholdVarCutoff;
    const int minClusterSize = cfgData.minClusterSize;
    string filePath = cfgData.filePath;
    string inFilename = cfgData.inFilename;
    string clusterFilename = cfgData.clusterFilename;

    //------------------Do analysis------------------------------------|

    // only do on root
    if (MPI::COMM_WORLD.Get_size() > 1)
    {
        cerr << "This is a serial program. Rerun using only one processor." << endl;
        exit(1);
    }

    string filename(filePath+inFilename);

    cout << "Reading file " << filename << endl;
    cout << "Calculating clusters for ignition times less than " << thresholdVarCutoff << " s" << endl;

    // std::vector<std::string> meshVars = {"igtm"};
    std::vector<std::string> meshVars = {"dens", "temp", "c12 "};

    Clustering::FloodFill floodClusterObj(filename, thresholdVarCutoff, &MaskCritIgTime, meshVars);

    const int nDim = floodClusterObj.GetNDim();

    for (int i=0; i<nDim; i++)
    {
        for (int bound=0; bound<2; bound++)
            floodClusterObj.setBoundaryType(i,bound,Clustering::FloodFill::PERIODIC);
    }

    floodClusterObj.SetMinClusterSize(minClusterSize);
    floodClusterObj.process();

    int nClusters = floodClusterObj.getNumberOfClusters();
    
    cout << "Total number of clusters: " << nClusters << endl;

    // std::cout << "Number of cells per cluster: " << '\n';
    // for (auto thisCluster : floodClusterObj.clusters)
    //     cout << thisCluster.size() << '\n';

    // cout << endl;

    // save the marked clusters as a mesh var in a new hdf5 file

    string outFilename(filePath+clusterFilename);
    
    std::vector<std::string> emptyMeshVars  = {"dens"};
    std::vector<std::string> auxVars  = {"clst"};

    FlashAmrMesh mesh(filename, outFilename, emptyMeshVars, auxVars);

    int clusterIDMeshIndex = mesh.findVarIndex("clst");

    int numLocalBlocks = mesh.getLocalNumBlocks();
    int numLeafBlocks = 0;
    int *blockList = new int[numLocalBlocks];
    // double* coords[8];
    mesh.getListOfBlocks(LEAF, blockList, numLeafBlocks);

    int startingPos[MESH_MDIM];
    startingPos[0] = 0;
    startingPos[1] = 0;
    startingPos[2] = 0;

    double domainBoundBox[2][MESH_MDIM];
    mesh.getDomainBoundBox(domainBoundBox);

    double cellSideLengths[MESH_MDIM];
    mesh.getCellSideLengths(blockList[0], cellSideLengths);

    // std::vector<int> numCellsWritten(20, 0);

    for (int lb=0; lb<numLeafBlocks; lb++)
    {
        int currentBlock = blockList[lb];
        double ****solnData;
        //Get pointer to the block's data
        mesh.getBlkPtr(currentBlock, solnData, CENTER);

        int blkLimits[2][MESH_MDIM];
        int blkLimitsGC[2][MESH_MDIM];

        mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

        for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++)
        {
            for (int j=blkLimits[LOWER][JAXIS]; j<= blkLimits[UPPER][JAXIS]; j++)
            {
                for (int k=blkLimits[LOWER][KAXIS]; k<= blkLimits[UPPER][KAXIS]; k++)
                {
                    // mark cells
                    // get the indices of the corresponding cell in the image array
                    double cellBoundBox[2][MESH_MDIM];
                    mesh.getCellBoundBox(currentBlock, INTERIOR, i, j, k, cellBoundBox);

                    int indicesImageArray[nDim];

                    for (int l=0; l<nDim; l++)
                    {
                        indicesImageArray[l] = lrint((cellBoundBox[LOWER][l] - domainBoundBox[LOWER][l]) / cellSideLengths[l]);
                    }

                    int iClusterMesh = indicesImageArray[0],
                        jClusterMesh = indicesImageArray[1],
                        kClusterMesh = indicesImageArray[2];


                    int massTracerID = floodClusterObj.getPixelValue(iClusterMesh, jClusterMesh, kClusterMesh);
                    solnData[clusterIDMeshIndex][i][j][k] = massTracerID;

                    // numCellsWritten[massTracerID]++;
                }
            }
        }

        mesh.putBlkData(currentBlock, CENTER, "clst", INTERIOR, startingPos, solnData);
    }

    delete[] blockList;

    cout << "Wrote cluster data to file " << outFilename << endl;

    // std::cout << "Number of cells copied to mesh in memory: " << '\n';
    // int totalCellsWritten = 0;
    // for (auto thisCluster : numCellsWritten)
    // {
    //     cout << thisCluster << '\n';
    //     totalCellsWritten += thisCluster;
    // }

    // cout << "Cells unwritten: " << 512*512*512 - totalCellsWritten;

    // cout << endl;

    mesh.writeOutFile();
    
} 



double GetIgnitionTime(double temp, double dens, double c12Frac)
{
    double T9 = temp / 1e9;
    double rho8 = dens / 1e8;

    double f = pow(T9 - 0.214, -7.566);
    double igTime = 1.15e-5*pow(c12Frac*rho8,-1.9) * f * (1+1193*f);

    return igTime;
}

bool MaskCritIgTime(double**** solnData, double thresholdVarCutoff, int varIndices[], int i, int j, int k)
{
    double dens = solnData[varIndices[0]][i][j][k];
    double temp = solnData[varIndices[1]][i][j][k];
    double c12 = solnData[varIndices[2]][i][j][k];
    // double ignTime = solnData[0][i][j][k];

    return GetIgnitionTime(temp, dens, c12) < thresholdVarCutoff;
    // return ignTime < thresholdVarCutoff && ignTime > 0;
}
