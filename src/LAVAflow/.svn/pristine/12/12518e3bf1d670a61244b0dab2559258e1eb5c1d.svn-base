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

using namespace std;

// int*** BuildMaskArray(string , const double , std::array<int, MESH_MDIM>& );
double GetIgnitionTime(double , double , double );
bool MaskCritIgTime(double**** , double, int[], int, int, int);


int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{

    LavaMPI mpiObj(argc, argv);

    //------------------Do analysis------------------------------------|

    // only do on root
    if (MPI::COMM_WORLD.Get_rank()==0)
    {
        string filename("/data1/df11c/data/data/tburn_hdf5_plt_cnt_0002");
        const double thresholdVarCutoff = 10.;
        const int minClusterSize = 10;

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


        // output the size of each cluster
        int numClusters = floodClusterObj.getNumberOfClusters();
        vector<int> clusterSizes;
        clusterSizes.reserve(numClusters);

        for (int i=0; i<numClusters; i++)
        {
            int clusterSize = floodClusterObj.getClusterSize(i);
            // if (clusterSize >minClusterSize)
            clusterSizes.push_back(clusterSize);
        }

        ofstream outClusterSize("clusterSizes.dat", ios::out);

        // cout << "Cluster sizes: " << endl;

        for (auto i : clusterSizes)
            outClusterSize << i << endl;

        outClusterSize.close();
    } // only do on root
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

    return GetIgnitionTime(temp, dens, c12) < thresholdVarCutoff;
}
