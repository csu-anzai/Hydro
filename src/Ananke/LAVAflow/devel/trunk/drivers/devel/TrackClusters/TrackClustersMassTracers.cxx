#include <cmath>
#include <unistd.h>
#include <vector>
#include <numeric>
#include "Driver_main.h"
#include "ArrayOperations.h"
#include "FlashAmrMesh.h"
#include "FloodFill.h"
#include "LavaMPI.h"
#include "FlashParticles.h"
#include "IndexOperations.h"
#include "TrackClusters.h"

using namespace std;

// const double unmixedClusterTol = 1e-4;
// const double clusterIDTol = 0.2;
// const double thresholdVarCutoff = 10.;
// const int minClusterSize = 100;
// const int backgroundClusterID = 1;

double GetIgnitionTime(double , double , double );
bool MaskCritIgTime(double**** , double, int[], int, int, int);
void SaveClusterMesh(FlashAmrMesh& , Clustering::FloodFill& );


void TrackClustersMassTracers(std::vector<string> &meshFilenames, std::string originalClusterFilename, ConfigData &cfgData)
{
    
    const double unmixedClusterTol = cfgData.unmixedClusterTol;
    const double clusterIDTol = cfgData.clusterIDTol;
    const double thresholdVarCutoff = cfgData.thresholdVarCutoff;
    const int minClusterSize = cfgData.minClusterSize;
    const int backgroundClusterID = cfgData.backgroundClusterID;
    const bool saveNewClusterMeshes = cfgData.saveNewClusterMeshes;
    std::string massTracerDataFilename = cfgData.massTracerDataFilename;
    std::string filePath = cfgData.filePath;

    std::vector<std::string> calcClustersMeshVars = {"dens", "temp", "c12 "};
    std::vector<std::string> clusterMeshVar = {"clst", "dens"};
    
    std::vector<std::string> auxVars;
    if (saveNewClusterMeshes)
        auxVars.push_back("nwcr");
        
    std::vector<double> timeSeps;
    
    // load the mesh that we're comparing against (the one marking the clusters)
    // use it to get the initial simulation time, as well as the number of clusters
    int nOldClusters;
    double origSimTime;
    std::vector<double> origClusterMasses;
    GetOrigMeshInfo(originalClusterFilename, nOldClusters, origSimTime, origClusterMasses);

    // cout << "Original cluster masses: " << '\n';
    // for (auto clusterMass : origClusterMasses)
    //     cout << clusterMass << '\n';

    // cout << endl;
    
    std::vector<std::vector<double>> clusterPrimaryCompVol;
    std::vector<std::vector<double>> clusterPrimaryCompSharedMass;
    std::vector<std::vector<double>> clusterPrimaryCompSharedMassFrac;
    std::vector<std::vector<int>> clusterPrimaryCompID;
    std::vector<std::vector<double>> clusterStdDevsTracerID;
    std::vector<std::vector<double>> clusterMeansTracerID;

    // for each output file
    int index = 0;
    for (auto thisFile : meshFilenames)
    {
        cout << "Processing file " << thisFile << "..." << endl;

        // calc clusters in current dataset
        Clustering::FloodFill floodClusterObj(thisFile, thresholdVarCutoff, &MaskCritIgTime, calcClustersMeshVars);
        const int nDim = floodClusterObj.GetNDim();

        for (int i=0; i<nDim; i++)
        {
            for (int bound=0; bound<2; bound++)
                floodClusterObj.setBoundaryType(i,bound,Clustering::FloodFill::PERIODIC);
        }

        floodClusterObj.SetMinClusterSize(minClusterSize);
        floodClusterObj.process();
        int nClusters = floodClusterObj.getNumberOfClusters();

        // now FLASH mesh is unloaded. the only mesh in memory is the cluster mesh        
        // load old cluster data (tracers) from mesh
        std::string outFilename = thisFile+"newCluster";
        FlashAmrMesh mesh(thisFile, outFilename, clusterMeshVar, auxVars);
        

        int clusterIDMeshIndex = mesh.findVarIndex("clst");
        int densMeshIndex = mesh.findVarIndex("dens");

        // some mesh info
        auto blockLinearIndexMap = BuildBlockLinearIndexMap(mesh);
        int nCellsBlockX = mesh.getInt("nxb"),
            nCellsBlockY = mesh.getInt("nyb"),
            nCellsBlockZ = mesh.getInt("nzb");

        int nBlocksI = mesh.getInt("nblockx");

        double deltaT = mesh.getTime() - origSimTime;
        timeSeps.push_back(deltaT);

        // only works on a uniform mesh!
        int firstCell[] = {0,0,0};
        double cellVol = mesh.getSingleCellVol(0, INTERIOR, firstCell);

        // double sideLengths[MESH_MDIM];
        // mesh.getCellSideLengths(0, sideLengths);

        // double blockLengths[MESH_MDIM];
        // mesh.getBlkPhysicalSize(0, blockLengths);

        // double tempSize = mesh.getSize(0,IAXIS);
        // double tempSize2 = mesh.getSize(0,JAXIS);
        // double tempSize3 = mesh.getSize(0,KAXIS);
            
        int nNewClusters = floodClusterObj.getNumberOfClusters() + 2;
        // int nOldClusters = numPartsByTag.size() - 2;

        // assume that the old and new meshes are exactly the same (ie--no refinement changes)
        // check each new cluster against groups of mass tracers from old runs
        // std::vector<std::vector<int>> uniqueIDsPerCluster(nNewClusters);
        std::vector<double> primaryTracerIDVol(nNewClusters);
        std::vector<double> primaryTracerIDMass(nNewClusters);
        std::vector<double> primaryTracerIDMassFrac(nNewClusters);
        std::vector<int> primaryTracerID(nNewClusters);
        std::vector<double> meanTracerID(nNewClusters);
        std::vector<double> stdDevTracerID(nNewClusters);

        // for each cluster in new mesh
        int clusterIDindex = 2;
        for (auto thisCluster : floodClusterObj.clusters)
        {
            // cout << "Cluster " << clusterIDindex << " has " << thisCluster.size() << " elements in new mesh." << endl;

            std::vector<double> massTracerVals;
            std::vector<double> cellMasses;
            std::vector<int> uniqueClusterIDs;

            // double thisClusterVolume = 0;

            // cout << "Here1" << endl;

            // preliminary pass--calculate average and std dev of mass tracers that are > 1.5 (to not count mixed background material)
            // identify cluster IDs by assuming that at least one cell will be essentially unmixed
            for (auto thisClusterCell : thisCluster)
            {
                // thisClusterVolume += cellVol;

                // get the row-major linear block index, assuming an equal number of cells per dimension in the block, and a cubic domain
                int linBlockIndex = GetLinearBlockIndex( thisClusterCell.i,  thisClusterCell.j,  thisClusterCell.k, nCellsBlockX, nBlocksI);
                int blockID = blockLinearIndexMap[linBlockIndex];

                int blockI = thisClusterCell.i % nCellsBlockX,
                    blockJ = thisClusterCell.j % nCellsBlockY,
                    blockK = thisClusterCell.k % nCellsBlockZ;

                double ****solnData;
                //Get pointer to the block's data
                mesh.getBlkPtr(blockID, solnData, CENTER);

                double clstMassTracer = solnData[clusterIDMeshIndex][blockI][blockJ][blockK];
                double cellMass = solnData[densMeshIndex][blockI][blockJ][blockK] * cellVol;

                if (clstMassTracer > 1)
                {
                    massTracerVals.push_back(clstMassTracer);
                    cellMasses.push_back(cellMass);

                    // check if this value is very close to an integer value
                    // if so, we assume that it is unmixed
                    if (abs(clstMassTracer - round(clstMassTracer)) < unmixedClusterTol)
                    {
                        int thisClusterID = round(clstMassTracer);

                        // don't add to the list if this cluster doesn't exist in the original file
                        if (thisClusterID <= nOldClusters+2)
                        {
                            // only add to the list if not already there
                            if ( std::find(uniqueClusterIDs.begin(), uniqueClusterIDs.end(), thisClusterID) ==  uniqueClusterIDs.end())
                                uniqueClusterIDs.push_back(thisClusterID);
                        }

                    }

                }

            } // end loop over cluster cells

            // cout << "Here2" << endl;

            // get some statistics for this cluster
            // double sum = std::accumulate(massTracerVals.begin(), massTracerVals.end(), 0.0);
            // double mean = sum / massTracerVals.size();

            // std::vector<double> diff(massTracerVals.size());
            // std::transform(massTracerVals.begin(), massTracerVals.end(), diff.begin(), [mean](double x) { return x - mean; });
            // double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            // double stdev = std::sqrt(sq_sum / massTracerVals.size());


            double sum = std::accumulate(std::begin(massTracerVals), std::end(massTracerVals), 0.0);
            double mean =  sum / massTracerVals.size();

            double accum = 0.0;
            std::for_each (std::begin(massTracerVals), std::end(massTracerVals), [&](const double d) {
                accum += (d - mean) * (d - mean);
            });

            double stdev = sqrt(accum / (massTracerVals.size()-1));


            meanTracerID[clusterIDindex] = mean;
            stdDevTracerID[clusterIDindex] = stdev;

            // cout << "Mean, std dev for mass tracer values in cluster " << clusterIDindex << ": " << mean << ", " << stdev << endl;
            // cout << "Unique mass tracer values in this cluster: ";
            // for (auto clusterID : uniqueClusterIDs)
            //     cout << clusterID << ", ";

            // cout << endl;

            // cout << "Number of old clusters: " << nOldClusters << endl;
            std::vector<double> tracerIDMass(nOldClusters+2, 0.);
            std::vector<double> tracerIDVol(nOldClusters+2, 0.);

            // cout << "Here3" << endl;

            // calculate the mass, vol comprising each prominent set of mass tracer values
            // loop over unique IDs in the cluster
            for (auto thisMassTracerID : uniqueClusterIDs)
            {
                // loop over all this cluster's cells, and compare with the ID
                int clusterCellIndex = 0;
                for (auto &thisClusterCellID : massTracerVals)
                {
                    if (abs(thisClusterCellID - thisMassTracerID) < clusterIDTol*thisMassTracerID)
                    {
                        // add the mass, vol to the total
                        if (clusterCellIndex > cellMasses.size() - 1)
                        {
                            cout << "Out of bounds error in cellMasses" << endl;
                        }

                        if (thisMassTracerID > tracerIDMass.size() - 1)
                        {
                            cout << "Out of bounds error in tracerIDMass" << endl;
                            cout << "Size of tracerIDMass: " << tracerIDMass.size() << endl;
                            cout << "Index: " << thisMassTracerID << endl;
                        }

                        tracerIDMass[thisMassTracerID] += cellMasses[clusterCellIndex];
                        tracerIDVol[thisMassTracerID] += cellVol;

                        // mark it so we don't double count it
                        thisClusterCellID = -1;
                    }
                    clusterCellIndex++;
                }
                
            } 
            
            // cout << "Size of massFracsThisCluster: " << tracerIDMass.size() << endl;
            std::vector<double> massFracsThisCluster(tracerIDMass.size());
            // cout << "Volumes of original clusters present in this cluster (normalized to original cluster volume): " << endl;
            
            
            double maxIDMassThisCluster = 0.;
            int tracerIDMaxMass = -1;

            // cout << "Number of unique clusters: " << uniqueClusterIDs.size() << endl;

            
            for (auto thisMassTracerID : uniqueClusterIDs)
            {
                // cout << "thisMassTracerID: " << thisMassTracerID << endl;
                double clusterFrac = tracerIDMass[thisMassTracerID] / origClusterMasses[thisMassTracerID];
                // cout << "Cluster #" << thisMassTracerID << ": " << clusterFrac << endl;
                if (thisMassTracerID > tracerIDMass.size() - 1)
                {
                    cout << "Error: out of bounds." << endl;
                }

                massFracsThisCluster[thisMassTracerID] = clusterFrac;

                // cout << "Mass: " << tracerIDMass[thisMassTracerID] << endl;

                if (tracerIDMass[thisMassTracerID] > maxIDMassThisCluster)
                {
                    maxIDMassThisCluster = tracerIDMass[thisMassTracerID];
                    tracerIDMaxMass = thisMassTracerID;
                }
            }
            
            // massFracsPerCluster[clusterIDindex] = massFracsThisCluster;

            // cout << "tracerIDMaxMass: " << tracerIDMaxMass << endl;
            
            if (tracerIDMaxMass > 0)
            {
                primaryTracerID[clusterIDindex] = tracerIDMaxMass;
                primaryTracerIDMassFrac[clusterIDindex] = massFracsThisCluster[tracerIDMaxMass];
                primaryTracerIDMass[clusterIDindex] = tracerIDMass[tracerIDMaxMass];
                primaryTracerIDVol[clusterIDindex] = tracerIDVol[tracerIDMaxMass];
            }
            else
            {
                primaryTracerID[clusterIDindex] = -1;
                primaryTracerIDMassFrac[clusterIDindex] = -1;
                primaryTracerIDMass[clusterIDindex] = -1;
                primaryTracerIDVol[clusterIDindex] = -1; 
            }
            
            clusterIDindex++;
            
        } // end loop over clusters

        clusterPrimaryCompVol.push_back(primaryTracerIDVol);
        clusterPrimaryCompSharedMass.push_back(primaryTracerIDMass);        
        clusterPrimaryCompSharedMassFrac.push_back(primaryTracerIDMassFrac);
        clusterPrimaryCompID.push_back(primaryTracerID);
        clusterStdDevsTracerID.push_back(stdDevTracerID);
        clusterMeansTracerID.push_back(meanTracerID);
        
        
        if (saveNewClusterMeshes)
            SaveClusterMesh(mesh, floodClusterObj);
        
        index++;
    } // end loop over files

    ofstream clusterMassFrac(filePath+massTracerDataFilename, ios::out);

    const int colWidth = 22;

    // column headers
    const char *vinit[] = {"deltaT", "origClusterID", "remainingMass", "remainingMassFrac", "meanClusterDens", "meanTracerID", "stdDevTracerID"};
    std::vector<std::string> colHeaders(vinit, end(vinit)); // definition

    clusterMassFrac << std::left;
    for (auto thisColHeader : colHeaders)
        clusterMassFrac << std::setw(colWidth) << thisColHeader.c_str();

    clusterMassFrac << '\n';

    for (int fileIndex=0; fileIndex < clusterPrimaryCompID.size(); fileIndex++)
    {
        for (int i=2; i<clusterPrimaryCompID[fileIndex].size(); i++)
        {
            if (clusterPrimaryCompID[fileIndex][i] > 1)
            clusterMassFrac <<   std::left << 
                                std::setw(colWidth) << timeSeps[fileIndex] <<   
                                std::setw(colWidth) << clusterPrimaryCompID[fileIndex][i] <<
                                std::setw(colWidth) << clusterPrimaryCompSharedMass[fileIndex][i] << 
                                std::setw(colWidth) << clusterPrimaryCompSharedMassFrac[fileIndex][i] << 
                                std::setw(colWidth) << clusterPrimaryCompSharedMass[fileIndex][i] / clusterPrimaryCompVol[fileIndex][i] <<
                                std::setw(colWidth) << clusterMeansTracerID[fileIndex][i] << 
                                std::setw(colWidth) << clusterStdDevsTracerID[fileIndex][i] <<
                                endl;
        }
    }

    clusterMassFrac.close();

    cout << "Wrote mass tracer analysis to " << filePath+massTracerDataFilename << '\n';
    
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

    double igTime = GetIgnitionTime(temp, dens, c12);

    return igTime < thresholdVarCutoff && igTime > 0;
}


void SaveClusterMesh(FlashAmrMesh &mesh, Clustering::FloodFill &floodClusterObj)
{
    const int nDim = floodClusterObj.GetNDim();

    int clusterIDMeshIndex = mesh.findVarIndex("nwcr");

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

    for (int lb=0; lb<numLeafBlocks; lb++)
    {
        int currentBlock = blockList[lb];
        double ****solnData;
        //Get pointer to the block's data
        mesh.getBlkPtr(currentBlock, solnData, CENTER);

        int blkLimits[2][MESH_MDIM];
        int blkLimitsGC[2][MESH_MDIM];

        mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

        for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++){
            for (int j=blkLimits[LOWER][JAXIS]; j<= blkLimits[UPPER][JAXIS]; j++){
                for (int k=blkLimits[LOWER][KAXIS]; k<= blkLimits[UPPER][KAXIS]; k++){
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


                    solnData[clusterIDMeshIndex][i][j][k] = floodClusterObj.getPixelValue(iClusterMesh, jClusterMesh, kClusterMesh);
                }
            }
        }

        mesh.putBlkData(currentBlock, CENTER, "nwcr", INTERIOR, startingPos, solnData);
    }

    delete[] blockList;

    mesh.writeOutFile();
}
