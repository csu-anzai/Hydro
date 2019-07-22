#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <deque>
#include <iomanip>
#include "Driver_main.h"
#include "ArrayOperations.h"
#include "FlashAmrMesh.h"
#include "LavaMPI.h"
#include "IO.h"

using namespace std;

struct clusterBoundaryParameters
{   
    bool isBoundaryCluster = false;
    bool crossesBoundary[MESH_MDIM] = {false, false, false};
    double clusterBounds[2][MESH_MDIM];
};

void GetClusterIndices(FlashAmrMesh&, std::vector<int>&, std::vector<std::deque<std::array<int,MESH_MDIM+1>>>&, std::vector<clusterBoundaryParameters>&);
double GetEnclosedMass(std::array<double,MESH_MDIM> , double , std::deque<std::array<double,MESH_MDIM>>&, std::deque<double>&, int&);

void GetClusterCentroid(FlashAmrMesh&, std::deque<std::array<int,MESH_MDIM+1>>&, std::deque<std::array<double,MESH_MDIM>>&, 
                        std::deque<double>&, std::array<double,MESH_MDIM>&, clusterBoundaryParameters&);

void FindBoundaryClusterBounds(FlashAmrMesh&, std::deque<std::array<int,MESH_MDIM+1>>&, clusterBoundaryParameters&);





int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{
    LavaMPI mpiObj(argc, argv);

    // only do on root
    if (MPI::COMM_WORLD.Get_size() > 1)
    {
        cerr << "This is a serial program. Rerun using only one processor." << endl;
        exit(1);
    }

    ConfigData cfgData;
    std::string cfgFile = "config";
    ReadConfigFile(cfgData, cfgFile);

    std::string filePath = cfgData.filePath;
    std::string inFilename = cfgData.inFilename;
    std::string outFilename = cfgData.outFilename;
    const double radiusStep = cfgData.radiusStep;

    string filename(filePath+inFilename);

    std::vector<std::string> meshVars = {"clst", "dens"};


    FlashAmrMesh mesh(filename, meshVars);

    int clusterIDMeshIndex = mesh.findVarIndex("clst");

    // only works on a uniform mesh!
    int firstCell[] = {0,0,0};
    double cellVol = mesh.getSingleCellVol(0, INTERIOR, firstCell);
    
    // only works for uniform meshes 
    double sideLengths[MESH_MDIM];
    mesh.getCellSideLengths(0, sideLengths);

    double domainBoundBox[2][MESH_MDIM];
    mesh.getDomainBoundBox(domainBoundBox);

    
    // check that the radius steps are larger than a cell length
    if (radiusStep < sideLengths[IAXIS])
    {
        cerr << "Radius steps are smaller than the cell side length of " << sideLengths[IAXIS] << ". Increase the step size. Exiting..." << endl;
        exit(-1);
    }

    std::vector<int> uniqueClusterIDs;
    std::vector<std::deque<std::array<int,MESH_MDIM+1>>> clusterIndices(2);
    std::vector<clusterBoundaryParameters> boundaryParameters;
    
    GetClusterIndices(mesh, uniqueClusterIDs, clusterIndices, boundaryParameters);

    // now we have a list of all clusters present, as well as indices for each cluster
    // cout << "Cluster IDs present: " ;
    // for (auto thisClusterID : uniqueClusterIDs)
    // {
    //     cout << thisClusterID << ", " << '\n';
    //     // cout << "Boundary cluster? " << boundaryParameters[thisClusterID][0] << '\n';
    //     // cout << "Bound dimension: " << boundaryParameters[thisClusterID][1] << '\n';
    // }

    // cout << endl; 

    int numClusters = uniqueClusterIDs.size();
    std::vector<std::array<double,MESH_MDIM>> clusterCentroids(numClusters+2);
    std::vector<std::deque<std::array<double,MESH_MDIM>>> clusterCoords(numClusters+2);
    std::vector<std::deque<double>> clusterCellMasses(numClusters+2);
    std::vector<std::vector<double>> clusterShellMasses(numClusters+2);
    std::vector<std::vector<int>> clusterNumCellsShell(numClusters+2);

    // loop over clusters
    for (auto thisClusterID : uniqueClusterIDs)
    {
        cout << "Binning mass in cluster " << thisClusterID << endl;

        if (boundaryParameters[thisClusterID].isBoundaryCluster)
        {
            // cout << "Finding boundaries for cluster " << thisClusterID << endl;
            FindBoundaryClusterBounds(mesh, clusterIndices[thisClusterID], boundaryParameters[thisClusterID]);

            // for (int dim=0; dim < MESH_MDIM; dim++)
            // {
            //     cout << "Bounds for dimension " << dim << endl;
            //     cout << "LOWER: " << boundaryParameters[thisClusterID].clusterBounds[LOWER][dim] << endl;
            //     cout << "UPPER: " << boundaryParameters[thisClusterID].clusterBounds[UPPER][dim] << endl << endl;
            // }
        }
        

        GetClusterCentroid(mesh, clusterIndices[thisClusterID], clusterCoords[thisClusterID], 
                           clusterCellMasses[thisClusterID], clusterCentroids[thisClusterID], boundaryParameters[thisClusterID]);

        // cout << "Mass this cluster: " << totalClusterMass << endl;

        // now we have the cluster's centroid
        // find the enclosed mass as a function of radius from the centroid for various radii
        std::vector<double> shellMasses;
        std::vector<int> numCellsShell;

        double shellRadius = radiusStep;
        int numCellsThisShell;
        double massThisShell = GetEnclosedMass(clusterCentroids[thisClusterID], shellRadius, 
                                               clusterCoords[thisClusterID], clusterCellMasses[thisClusterID], numCellsThisShell);

        shellMasses.push_back(massThisShell);
        numCellsShell.push_back(numCellsThisShell);

        double cumMass = massThisShell;
        double cumNumCells = numCellsThisShell;
        
        while (massThisShell > 0)
        {
            shellRadius += radiusStep;
            massThisShell = GetEnclosedMass(clusterCentroids[thisClusterID], shellRadius, 
                                            clusterCoords[thisClusterID], clusterCellMasses[thisClusterID], numCellsThisShell) - cumMass;

            numCellsThisShell -= cumNumCells;

            shellMasses.push_back(massThisShell);
            numCellsShell.push_back(numCellsThisShell);

            cumMass += massThisShell;
            cumNumCells += numCellsThisShell;

            // cout << massThisShell << endl;
        }

        clusterShellMasses[thisClusterID] = shellMasses;
        clusterNumCellsShell[thisClusterID] = numCellsShell;

    } // end loop over clusters

    ofstream clusterMassDist(filePath+outFilename, ios::out);

    const int colWidth = 22;

    // column headers
    const char *vinit[] = {"clusterID", "centroidX", "centroidY", "centroidZ", "binLB", "binUB", "binMass", "totalMass", "meanDens", "occupiedVolume"};
    std::vector<std::string> colHeaders(vinit, end(vinit));

    clusterMassDist << std::left;
    for (auto thisColHeader : colHeaders)
        clusterMassDist << std::setw(colWidth) << thisColHeader.c_str();

    clusterMassDist << '\n';

    for (int i=2; i<numClusters+2; i++)
    {
        double displayCentroid[MESH_MDIM] = {clusterCentroids[i][IAXIS], clusterCentroids[i][JAXIS], clusterCentroids[i][KAXIS]};
        
        for (int dim=0; dim<MESH_MDIM; dim++)
        {
            if (displayCentroid[dim] < domainBoundBox[LOWER][dim])
                displayCentroid[dim] = domainBoundBox[UPPER][dim] - (domainBoundBox[LOWER][dim] - displayCentroid[dim]);
        }

        int numCellsThisCluster = clusterIndices[i].size();
        double totalClusterMass = 0;
        for (auto thisClusterCellMass : clusterCellMasses[i])
            totalClusterMass += thisClusterCellMass;
        double meanClusterDens = totalClusterMass / (numCellsThisCluster*cellVol);
        

        int numBinsThisCluster = clusterShellMasses[i].size();
        for (int bin=0; bin<numBinsThisCluster; bin++)
        {
            double occupiedVolume = clusterNumCellsShell[i][bin] * cellVol;

            double binLB = bin*radiusStep,
                   binUB = binLB + radiusStep,
                   binCenter = (binLB+binUB) / 2.;

            clusterMassDist << std::left <<
                            std::setw(colWidth) << i << 
                            std::setw(colWidth) << displayCentroid[IAXIS] <<
                            std::setw(colWidth) << displayCentroid[JAXIS] <<
                            std::setw(colWidth) << displayCentroid[KAXIS] <<
                            std::setw(colWidth) << binLB <<
                            std::setw(colWidth) << binUB <<
                            std::setw(colWidth) << clusterShellMasses[i][bin] <<
                            std::setw(colWidth) << totalClusterMass <<
                            std::setw(colWidth) << meanClusterDens <<
                            std::setw(colWidth) << occupiedVolume <<
                            endl;
        }
    }

    clusterMassDist.close();

    cout << "Wrote mass distribution analysis to " << filePath+outFilename << '\n';
    cout << endl;
}

void GetClusterIndices(FlashAmrMesh &mesh, std::vector<int> &uniqueClusterIDs, 
                       std::vector<std::deque<std::array<int,MESH_MDIM+1>>> &clusterIndices,
                       std::vector<clusterBoundaryParameters> &boundaryParameters)
{
    // only works on a uniform mesh!
    int firstCell[] = {0,0,0};
    double cellVol = mesh.getSingleCellVol(0, INTERIOR, firstCell);
    
    int clusterIDMeshIndex = mesh.findVarIndex("clst");
    int densMeshIndex = mesh.findVarIndex("dens");

    // we need to check each block to see if it's on a boundary
    // this is probably better accomplished using the gid array,
    // but there's still some functionality missing from FlashAMRMesh with regard
    // to the gid array, and I don't want to write it because I don't have time
    double domainBoundBox[2][MESH_MDIM];
    mesh.getDomainBoundBox(domainBoundBox);

    int numLocalBlocks = mesh.getLocalNumBlocks();
    int numLeafBlocks = 0;
    int *blockList = new int[numLocalBlocks];
    // double* coords[8];
    mesh.getListOfBlocks(LEAF, blockList, numLeafBlocks);


    for (int lb=0; lb<numLeafBlocks; lb++)
    {
        int currentBlock = blockList[lb];
        double ****solnData;
        //Get pointer to the block's data
        mesh.getBlkPtr(currentBlock, solnData, CENTER);


        // check if any of the block's bounds correspond to the domain's
        bool isBoundaryBlock = false;
        bool isBoundaryBlockDimension[MESH_MDIM] = {false, false, false};
        double blockBoundBox[2][MESH_MDIM];
        mesh.getBlkBoundBox(currentBlock, blockBoundBox);
        // for (int bound=LOWER; bound<UPPER+1; bound++)
        for (int dimension=0; dimension<MESH_MDIM; dimension++)
        {
            if (domainBoundBox[LOWER][dimension] == blockBoundBox[LOWER][dimension])
            {
                isBoundaryBlock = true;
                isBoundaryBlockDimension[dimension] = true;
            }
        }

        // std::vector<double[2][MESH_MDIM]> clusterBounds;

        int blkLimits[2][MESH_MDIM];
        int blkLimitsGC[2][MESH_MDIM];

        mesh.getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

        for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++)
        {
            for (int j=blkLimits[LOWER][JAXIS]; j<= blkLimits[UPPER][JAXIS]; j++)
            {
                for (int k=blkLimits[LOWER][KAXIS]; k<= blkLimits[UPPER][KAXIS]; k++)
                {
                    int thisClusterID = int(round(solnData[clusterIDMeshIndex][i][j][k]));
                    double cellMass = solnData[densMeshIndex][i][j][k] * cellVol;

                    if (thisClusterID > 1)
                    {
                        // only add to the list if not already there
                        if ( std::find(uniqueClusterIDs.begin(), uniqueClusterIDs.end(), thisClusterID) ==  uniqueClusterIDs.end())
                        {
                            uniqueClusterIDs.push_back(thisClusterID);
                            
                            // only resize if this ID makes the array bigger
                            int currentMaxIndex = clusterIndices.size()-1;
                            if (thisClusterID > currentMaxIndex)
                            {
                                // clusterMasses.resize(thisClusterID+1, 0.);
                                // numCells.resize(thisClusterID+1, 0.);
                                clusterIndices.resize(thisClusterID+1);
                                boundaryParameters.resize(thisClusterID+1);
                            }

                        }

                        // clusterMasses[thisClusterID] += cellMass;
                        // numCells[thisClusterID]++;
                        clusterIndices[thisClusterID].push_back({currentBlock,i,j,k});

                        if (isBoundaryBlock)
                        {
                            if (i == 0 && isBoundaryBlockDimension[IAXIS])
                            {
                                boundaryParameters[thisClusterID].isBoundaryCluster = true;
                                boundaryParameters[thisClusterID].crossesBoundary[IAXIS] = true;
                            }

                            if (j == 0 && isBoundaryBlockDimension[JAXIS])
                            {
                                boundaryParameters[thisClusterID].isBoundaryCluster = true;
                                boundaryParameters[thisClusterID].crossesBoundary[JAXIS] = true;
                            }

                            if (k == 0 && isBoundaryBlockDimension[KAXIS])
                            {
                                boundaryParameters[thisClusterID].isBoundaryCluster = true;
                                boundaryParameters[thisClusterID].crossesBoundary[KAXIS] = true;
                            }
                                
                        }

                    }
                }
            }
        }
    }    
}

void FindBoundaryClusterBounds(FlashAmrMesh &mesh, std::deque<std::array<int,MESH_MDIM+1>> &clusterIndices, clusterBoundaryParameters &boundaryParameters)
{
    double sideLengths[MESH_MDIM];
    mesh.getCellSideLengths(0, sideLengths);

    double domainBoundBox[2][MESH_MDIM];
    mesh.getDomainBoundBox(domainBoundBox);

    for (int dim=0; dim<MESH_MDIM; dim++)
    {
        boundaryParameters.clusterBounds[LOWER][dim] = domainBoundBox[UPPER][dim] + sideLengths[dim] / 2.;
        boundaryParameters.clusterBounds[UPPER][dim] = domainBoundBox[LOWER][dim] - sideLengths[dim] / 2.;
    }

    for (auto thisClusterCell : clusterIndices)
    {
        int thisCoordBlock = thisClusterCell[0],
        thisCoordI = thisClusterCell[1],
        thisCoordJ = thisClusterCell[2],
        thisCoordK = thisClusterCell[3];

        double thisCellCoords[MESH_MDIM];
        mesh.getCellCoords(thisCoordBlock, thisCoordI, thisCoordJ, thisCoordK, thisCellCoords);

        
        for (int dim=0; dim < MESH_MDIM; dim++)
        {
            int effectiveCellIndex = thisCellCoords[dim] / sideLengths[dim];
            if (boundaryParameters.crossesBoundary[dim])
            {
                
                if (thisCellCoords[dim] > boundaryParameters.clusterBounds[UPPER][dim])
                {
                    bool connectedAbove = false;
                    // now search the entire cluster to see if there are connecting cells one index higher
                    for (auto searchNeighborClusterCell : clusterIndices)
                    {
                        
                        double searchNeighborCoords[MESH_MDIM];
                        mesh.getCellCoords(searchNeighborClusterCell[0], searchNeighborClusterCell[1], 
                                           searchNeighborClusterCell[2], searchNeighborClusterCell[3], searchNeighborCoords);

                        int effectiveNeighborIndex = searchNeighborCoords[dim] / sideLengths[dim];

                        if (effectiveNeighborIndex == effectiveCellIndex+1)
                        {
                            connectedAbove = true;                            
                            break;
                        }

                    }
                    if (!connectedAbove && abs(thisCellCoords[dim] - domainBoundBox[UPPER][dim]) > sideLengths[dim] / 2.)
                    {
                        boundaryParameters.clusterBounds[UPPER][dim] = thisCellCoords[dim];
                    }
                }

                if (thisCellCoords[dim] < boundaryParameters.clusterBounds[LOWER][dim])
                {
                    bool connectedBelow = false;
                    // now search the entire cluster to see if there are connecting cells one index higher
                    for (auto searchNeighborClusterCell : clusterIndices)
                    {
                        
                        double searchNeighborCoords[MESH_MDIM];
                        mesh.getCellCoords(searchNeighborClusterCell[0], searchNeighborClusterCell[1], 
                                           searchNeighborClusterCell[2], searchNeighborClusterCell[3], searchNeighborCoords);

                        int effectiveNeighborIndex = searchNeighborCoords[dim] / sideLengths[dim];

                        if (effectiveNeighborIndex == effectiveCellIndex-1)
                        {
                            connectedBelow = true;                            
                            break;
                        }

                    }
                    if (!connectedBelow && abs(thisCellCoords[dim] - domainBoundBox[LOWER][dim]) > sideLengths[dim] / 2.)
                    {
                        boundaryParameters.clusterBounds[LOWER][dim] = thisCellCoords[dim];
                    }
                }
            }
        }
    }
}

void GetClusterCentroid(FlashAmrMesh &mesh,
                        std::deque<std::array<int,MESH_MDIM+1>> &clusterIndices, 
                        std::deque<std::array<double,MESH_MDIM>> &clusterCoords, 
                        std::deque<double> &clusterCellMasses,
                        std::array<double,MESH_MDIM> &clusterCentroids,
                        clusterBoundaryParameters &boundaryParameters)
{
    
        // find the cluster's centroid

        double domainBoundBox[2][MESH_MDIM];
        mesh.getDomainBoundBox(domainBoundBox);

        int densMeshIndex = mesh.findVarIndex("dens");
        // only works for uniform meshes 
        // double sideLengths[MESH_MDIM];
        // mesh.getCellSideLengths(0, sideLengths);

        // only works on a uniform mesh!
        int firstCell[] = {0,0,0};
        double cellVol = mesh.getSingleCellVol(0, INTERIOR, firstCell);


        double cumWeightedXCoord = 0,
               cumWeightedYCoord = 0, 
               cumWeightedZCoord = 0;

        double totalClusterMass = 0;

        // loop over cells in this cluster
        // int cellIndex = 0;
        for (auto thisClusterCell : clusterIndices)
        {
            int thisCoordBlock = thisClusterCell[0],
                thisCoordI = thisClusterCell[1],
                thisCoordJ = thisClusterCell[2],
                thisCoordK = thisClusterCell[3];

            // get the density of this cell
            int cellIndices[MESH_MDIM] = {thisCoordI,thisCoordJ,thisCoordK};
            double density = mesh.getPointData(thisCoordBlock, CENTER, densMeshIndex, INTERIOR, cellIndices);

            double cellCoords[MESH_MDIM];
            mesh.getCellCoords(thisCoordBlock, thisCoordI, thisCoordJ, thisCoordK, cellCoords);
            std::array<double,MESH_MDIM> thisCellCoords = {cellCoords[IAXIS], cellCoords[JAXIS], cellCoords[KAXIS]};


            // if this is a boundary cluster, shift the coordinates of the appropriate dimension so that the cluster appears contiguous
            // the recorded cell masses still represent the correct mesh cells, so the weighting and binning will still be correct
            // transform the coordinates to the side of the lower boundary
            if (boundaryParameters.isBoundaryCluster)
            {
                for (int dim=0; dim<MESH_MDIM; dim++)
                {
                    if (boundaryParameters.crossesBoundary[dim])
                    {
                        if (thisCellCoords[dim] >= boundaryParameters.clusterBounds[LOWER][dim])
                        {
                            thisCellCoords[dim] = thisCellCoords[dim] - (domainBoundBox[UPPER][dim] - domainBoundBox[LOWER][dim]);
                        }
                    }
                }
                
            }
            
            clusterCoords.push_back(thisCellCoords);

            double cellMass = density * cellVol;
            clusterCellMasses.push_back(cellMass);
            // cellMass = 1;
            cumWeightedXCoord += cellMass * thisCellCoords[IAXIS];
            cumWeightedYCoord += cellMass * thisCellCoords[JAXIS];
            cumWeightedZCoord += cellMass * thisCellCoords[KAXIS];

            totalClusterMass += cellMass;
            // cellIndex++;

        } // end loop over cells in this cluster

        // int numCellsThisCluster = clusterIndices.size();
        double normFactor = totalClusterMass;
        clusterCentroids = {cumWeightedXCoord/normFactor, cumWeightedYCoord/normFactor, cumWeightedZCoord/normFactor};

        // for (int dim=0; dim<MESH_MDIM; dim++)
}

double GetEnclosedMass(std::array<double,MESH_MDIM> origin, double radius, 
                       std::deque<std::array<double,MESH_MDIM>> &clusterCoords, std::deque<double> &clusterCellMasses, int &numCells)
{
    double radiusSquared = radius*radius;

    double cumEnclosedMass = 0.;

    // loop over all cluster cells
    int cellIndex = 0;
    numCells = 0;
    for (auto thisClusterCell : clusterCoords)
    {
        // calculate the distance from the origin point
        double diffX = thisClusterCell[IAXIS] - origin[IAXIS],
               diffY = thisClusterCell[JAXIS] - origin[JAXIS],
               diffZ = thisClusterCell[KAXIS] - origin[KAXIS];

        double distanceFromOriginSq = diffX*diffX + diffY*diffY + diffZ*diffZ;

        if (distanceFromOriginSq <= radiusSquared)
        {
            cumEnclosedMass += clusterCellMasses[cellIndex];
            numCells++;
        }
        

        cellIndex++;
    }

    return cumEnclosedMass;
}
