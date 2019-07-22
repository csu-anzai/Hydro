// DF NOTE: this doesn't seem to be working in parallel. it takes an unreasonable amount of time to write the file
// not sure what's going on, but it must be in the underlying UniformRefinement code
// for now, running in serial


#include <cstdlib>
#include <iostream>
#include "FlashAmrMesh.h"
#include "Driver_main.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include "IO.h"

using namespace std;

void SmoothFace(double**, double***, int*, int , int, int, int );
void CalculateNeighborSet(int[][MESH_MDIM], int , int , int );
void WriteFaceData(double**, double***, int*, int , int, int );
void AdjustCorners(double***, int*, int, int );
void SmoothPlaneBorder(double ***, int *, int );

int LAVAFLOW_DRIVER_BASENM(int argc, char **argv) {

    ConfigData cfgData;
    std::string cfgFile = "config";
    ReadConfigFile(cfgData, cfgFile);

    std::string inFilePath = cfgData.inFilePath;
    std::string outFilePath = cfgData.outFilePath;
    std::string inFileName = cfgData.inFilename;
    std::string outFileName = cfgData.outFilename;

    /* Parameters */
    // std::string inFileName("/data1/df11c/data/data/wdm/wdm_hdf5_chk_0050");
    // std::string outFileName("/data1/df11c/Simulations/MergerEneryDecomp/wdm_hdf5_chk_0050");

    /* A uniform mesh will be output at this refinement level. 
     * If the refinement level is too high, the highest refinement
     * level will automatically be used.
     */
    int  refineLevel = -1;

    std::vector<std::string> meshVars = {"velx", "vely", "velz", "dens", "pres", "eint", "temp", "c12 "};
    //std::vector<std::string> meshVars = {"velx", "vely", "velz", "dens"};
    

    FlashAmrMesh mesh(inFilePath+inFileName, outFilePath+outFileName, meshVars);

    
    /* Only need to declare the 4d pointer here. 
     * Memory allocation is done within the refineToFinest function. 
     */
    double**** fineArray;
    
    /* Subdomain specification -- Must be within the original mesh
     * and must be large enough to contain at least one cell
     * per processor. 
     *
     * If no subdomain is required, set all values to 0.
     *
     * Need to specify all dimensions -- Ex: If you're using a 2D mesh,
     * Leave the KAXIS values 0 and set the IAXIS and JAXIS values.
     *
     * Due to discretization, the output coordinates will be slightly offset
     * from the original input coordinates.
     */
 
    double subdomainCoords[2][MESH_MDIM];

    // TODO: adjust coords to give an even number of cells per block


    subdomainCoords[LOWER][IAXIS] = cfgData.subdomainBoundsX[0];
    subdomainCoords[UPPER][IAXIS] = cfgData.subdomainBoundsX[1];
    subdomainCoords[LOWER][JAXIS] = cfgData.subdomainBoundsY[0];
    subdomainCoords[UPPER][JAXIS] = cfgData.subdomainBoundsY[1];
    subdomainCoords[LOWER][KAXIS] = cfgData.subdomainBoundsZ[0];
    subdomainCoords[UPPER][KAXIS] = cfgData.subdomainBoundsZ[1];

   
    mesh.refineToFinest(fineArray, subdomainCoords, refineLevel);

    // std::cout << "Initial value in fine array: " << fineArray[3][23][0][5] << std::endl;

#if 0   
 // smooth it around the edges so we can use periodic BCs
    
    // for a given cell, this gives how many neighbors should be used for the averaging
    const int smoothWindow = 3;
    
    // this specifies which cells should be smoothed
    // a smoothing border of 1 means just the border cells.
    const int smoothBorder = 5;

    // THIS SHOULD ONLY BE RUN IN SERIAL!!
    // I don't want to have to pass data between processors for the averaging
    // note that this is different than the issue mentioned at the top of this file concerning running in parallel
    int* totalCells = mesh.getTotalCellsUniform();

    double ***averagedData;

    averagedData = new double**[MESH_MDIM*2];
    for (int i=0; i<MESH_MDIM*2; i++)
    {
        averagedData[i] = new double*[totalCells[JAXIS]];
        for (int j=0;j<totalCells[JAXIS]; j++)
        {
            averagedData[i][j] = new double[totalCells[KAXIS]];
        }
    }
   
    // loop over mesh variables
    for (int meshVar = 0; meshVar < meshVars.size(); meshVar++)
    {

        // zero out the intersections
        // this only works on a cubic domain right now, because I need to graduate and I don't really feel like doing it right
        // in fact, that statement applies to most of this code

        
    
        // loop over the planes defining the border
        for (int planeOffset=0; planeOffset<smoothBorder; planeOffset++)
        {
            // loop over the faces in each dimension
            int faceCount = 0;
            for (int dimension=0; dimension<MESH_MDIM; dimension++)
            {
                // average over the lower face
                SmoothFace(averagedData[faceCount], fineArray[meshVar], totalCells, planeOffset, dimension, smoothWindow, planeOffset);
                SmoothFace(averagedData[faceCount+1], fineArray[meshVar], totalCells, totalCells[dimension]-1-planeOffset, dimension, smoothWindow, planeOffset);

                faceCount+=2; 
            }
        // }

            
            //AdjustCorners(fineArray[meshVar], totalCells, 0, planeOffset);
            
            faceCount = 0;
        // for (int planeOffset=0; planeOffset<smoothBorder; planeOffset++)
        // {
            // write the values from the averagedData array into the original mesh
            for (int dimension=0; dimension<MESH_MDIM; dimension++)
            {
              //  WriteFaceData(averagedData[faceCount], fineArray[meshVar], totalCells, planeOffset, dimension, planeOffset);
              //  WriteFaceData(averagedData[faceCount+1], fineArray[meshVar], totalCells, totalCells[dimension]-1-planeOffset, dimension, planeOffset);

                faceCount+=2; 
            }

            //AdjustCorners(fineArray[meshVar], totalCells, -1, planeOffset);
        }

        
    }
#endif
    
    mesh.writeOutFileUniform(fineArray);

    mesh.deallocateFinest(fineArray);

    
#if 0
    for (int i=0; i<MESH_MDIM*2; i++)
    {
        for (int j=0;j<totalCells[JAXIS]; j++)
        {
            delete[] averagedData[i][j];
        }
        delete[] averagedData[i];
    }

    delete[] averagedData;
#endif

    return 0;
}



void SmoothFace(double **faceData, double ***origData, int *nCells, int faceIndex, int normalDir, int smoothWindow, int planeOffset)
{
    
    int neighborIndices[2][MESH_MDIM];
    CalculateNeighborSet(neighborIndices, faceIndex, normalDir, nCells[normalDir]);
    
    if (normalDir == IAXIS)
    {       
        for (int i=planeOffset; i<nCells[JAXIS]-planeOffset; i++)
        {
            CalculateNeighborSet(neighborIndices, i, JAXIS, nCells[JAXIS]);

            for (int j=planeOffset; j<nCells[KAXIS]-planeOffset; j++)
            {
                CalculateNeighborSet(neighborIndices, j, KAXIS, nCells[KAXIS]);

                double cumNeighborSum = 0;

                for (int windowPosition=0; windowPosition<smoothWindow; windowPosition++)
                {
                    int lowerBound = neighborIndices[LOWER][IAXIS] - windowPosition;
                    int upperBound = neighborIndices[UPPER][IAXIS] + windowPosition;
                    
                    // lowerBound = lowerBound ;
                    if (lowerBound < 0)
                        lowerBound = nCells[IAXIS] + lowerBound;

                    // upperBound = upperBound ;
                    if (upperBound > nCells[IAXIS]-1)
                        upperBound = upperBound - nCells[IAXIS];
                    cumNeighborSum += origData[lowerBound][i][j] + origData[upperBound][i][j];
                }

                faceData[i][j] = (cumNeighborSum + origData[faceIndex][i][j]) / (1+1*smoothWindow*2.);

                // if (i == 23 && j == 5 && faceIndex == 23)
                // {
                //     std::cout << "x plane corner data: " << faceData[i][j] << std::endl;
                //     std::cout << "lower neighbor: " << neighborIndices[LOWER][IAXIS] << std::endl;
                //     std::cout << "upper neighbor: " << neighborIndices[UPPER][IAXIS] << std::endl;
                // }

                // faceData[i][j] = 10;

            }
        }
    }

    if (normalDir == JAXIS)
    {       
        for (int i=planeOffset; i<nCells[IAXIS]-planeOffset; i++)
        {
            CalculateNeighborSet(neighborIndices, i, IAXIS, nCells[IAXIS]);

            for (int j=planeOffset; j<nCells[KAXIS]-planeOffset; j++)
            {
                CalculateNeighborSet(neighborIndices, j, KAXIS, nCells[KAXIS]);

                double cumNeighborSum = 0;

                
                // int lowerBound = neighborIndices[LOWER][JAXIS];
                // int upperBound = neighborIndices[UPPER][JAXIS];
                for (int windowPosition=0; windowPosition<smoothWindow; windowPosition++)
                {
                    int lowerBound = neighborIndices[LOWER][JAXIS] - windowPosition;
                    int upperBound = neighborIndices[UPPER][JAXIS] + windowPosition;
                    
                    // lowerBound = lowerBound - windowPosition;
                    if (lowerBound < 0)
                        lowerBound = nCells[JAXIS] + lowerBound;

                    // upperBound = upperBound + windowPosition;
                    if (upperBound > nCells[JAXIS]-1)
                        upperBound = upperBound - nCells[JAXIS];
                    cumNeighborSum += origData[i][lowerBound][j] + origData[i][upperBound][j];

                    // if (i == 23 && j == 5 && faceIndex == 23)
                    // {
                        // std::cout << "bounds (lower, upper): " << lowerBound << ", " << upperBound << std::endl;
                    // }
                }
                

                faceData[i][j] = (cumNeighborSum + origData[i][faceIndex][j]) / (1+1*smoothWindow*2.);

                // if (i == 23 && j == 5 && faceIndex == 23)
                // {
                //     std::cout << "y plane corner data: " << faceData[i][j] << std::endl;
                //     std::cout << "lower neighbor: " << neighborIndices[LOWER][JAXIS] << std::endl;
                //     std::cout << "upper neighbor: " << neighborIndices[UPPER][JAXIS] << std::endl;
                //     std::cout << "neighbor values: " << origData[23][2][5] << std::endl;
                //     std::cout << "neighbor values: " << origData[23][1][5] << std::endl;
                //     std::cout << "neighbor values: " << origData[23][0][5] << std::endl;
                //     std::cout << "neighbor values: " << origData[23][22][5] << std::endl;
                //     std::cout << "neighbor values: " << origData[23][21][5] << std::endl;
                //     std::cout << "neighbor values: " << origData[23][20][5] << std::endl;
                // }

                // faceData[i][j] = 7;

            }
        }
    }

    if (normalDir == KAXIS)
    {       
        for (int i=planeOffset; i<nCells[IAXIS]-planeOffset; i++)
        {
            CalculateNeighborSet(neighborIndices, i, IAXIS, nCells[IAXIS]);

            for (int j=planeOffset; j<nCells[JAXIS]-planeOffset; j++)
            {
                CalculateNeighborSet(neighborIndices, j, JAXIS, nCells[JAXIS]);

                double cumNeighborSum = 0;

                
                for (int windowPosition=0; windowPosition<smoothWindow; windowPosition++)
                {
                    int lowerBound = neighborIndices[LOWER][KAXIS] - windowPosition;
                    int upperBound = neighborIndices[UPPER][KAXIS] + windowPosition;
                    
                    // lowerBound = lowerBound;
                    if (lowerBound < 0)
                        lowerBound = nCells[KAXIS] + lowerBound;

                    // upperBound = upperBound + windowPosition;
                    if (upperBound > nCells[KAXIS]-1)
                        upperBound = upperBound - nCells[KAXIS];
                    cumNeighborSum += origData[i][j][lowerBound] + origData[i][j][upperBound];
                }
                

                faceData[i][j] = (cumNeighborSum + origData[i][j][faceIndex]) / (1+1*smoothWindow*2.);

                // faceData[i][j] = 5;

            }
        }
    }


}

void CalculateNeighborSet(int neighborIndices[][MESH_MDIM], int scanIndex, int neighborIndex, int nCells)
{
    if (scanIndex == 0)
        neighborIndices[LOWER][neighborIndex] = nCells - 1;
    else
        neighborIndices[LOWER][neighborIndex] = scanIndex-1;

    if (scanIndex == nCells - 1)
        neighborIndices[UPPER][neighborIndex] = 0;
    else
        neighborIndices[UPPER][neighborIndex] = scanIndex+1;   
}

void WriteFaceData(double **faceData, double ***origData, int *nCells, int faceIndex, int normalDir, int planeOffset)
{
    if (normalDir == IAXIS)
    {
        for (int i=planeOffset; i<nCells[JAXIS]-planeOffset; i++)
            for (int j=planeOffset; j<nCells[KAXIS]-planeOffset; j++)
            {
                if (i == planeOffset || j == planeOffset || i == nCells[JAXIS]-1-planeOffset || j == nCells[KAXIS]-1-planeOffset)
                {
                    // origData[faceIndex][i][j] += 7.;
                    origData[faceIndex][i][j] += faceData[i][j];
                }
                else
                    origData[faceIndex][i][j] = faceData[i][j];
            }

    }
    else if (normalDir == JAXIS)
    {
        
        for (int i=planeOffset; i<nCells[IAXIS]-planeOffset; i++)
            for (int j=planeOffset; j<nCells[KAXIS]-planeOffset; j++)
            {
                if (i == planeOffset || j == planeOffset || i == nCells[IAXIS]-1-planeOffset || j == nCells[KAXIS]-1-planeOffset)
                {
                    // origData[i][faceIndex][j] += 7.;
                    origData[i][faceIndex][j] += faceData[i][j];
                }
                else
                    origData[i][faceIndex][j] = faceData[i][j];
            }

    }
    else if (normalDir == KAXIS)
    {
        for (int i=planeOffset; i<nCells[IAXIS]-planeOffset; i++)
            for (int j=planeOffset; j<nCells[JAXIS]-planeOffset; j++)
            {
                if (i == planeOffset || j == planeOffset || i == nCells[JAXIS]-1-planeOffset || j == nCells[KAXIS]-1-planeOffset)
                {
                    // origData[i][j][faceIndex] += 7.;
                    origData[i][j][faceIndex] += faceData[i][j];
                }
                else
                    origData[i][j][faceIndex] = faceData[i][j];
            }
    }
}

void AdjustCorners(double ***data, int *nCells, int value, int planeOffset)
{
    //  value is a flag
    // if it's negative, average the existing values
    // otherwise, zero it out
    if (value >= 0)
    {
        // all x-running edges
        for (int i=planeOffset; i<nCells[IAXIS]-planeOffset; i++)
        {
            data[i][planeOffset][planeOffset] = 0.;
            data[i][planeOffset][nCells[KAXIS]-1-planeOffset] = 0.;
            data[i][nCells[JAXIS]-1-planeOffset][planeOffset] = 0.;
            data[i][nCells[JAXIS]-1-planeOffset][nCells[KAXIS]-1-planeOffset] = 0.;
        }

        // all y-running edges
        for (int j=planeOffset; j<nCells[JAXIS]-planeOffset; j++)
        {
            data[planeOffset][j][planeOffset] = 0.;
            data[planeOffset][j][nCells[KAXIS]-1-planeOffset] = 0.;
            data[nCells[IAXIS]-1-planeOffset][j][planeOffset] = 0.;
            data[nCells[IAXIS]-1-planeOffset][j][nCells[KAXIS]-1-planeOffset] = 0.;
        }

        // all z-running edges
        for (int k=planeOffset; k<nCells[KAXIS]-planeOffset; k++)
        {
            data[planeOffset][planeOffset][k] = 0.;
            data[planeOffset][nCells[JAXIS]-1-planeOffset][k] = 0.;
            data[nCells[IAXIS]-1-planeOffset][planeOffset][k] = 0.;
            data[nCells[IAXIS]-1-planeOffset][nCells[JAXIS]-1-planeOffset][k] = 0.;
        }        
    }
    else
    {
        // all x-running edges
        for (int i=planeOffset+1; i<nCells[IAXIS]-1-planeOffset; i++)
        {
            data[i][planeOffset][planeOffset] = data[i][planeOffset][planeOffset] / 2. ;
            data[i][planeOffset][nCells[KAXIS]-1-planeOffset] = data[i][planeOffset][nCells[KAXIS]-1-planeOffset] / 2. ;
            data[i][nCells[JAXIS]-1-planeOffset][planeOffset] = data[i][nCells[JAXIS]-1-planeOffset][planeOffset] / 2. ;
            data[i][nCells[JAXIS]-1-planeOffset][nCells[KAXIS]-1-planeOffset] = data[i][nCells[JAXIS]-1-planeOffset][nCells[KAXIS]-1-planeOffset] / 2. ;
        }

        // all y-running edges
        for (int j=planeOffset+1; j<nCells[JAXIS]-1-planeOffset; j++)
        {
            data[planeOffset][j][planeOffset] = data[planeOffset][j][planeOffset] / 2. ;
            data[planeOffset][j][nCells[KAXIS]-1-planeOffset] = data[planeOffset][j][nCells[KAXIS]-1-planeOffset] / 2. ;
            data[nCells[IAXIS]-1-planeOffset][j][planeOffset] = data[nCells[IAXIS]-1-planeOffset][j][planeOffset] / 2. ;
            data[nCells[IAXIS]-1-planeOffset][j][nCells[KAXIS]-1-planeOffset] = data[nCells[IAXIS]-1-planeOffset][j][nCells[KAXIS]-1-planeOffset] / 2. ;
        }

        // all z-running edges
        for (int k=planeOffset+1; k<nCells[KAXIS]-1-planeOffset; k++)
        {
            data[planeOffset][planeOffset][k] = data[planeOffset][planeOffset][k] / 2. ;
            data[planeOffset][nCells[JAXIS]-1-planeOffset][k] = data[planeOffset][nCells[JAXIS]-1-planeOffset][k] / 2. ;
            data[nCells[IAXIS]-1-planeOffset][planeOffset][k] = data[nCells[IAXIS]-1-planeOffset][planeOffset][k] / 2. ;
            data[nCells[IAXIS]-1-planeOffset][nCells[JAXIS]-1-planeOffset][k] = data[nCells[IAXIS]-1-planeOffset][nCells[JAXIS]-1-planeOffset][k] / 2. ;
        }

        // do the corners
        for (int i=planeOffset; i<nCells[IAXIS]-planeOffset; i+=nCells[IAXIS]-1-(planeOffset*2))
            for (int j=planeOffset; j<nCells[JAXIS]-planeOffset; j+=nCells[JAXIS]-1-(planeOffset*2))
                for (int k=planeOffset; k<nCells[KAXIS]-planeOffset; k+=nCells[KAXIS]-1-(planeOffset*2))
                    data[i][j][k] = data[i][j][k] / 3.;


        SmoothPlaneBorder(data, nCells, planeOffset);

    }
}

void SmoothPlaneBorder(double ***data, int *nCells, int planeOffset)
{
    
    for (int planeNormal=0; planeNormal<MESH_MDIM; planeNormal++)
    {
        if (planeNormal == IAXIS)
        {
            for (int k=planeOffset+1; k<nCells[KAXIS]-planeOffset-1; k+=nCells[KAXIS]-3-(planeOffset*2))
                for (int j=planeOffset+1; j<nCells[JAXIS]-1-planeOffset; j++)
                {
                    data[planeOffset][j][k] = (data[planeOffset][j][k+1] + data[planeOffset][j][k-1]) / 2.;
                    data[nCells[IAXIS]-1-planeOffset][j][k] = (data[nCells[IAXIS]-1-planeOffset][j][k+1] + data[nCells[IAXIS]-1-planeOffset][j][k-1]) / 2.;
                }

            for (int j=planeOffset+1; j<nCells[KAXIS]-planeOffset-1; j+=nCells[JAXIS]-3-(planeOffset*2))
                for (int k=planeOffset+1; k<nCells[KAXIS]-1-planeOffset; k++)
                {
                    data[planeOffset][j][k] = (data[planeOffset][j+1][k] + data[planeOffset][j-1][k]) / 2.;
                    data[nCells[IAXIS]-1-planeOffset][j][k] = (data[nCells[IAXIS]-1-planeOffset][j+1][k] + data[nCells[IAXIS]-1-planeOffset][j-1][k]) / 2.;
                }
        }

        if (planeNormal == JAXIS)
        {
            for (int k=planeOffset+1; k<nCells[KAXIS]-planeOffset-1; k+=nCells[KAXIS]-3-(planeOffset*2))
                for (int i=planeOffset+1; i<nCells[IAXIS]-1-planeOffset; i++)
                {
                    data[i][planeOffset][k] = (data[i][planeOffset][k+1] + data[i][planeOffset][k-1]) / 2.;
                    data[i][nCells[JAXIS]-1-planeOffset][k] = (data[i][nCells[JAXIS]-1-planeOffset][k+1] + data[i][nCells[JAXIS]-1-planeOffset][k-1]) / 2.;
                }

            for (int i=planeOffset+1; i<nCells[IAXIS]-planeOffset-1; i+=nCells[KAXIS]-3-(planeOffset*2))
                for (int k=planeOffset+1; k<nCells[KAXIS]-1-planeOffset; k++)
                {
                    data[i][planeOffset][k] = (data[i+1][planeOffset][k] + data[i-1][planeOffset][k]) / 2.;
                    data[i][nCells[JAXIS]-1-planeOffset][k] = (data[i+1][nCells[JAXIS]-1-planeOffset][k] + data[i-1][nCells[JAXIS]-1-planeOffset][k]) / 2.;
                }
        }

        if (planeNormal == KAXIS)
        {
            for (int j=planeOffset+1; j<nCells[JAXIS]-planeOffset-1; j+=nCells[JAXIS]-3-(planeOffset*2))
                for (int i=planeOffset+1; i<nCells[IAXIS]-1-planeOffset; i++)
                {
                    data[i][j][planeOffset] = (data[i][j+1][planeOffset] + data[i][j-1][planeOffset]) / 2.;
                    data[i][j][nCells[KAXIS]-1-planeOffset] = (data[i][j+1][nCells[KAXIS]-1-planeOffset] + data[i][j-1][nCells[KAXIS]-1-planeOffset]) / 2.;
                }

            for (int i=planeOffset+1; i<nCells[IAXIS]-planeOffset-1; i+=nCells[IAXIS]-3-(planeOffset*2))
                for (int j=planeOffset+1; j<nCells[JAXIS]-1-planeOffset; j++)
                {
                    data[i][j][planeOffset] = (data[i+1][j][planeOffset] + data[i-1][j][planeOffset]) / 2.;
                    data[i][j][nCells[KAXIS]-1-planeOffset] = (data[i+1][j][nCells[KAXIS]-1-planeOffset] + data[i-1][j][nCells[KAXIS]-1-planeOffset]) / 2.;
                }
        }
    }
}
