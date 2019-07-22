
#include "ShellAverage.h"
#include <vector>
#include "FlashAmrMesh.h"
#include <iostream>




ShellAverage::ShellAverage(void)
{
    //classData averageData

}


ShellAverage::~ShellAverage(void)
{
    

}

void ShellAverage::calcShellAverage(FlashAmrMesh* mesh,int sweepDir, int nBins, int varIndex )

{
 
    int secondDir, thirdDir, nDim=2, sweepInd;
    double sweepMin, sweepMax, dbin;
    double domainBounds [2][MESH_MDIM]; 
    averageData.binFaces.resize(nBins+1); 
    averageData.binAvg.resize(nBins);
    averageData.binVol.resize(nBins);  
    mesh->getDomainBoundBox(domainBounds);
    sweepMin = domainBounds[LOWER][sweepDir];
    sweepMax = domainBounds[UPPER][sweepDir];

    switch(sweepDir)
    {
        case IAXIS:
        {
            secondDir = JAXIS;
            thirdDir = KAXIS;
            break;
        }
        case JAXIS:
        {
            secondDir = IAXIS;
            thirdDir = KAXIS;
            break;
        }
        case KAXIS:
        {
            secondDir = IAXIS;
            thirdDir = JAXIS;
            break;
        }

        default:
        {
            cout << "calcShellAverage : invalid sweepDir, must be IAXIS, JAXIS, or KAXIS" << endl;
		    exit(EXIT_FAILURE);
        }
    }

    dbin = (sweepMax - sweepMin)/double(nBins);
    
    for(int i=0; i<=nBins; i++)
    {
       averageData.binFaces[i] = sweepMin + i*dbin; 
    }
    vector <double> binMin(this->averageData.binFaces.begin(), this->averageData.binFaces.end()-1);
    vector <double> binMax(this->averageData.binFaces.begin()+1, this->averageData.binFaces.end());

   	double cumDens = 0;
	int numLocalBlocks = mesh->getNBlocks();
	int numLeafBlocks = 0; 
	int *blockList = new int[numLocalBlocks];
    double * sweepDirMinFace;
    double * sweepDirMaxFace;
    double * secondDirCoords;
    double * thirdDirCoords;
    sweepDirMinFace = new double [16]; 
    sweepDirMaxFace = new double [16]; 
    secondDirCoords = new double [16]; 
    thirdDirCoords  = new double [16]; 
    nDim = 3;
	mesh->getListOfBlocks(LEAF, blockList, numLeafBlocks);

	int cellCount = 0;
	double totalMass = 0;
    double lowerBound, upperBound, dx1, dx2, dx3=1.0, dVol;
	for (int lb=0; lb<numLeafBlocks; lb++)
	{
		int currentBlock = blockList[lb];
		double ****solnData;
		//Get pointer to the block's data
		mesh->getBlkPtr(currentBlock, solnData, CENTER);
		double cellVol = mesh->getSingleCellVol(currentBlock);

		int blkLimits[2][MESH_MDIM];
		int blkLimitsGC[2][MESH_MDIM];
    
		mesh->getBlkIndexLimits(currentBlock, blkLimits, blkLimitsGC);

        mesh->getCellCoords(sweepDir, currentBlock, LEFT_EDGE, 0, sweepDirMinFace);
        mesh->getCellCoords(sweepDir, currentBlock, RIGHT_EDGE, 0, sweepDirMaxFace);
        mesh->getCellCoords(secondDir, currentBlock, CENTER, 0, secondDirCoords);
        mesh->getCellCoords(thirdDir, currentBlock, LEFT_EDGE, 0, thirdDirCoords);
        dx2 = secondDirCoords[1]-secondDirCoords[0]; 
        if (nDim > 2) 
        {
            
            dx3 = thirdDirCoords[1] - thirdDirCoords[0]; 

        }
        else dx3 = 1;
            for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++)
            {
                for (int j=blkLimits[LOWER][JAXIS]; j<= blkLimits[UPPER][JAXIS]; j++)
                {
				
                     for (int k=blkLimits[LOWER][KAXIS]; k<= blkLimits[UPPER][KAXIS]; k++)
                     {	
                        
	                        switch(sweepDir)
                            {
                                case IAXIS:
                                {
                                    sweepInd = i;
                                    break;
                                }
                                case JAXIS:
                                {
                                    sweepInd = j;
                                    break;
                                }
                                case KAXIS:
                                {
                                    sweepInd = k;
                                    break;
                                }

                                default:
                                {
                                    cout << "calcShellAverage : invalid sweepDir, must be IAXIS, JAXIS, or KAXIS" << endl;
		                            exit(EXIT_FAILURE);
                                }
                            }




                            for (int b = 0; b <= nBins-1; b++)
                            {                            
                                if(binMin[b] <= sweepDirMaxFace[sweepInd] && sweepDirMinFace[sweepInd] < binMax[b]) 
                                {
                                   // cout << binMin[b] << " " << (sweepDirMaxFace[sweepInd] + sweepDirMinFace[sweepInd]) /2.0 <<" "<< binMax[b] << " " <<secondDirCoords[j]  << endl;
                                    lowerBound = max(sweepDirMinFace[sweepInd], binMin[b]);
                                    upperBound = min(sweepDirMaxFace[sweepInd], binMax[b]); 
                                    dx1 =  upperBound - lowerBound;
                                    dVol = dx1 * dx2 * dx3;

                                    averageData.binAvg[b] += solnData[varIndex][i][j][k] * dVol; 
                                    averageData.binVol[b] += dVol;
                                   // cout << averageData.binVol[b] << " " << dx1 << " " << b  << endl;

                                }
 
                            }

                     }
                }
            }
                         
 
    




    }
	
    int binAvetest;
    for (int i=0; i <= nBins-1; i++)
    {
                binAvetest = (binMax[i]-binMin[i]) * 2 * 2;
                cout << binMax[i] << " " << binMin[i] << endl;
                averageData.binAvg[i] = averageData.binAvg[i] / averageData.binVol[i]; 
             	cout << "Average temperature in bin is " << averageData.binAvg[i] << " " << averageData.binVol[i] << " " << binAvetest << endl;
    }


	delete blockList;

}


        

    








