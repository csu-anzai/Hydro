for (int b = 0; b <= nBins; b++)
{
    for (int i=blkLimits[LOWER][IAXIS]; i<= blkLimits[UPPER][IAXIS]; i++)
    {
        for (int j=blkLimits[LOWER][JAXIS]; j<= blkLimits[UPPER][JAXIS]; j++)
        {
				
             for (int k=blkLimits[LOWER][KAXIS]; k<= blkLimits[UPPER][KAXIS]; k++)
             {		
                  if(BinMin[b] <= SweepDirMaxFace[k] && SweepDirMinFace[k] <= BinMax[b]) 
                  {
                      lowerBound = max(SweepDirMinFace[k], BinMin[b]);
                      upperBound = min(SweepDirMaxFace[k], BinMax[b]); 
                      dx3 =  upperBound - lowerBound;
                      dVol = dx1 * dx2 * dx3;
                      BinAvg[b] += solnData[varIndex][i][j][k] * dVol; 
                  }
             }
        }
    }
}                  
