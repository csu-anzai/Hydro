#include "TimeAutoCorr.h"
#include <vector>


void CollectAndOutput(TimeAutoCorr* autoCorr, std::string outFilename)
{
    int numVars = autoCorr->numVars;
    int numTimeSeps = autoCorr->numTimeSeps;
    int packedVectorLength = numTimeSeps*numVars;

    std::vector<double> globalAvgTimeSeps(packedVectorLength),
    globalAvgAutoCorr(packedVectorLength), localAvgTimeSeps(packedVectorLength),
    localAvgAutoCorr(packedVectorLength);

    // pack the vectors
    int packedIndex=0;
    for (int i=0; i<numTimeSeps; i++)
    {
        localAvgTimeSeps[packedIndex] = autoCorr->timeSeps[i];
        for (int j=0; j<numVars; j++)
        {            
            localAvgAutoCorr[packedIndex] = autoCorr->autoCorrelations[j][i];
            packedIndex++;
        }
    }

    MPI::COMM_WORLD.Reduce(&localAvgTimeSeps.front(), &globalAvgTimeSeps.front(),
                            localAvgTimeSeps.size(), MPI::DOUBLE, MPI_SUM, 0);

    MPI::COMM_WORLD.Reduce(&localAvgAutoCorr.front(), &globalAvgAutoCorr.front(),
                            localAvgAutoCorr.size(), MPI::DOUBLE, MPI_SUM, 0);


    //---------------------Output data on screen and/or file----------------|

    if (MPI::COMM_WORLD.Get_rank() ==0)
    {
        int commSize = MPI::COMM_WORLD.Get_size();

        std::vector<double> timeSeps(numTimeSeps);
        std::vector<std::vector<double>> autoCorrelations(numTimeSeps);

        for (auto &acThisTime:autoCorrelations)
          acThisTime.resize(numVars);

        // unpack and average
        int packedIndex=0;
        for (int i=0; i<numTimeSeps; i++)
        {
            timeSeps[i] = globalAvgTimeSeps[packedIndex] / commSize;
            for (int j=0; j<numVars; j++)
            {                
                autoCorrelations[i][j] = globalAvgAutoCorr[packedIndex] / commSize;
                packedIndex++;
            }
        }

        // std::string outFilename(baseFilename+"results");
        std::ofstream outFile;
        outFile.open (outFilename);

        std::cout << "delta t: " ;
        outFile << "delta t: " ;
        for (auto i : timeSeps)
        {
            std::cout << i << ' ';
            outFile << i << ' ';
        }
        std::cout << std::endl;
        outFile << std::endl;

        std::cout << "autocorrelation: " ;
        outFile << "autocorrelation: " ;
        for (int dim=0; dim<numVars; dim++)
        {
            for (auto i : autoCorrelations)
            {
                std::cout << i[dim] << ' ';
                outFile << i[dim] << ' ';
            }

            std::cout << std::endl;
            outFile << std::endl;
        }

        outFile.close();
    }

    MPI::COMM_WORLD.Barrier();
}