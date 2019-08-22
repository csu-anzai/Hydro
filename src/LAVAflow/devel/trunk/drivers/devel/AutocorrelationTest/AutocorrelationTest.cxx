#include "EulTimeAutoCorr.h"
#include "LavaMPI.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include "Driver_main.h"
using namespace std;

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{

	// LavaMPI mpiObj(argc, argv);

	// string baseFilename("/data1/data/tburn_hdf5_plt_cnt_");
	// int fileNumBounds[2] = {1, 3};
	// int numPoints = 1000;

	// int aflag = 0;
	// int bflag = 0;
	// int c;

	// /* options are: b: base filename
	// 				i: initial file number
	// 				f: final file number
	// 				n: number of sample points */
	// while ((c = getopt (argc, argv, "b:i:f:n:")) != -1)
	//   switch (c)
	//     {
	//     case 'b':
	//       baseFilename = optarg;
	//       break;
	//     case 'i':
	//       fileNumBounds[0] = atoi(optarg);
	//       break;
	//     case 'f':
	//       fileNumBounds[1] = atoi(optarg);
	//       break;
	//     case 'n':
	//       numPoints = atoi(optarg);
	//       break;
	//     case '?':
	//       if (isprint (optopt))
	//         fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	//       else
	//         fprintf (stderr,
	//                  "Unknown option character `\\x%x'.\n",
	//                  optopt);
	//       return 1;
	//     default:
	//       abort ();
	//     }
	// // cout << "base filename: " << baseFilename << ", file bounds: " << fileNumBounds[0] <<
	// // " " << fileNumBounds[1] << ", number of points: " << numPoints << endl;

 //  	for (int index = optind; index < argc; index++)
 //    	printf ("Non-option argument %s\n", argv[index]);

	// EulTimeAutoCorr autoCorr(baseFilename, fileNumBounds, numPoints);
	// autoCorr.GetAutoCorrelations();

	// int numTimeSeps = autoCorr.numTimeSeps;

	// int packedVectorLength = numTimeSeps*MESH_MDIM;

	// vector<double> globalAvgTimeSeps(packedVectorLength),
	// globalAvgAutoCorr(packedVectorLength), localAvgTimeSeps(packedVectorLength),
	// localAvgAutoCorr(packedVectorLength);

	// // pack the vectors
	// int packedIndex=0;
	// for (int i=0; i<autoCorr.timeSeps.size(); i++)
	// 	for (int j=0; j<MESH_MDIM; j++)
	// 	{
	// 		localAvgTimeSeps[packedIndex] = autoCorr.timeSeps[i][j];
	// 		localAvgAutoCorr[packedIndex] = autoCorr.autoCorrelations[i][j];
	// 		packedIndex++;
	// 	}

	// MPI::COMM_WORLD.Reduce(&localAvgTimeSeps.front(), &globalAvgTimeSeps.front(),
	// 						localAvgTimeSeps.size(), MPI::DOUBLE, MPI_SUM, 0);

	// MPI::COMM_WORLD.Reduce(&localAvgAutoCorr.front(), &globalAvgAutoCorr.front(),
	// 						localAvgAutoCorr.size(), MPI::DOUBLE, MPI_SUM, 0);

	// if (MPI::COMM_WORLD.Get_rank() ==0)
	// {
	// 	int commSize = MPI::COMM_WORLD.Get_size();

	// 	vector<array<double, MESH_MDIM>> timeSeps(numTimeSeps);
	// 	vector<array<double, MESH_MDIM>> autoCorrelations(numTimeSeps);

	// 	// unpack and average
	// 	int packedIndex=0;
	// 	for (int i=0; i<numTimeSeps; i++)
	// 		for (int j=0; j<MESH_MDIM; j++)
	// 		{
	// 			timeSeps[i][j] = globalAvgTimeSeps[packedIndex] / commSize;
	// 			autoCorrelations[i][j] = globalAvgAutoCorr[packedIndex] / commSize;
	// 			packedIndex++;
	// 		}

	// 	string outFilename(baseFilename+"results");
	// 	ofstream outFile;
	// 	outFile.open (outFilename);

	// 	cout << "delta t: " ;
	// 	outFile << "delta t: " ;
	// 	for (auto i : timeSeps)
	// 	{
	//     	std::cout << i[IAXIS] << ' ';
	//     	outFile << i[IAXIS] << ' ';
	// 	}
	//     cout << endl;
	//     outFile << endl;

	//     cout << "autocorrelation: " ;
	//     outFile << "autocorrelation: " ;
	// 	for (int dim=0; dim<MESH_MDIM; dim++)
	// 	{
	// 		for (auto i : autoCorrelations)
	// 		{
	// 	    	std::cout << i[dim] << ' ';
	// 	    	outFile << i[dim] << ' ';
	// 		}

	//     	cout << endl;
	//     	outFile << endl;
	//     }

	//     outFile.close();
	// }
}
