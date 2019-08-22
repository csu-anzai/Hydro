#ifndef LAVA_TURBULENCE_TIMEAUTOCORR_H
#define LAVA_TURBULENCE_TIMEAUTOCORR_H

#include <string>
#include <vector>
#include <array>
#include "../Mesh/FlashAmrMesh.h"

typedef std::tuple<int,int,int> i3tuple;
typedef std::tuple<double,double,double> d3tuple;

class TimeAutoCorr
{
public:
	int numTimeSeps;
	std::vector<std::array<double, MESH_MDIM>> timeSeps;
	std::vector<std::array<double, MESH_MDIM>> autoCorrelations;

	std::string baseFilename;
	int fileIndexBounds[2];
	std::string *filenames;
	
	int numPoints;
	std::vector<std::array<double, MESH_MDIM>> points;
	int *pointBlockIDs;
	std::vector<std::array<int, MESH_MDIM>> pointZoneIDs;

	MPI::Intracomm meshComm;

	// TimeAutoCorr();
	TimeAutoCorr(std::string, int[2], int, MPI::Intracomm = MPI::COMM_WORLD);
	~TimeAutoCorr();

	void GeneratePoints(FlashAmrMesh&);
	void GetAutoCorrelations();
	void GetQuantityAtPoints(int , FlashAmrMesh &, double[] );
};


#endif
