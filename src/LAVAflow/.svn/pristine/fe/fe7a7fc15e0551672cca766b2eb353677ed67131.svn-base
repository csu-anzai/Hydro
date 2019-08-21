#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include "FlashAmrMesh.h"
#include "ArrayOperations.h"
#include "FlashParticles.h"
#include "IndexOperations.h"
#include "IO.h"


void TrackClustersMassTracers(std::vector<std::string> &, std::string, ConfigData&);
void TrackClustersParticlesMT(std::vector<std::string> &, std::vector<std::string> &, std::string, ConfigData& );
void TrackClustersParticles(std::string, ConfigData& );
std::vector<int> BuildBlockLinearIndexMap(FlashAmrMesh&);
int GetLinearBlockIndex(int , int , int , int , int );
void GetOrigMeshInfo(std::string , int &, double &, std::vector<double> &);
