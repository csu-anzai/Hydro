#include "Driver_main.h"
#include <iostream>
#include "FlashParticles.h"
#include <vector>
#include <algorithm>
#include "ArrayOperations.h"

using namespace std;

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{
	string filename = "/data1/data/tburn_hdf5_part_0001";
	vector<string> vars = {"blk", "velx"};

	std::vector<const char*> varsChar = VecStringToVecChar(vars);
    FlashParticles particles(filename.c_str(), &varsChar[0], vars.size(), PART_SORT_TAG);

	int nParts = particles.nParticlesLocal;

	int blkInd = particles.findVarIndex("blk");

	double** partData = particles.getDataPtr();

	cout << "Rank " << MPI::COMM_WORLD.Get_rank() << " block IDs" << endl;
	for (int i=0; i<5; i++)
		cout << partData[blkInd][i] << endl;

    int nPartsLocl = particles.nParticlesLocal;
    int nPartsGlobl = particles.nParticlesGlobal;

    cout << "Number of local particles is " << nPartsLocl << endl;
    cout << "Number of global particles is " << nPartsGlobl << endl;
}

