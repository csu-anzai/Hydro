#ifndef _LAVA_MPI_H
#define _LAVA_MPI_H

#include <mpi.h>

class LavaMPI : MPI::Intracomm
{
	int mpiID;
	int nProcs;

public:

	LavaMPI();
	LavaMPI(int argc, char** argv);
	~LavaMPI();

	void init( int, char** );
	void finalize();
	int getID();
	int getNumProcs();
};

#endif
