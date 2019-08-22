#ifndef _LF_MPI_H
#define _LF_MPI_H

#include <mpi.h>


class LavaFlowMPI : MPI::Intracomm
{
	int mpiID;
	int nProcs;

public:

	LavaFlowMPI();
	LavaFlowMPI(int argc, char** argv);
	~LavaFlowMPI();

	void Init( int, char** );
	void Finalize();
	int GetID();
	int GetNumProcs();
};

#endif