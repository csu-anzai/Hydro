#include "LavaFlowMPI.h"


LavaFlowMPI::LavaFlowMPI()
{
	if (!MPI::Is_initialized()) MPI::Init();
	mpiID = MPI::COMM_WORLD.Get_rank();
	nProcs = MPI::COMM_WORLD.Get_size();
}

LavaFlowMPI::LavaFlowMPI(int argc, char** argv)
{
	if (!MPI::Is_initialized()) MPI::Init( argc, argv );
	mpiID = MPI::COMM_WORLD.Get_rank();
	nProcs = MPI::COMM_WORLD.Get_size();
}

LavaFlowMPI::~LavaFlowMPI()
{
	if (!MPI::Is_finalized()) MPI::Finalize();
}

void LavaFlowMPI::Init( int argc, char** argv ) {
	if (!MPI::Is_initialized()) MPI::Init( argc, argv );
	mpiID = MPI::COMM_WORLD.Get_rank();
	nProcs = MPI::COMM_WORLD.Get_size();
}

void LavaFlowMPI::Finalize() {
	if (!MPI::Is_finalized()) MPI::Finalize();
}

int LavaFlowMPI::GetID() {
	return mpiID;
}

int LavaFlowMPI::GetNumProcs() {
	return nProcs;
}