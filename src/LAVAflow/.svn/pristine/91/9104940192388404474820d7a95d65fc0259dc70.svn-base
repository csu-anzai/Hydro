#include "LavaMPI.h"


LavaMPI::LavaMPI()
{
	if (!MPI::Is_initialized()) MPI::Init();
	mpiID = MPI::COMM_WORLD.Get_rank();
	nProcs = MPI::COMM_WORLD.Get_size();
}

LavaMPI::LavaMPI(int argc, char** argv)
{
	if (!MPI::Is_initialized()) MPI::Init( argc, argv );
	mpiID = MPI::COMM_WORLD.Get_rank();
	nProcs = MPI::COMM_WORLD.Get_size();
}

LavaMPI::~LavaMPI()
{
	if (!MPI::Is_finalized()) MPI::Finalize();
}

void LavaMPI::init( int argc, char** argv )
{
	if (!MPI::Is_initialized()) MPI::Init( argc, argv );
	mpiID = MPI::COMM_WORLD.Get_rank();
	nProcs = MPI::COMM_WORLD.Get_size();
}

void LavaMPI::finalize()
{
	if (!MPI::Is_finalized()) MPI::Finalize();
}

int LavaMPI::getID()
{
	return mpiID;
}

int LavaMPI::getNumProcs()
{
	return nProcs;
}
