//****f* LAVAflow/Driver::finalize
//  NAME
//    Driver::finalize
//
//  SYNOPSIS
//    Driver::finalize()
//
//  DESCRIPTION
//    Finalizes any operations for the Driver class. This includes
//    closing logging files, shutting down MPI, etc.
//
//  NOTES
//    This function is currently empty
//
//  SEE ALSO
//    Driver::initialize
//****

#include "Driver.h"
void Driver::finalize()
{
	// right now this probably isn't necessary, because the 
	// LavaMPI destructor already finalizes MPI
	driverMPI.finalize();
}
