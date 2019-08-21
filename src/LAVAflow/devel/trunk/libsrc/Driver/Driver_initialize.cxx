//****f* LAVAflow/Driver::initialize
//  NAME
//    Driver::initialize
//
//  SYNOPSIS
//    Driver::initialize()
//
//  DESCRIPTION
//    Initializes any operations for the Driver class. This includes
//    opening logging files, starting MPI, etc.
//
//  NOTES
//    This function is currently empty
//
//  SEE ALSO
//    Driver::finalize
//****

#include "Driver.h"
void Driver::initialize(int argc, char** argv)
{
	// right now this probably doesn't do anything, because the 
	// LavaMPI default constructor already initializes MPI
	driverMPI.init(argc, argv);
}
