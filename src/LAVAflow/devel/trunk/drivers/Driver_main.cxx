//****f* drivers/Driver_main
//  NAME
//    Driver_main
//
//  SYNOPSIS
//    int main(int argc, char** argv)
//
//  DESCRIPTION
//    This function provides MAIN for drivers which include it in their
//    CMakeLists.txt sources. It is meant to wrap the analysis function
//    for sharable functionality across all drivers. The actual analysis
//    is executed by providing a function pointer to the analysis routine
//    to a Driver object, which then calls Driver::execute.
//
//  NOTES
//    For this to function properly, the analysis function must be declared as
//
//       int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
//
//    where LAVAFLOW_DRIVER_BASENM is defined in Driver_main.h. This macro
//    will be expanded at compilation time to provide a true function name.
//
//    The rationale behind defining the analysis function name as a macro
//    is to avoid having problems during the linking step. By doing it this way,
//    we can ensure that the user-implemented function adheres to the same
//    naming expected here to create the function pointer. Note that
//    instead of declaring the LAVAFLOW_DRIVER_BASENM routine in an additional
//    header file, we can declare it instead inside this file assuming the
//    user used the LAVAFLOW_DRIVER_BASENM as their function name.
//
//  TODO
//    Add more features.
//
//  SEE ALSO
//    Driver
//
//****

#include "Driver.h"
#include "Driver_main.h"

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv);

int main(int argc, char** argv)
{
    // Setup the Driver object
    Driver driver(argc, argv);
    driver.initialize(argc,argv);
    driver.setExecuteFunction(LAVAFLOW_DRIVER_BASENM);

    // Run the LAVAflow driver
    driver.execute(argc, argv);

    // Finalize Driver operations
    driver.finalize();
}
