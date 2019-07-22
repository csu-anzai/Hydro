//****f* LAVAflow/Driver::execute
//  NAME
//    Driver::execute
//
//  SYNOPSIS
//    void Driver::execute(int argc, char** argv)
//
//  DESCRIPTION
//    Provided a user-defined function (with an interface mimicing "main")
//    has been set using Driver::setExecuteFunction, this function:
//    - Performs any high-level, pre-execution tasks
//    - Executes the provided user-defined function
//    - Performs any high-level, post-execution tasks (i.e. error handling)
//
//  EXAMPLE
//
//    int LAVAFLOW_DRIVER_BASENM(int argc, char** argv);
//    ...
//    Driver driver(...)
//    driver.setExecuteFunction(LAVAFLOW_DRIVER_BASENM);
//    driver.execute(argc, argv);
//
//  TODO
//    - Provide abort functionality if the execution function has not been set
//    - Handle return values from the execution function
//
//  SEE ALSO
//    Driver::setExecuteFunction
//
//****
#include "Driver.h"
#include <iostream>

void Driver::execute(int argc, char** argv)
{
    // Confirm that the execution function is set
    if(!this->isExecuteFunctionSet)
    {
        std::cerr<<"[Driver::execute] Execution function is not set."<<std::endl;
        // TODO: Abort
    }

    // Call the execution function (i.e. call the primary analysis driver)
    int res = this->executeFunction(argc, argv);

    // TODO: Handle executeFunction return values

}
