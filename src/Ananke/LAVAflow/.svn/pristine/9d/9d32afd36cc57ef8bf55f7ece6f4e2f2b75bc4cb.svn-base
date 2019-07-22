//****f* LAVAflow/Driver::setExecuteFunction
//  NAME
//    Driver::setExecuteFunction
//
//  SYNOPSIS
//    Driver::setExecuteFunction(executeFunctionType exefunc)
//
//  DESCRIPTION
//    This accessor stores a function pointer which the Driver object executes.
//    The provided function is a user-defined routine whose interface
//    mimics that of "main".
//
//  EXAMPLE
//
//    int LAVAFLOW_DRIVER_BASENM(int argc, char** argv);
//    ...
//    Driver driver(...)
//    driver.setExecuteFunction(LAVAFLOW_DRIVER_BASENM);
//
//  NOTES
//    The executeFunctionType is defined in Driver and is a pointer to a
//    function whose declaration is identifcal to that of "main":
//    int (*f)(int, char**).
//
//  SEE ALSO
//    Driver
//    Driver::execute
//
//****

#include "Driver.h"
void Driver::setExecuteFunction(executeFunctionType exefunc)
{
    // Set the function to execute in Driver::execute and
    // flag that it has been set.
    this->executeFunction = exefunc;
    this->isExecuteFunctionSet = true;
}
