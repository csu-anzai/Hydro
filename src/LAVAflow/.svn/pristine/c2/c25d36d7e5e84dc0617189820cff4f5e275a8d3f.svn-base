//****h* LAVAflow/Driver
//  NAME
//    Driver
//
//  DESCRIPTION
//    Driver class which is meant to serve as a wrapper for user-written
//    analyses. This may initialize/finalize things such as high-level IO
//    (logging), MPI, etc. It is meant to be used inside of a "main" function,
//    from which a pointer to the user-defined analysis function
//    (mimicing the "main") is provided and executed by the Driver object.
//
//  METHODS
//     Driver::Driver
//     Driver::~Driver
//     Driver::initialize
//     Driver::finalize
//     Driver::setExecuteFunction
//     Driver::execute
//****

#include "LavaMPI.h"

#ifndef LAVA_DRIVER_H
#define LAVA_DRIVER_H

// Definition of a function pointer which mimics "main".
// All user-defined functions wishing to use the Driver
// class must adhere to this function definition.
typedef int (*executeFunctionType)(int, char**);

class Driver
{
    private:

        //===========================
        // Executor function
        //===========================
        bool isExecuteFunctionSet;
        executeFunctionType executeFunction;

    public:

        LavaMPI driverMPI;

        //===========================
        // Constructors/Destructors
        //===========================
        Driver();
        Driver(int argc, char** argv);
        ~Driver();

        //===========================
        // Initializer/Finalizer
        //===========================
        void initialize(int argc, char** argv);
        void finalize();

        //===========================
        // Executor
        //===========================
        void setExecuteFunction(executeFunctionType ef);
        void execute(int argc, char** argv);

};

#endif

