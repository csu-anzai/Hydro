#============================================================
#
# Top-level CMakeLists for the LAVAflow library + drivers
#
#
#============================================================
cmake_minimum_required(VERSION 2.8)

#========================================
# Set compiler locations
#========================================
set(CMAKE_CXX_COMPILER     $ENV{MPIHOME}/bin/mpicxx)
set(CMAKE_C_COMPILER       $ENV{MPIHOME}/bin/mpicc)
set(CMAKE_Fortran_COMPILER $ENV{MPIHOME}/bin/mpif90)


#========================================
# Project information
#========================================
PROJECT(LAVAflow)
ENABLE_LANGUAGE(CXX)

#========================================
# Set compiler flags
#========================================
SET( CMAKE_CXX_FLAGS_RELEASE "-O2 -Wuninitialized -std=c++11")
SET( CMAKE_CXX_FLAGS_DEBUG   "-ggdb -pg -Wno-div-by-zero -Wundef -Wconversion -Wstrict-prototypes -Wunreachable-code -pedantic -Wall -Wextra -Winit-self -ftree-vrp -Wfloat-equal -Wunsafe-loop-optimizations -Wpadded -fstack-check -fstack-protector-all -mcmodel=large -std=c++11")
SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})

SET( STOP_ON_FIRST_ERROR "-Wfatal-errors")
SET( CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} ${STOP_ON_FIRST_ERROR}" )
SET( CMAKE_CXX_FLAGS_DEBUG    "${CMAKE_CXX_FLAGS_DEBUG}   ${STOP_ON_FIRST_ERROR}" )
SET( CMAKE_CXX_FLAGS          "${CMAKE_CXX_FLAGS}         ${STOP_ON_FIRST_ERROR}" )


#========================================
# Set linker flags
#========================================
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -ldl")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG   "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -ldl")
set(CMAKE_EXE_LINKER_FLAGS         "${CMAKE_EXE_LINKER_FLAGS} -ldl")

#========================================
# External libraries
#========================================

#-------- HDF5 --------
set(HDF5_DIR $ENV{HDF5HOME})
include_directories(${HDF5_DIR}/include)
link_directories(${HDF5_DIR}/lib)
set(HDF5_LIBRARIES ${HDF5_DIR}/lib/libhdf5.a ${HDF5_DIR}/lib/libhdf5_cpp.a)

#-------- SZIP --------
set(SZIP_DIR $ENV{SZIPHOME})
include_directories(${SZIP_DIR}/include)
link_directories(${SZIP_DIR}/lib)
list(APPEND HDF5_LIBRARIES ${SZIP_DIR}/lib/libsz.a)

#-------- ZLIB --------
set(ZLIB_DIR $ENV{ZLIBHOME})
include_directories(${ZLIB_DIR}/include)
link_directories(${ZLIB_DIR}/lib)
list(APPEND HDF5_LIBRARIES ${ZLIB_DIR}/lib/libz.a)

#-------- MPI --------
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

#========================================
# Diagnostics
#========================================
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)

message("")
message("#===============#")
message("# Build Info    #")
message("#===============#")
message("Build type:        ${CMAKE_BUILD_TYPE}")
message("CXX compiler:      ${CMAKE_CXX_COMPILER}")
message("CXX flags NONE:    ${CMAKE_CXX_FLAGS}")
message("CXX flags RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message("CXX flags DEBUG:   ${CMAKE_CXX_FLAGS_DEBUG}")
message("Link flags:        ${CMAKE_EXE_LINKER_FLAGS}")
message("")

message("#===============#")
message("# MPI Info      #")
message("#===============#")
message("MPI includes:      ${MPI_INCLUDE_PATH}")
message("MPI libraries:     ${MPI_LIBRARIES}")
message("")

message("#===============#")
message("# HDF5 Info     #")
message("#===============#")
message("HDF5 libraries:      ${HDF5_LIBRARIES}")
message("")

#========================================
# Initialize include directories
#========================================
list(APPEND LAVAflow_INCLUDES "${LAVAflow_SOURCE_DIR}/libsrc/includes")

#========================================
# Add subdirectories to traverse
#========================================
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/drivers)
