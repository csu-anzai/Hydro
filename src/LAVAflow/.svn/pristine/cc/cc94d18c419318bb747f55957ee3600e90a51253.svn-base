

# Prepare list of include directories
#list(APPEND LAVAflow_INCLUDES "${LAVAflow_SOURCE_DIR}/libsrc/Turbulence")

ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Mesh)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Math)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Geometry)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Utilities)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Clustering)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Graph)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Particles)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Statistics)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/MPI)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Driver)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Turbulence)


set(LAVAflow_INCLUDES ${LAVAflow_INCLUDES} PARENT_SCOPE)

MESSAGE( STATUS "Include directories:  ${LAVAflow_INCLUDES}" )

set(LAVAflow_SOURCES
	${LAVAflow_SOURCE_DIR}/libsrc/Turbulence/TimeAutoCorr
)

INCLUDE_DIRECTORIES(${LAVAflow_INCLUDES})
ADD_LIBRARY(LAVAflow STATIC ${LAVAflow_SOURCES})

TARGET_LINK_LIBRARIES(LAVAflow
	LAVAMATH
	LAVAGEOM
	LAVAUTIL
	LAVACLUST
	LAVAGRAPH
	LAVAPARTICLES
	LAVASTAT
	LAVAMESH
	LAVAMPI
	LAVADRIVER
	LAVATURB
	LAVAUNIFORMREFINEMENT
)

MESSAGE( STATUS "inc_list:  ${inc_list}" )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET(GCC_C++11_COMPILE_FLAGS "-std=gnu++11")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_C++11_COMPILE_FLAGS}" )
else()
	MESSAGE("Not using gnu compiler. Compilation may fail if compiler does not support c++11 without additional flags.")
endif()
