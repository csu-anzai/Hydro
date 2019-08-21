

ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Math) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Geometry) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Utilities) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Clustering) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Graph) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Particles) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Statistics)
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Mesh) 

set(LAVAflow_INCLUDES
	${MPI_INCLUDE_PATH}
	${VTK_INCLUDE_DIR}
	${LAVAflow_SOURCE_DIR}/libsrc/
	${LAVAflow_SOURCE_DIR}/libsrc/includes
	${LAVAflow_SOURCE_DIR}/libsrc/Operators
	${LAVAflow_SOURCE_DIR}/libsrc/Selectors
#	${LAVAflow_SOURCE_DIR}/libsrc/Renderer
	${LAVAflow_SOURCE_DIR}/libsrc/Readers
#	${LAVAflow_SOURCE_DIR}/libsrc/Geometry
	${LAVAflow_SOURCE_DIR}/libsrc/Math
	${LAVAflow_SOURCE_DIR}/libsrc/Utilities
	${LAVAflow_SOURCE_DIR}/libsrc/Mesh
	CACHE INTERNAL "Global includes for all files"
)

MESSAGE( STATUS "Include directories: " ${LAVAflow_INCLUDES} )

set(LAVAflow_SOURCES
#	${LAVAflow_SOURCE_DIR}/libsrc/Operators/MultiplyOperator
#	${LAVAflow_SOURCE_DIR}/libsrc/Operators/ContourOperator
#	${LAVAflow_SOURCE_DIR}/libsrc/Operators/SurfaceAreaOperator
#	${LAVAflow_SOURCE_DIR}/libsrc/Operators/CurvatureOperator
	${LAVAflow_SOURCE_DIR}/libsrc/Operators/shellAveraging/ShellAverageOperator
	${LAVAflow_SOURCE_DIR}/libsrc/Operators/shellAveraging/planeShells/ShellAveragePlaneOperator
	${LAVAflow_SOURCE_DIR}/libsrc/Operators/surfaceAveraging/SurfaceAverageOperator
	${LAVAflow_SOURCE_DIR}/libsrc/Operators/surfaceAveraging/planarSurface/SurfaceAveragePlaneOperator
#	${LAVAflow_SOURCE_DIR}/libsrc/Selectors/NullSelector
#	${LAVAflow_SOURCE_DIR}/libsrc/Rendering/Renderer
	${LAVAflow_SOURCE_DIR}/libsrc/Readers/VTKReader
	${LAVAflow_SOURCE_DIR}/libsrc/MPI/LavaFlowMPI
	${LAVAflow_SOURCE_DIR}/libsrc/Turbulence/TimeAutoCorr
)

INCLUDE_DIRECTORIES(${LAVAflow_INCLUDES})
ADD_LIBRARY(LAVAflow STATIC ${LAVAflow_SOURCES})

MESSAGE( STATUS "LAVA_LIBS: " ${LAVA_LIBS})
TARGET_LINK_LIBRARIES(LAVAflow
	LAVAMATH
	LAVAGEOM
	LAVAUTIL
	LAVACLUST
	LAVAGRAPH
	LAVAPARTICLES
	LAVASTAT
	LAVAMESH
	${VTK_LIBRARIES}
)


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET(GCC_C++11_COMPILE_FLAGS "-std=gnu++11")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_C++11_COMPILE_FLAGS}" )
else()
	MESSAGE("Not using gnu compiler. Compilation may fail if compiler does not support c++11 without additional flags.")
endif()
