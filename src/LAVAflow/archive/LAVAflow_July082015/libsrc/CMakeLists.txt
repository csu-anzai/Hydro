

ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Math) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Geometry) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Utilities) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Clustering) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Graph) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Particles) 
ADD_SUBDIRECTORY(${LAVAflow_SOURCE_DIR}/libsrc/Statistics) 


set(LAVAflow_INCLUDES
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
	${LAVAflow_SOURCE_DIR}/libsrc/Mesh/Mesh
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
	${VTK_LIBRARIES}
)
