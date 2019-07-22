set(MODULENAME LAVAGEOM)

INCLUDE_DIRECTORIES(${LAVAflow_INCLUDES})

set(${MODULENAME}_SOURCES
	Geometry.cxx
	Primitives/VecTN.cxx
	Primitives/Vec3.cxx
	Primitives/Edge3.cxx

	intersectionLinePlane.cxx
	clipPolygonUsingPlane.cxx
	isPointInfrontOfPlane.cxx

	HEMesh/HEDataStructure.cxx
	HEMesh/HEPolytope.cxx
	HEMesh/HEFace.cxx
	HEMesh/HEEdge.cxx
	HEMesh/HEVertex.cxx

)

ADD_LIBRARY(${MODULENAME} STATIC ${${MODULENAME}_SOURCES})
TARGET_LINK_LIBRARIES(${MODULENAME})

list(APPEND LAVAflow_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND LAVAflow_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/HEMesh)
list(APPEND LAVAflow_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/Primitives)
set(LAVAflow_INCLUDES ${LAVAflow_INCLUDES} PARENT_SCOPE)
