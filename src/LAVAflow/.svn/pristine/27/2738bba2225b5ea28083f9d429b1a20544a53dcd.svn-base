set(MODULENAME LAVAGEOM)

INCLUDE_DIRECTORIES(${LAVAflow_INCLUDES})

MESSAGE( STATUS "LAVA_LIBS: " ${LAVA_LIBS})

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
TARGET_LINK_LIBRARIES(${MODULENAME} ${VTK_LIBRARIES})
