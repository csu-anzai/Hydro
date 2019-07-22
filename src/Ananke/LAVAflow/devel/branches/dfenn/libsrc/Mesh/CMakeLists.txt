set(MODULENAME LAVAMESH)

INCLUDE_DIRECTORIES(${LAVAflow_INCLUDES})

set(${MODULENAME}_SOURCES
	Mesh.cxx
	FlashAmrMesh.cxx
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET(GCC_C++11_COMPILE_FLAGS "-std=gnu++11")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_C++11_COMPILE_FLAGS}" )
else()
	MESSAGE("Not using gnu compiler. Compilation may fail if compiler does not support c++11 without additional flags.")
endif()

ADD_LIBRARY(${MODULENAME} STATIC ${${MODULENAME}_SOURCES})
TARGET_LINK_LIBRARIES(${MODULENAME} ${VTK_LIBRARIES})
