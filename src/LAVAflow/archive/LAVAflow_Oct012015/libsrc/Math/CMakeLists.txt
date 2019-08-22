INCLUDE_DIRECTORIES(${LAVAflow_INCLUDES})
MESSAGE( STATUS "Include directories: " ${LAVAflow_INCLUDES} )

set(LAVAMATH_SOURCES
	dot.cxx
	cross.cxx
	legendre.cxx
	sphericalHarmonic.cxx
	factorial.cxx
)

ADD_LIBRARY(LAVAMATH STATIC ${LAVAMATH_SOURCES})
TARGET_LINK_LIBRARIES(LAVAMATH ${VTK_LIBRARIES})
