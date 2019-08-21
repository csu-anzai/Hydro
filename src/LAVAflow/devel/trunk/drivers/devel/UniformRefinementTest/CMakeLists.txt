SET(DRIVER_NAME "UniformRefinementTest")

SET(DRIVER_SOURCES
    ${DRIVER_NAME}.cxx ${DRIVER_SOURCES}
    ${LAVAflow_SOURCE_DIR}/drivers/Driver_main.cxx
    ${LAVAflow_SOURCE_DIR}/libsrc/Mesh/UniformRefinement/UniformRefinement.cxx)

MESSAGE("LAVAflow source dir is ${LAVAflow_SOURCE_DIR}")

ADD_EXECUTABLE(${DRIVER_NAME} ${DRIVER_SOURCES})
TARGET_LINK_LIBRARIES(${DRIVER_NAME} LAVAflow ${HDF5_LIBRARIES})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET(GCC_C++11_COMPILE_FLAGS "-std=gnu++11")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_C++11_COMPILE_FLAGS}" )
else()
	MESSAGE("Not using gnu compiler. Compilation may fail if compiler does not support c++11 without additional flags.")
endif()

