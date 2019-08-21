# Set the name of the driver program
SET(DRIVER_NAME "FlashParticlesTest")

# Add additional source files that will live in the driver directory.
# This should be a space or newline separated list.
# e.g. SET(DRIVER_SOURCES file1.cxx file2.cxx)
# e.g. SET(DRIVER_SOURCES
#              file1.cxx
#              file2.cxx
#      )
SET(DRIVER_SOURCES ${LAVAflow_SOURCE_DIR}/drivers/Driver_main.cxx)

# Add the executable and link targets
ADD_EXECUTABLE(${DRIVER_NAME} ${DRIVER_NAME}.cxx ${DRIVER_SOURCES})
TARGET_LINK_LIBRARIES(${DRIVER_NAME} LAVAflow ${MPI_LIBRARIES})


