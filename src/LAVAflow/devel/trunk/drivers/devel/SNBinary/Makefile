
include ../../sites/site.mk

LIB_DIR = ../../src
INCLUDE_DIR = ../../include

EXTLIB = -L$(HDF_INSTALL)/lib -L$(MPI_INSTALL)/lib -L$(SZIP_INSTALL)/lib


CFLAGS = -Wall -g
LFLAGS = -lz -lm $(HOST_LIB_FLAGS) -lsz -ldl

INCLUDE = -I$(HDF_INSTALL)/include -I$(MPI_INSTALL)/include -I$(SZIP_INSTALL)/include -I$(INCLUDE_DIR) -I.


SOURCES=analysisCommon.cpp $(LIB_DIR)/flash_reader.cpp $(LIB_DIR)/outside_functions.cpp

# OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS=analysisCommon.o flash_reader.o outside_functions.o

EXECUTABLE1=colDensSNBinary
EXECUTABLE2=strippedMassSNBinary
EXECUTABLE3=eintTempSNBinary
EXECUTABLE4=eintTimeSNBinary
EXECUTABLE5=velDistSNBinary
EXECUTABLE6=contaminationSNBinary
EXECUTABLE7=radialDistributionSNBinary

EXECUTABLES=$(EXECUTABLE1) $(EXECUTABLE2) $(EXECUTABLE3) $(EXECUTABLE4) $(EXECUTABLE5) $(EXECUTABLE6) $(EXECUTABLE7)

all: $(EXECUTABLES)

$(EXECUTABLE1): colDensSNBinary.o $(OBJECTS)
	$(CC) colDensSNBinary.o $(OBJECTS) -o $@ $(EXTLIB) $(LFLAGS)

$(EXECUTABLE2): strippedMassSNBinary.o $(OBJECTS)
	$(CC) strippedMassSNBinary.o $(OBJECTS) -o $@ $(EXTLIB) $(LFLAGS)

$(EXECUTABLE5): eintTimeSNBinary.o $(OBJECTS)
	$(CC) eintTimeSNBinary.o $(OBJECTS) -o $@ $(EXTLIB) $(LFLAGS)

$(EXECUTABLE3): eintTempSNBinary.o $(OBJECTS)
	$(CC) eintTempSNBinary.o $(OBJECTS) -o $@ $(EXTLIB) $(LFLAGS)

$(EXECUTABLE4): velDistSNBinary.o $(OBJECTS)
	$(CC) velDistSNBinary.o $(OBJECTS) -o $@ $(EXTLIB) $(LFLAGS)

$(EXECUTABLE6): contaminationSNBinary.o $(OBJECTS)
	$(CC) contaminationSNBinary.o $(OBJECTS) -o $@ $(EXTLIB) $(LFLAGS)

$(EXECUTABLE7): radialDistributionSNBinary.o $(OBJECTS)
	$(CC) radialDistributionSNBinary.o $(OBJECTS) -o $@ $(EXTLIB) $(LFLAGS)

# .cpp.o:
# 	$(CC) $(CFLAGS) $< -o $@ $(INCLUDE)

analysisCommon.o: analysisCommon.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $<

colDensSNBinary.o: colDensSNBinary.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $< 

strippedMassSNBinary.o: strippedMassSNBinary.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $<

eintTempSNBinary.o: eintTempSNBinary.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $<

eintTimeSNBinary.o: eintTimeSNBinary.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $<

velDistSNBinary.o: velDistSNBinary.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $<

contaminationSNBinary.o: contaminationSNBinary.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $<

radialDistributionSNBinary.o: radialDistributionSNBinary.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $<

flash_reader.o: $(LIB_DIR)/flash_reader.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $< 

outside_functions.o: $(LIB_DIR)/outside_functions.cpp
	$(CC) -o $@ $(INCLUDE) $(CFLAGS) -c $<

clean:
	rm -rf *o $(EXECUTABLES)
