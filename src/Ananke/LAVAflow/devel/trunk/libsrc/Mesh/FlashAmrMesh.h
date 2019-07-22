#include <hdf5.h>
#include <mpi.h>
#include "Mesh.h"
#include "FlashIO.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <cstddef>        // std::size_t

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif

// #ifndef H5_NO_NAMESPACE
	// using namespace H5;
// #endif

#ifndef _FLASH_AMR_MESH_H
#define _FLASH_AMR_MESH_H

#define MESH_VAR_STRING_SIZE 4
#define FR_LEAF_NODE 1
#define FR_MAXVARS 100
// #define FR_MAXBLOCKS 

//Only using HDF5, at least for now.
#define FR_HDF5 1

//Make this boolean true if you want the program to crash when there's a failure
// with deallocation.
#define FORCEDEALLOCATECRASH true
#define USEGUARDCELLS false
#define NGUARD 0



#define SUPPRESSOUTPUT true

enum BLOCKTYPE      { ALL_BLKS, IBDRY_BLKS, JBDRY_BLKS, KBDRY_BLKS, ANY_BDRY_BLKS, 
                      ACTIVE_BLKS, LEAF, PARENT_BLK, ANCESTOR, REFINEMENT, INREGION};


typedef struct block {
	int blockNumber;			  /* Local processor's block number. Included for completeness... */
	int cpuNumber;                /* CPU number on which the block resided during the simulation. */

	int nodeType;				  /* Location in AMR tree. =1 means this is a leaf block */
	int refineLevel;              /* Refinement of the block. =1 means this is the most coarse block */

	double size[MESH_MDIM];		  /* Physical size of block */
	double coord[MESH_MDIM];	  /* Center of block */
	double lowerBound[MESH_MDIM]; /* coord - size/2 */
	double upperBound[MESH_MDIM]; /* coord + size/2 */

	double ****data;              /* 4D array containing all of the mesh variable data.
	                                 Indexed as data[variable][x][y][z] */

//	double *x;   				  /* Cell coordinates (not yet implemented in back end) */
//	double *y;
//	double *z;
//	double dx[MESH_MDIM];	      /* Cell lengths per dimension (not yet implemented) */

//	double vol;					  /* Cell volumes (not yet implemented) */

	int parent;                   /* Parent, children, neighbor information from GID array information */
	int *children;
	int *neighbor;                /* Relative positions in neighbor array: <x, >x, <y, >y, <z, >z */

//	block *allNeighbors[MESH_MDIM][MESH_MDIM][MESH_MDIM]; //All neighbors

} block;


class FlashAmrMesh:Mesh, public FlashIO {

private:
	block  *blocks;

	int **gid;                        /* Global ID array used to find parent, children, and neighbors for blocks. */
	int *procGidMap;

	int    nz[MESH_MDIM+1];     	  /* Dimensions + nvar ( auxiliary ) */
	int    nzGuardLow [MESH_MDIM+1];  /* Dimensions with guard cells + nvar (auxiliary) */
	int    nzGuardHigh[MESH_MDIM+1];
	int    nguard;  				  /* Number of guard cells in each direction */

	int    nFileVars; 				  /* Number of field variables found in file */
	int    nUserVars;	 			  /* Number of field variables requested by the user (file and derived vars) */
	int    nReqFileVars;              /* Number of field variables requested by the user (file variables only) */
	int    nReqAuxVars;               /* Number of field variables requested by the user (derived variables only) */
	char   k2d; 					  /* Set to 1 if at least 2d is present */
	char   k3d; 					  /* Set to 1 if at least 3d is present */

	int    nCellsVec[MESH_MDIM];	  /* Number of cells in each direction */
	int    nCells;					  /* Number of cells in a block: Product of nCellsVec's components. */

  	int    iteration;          		  /* Current cycle (iteration) */

	int    format;					  /* Format of the file (Will just be HDF5) */
	int    nBlocks;                   /* Total number of blocks on the mesh */
	
	hid_t inFile;                     /* Input HDF5 file */
	hid_t outFile;                    /* Output HDF5 file (if used) */

	char   variableNames[FR_MAXVARS][MESH_VAR_STRING_SIZE+1]; /* Names of all stored variables */
	char   variablesInFile[FR_MAXVARS][MESH_VAR_STRING_SIZE+1]; /* Names of stored variables found in read-in file */
	char   variablesAux[FR_MAXVARS][MESH_VAR_STRING_SIZE+1]; /* Names of derived (auxiliary) stored variables */

	int fileType;                     /* FLASH file type: Used to see if it is type '7' or not. 
	                                     Some different functions for reading of scalars and parameters 
	                                     must be used if it is. */

	int nBlocksAssigned;              /* Total number of blocks assigned to an MPI processor. */
	MPI::Intracomm comm;

	std::map<std::string,bool> writeVarMap; /* Whether or not the specified variable (read in or derived) will be written.
	                                          Default is true. */


	// /* Variables used only in the uniform refinement function refineToFinest() */
	int    fineCells[MESH_MDIM];
	int    validCells[MESH_MDIM];
	int    totalCells[MESH_MDIM];
	int    numProcs[MESH_MDIM];
	double gridDelta[MESH_MDIM];
	double refinedDomainBoundBox[2][MESH_MDIM];
	int    **procIndexMin;
	bool   variablesAllocated;

public:

	FlashAmrMesh();
	FlashAmrMesh(std::string, std::vector<std::string>);
	FlashAmrMesh(std::string, std::vector<std::string>, std::vector<std::string>);
	FlashAmrMesh(std::string, std::string, std::vector<std::string>);
	FlashAmrMesh(std::string, std::string, std::vector<std::string>, std::vector<std::string>);
	~FlashAmrMesh();

	//These two should be removed once I figure out how to initialize the file correctly
	//in the constructor.
	void setGridFileName(std::string);
	std::string getGridInFileName(void);

	void setGridOutFileName(std::string);
    std::string getGridOutFileName(void);
	// void    readFile(std::string);

	// void LoadMesh(std::string , std::vector<std::string> );

	void    generalRead();

	void    readInts(const char*);
	void    readReals(const char*);
	void    readStrings(const char*);
	void    readLogicals(const char*);
	void    readSimParameters(void);

	int     getInt(const char*);
	double  getReal(const char*);
	int     getFileType();
	void	setLogical(const char*, const bool);

	void    allocateBlocks();
	void    deallocateBlocks();

	void    allocateGid();
	void    deallocateGid();

	void    readBlockNumber();
	int     getBlockNumber(block);

	void    readNcells(void);
	int     getNcells(int); 
	int*    getNcells(void);
	int     getTotalCells(void);

	void    readDim(void);
	int     getDim(void);
	
	void    readNBlocks(void);
	int     getNBlocks(void);

	void    readGeometryType(void);
	int     getGeometryType(void);

	void    readCoordinates(void);
	double* getCoordinates(int);
	double  getCoordinate(int,int);

	void    readSizes(void);
	double* getSizes(int);
	double  getSize(int,int);

	void    readBounds(void);
	double* getBounds(int, int);
	double  getBound(int, int, int);

	bool    datasetAvailable(const char*, const char*);
	void    readCPUNumbers(void);
	int     getCPUNumber(int);

	void    readNodeType(void);
	int     getNodeType(int);

	void    readRefineLevel(void);
	int     getRefineLevel(int);

	void    readTime(void);
	double  getTime(void);

	void    readGid(void);
	void 	assignGidArrays(void);
	void    readFileType(void);

	//std::vector<std::vector<int> >&    getGidVector();

	void    addVariableNames(std::vector<std::string>,std::vector<std::string>);
	void    readNFileVarsList(void);
	int     getNFileVars(void);
	int     getNUserVars(void);

	void    allocateVarArray(void);
	void    deallocateVarArray(void);
	int     findVarIndex(const char*);
	std::string  findVarName(const int);
	void    readVarData(const int, const char*);
	void    initAuxVarData(const int, const char*);

	void    readVariables(void);
    void    noWrite(const char*);
	void    addUserVariable(const char*, const char*);
	double  getVariable(int, int, int, int, int);

	void    deallocateVariables(void);

	void 	openInFile(void);
	void    closeInFile(void);

	void    readInFile(std::vector<std::string>,std::vector<std::string>);

    // Writing functions
    void    createOutFile();
    void    generalWrite();
    void    closeOutFile();
    void    writeOutFile();

	void    writeFileType();
	void    writeSimParameters();
	void    writeSimInfo();
	void    writeInts(const char*);
	void    writeLogicals(const char*);
	void    writeReals(const char*);
	void    writeStrings(const char*);
	void    writeCoordinates();
	void    writeCPUNumbers();
	void    writeSizes();
	void    writeBounds();
	void    writeNodeType();
	void    writeRefineLevel();
	void    writeGid();
	void    writeWhichChild();
	void    writeBFlags();
	void    writeNFileVarsList();
	void    writeVarData(const char*);
	void    writeVariables();


	//MPI Functions
	int mpiGetID();
	int mpiGetProcs();
	void mpiAssignBlocks();
	int getLocalNumBlocks();
	int mpiGetGlobalBlockID(int);
	int mpiGetLocalBlockID(int);

	//Flash equivalents
	void fillGuardCells(void); //done
	void getBlkBoundBox(int, double[2][MESH_MDIM]); //done
	void getBlkData(const int, GRIDDATASTRUCT, const char[MESH_VAR_STRING_SIZE+1], 
					BEGINCOUNT, int[MESH_MDIM], double****& /*, int* */); //done
	void getBlkData(const int, GRIDDATASTRUCT, int, 
				BEGINCOUNT, int[MESH_MDIM], double****& /*, int* */); //done
	void getBlkIndexLimits(int, int[2][MESH_MDIM], int[2][MESH_MDIM], GRIDDATASTRUCT=CENTER); //not done
	void getBlkPhysicalSize(int, double[MESH_MDIM]); //Done
	void getBlkPtr(int, double****&, GRIDDATASTRUCT); //Done
	void getDomainBoundBox(double[2][MESH_MDIM]); //Done
	void getListOfBlocks(BLOCKTYPE, int*, int&, int=-1, double[2][MESH_MDIM]=NULL, bool=false);
	void putBlkData(int, GRIDDATASTRUCT, const char[MESH_VAR_STRING_SIZE+1], 
					BEGINCOUNT, int[MESH_MDIM], double****& /*, int* */);
	void releaseBlkPtr(int, double****&, GRIDDATASTRUCT); //Done
	void getCellBoundBox(int, BEGINCOUNT, int, int, int, double[2][MESH_MDIM]);
	double getSingleCellVol(int, BEGINCOUNT=EXTERIOR, int[MESH_MDIM]=NULL);
	void getCellCoords(int, int, int, bool, double*);
	void getCellSideLengths(int , double[MESH_MDIM] );
	double getPointData(int, GRIDDATASTRUCT, int, BEGINCOUNT, int[MESH_MDIM]);
	void getCellCoords(int , int , int , int , double[MESH_MDIM]);


	void refineToFinest(double****&, double[2][MESH_MDIM], int=-1);
	void printRefinedData(double****, const char*, const char*, const char*);
	void printRefinedDataHDF5(double****, const char*);
	void deallocateFinest(double****&);
	void writeOutFileUniform(double****&);
	void generalWriteUniform();
	void writeVariablesUniform(double****&);
	void writeCoordinatesUniform();
	void writeBoundsUniform();
	void writeNodeTypeUniform();
	void writeRefineLevelUniform();
	void writeGidUniform();
	void writeVarDataUniform(const char[], double****&);
	void writeVariablesUniform();
	void writeIntsUniform(const char[]);
	void writeRealsUniform(const char[]);
	void writeSimParametersUniform();
	void writeWhichChildUniform();
	void writeBFlagsUniform();
	void writeSizesUniform();
	int* getTotalCellsUniform();


	void setVariable(int , std::string , int , int , int , double );
};


#endif
