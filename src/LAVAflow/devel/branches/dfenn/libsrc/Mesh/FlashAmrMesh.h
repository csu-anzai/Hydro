#include <hdf5.h>
#include <mpi.h>
#include "Mesh.h"

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
#define FORCEDEALLOCATECRASH false
#define USEGUARDCELLS false
#define NUMBEROFGUARDCELLS 0

#define SUPPRESSOUTPUT true

enum BLOCKTYPE      { ALL_BLKS, IBDRY_BLKS, JBDRY_BLKS, KBDRY_BLKS, ANY_BDRY_BLKS, 
                      ACTIVE_BLKS, LEAF, PARENT_BLK, ANCESTOR, REFINEMENT, INREGION};


typedef struct block {
	int blockNumber;			// Local processor's block number. Included for completeness...
	int cpuNumber;              // CPU number on which the block resided during simulation. NOT NECESSARILY 
								// THE MPI PROCESSOR CURRENTLY ACCESSING THE BLOCK.

	int nodeType;				// ==1 means this is a leaf block
	int refineLevel;

	double size[MESH_MDIM];		// Physical size of block
	double coord[MESH_MDIM];		// Center of block
	double lowerBound[MESH_MDIM]; /* coord - size/2 */
	double upperBound[MESH_MDIM]; /* coord + size/2 */

	double ****data;

//	double *x;   					  //Cell coordinates
//	double *y;
//	double *z;
//	double dx[MESH_MDIM];				  // Cell sizes

//	double vol;					
//	double ****data;
//	double ****udata;

	int parent;
	int *children;
	int *neighbor;       // <x, >x, <y, >y, <z, >z 
//	block *allNeighbors[MESH_MDIM][MESH_MDIM][MESH_MDIM]; //All neighbors

} block;


class FlashAmrMesh:Mesh {

private:
	block  *blocks;

	int **gid;
	int *procGidMap;

	// std::vector<std::vector<int> >    gidVector; 					/* Topological data of grid */

	int    nz[MESH_MDIM+1];     		/* Dimensions + nvar ( auxiliary ) */
	int    nzGuardLow [MESH_MDIM+1];  /* Dimensions with guard cells + nvar (auxiliary) */
	int    nzGuardHigh[MESH_MDIM+1];
	int    nguard;  				/* Number of guard cells in each direction */

	int    nFileVars; 				/* Number of field variables found in file */
	int    nUserVars;	 			/* Number of field variables requested by the user */
	char   k2d; 					/* Set to 1 if at least 2d is present */
	char   k3d; 					/* Set to 1 if at least 3d is present */

	int    nCellsVec[MESH_MDIM];		/* Number of cells in each direction */
	int    nCells;					/* Number of cells in a block: Product of nCellsVec's components. */

  	int    iteration;          		/* Current cycle (iteration) */

	int    format;					/* Format of the file (Will just be HDF5) */
	int    nBlocks;
	
	hid_t file;
	// H5File *file;

	char   variableNames[FR_MAXVARS][MESH_VAR_STRING_SIZE+1]; /* Names of stored variables */
	//char   userVariableNames[FR_MAXVARS][MESH_VAR_STRING_SIZE+1];
	//char   fileName[MAX_FILENAME_LENGTH];

	//Use map for variableNames. 

	// H5std_string fileName; //Done
	//H5File gridFile; //Can't get this working correctly. For now each method opens and 
	                   //closes a file.
	int fileType;

	int nBlocksAssigned;
	MPI::Intracomm comm;

	/* Variables used only in the uniform refinement function refineToFinest() */
	int    fineCells[MESH_MDIM];
	int    validCells[MESH_MDIM];
	int    totalCells[MESH_MDIM];
	int    numProcs[MESH_MDIM];
	double gridDelta[MESH_MDIM];
	double refinedDomainBoundBox[2][MESH_MDIM];
	int    **procIndexMin;

public:

	FlashAmrMesh();
	FlashAmrMesh(std::string , std::vector<std::string> );
	FlashAmrMesh(MPI::Intracomm);
	~FlashAmrMesh();

	//These two should be removed once I figure out how to initialize the file correctly
	//in the constructor.
	void setGridFileName(std::string);
	std::string getGridFileName(void);

	// void    readFile(std::string);

	void    generalRead();

	void    readInts(const char*);
	void    readReals(const char*);
	void    readStrings(const char*);
	void    readSimParameters(void);

	int     getInt(const char*);
	double  getReal(const char*);
	int     getFileType();

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

	void    readNFileVarsList(void);
	int     getNFileVars(void);
	int     getNUserVars(void);

	void    allocateVarArray(void);
	void    deallocateVarArray(void);
	int     findVarIndex(const char*);
	std::string  findVarName(const int);
	void    readVarData(const int, const char*);

	void    readVariables(void);

	void    addUserVariable(const char*);
	void    addAllVariables(void);
	double  getVariable(int, int, int, int, int);

	void    deallocateVariables(void);

	void 	openFile(void);
	void    closeFile(void);

	//MPI Functions
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

	void getCellSideLengths(int , double[MESH_MDIM] );
	void refineToFinest(double****&, double[2][MESH_MDIM], int=-1);
	void printRefinedData(double****, const char*, const char*, const char*);
	void printRefinedDataHDF5(double****, const char*);
	void deallocateFinest(double****&);

	
	
	void setVariable(int , std::string , int , int , int , double );
};


#endif
