 // #include "H5Cpp.h"
//#include <hdf5.h>
#include <mpi.h>

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

#ifndef LAVA_MESH_H
#define LAVA_MESH_H

#define MESH_MDIM 3
#define MESH_VAR_STRING_SIZE 4
#define MAX_STRING_LENGTH 80
#define MAX_FILENAME_LENGTH 1024


enum GEOMETRY       { UNKNOWN=-1, CARTESIAN, POLAR, CYLINDRICAL, SPHERICAL };
enum AXIS           { IAXIS, JAXIS, KAXIS };
enum BOUND          { LOWER, UPPER };
enum GRIDDATASTRUCT { CENTER, FACEX, FACEY, FACEZ, WORK, SCRATCH, SCRATCH_CTR, SCRATCH_FACEX, 
                      SCRATCH_FACEY, SCRATCH_FACEZ, CELL_VOLUME, CELL_FACEAREA };
enum BEGINCOUNT     { INTERIOR, EXTERIOR };

      
      //  CENTER cell centered variables
      //  FACEX  face centered variable on faces along IAXIS
      //  FACEY  face centered variable on faces along JAXIS
      //  FACEZ  face centered variable on faces along IAXIS
      //  WORK   single, cell centered variable, valid only for paramesh.
      //  SCRATCH scratch space that can fit cell and face centered variables
      //  SCRATCH_CTR scratch space for cell centered variables
      //  SCRATCH_FACEX scratch space facex variables
      //  SCRATCH_FACEY scratch space facey variables
      //  SCRATCH_FACEZ scratch space facez variables
      //  CELL_VOLUME volumes of specified cells 
      //  CELL_FACEAREA face area for the specified cells on specified face 

      // ALL_BLKS    all local blocks on a processor

      // IBDRY_BLKS  blocks that are on physical boundary along IAXIS
      // JBDRY_BLKS  blocks that are on physical boundary along JAXIS
      // KBDRY_BLKS  blocks that are on physical boundary along KAXIS
      // ANY_BDRY_BLKS  blocks that have any of their faces on a physical
      //                 boundary.
      // ACTIVE_BLKS all currently active blocks, in paramesh
      // context that means parent and leaf blocks
              
      // values that have meaning only for paramesh are :
      // LEAF, PARENT_BLK or ANCESTOR  representing
      // the type of node in the Oct-tree managing the blocks.
      // REFINEMENT the refinement level
      // INREGION    All blocks within the region defined by
      //             the accompanying optional argument
      //             region_bndBox. If the optional argument
      //             refinementLevel is also present then the blocks
      //             are also checked to see if they are at the specified
      //             refinement level, and only those that are get included
      //             in the list

class Mesh {

protected:
	int    dim;
  	double time; 	
	int    geometryType;

	char   runComment[MAX_STRING_LENGTH];
	char   geometryName[MAX_STRING_LENGTH];
	
	//Use map for variableNames. 

	std::string   fileName;
	//H5File gridFile; //Can't get this working correctly. For now each method opens and 
	                   //closes a file.
	int fileType;

	bool forceDeallocateCrash;

	std::map<int,std::string>     		uVarsMap_iToS;
	std::map<std::string,int>			uVarsMap_sToI;
	std::map<std::string,int>			fileVarsMap_sToI;
	std::map<int,std::string>			fileVarsMap_iToS;

	std::map<std::string,std::string>   stringScalarsMap;
	std::map<std::string,int>     		intScalarsMap;
	std::map<std::string,double>  		realScalarsMap;

public:

	Mesh();
	~Mesh();

	void setGridFileName(std::string);
	std::string getGridFileName(void);
};


#endif
