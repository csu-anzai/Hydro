/* C++ code used to read from FLASH files.
   Written by Philip Boehner
   July 3rd, 2013.

   Main function should look something like this:

#include "../include/flash_reader.h"
#include "../include/outside_functions.h"
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
	using namespace H5;
#endif

const H5std_string FILE_NAME( "snb_hdf5_chk_0000" );

int main() {

	FlashAmrMesh newGrid;
	
	newGrid.setGridFileName( FILE_NAME );

	newGrid.generalRead();

		//Read all variables
	//newGrid.readAllVariables();

		//Read some variables
	//newGrid.addUserVariable("dens");
	//newGrid.addUserVariable("pres");
	//newGrid.readPartialVariables();


		//User code goes here


	newGrid.deallocateVariables();

	return 0;

}

*/

/*	TODO:
		  [done] Suppress output in general read function
		  Check for memory leaks in refine (it was on the board)

		  MPI:

		  	[DONE] Basic MPI stuff (init, exit)
		  	[DONE] Broadcast GID information to all processors
		  	[DONE]	-Each processor can just have a copy of the Grid. Maybe.
		  	[done] Make GID function(s) for connectivity.
		  		-Populated children and neighbor arrays for each block, parent int.
		  	Split up blocks evenly across all processors. This will
		  	 cause problems for connectivity, but this is just the
		  	 first iteration.
		  	Keep only one grid.
		  	Maybe make a secondary gid array that keeps track of where 
			  	 everything is located in processor space.
			HDF reader has support for MPI stuff. Use it when reading in.
				There's a serial and parallel version for reading things in.
				Serial means everything is piped through processor 0 to begin with, I think.
				Parallel reads things in in parallel.
			Extra block information (nodeType, refineLevel, etc) needs to be read in
			on each block. Take out of general read function.
	
			Parallel:
				Fix reading cpu numbers
				Fix reading block numbers
				Fix reading node type
				Fix reading refineLevel -- Fix refineToFinest to look for all refine levels.
				Fix reading size, coord, lower, upper bounds.
				Fix variable reads.
				Fix refineToFinest 
*/

#include "FlashAmrMesh.h"

using namespace std;
	
FlashAmrMesh::FlashAmrMesh() {
	/* If using this constructor, these functions must be called in Main in this order:
	      setGridFileName(name of input HDF5 file as char array)
	      (optional) setGridOutFileName (name of ourput HDF5 file as char array) 
          openInFile()
          generalRead()
          addVariableNames(mesh variables in file as string array, auxiliary variables as string array)
          allocateVarArray()
		  readVariables()

		  NOTE: Writing of variables is always done using writeOutFile() in main.
	*/

	nBlocks = 0;
	k2d = 0;
	k3d = 0;
	format = 1;
	nFileVars = 0;
	nUserVars = 0;
	nReqFileVars = 0;
	nReqAuxVars = 0;
	blocks = NULL;
	gid = NULL;
	inFile = NULL;
	outFile = NULL;
	variablesAllocated = false;
	comm = MPI::COMM_WORLD;


}

FlashAmrMesh::FlashAmrMesh(string inFileStr, 
						   vector<string> fileVars) {
	/* Constructor with no support for writing files or creating new derived variables. */

	nBlocks = 0;
	k2d = 0;
	k3d = 0;
	format = 1;
	nFileVars = 0;
	nUserVars = 0;
	nReqFileVars = 0;
	nReqAuxVars = 0;
	blocks = NULL;
	gid = NULL;
	inFile = NULL;
	outFile = NULL;
	variablesAllocated = false;
	comm = MPI::COMM_WORLD;
	vector<string> auxVars;

	setGridFileName(inFileStr);
	readInFile(fileVars, auxVars);
}

FlashAmrMesh::FlashAmrMesh(string inFileStr, 
						   vector<string> fileVars, 
						   vector<string> auxVars) {
	/* Constructor with capabilities for creating new derived variables. */

	nBlocks = 0;
	k2d = 0;
	k3d = 0;
	format = 1;
	nFileVars = 0;
	nUserVars = 0;
	nReqFileVars = 0;
	nReqAuxVars = 0;
	blocks = NULL;
	gid = NULL;
	inFile = NULL;
	outFile = NULL;
	variablesAllocated = false;
	comm = MPI::COMM_WORLD;

	setGridFileName(inFileStr);
	readInFile(fileVars, auxVars);
}

FlashAmrMesh::FlashAmrMesh(string inFileStr, 
	                       string outFileStr, 
	                       vector<string> fileVars) {
	/* Constructor allowing for writing variables but no creation of new derived variables.
	   NOTE: Writing of variables is always done using writeOutFile() in main. 
	*/

	nBlocks = 0;
	k2d = 0;
	k3d = 0;
	format = 1;
	nFileVars = 0;
	nUserVars = 0;
	nReqFileVars = 0;
	nReqAuxVars = 0;
	blocks = NULL;
	gid = NULL;
	inFile = NULL;
	outFile = NULL;
	variablesAllocated = false;
	comm = MPI::COMM_WORLD;
	vector<string> auxVars;

	setGridFileName(inFileStr);
	setGridOutFileName(outFileStr);
	readInFile(fileVars, auxVars);
}

FlashAmrMesh::FlashAmrMesh(string inFileStr, 
	                       string outFileStr, 
	                       vector<string> fileVars, 
	                       vector<string> auxVars) {
	/* Constructor allowing for writing variables and adding new derived variables.
	   NOTE: Writing of variables is always done using writeOutFile() in main. 
	*/

	nBlocks = 0;
	k2d = 0;
	k3d = 0;
	format = 1;
	nFileVars = 0;
	nUserVars = 0;
	nReqFileVars = 0;
	nReqAuxVars = 0;
	blocks = NULL;
	gid = NULL;
	inFile = NULL;
	outFile = NULL;
	variablesAllocated = false;
	comm = MPI::COMM_WORLD;

	setGridFileName(inFileStr);
	setGridOutFileName(outFileStr);
	readInFile(fileVars, auxVars);
}

FlashAmrMesh::~FlashAmrMesh() {
	deallocateVariables();
	closeInFile();
}

// void FlashAmrMesh::LoadMesh(std::string inFileStr, 
// 	                        std::vector<std::string> fileVars)
// {
// 	setGridFileName(inFileStr);
// 	openInFile();
// 	generalRead();
//  vector<string> auxVars;
//  addVariableNames(fileVars, auxVars);
// 	allocateVarArray();
// 	readVariables();
// }

void FlashAmrMesh::setGridFileName(std::string inFileStr){
	inFileName = inFileStr;
}

void FlashAmrMesh::setGridOutFileName(std::string outFileStr){
	outFileName = outFileStr;
}

std::string FlashAmrMesh::getGridInFileName() {
	return inFileName;
}

std::string FlashAmrMesh::getGridOutFileName() {
	return outFileName;
}

void FlashAmrMesh::readInts(const char *which){
	//This routine returns an integer value corresponding to the 
	// requested scalar char array argument from either the
	// integer scalar list or the integer runtime parameters list. 
	// These values go into the map object intScalarsMap.

	typedef struct intScalarStruct {
		char name[MAX_STRING_LENGTH+1];
		int  value;
	} intScalarStruct;

	char datasetName[50] = {};
	if (string(which) == "scalars") {
		strcpy(datasetName,"integer scalars");
	}
	else if (string(which) == "parameters") {
		strcpy(datasetName,"integer runtime parameters");
	}

	if (!datasetAvailable(datasetName,"required")) {
		return;
	}

	//Get scalars list.
	hid_t dset = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	
	//Get dimensions of dataset.
	hsize_t dims_out[1];
	hid_t space = H5Dget_space( dset );
	// hid_t ndims = H5Sget_simple_extent_dims( space, dims_out, NULL );
	H5Sget_simple_extent_dims( space, dims_out, NULL );

	int sizeOfScalarArray = (int)dims_out[0];

	if (!SUPPRESSOUTPUT) cout << "Size of " << datasetName << " array: " << sizeOfScalarArray << endl;

	//Make specific type for the 80 character string.
	hid_t strtype = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype, MAX_STRING_LENGTH+1);

	//"integer scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(intScalarStruct) );
	H5Tinsert( memtype, "name", HOFFSET(intScalarStruct, name), strtype);
	H5Tinsert( memtype, "value", HOFFSET(intScalarStruct, value), H5T_NATIVE_INT);

	intScalarStruct* out = new intScalarStruct[sizeOfScalarArray];
	H5Dread( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);

	//Make array of strings for easy comparisons
	string* strArray = new string[sizeOfScalarArray];

	//Erase white space after name.
	string whiteSpaces(" \t\f\v\n\r");
	for (int i = 0; i < sizeOfScalarArray; i++) {
		strArray[i] = string(out[i].name);
		size_t firstWhiteSpace = strArray[i].find_last_not_of(whiteSpaces);
		strArray[i].erase(firstWhiteSpace+1);
	}

	for (int i = 0; i < sizeOfScalarArray; i++) {
		intScalarsMap.insert(pair<string,int>(strArray[i],out[i].value));
	}

	H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, out);
	H5Dclose(dset);
	H5Sclose(space);
	H5Tclose(memtype);
	H5Tclose(strtype);

	delete[] out;
	delete[] strArray;

}

void FlashAmrMesh::readLogicals(const char *which) {
	/* Similar to the readInts function, but not in every flash file.
	   Suppress the error output temporarily. We want to see if the locals list even exists.
	   If it does not, exit out of this function. If it does, read it.
	*/
	
	typedef struct logicalScalarStruct {
		char name[MAX_STRING_LENGTH+1];
		int  value;
	} logicalScalarStruct;

	char datasetName[50] = {};
	if (string(which) == "scalars") {
		strcpy(datasetName,"logical scalars");
	}
	else if (string(which) == "parameters") {
		strcpy(datasetName,"logical runtime parameters");
	}

	if (!datasetAvailable(datasetName,"optional")) {
		return;
	}
	
	//Get scalars list.
	hid_t dset = H5Dopen(inFile, datasetName, H5P_DEFAULT);

	//Get dimensions of dataset.
	hsize_t dims_out[1];
	hid_t space = H5Dget_space( dset );
	// hid_t ndims = H5Sget_simple_extent_dims( space, dims_out, NULL );
	H5Sget_simple_extent_dims( space, dims_out, NULL );

	int sizeOfScalarArray = (int)dims_out[0];

	if (!SUPPRESSOUTPUT) cout << "Size of " << datasetName << " array: " << sizeOfScalarArray << endl;

	//Make specific type for the 80 character string.
	hid_t strtype = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype, MAX_STRING_LENGTH+1);

	//"integer scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(logicalScalarStruct) );
	H5Tinsert( memtype, "name", HOFFSET(logicalScalarStruct, name), strtype);
	H5Tinsert( memtype, "value", HOFFSET(logicalScalarStruct, value), H5T_NATIVE_INT);

	logicalScalarStruct* out = new logicalScalarStruct[sizeOfScalarArray];
	H5Dread( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);

	//Make array of strings for easy comparisons
	string* strArray = new string[sizeOfScalarArray];

	//Erase white space after name.
	string whiteSpaces(" \t\f\v\n\r");
	for (int i = 0; i < sizeOfScalarArray; i++) {
		strArray[i] = string(out[i].name);
		size_t firstWhiteSpace = strArray[i].find_last_not_of(whiteSpaces);
		strArray[i].erase(firstWhiteSpace+1);
	}

	for (int i = 0; i < sizeOfScalarArray; i++) {
		logicalScalarsMap.insert(pair<string,int>(strArray[i],out[i].value));
	}
    
	H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, out);
	H5Dclose(dset);
	H5Sclose(space);
	H5Tclose(memtype);
	H5Tclose(strtype);

	delete[] out;
	delete[] strArray;

}

void FlashAmrMesh::readSimParameters() {
	//For use with file type 7 FLASH files only.

	typedef struct simParamStruct {
		//Instead of this file being formatted like the rest where there is a name
		// column and a value column, flie format 7 has 8 columns with no names.
		// The values have to be read directly into the variables here.
		int totalBlocks;
		float time;
		float timeStep;
		float redShift;
		int numberOfSteps;
		int nxb;
		int nyb;
		int nzb;
	} realParamStruct;

	if (!datasetAvailable("simulation parameters","required")) {
		return;
	}

	//Get scalars list.
	hid_t dset = H5Dopen(inFile, "simulation parameters", H5P_DEFAULT);
	
	//Get dimensions of dataset.
	hsize_t dims_out[1];
	hid_t space = H5Dget_space( dset );
	// hid_t ndims = H5Sget_simple_extent_dims( space, dims_out, NULL );
	H5Sget_simple_extent_dims( space, dims_out, NULL );

	int sizeOfScalarArray = (int)dims_out[0];

	if (!SUPPRESSOUTPUT) cout << "Size of simulation parameter array: " << sizeOfScalarArray << endl;

	//"simulation parameters" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(simParamStruct));

	//NOTE: IF THERE ARE ANY EXTRA PARAMETERS HERE, THEY WILL NOT BE READ.
	//      THIS IS DUE TO THE WEIRD WAY THAT THE PARAMETERS WERE WRITTEN TO THE HDF5 FILE.
	H5Tinsert( memtype, "total blocks", 	HOFFSET(simParamStruct, totalBlocks), 	H5T_NATIVE_INT);
	H5Tinsert( memtype, "time", 			HOFFSET(simParamStruct, time), 			H5T_NATIVE_FLOAT);
	H5Tinsert( memtype, "timestep", 		HOFFSET(simParamStruct, timeStep), 		H5T_NATIVE_FLOAT);
	H5Tinsert( memtype, "redshift", 		HOFFSET(simParamStruct, redShift), 		H5T_NATIVE_FLOAT);
	H5Tinsert( memtype, "number of steps", 	HOFFSET(simParamStruct, numberOfSteps), H5T_NATIVE_INT);
	H5Tinsert( memtype, "nxb", 				HOFFSET(simParamStruct, nxb), 			H5T_NATIVE_INT);
	H5Tinsert( memtype, "nyb", 				HOFFSET(simParamStruct, nyb), 			H5T_NATIVE_INT);
	H5Tinsert( memtype, "nzb", 				HOFFSET(simParamStruct, nzb), 			H5T_NATIVE_INT);
	
	simParamStruct* out = new simParamStruct[sizeOfScalarArray];
	H5Dread( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);

	intScalarsMap.insert(pair<std::string,int>    ("total blocks",out[0].totalBlocks));
	realScalarsMap.insert(pair<std::string,double>("time", out[0].time));
	realScalarsMap.insert(pair<std::string,double>("timestep", out[0].timeStep));
	realScalarsMap.insert(pair<std::string,double>("redshift", out[0].redShift));
	intScalarsMap.insert(pair<std::string,int>    ("number of steps",out[0].numberOfSteps));
	intScalarsMap.insert(pair<std::string,int>    ("nxb",out[0].nxb));
	intScalarsMap.insert(pair<std::string,int>    ("nyb",out[0].nyb));
	intScalarsMap.insert(pair<std::string,int>    ("nzb",out[0].nzb));

	H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, out);
	H5Dclose(dset);
	H5Sclose(space);
	H5Tclose(memtype);

	delete[] out;
}

void FlashAmrMesh::readReals(const char *which){
	//This routine returns a double corresponding to the 
	// requested scalar char array argument from either the 
	// real scalar list or the real runtime parameters list. 
	// These values go into the map object realScalarsMap.

	typedef struct realScalarStruct {
		char 	name[MAX_STRING_LENGTH+1];
		double  value;
	} realScalarStruct;


	char datasetName[50] = {};
	if (string(which) == "scalars") {
		strcpy(datasetName,"real scalars");
	}
	else if (string(which) == "parameters") {
		strcpy(datasetName,"real runtime parameters");
	}

	if (!datasetAvailable(datasetName,"required")) {
		return;
	}

	//Get scalars list.
	hid_t dset = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	
	//Get dimensions of dataset.
	hsize_t dims_out[1];
	hid_t space = H5Dget_space( dset );
	// hid_t ndims = H5Sget_simple_extent_dims( space, dims_out, NULL );
	H5Sget_simple_extent_dims( space, dims_out, NULL );

	int sizeOfScalarArray = (int)dims_out[0];

	if (!SUPPRESSOUTPUT) cout << "Size of " << datasetName << " array: " << sizeOfScalarArray << endl;

	//Make specific type for the 80 character string.
	hid_t strtype = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype, MAX_STRING_LENGTH+1);

	//"integer scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(realScalarStruct) );
	H5Tinsert( memtype, "name", HOFFSET(realScalarStruct, name), strtype);
	H5Tinsert( memtype, "value", HOFFSET(realScalarStruct, value), H5T_NATIVE_DOUBLE);

	realScalarStruct* out = new realScalarStruct[sizeOfScalarArray];
	H5Dread( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);

	//Make array of strings for easy comparisons
	string* strArray = new string[sizeOfScalarArray];

	//Erase white space after name.
	string whiteSpaces(" \t\f\v\n\r");
	for (int i = 0; i < sizeOfScalarArray; i++) {
		strArray[i] = string(out[i].name);
		size_t firstWhiteSpace = strArray[i].find_last_not_of(whiteSpaces);
		strArray[i].erase(firstWhiteSpace+1);
	}

	for (int i = 0; i < sizeOfScalarArray; i++) {
		realScalarsMap.insert(pair<std::string,double>(strArray[i],out[i].value));
	}

	H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, out);
	H5Dclose(dset);
	H5Sclose(space);
	H5Tclose(memtype);
	H5Tclose(strtype);

	delete[] out;
	delete[] strArray;
}

void FlashAmrMesh::readStrings(const char *which){
	//This routine returns a string corresponding to the 
	// requested scalar char array argument from either the
	// string scalar list or the string runtime parameters list. 
	// These values go into the map object stringScalarsMap.

	typedef struct stringScalarStruct {
		char name[MAX_STRING_LENGTH+1];
		char value[MAX_STRING_LENGTH+1];
	} stringScalarStruct;

	char datasetName[50] = {};
	if (string(which) == "scalars") {
		strcpy(datasetName,"string scalars");
	}
	else if (string(which) == "parameters") {
		strcpy(datasetName,"string runtime parameters");
	}

	if (!datasetAvailable(datasetName,"required")) {
		return;
	}

	//Get scalars list.
	hid_t dset = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	
	//Get dimensions of dataset.
	hsize_t dims_out[1];
	hid_t space = H5Dget_space( dset );
	// hid_t ndims = H5Sget_simple_extent_dims( space, dims_out, NULL );
	H5Sget_simple_extent_dims( space, dims_out, NULL );

	int sizeOfScalarArray = (int)dims_out[0];

	if (!SUPPRESSOUTPUT) cout << "Size of string scalar array: " << sizeOfScalarArray << endl;

	//Make specific type for the 80 character string.
	hid_t strtype = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype, MAX_STRING_LENGTH+1);
	hid_t strval = H5Tcopy(H5T_C_S1);
	H5Tset_size (strval, MAX_STRING_LENGTH+1);

	//"integer scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(stringScalarStruct) );
	H5Tinsert( memtype, "name", HOFFSET(stringScalarStruct, name), strtype);
	H5Tinsert( memtype, "value", HOFFSET(stringScalarStruct, value), strval);

	stringScalarStruct* out = new stringScalarStruct[sizeOfScalarArray];
	H5Dread( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);

	//Make array of strings for easy comparisons
	string* nameArray = new string[sizeOfScalarArray];
	string* valueArray = new string[sizeOfScalarArray];

	//Erase white space after name and value.
	string whiteSpaces(" \t\f\v\n\r");
	for (int i = 0; i < sizeOfScalarArray; i++) {
		nameArray[i] = string(out[i].name);
		size_t firstWhiteSpace = nameArray[i].find_last_not_of(whiteSpaces);
		nameArray[i].erase(firstWhiteSpace+1);

		valueArray[i] = string(out[i].value);
		firstWhiteSpace = valueArray[i].find_last_not_of(whiteSpaces);
		valueArray[i].erase(firstWhiteSpace+1);
	}

	for (int i = 0; i < sizeOfScalarArray; i++) {
		stringScalarsMap.insert(pair<std::string,std::string>(nameArray[i],valueArray[i]));
	}

	H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, out);
	H5Dclose(dset);
	H5Sclose(space);
	H5Tclose(memtype);
	H5Tclose(strtype);
	H5Tclose(strval);

	delete[] out;
	delete[] nameArray;
	delete[] valueArray;
}

int FlashAmrMesh::getInt(const char *whichValue) {
	return intScalarsMap[whichValue];
}

double FlashAmrMesh::getReal(const char *whichValue) {
	return realScalarsMap[whichValue];
}

int FlashAmrMesh::getFileType() {
	return fileType;
}

void FlashAmrMesh::readDim() {
	//Reads the dimensionality of the grid and sets the
	// dim variables.
	if (fileType != 7) {
		dim = intScalarsMap["dimensionality"];
		if (dim == 2) {
			k2d = 1;
		}
		if (dim == 3) {
			k2d = 1;
			k3d = 1;
		}
	
	}
}

int FlashAmrMesh::getDim() {
	return dim; 
}

void FlashAmrMesh::readNBlocks() {
	if (fileType == 7) {
		nBlocks = intScalarsMap["total blocks"];		
	}
	else {
		nBlocks = intScalarsMap["globalnumblocks"];
	}
}

int FlashAmrMesh::getNBlocks() {
	return nBlocks;
}

void FlashAmrMesh::allocateBlocks() {
	//Allocates the blocks array in the grid. Can only be called
	// after the setNBlocks function has been called.

	if (nBlocks == 0) {
		cout << "[flash_reader] Must read in number of blocks before allocating blocks array." << endl;
		cout << "[flash_reader] Use function setNBlocks()." << endl;
	}
	else {
		try
		{
			//MPI test
			blocks = new block[nBlocksAssigned];
  			//blocks = new block [nBlocks];
		}
		catch (bad_alloc&)
		{
  			cout << "[flash_reader] Problem with allocation in allocateBlocks function." << endl;
			exit(EXIT_FAILURE);
		}
	}
}

void FlashAmrMesh::deallocateBlocks() {

	// for (int i = 0; i < nBlocksAssigned; i++) {
	// 	delete[] blocks[i].children;
	// 	delete[] blocks[i].neighbor;
	// }

	if (!SUPPRESSOUTPUT) std::cout << "Processor " << comm.Get_rank() << " blocks: " << blocks << std::endl;
	if (blocks != NULL) {
		delete[] blocks;
		blocks = NULL;
	}
	else {
		if (FORCEDEALLOCATECRASH) {
			cout << "[flash_reader] Failure in deallocateBlocks function." << endl;
			exit(EXIT_FAILURE);	
		}
	}
}

void FlashAmrMesh::allocateGid() {
	//Allocates the gid array in the grid. Can only be called 
	// after the setNBlocks function has been called.

	if (nBlocks == 0) {
		cout << "[flash_reader] Must read in number of blocks before allocating gid array." << endl;
		cout << "[flash_reader] Use function setNBlocks()." << endl;
	}
	else {
		try
		{
  			//gid = new int [nBlocks][MESH_MDIM*2+1+(int)pow(2,MESH_MDIM)];

  			int nGidParams = MESH_MDIM*2 + 1 + (int)pow(2,MESH_MDIM);

			// gid for general codes
			gid = new int*[nBlocks];
			for (int i = 0; i < nBlocks; i++) {
								//neighbors + parent + children
				gid[i] = new int[nGidParams];
			}

			int maxNumChildren = (int)pow(2,MESH_MDIM);
			int maxNeighbors = 2*MESH_MDIM;

			for (int i = 0; i < nBlocksAssigned; i++) {
				blocks[i].children = new int[maxNumChildren];
				blocks[i].neighbor = new int[maxNeighbors];
			}

  			// // gidVector for Dan's code
  			// int GID_NPARAMS = MESH_MDIM*2 + 1 + (int)pow(2,MESH_MDIM);
  			// cout << "GID_NPARAMS is " << GID_NPARAMS << endl;
  			
  			// gidVector.resize(nBlocks);
  			// for (auto& i : gidVector)
  			// 	i.reserve(GID_NPARAMS);
		}
		catch (bad_alloc&)
		{
  			cout << "[flash_reader] Problem with allocation in allocateGid function." << endl;
			exit(EXIT_FAILURE);
		}
	}	
}

void FlashAmrMesh::assignGidArrays() {

	for (int i = 0; i < nBlocksAssigned; i++) {
		int globalBlockI = mpiGetGlobalBlockID(i);
		blocks[i].parent = gid[globalBlockI][getDim()*2];
		for (int j = 0; j < getDim()*2; j++) {
			blocks[i].neighbor[j] = gid[globalBlockI][j];
		}
		for (int j = 0; j < (int)pow(2,getDim()); j++) {
			blocks[i].children[j] = gid[globalBlockI][getDim()*2+1 + j];
		}
	}

	// std::cout << "Dimensions: " << getDim() << std::endl;
	// std::cout << "Assigning Block 1" << std::endl;
	// for (int i = 0; i < getDim()*2; i++) {
	// 	std::cout << blocks[0].neighbor[i] << std::endl;
	// }
	// std::cout << blocks[0].parent << std::endl;
	// for (int i = 0; i < (int)pow(2,getDim()); i++) {
	// 	std::cout << blocks[0].children[i] << std::endl;
	// }

	// enum GIDDATA 		{ GID_NEIGHBOR_END = MESH_MDIM*2-1, GID_PARENT = GID_NEIGHBOR_END+1, 
	// 						GID_CHILDREN_BEGIN = GID_PARENT + 1, GID_CHILDREN_END = GID_CHILDREN_BEGIN+(int)pow(2,MESH_MDIM),
	// 						GID_NPARAMS = GID_CHILDREN_END};

}

void FlashAmrMesh::deallocateGid() {
	if (gid != NULL) {
		for (int i = 0; i < nBlocks; i++) {
			delete [] gid[i];
		}
		delete[] gid;

		for (int i = 0; i < nBlocksAssigned; i++) {
			delete[] blocks[i].children;
			delete[] blocks[i].neighbor;

			blocks[i].children = NULL;
			blocks[i].neighbor = NULL;
		}

		gid = NULL;

	}
	else {
		if (FORCEDEALLOCATECRASH) {
			cout << "[flash_reader] Failure in deallocateGid function." << endl;
			exit(EXIT_FAILURE);	
		}
	}
}

void FlashAmrMesh::readBlockNumber() {
	/* Local block number. Mostly useless. */
	for (int i = 0; i < nBlocksAssigned; i++) {
		blocks[i].blockNumber = i;
	}

	// std::cout << "BLOCK NUMBER TEST PROC: " << comm.Get_rank() << " BLOCK: " << mpiGetGlobalBlockID(nBlocksAssigned-1) << std::endl;
	// std::cout << "Local block number: " << nBlocksAssigned-1 << " Assigned block number: " << blocks[nBlocksAssigned-1].blockNumber << std::endl;
}

int FlashAmrMesh::getBlockNumber(block inBlock) {
	return inBlock.blockNumber;
}

void FlashAmrMesh::readNcells() {
	//Reads the size of the cells for each block. 
	// Note that this is constant across all blocks, so it
	// is set inside the grid rather than each block.

	nCellsVec[0] = intScalarsMap["nxb"];
	nCellsVec[1] = intScalarsMap["nyb"];
	nCellsVec[2] = intScalarsMap["nzb"];
	nCells = nCellsVec[0] * nCellsVec[1] * nCellsVec[2];

}

int FlashAmrMesh::getNcells(int axis) {
	return nCellsVec[axis];
}

int* FlashAmrMesh::getNcells() {
	return nCellsVec;
}

int FlashAmrMesh::getTotalCells() {
	return nCells;
}

void FlashAmrMesh::readTime() {
	time = realScalarsMap["time"];
}

double FlashAmrMesh::getTime() {
	return time;
}

void FlashAmrMesh::readGeometryType() {
	//Retrieves the geometry from the stringScalars map and
	// returns both the geometryType int value as well as 
	// the name of the geometry type in the geometryName array.

	string tempString = stringScalarsMap["geometry"];
	
	if (tempString == "cartesian") {
		geometryType = CARTESIAN;
		strcpy(geometryName,"cartesian");
	}
	else if (tempString == "polar") {
		geometryType = POLAR;
		strcpy(geometryName,"polar");
	}
	else if (tempString == "cylindrical") {
		geometryType = CYLINDRICAL;
		strcpy(geometryName,"cylindrical");
	}
	else if (tempString == "spherical") {
		geometryType = SPHERICAL;
		strcpy(geometryName,"spherical");
	}
	else {
		cout << "[flash_reader] Possible error: Unknown geometry type." << endl;
		geometryType = UNKNOWN;
		strcpy(geometryName,"spherical");
	}
}

int FlashAmrMesh::getGeometryType() {
	return geometryType;
}

void FlashAmrMesh::readCoordinates() {
	//Reads the coordinates data from the HDF5 file 
	// and places it in the coord array in each block.
	// These coordinates correspond to the centers of each block.
	//
	//Special functionality for file format 7 at the bottom of this function.

	if (!datasetAvailable("coordinates","required")) {
		return;
	}

	hid_t dset = H5Dopen(inFile, "coordinates", H5P_DEFAULT);

	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[2];

	int rank = H5Sget_simple_extent_dims( space, dims_out, NULL );

	//MPI stuff
	hsize_t offset[2]; //hyperslab offset in file
	hsize_t count[2];  //size of hyperslab
	offset[0] = mpiGetGlobalBlockID(0);
	offset[1] = 0;
	count[0]  = nBlocksAssigned;
	count[1]  = dims_out[1];

	H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

	hid_t memspace = H5Screate_simple(2, dims_out, NULL);

    hsize_t      offset_out[2];   // hyperslab offset in memory
    hsize_t      count_out[2];    // size of the hyperslab in memory
    offset_out[0] = 0;
    offset_out[1] = 0;
    count_out[0]  = nBlocksAssigned;
    count_out[1]  = dims_out[1];

	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);

	//Create temp 2d array.
	double *tempCoordsOut = new double[nBlocksAssigned*dims_out[1]];

	H5Dread( dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, tempCoordsOut);

	//Read temp data into corresponding block coord arrays.
	for (int i = 0; i < nBlocksAssigned; i++) {
		for (int j = 0; j < (int)dims_out[1]; j++) {
			blocks[i].coord[j] = tempCoordsOut[MESH_MDIM*i + j];
		}
	}	
	
	//This is only specific to file format 7 FLASH output files. 
	//There is no dim variable in the scalar/parameter lists (simulation, string, real, etc)
	//like there is for the other formats, so we're using coords to determine what dim is.
	if (fileType == 7) {
		if (rank == 3) {
			dim = 3;
			k2d = 1;
			k3d = 1;
		}
		else if (rank == 2) {
			dim = 2;
			k2d = 1;
			k3d = 0;
		}
		else {
			dim = 1;
			k2d = 0;
			k3d = 0;
		}
	}

	delete[] tempCoordsOut;

	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}

void FlashAmrMesh::readSizes() {
	//Reads the physical sizes of each block from the HDF5 file 
	// and places them in the size array in each block.

	// This dataset is not required for writing, but may be used in internal
	// functions to get the block sizes. Setting this to required, even though it's
	// optional in the corresponding writeSizes() function.
	if (!datasetAvailable("block size","required")) {
		return;
	}

	hid_t dset = H5Dopen(inFile, "block size", H5P_DEFAULT);

	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[2];

	H5Sget_simple_extent_dims( space, dims_out, NULL );

	//MPI stuff
	hsize_t offset[2]; //hyperslab offset in file
	hsize_t count[2];  //size of hyperslab
	offset[0] = mpiGetGlobalBlockID(0);
	offset[1] = 0;
	count[0]  = nBlocksAssigned;
	count[1]  = dims_out[1];

	H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

	hid_t memspace = H5Screate_simple(2, dims_out, NULL);

    hsize_t      offset_out[2];   // hyperslab offset in memory
    hsize_t      count_out[2];    // size of the hyperslab in memory
    offset_out[0] = 0;
    offset_out[1] = 0;
    count_out[0]  = nBlocksAssigned;
    count_out[1]  = dims_out[1];

	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
	
	//Create temp 2d array.
	double *tempSizesOutD = new double[nBlocksAssigned*dims_out[1]];
	float *tempSizesOutF = new float[nBlocksAssigned*dims_out[1]];

	if (fileType == 7) {
		H5Dread( dset, H5T_NATIVE_FLOAT, memspace, space, H5P_DEFAULT, tempSizesOutF);
		// sizesData.read( tempSizesOutF, PredType::IEEE_F32LE, memspace, sizesDataSpace );
		//Read temp data into corresponding block coord arrays.
		for (int i = 0; i < nBlocksAssigned; i++) {
			for (int j = 0; j < (int)dims_out[1]; j++) {
				blocks[i].size[j] = tempSizesOutF[(int)dims_out[1]*i + j];
			}	
		}
	}
	else {
		H5Dread( dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, tempSizesOutD);
		// sizesData.read( tempSizesOutD, PredType::IEEE_F64LE, memspace, sizesDataSpace );	
		
		//First set all physical sizes to 0. This way we don't need special exceptions for 1D or 2D spaces.
		for (int i = 0; i < nBlocksAssigned; i++) {
			for (int j = 0; j < MESH_MDIM; j++) {
				blocks[i].size[j] = 0;
			}
		}

		//Read temp data into corresponding block coord arrays.
		for (int i = 0; i < nBlocksAssigned; i++) {
			for (int j = 0; j < (int)dims_out[1]; j++) {
				blocks[i].size[j] = tempSizesOutD[MESH_MDIM*i + j];
			}
		}

	}

	// std::cout << "SIZES TEST PROC: " << comm.Get_rank() << " BLOCK: " << mpiGetGlobalBlockID(0) << std::endl;
	// std::cout << blocks[1].size[1] << std::endl;

	delete[] tempSizesOutD;
	delete[] tempSizesOutF;

	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}


void FlashAmrMesh::readBounds() {
	//Reads the lower and upper bounds for each block from the HDF5 file 
	// and places them in the lowerBound and upperBound arrays in each block.

	if (!datasetAvailable("bounding box","required")) {
		return;
	}

	hid_t dset = H5Dopen(inFile, "bounding box", H5P_DEFAULT);

	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[3];

	H5Sget_simple_extent_dims( space, dims_out, NULL );


	//MPI stuff
	hsize_t offset[3]; //hyperslab offset in file
	hsize_t count[3];  //size of hyperslab

	offset[0] = mpiGetGlobalBlockID(0);
	offset[1] = 0;
	offset[2] = 0;

	count[0]  = nBlocksAssigned;
	count[1]  = dims_out[1];
	count[2]  = dims_out[2];

	H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

	hid_t memspace = H5Screate_simple(3, dims_out, NULL);

    hsize_t      offset_out[3];   // hyperslab offset in memory
    hsize_t      count_out[3];    // size of the hyperslab in memory
    offset_out[0] = 0;
    offset_out[1] = 0;
    offset_out[2] = 0;
    count_out[0]  = nBlocksAssigned;
    count_out[1]  = dims_out[1];
    count_out[2]  = dims_out[2];

	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
	
	//Create temp 3d array.
	double *tempBoundsOut = new double[nBlocksAssigned*dims_out[1]*dims_out[2]];

	H5Dread( dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, tempBoundsOut);

	//Read temp data into corresponding block coord arrays.
	for (int i = 0; i < nBlocksAssigned; i++) {
		for (int j = 0; j < (int)dims_out[1]; j++) {
			blocks[i].lowerBound[j] = tempBoundsOut[(i*dims_out[1]+j)*dims_out[2] ];
			blocks[i].upperBound[j] = tempBoundsOut[(i*dims_out[1]+j)*dims_out[2] + 1];
		}
	}

	// std::cout << "LOWER AND UPPER BOUND TEST PROC: " << comm.Get_rank() << " BLOCK: " << mpiGetGlobalBlockID(5) << std::endl;
	// std::cout << "LOWER: " << blocks[5].lowerBound[1] << std::endl;
	// std::cout << "UPPER: " << blocks[5].upperBound[1] << std::endl;

	delete[] tempBoundsOut;
	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}


void FlashAmrMesh::readCPUNumbers(){
	//Reads which CPU each block is assigned to from the HDF5 file 
	// and places it in the cpuNumber variable in each block.

	if (!datasetAvailable("processor number","optional")) {	
		return;
	}

	hid_t dset = H5Dopen(inFile, "processor number", H5P_DEFAULT);

	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[1];

	H5Sget_simple_extent_dims( space, dims_out, NULL );

	//MPI stuff
	hsize_t offset[1]; //hyperslab offset in file
	hsize_t count[1];  //size of hyperslab

	offset[0] = mpiGetGlobalBlockID(0);
	count[0]  = nBlocksAssigned;

	H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

	hid_t memspace = H5Screate_simple(1, dims_out, NULL);

    hsize_t      offset_out[1];   // hyperslab offset in memory
    hsize_t      count_out[1];    // size of the hyperslab in memory
    offset_out[0] = 0;
    count_out[0]  = nBlocksAssigned;

	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
	
	//Create temp 1d array.
	int *tempProcessorOut = new int[nBlocksAssigned];

	H5Dread( dset, H5T_NATIVE_INT, memspace, space, H5P_DEFAULT, tempProcessorOut);

	//Read from temp into blocks.
	for (int i = 0; i < nBlocksAssigned; i++) {
		blocks[i].cpuNumber = tempProcessorOut[i];
	}

	// std::cout << "CPU READ TEST PROC: " << comm.Get_rank() << " BLOCK: " << mpiGetGlobalBlockID(0) << std::endl;
	// std::cout << blocks[0].cpuNumber << std::endl;

	delete[] tempProcessorOut;
	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}

int FlashAmrMesh::getCPUNumber(int blockIndex) {
	return blocks[blockIndex].cpuNumber;
}


void FlashAmrMesh::readNodeType(){
	//Reads the note type data from the HDF5 file 
	// and places it in the nodeType variable in each block.
	// These correspond to the location on the tree. For example, 
	// a nodeType equal to 1 would represent a    node.

	if (!datasetAvailable("node type","required")) {
		if (!SUPPRESSOUTPUT) std::cout << "[flash_reader] Warning: No 'node type' dataset found in the input file." << std::endl;
		return;
	}

	hid_t dset = H5Dopen(inFile, "node type", H5P_DEFAULT);

	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[1];

	H5Sget_simple_extent_dims( space, dims_out, NULL );

	//MPI stuff
	hsize_t offset[1]; //hyperslab offset in file
	hsize_t count[1];  //size of hyperslab

	offset[0] = mpiGetGlobalBlockID(0);
	count[0]  = nBlocksAssigned;

	H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

	hid_t memspace = H5Screate_simple(1, dims_out, NULL);

    hsize_t      offset_out[1];   // hyperslab offset in memory
    hsize_t      count_out[1];    // size of the hyperslab in memory
    offset_out[0] = 0;
    count_out[0]  = nBlocksAssigned;

	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
	
	//Create temp 1d array.
	int *tempNodeOut = new int[nBlocksAssigned];

	H5Dread( dset, H5T_NATIVE_INT, memspace, space, H5P_DEFAULT, tempNodeOut);

	//Read from temp into blocks.
	for (int i = 0; i < nBlocksAssigned; i++) {
		blocks[i].nodeType = tempNodeOut[i];
	}
	
	// std::cout << "NODE READ TEST PROC: " << comm.Get_rank() << " BLOCK: " << mpiGetGlobalBlockID(0) << std::endl;
	// std::cout << blocks[0].nodeType << std::endl;

	delete[] tempNodeOut;
	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}

int FlashAmrMesh::getNodeType(int blockIndex) {
	return blocks[blockIndex].nodeType;
}


void FlashAmrMesh::readRefineLevel(){
	//Reads the refinement level data from the HDF5 file 
	// and places it in the refineLevel array in each block.
	// The refinement level just keeps track of how refined the
	// block is. A low number represents a coarse refinement, and
	// a high number represents a finer refinement.

	if (!datasetAvailable("refine level","required")) {
		return;
	}

	hid_t dset = H5Dopen(inFile, "refine level", H5P_DEFAULT);

	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[1];

	H5Sget_simple_extent_dims( space, dims_out, NULL );

	//MPI stuff
	hsize_t offset[1]; //hyperslab offset in file
	hsize_t count[1];  //size of hyperslab

	offset[0] = mpiGetGlobalBlockID(0);
	count[0]  = nBlocksAssigned;

	H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

	hid_t memspace = H5Screate_simple(1, dims_out, NULL);

    hsize_t      offset_out[1];   // hyperslab offset in memory
    hsize_t      count_out[1];    // size of the hyperslab in memory
    offset_out[0] = 0;
    count_out[0]  = nBlocksAssigned;

	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
	
	//Create temp 1d array.
	int *tempRefineOut = new int[nBlocksAssigned];

	H5Dread( dset, H5T_NATIVE_INT, memspace, space, H5P_DEFAULT, tempRefineOut);

	//Read from temp into blocks.
	for (int i = 0; i < nBlocksAssigned; i++) {
		blocks[i].refineLevel = tempRefineOut[i];
	}

	delete[] tempRefineOut;
	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}

int FlashAmrMesh::getRefineLevel(int blockIndex) {
	return blocks[blockIndex].refineLevel;
}

void FlashAmrMesh::readGid() {
	//Reads the gid data from the HDF5 file 
	// and places it in the gid array in the grid.
	// Must be read after the gid array has been allocated.

	// Note: This read is not parallelized because every processor needs
	// to have the full gid array in memory.
	
	if (!datasetAvailable("gid","required")) {
		return;
	}

	hid_t dset = H5Dopen(inFile, "gid", H5P_DEFAULT);

	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[2];

	H5Sget_simple_extent_dims( space, dims_out, NULL );

	hid_t memspace = H5Screate_simple(2, dims_out, NULL);

	//Create temp 2d array.
	int *tempGidOut = new int[dims_out[0]*dims_out[1]];

	H5Dread( dset, H5T_NATIVE_INT, memspace, space, H5P_DEFAULT, tempGidOut);

	//Read temp data into corresponding block coord arrays.
	for (int i = 0; i < (int)dims_out[0]; i++) {
		for (int j = 0; j < (int)dims_out[1]; j++) {
			gid[i][j] = tempGidOut[dims_out[1]*i + j];
		}
	}

	delete[] tempGidOut;
	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);

	// std::cout << "GID READ TEST PROC: " << comm.Get_rank() << std::endl;
	// std::cout << gid[10][2] << std::endl;
}

void FlashAmrMesh::readFileType() {
	/* Suppress the error output temporarily. We want to see if we're using an HDF5
	   file with file format 7. If it exists, read the file format version. If it doesn't,
	   then exit with file type 0.
	*/
	
	if (!datasetAvailable("file format version","optional")) {
		fileType = 0;
		return;
	}
	
	hid_t dset = H5Dopen(inFile, "file format version", H5P_DEFAULT);

	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[1];

	H5Sget_simple_extent_dims( space, dims_out, NULL );

	hid_t memspace = H5Screate_simple(1, dims_out, NULL);
	
	//Create temp 1d array. Should be of size 1.
	int *fileTypeOut = new int[dims_out[0]];

	H5Dread( dset, H5T_NATIVE_INT, memspace, space, H5P_DEFAULT, fileTypeOut);

	fileType = fileTypeOut[0];

	delete[] fileTypeOut;
	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}


double* FlashAmrMesh::getCoordinates(int blockIndex) {
	return blocks[blockIndex].coord;
}

double FlashAmrMesh::getCoordinate(int blockIndex, int axis) {
	return blocks[blockIndex].coord[axis];
}

double* FlashAmrMesh::getSizes(int blockIndex) {
	return blocks[blockIndex].size;
}

double FlashAmrMesh::getSize(int blockIndex, int axis) {
	return blocks[blockIndex].size[axis];
}

double* FlashAmrMesh::getBounds(int blockIndex, int upperOrLower) {
	if (upperOrLower == LOWER) {
		return blocks[blockIndex].lowerBound;
	}
	else if (upperOrLower == UPPER) {
		return blocks[blockIndex].upperBound;
	}
	else {
		cout << "[flash_reader] Must use UPPER or LOWER in getBounds." << endl;
		exit(EXIT_FAILURE);
	}
}

double FlashAmrMesh::getBound(int blockIndex, int upperOrLower, int axis) {
	if (upperOrLower == LOWER) {
		return blocks[blockIndex].lowerBound[axis];
	}
	else if (upperOrLower == UPPER) {
		return blocks[blockIndex].upperBound[axis];
	}
	else {
		cout << "[flash_reader] Must use UPPER or LOWER in getBounds." << endl;
		exit(EXIT_FAILURE);
	}	
}

// std::vector<std::vector<int>>& FlashAmrMesh::getGidVector() {
// 	return gidVector;
// }

bool FlashAmrMesh::datasetAvailable(const char* dataset, const char* need) {
	// Checks the opened file to see if the dataset is present.

	// Though this checks for links, datasets can be checked for too.
	// Links are apparently just names in HDF5 files, where those
	// names can point to a number of hdf5 containers.
	htri_t returnVal = H5Lexists(inFile, dataset, H5P_DEFAULT );

	// Crash the code if this variable is absolutely required for VisIt.
	if (returnVal == false) {
		if (string(need) == "required") {
			std::cout << "[flash_reader] Critical error: No '" << dataset << "' dataset found in the input file." << std::endl;
			std::cout << "[flash_reader] This variable is required in order to read the HDF5 file into VisIt. Stopping execution."  << std::endl;
			exit(EXIT_FAILURE);
		}
		else if (string(need) == "optional") {
			if (!SUPPRESSOUTPUT) std::cout << "[flash_reader] Warning: No '" << dataset << "' dataset found in the input file." << std::endl;
			if (!SUPPRESSOUTPUT) std::cout << "[flash_reader] This dataset will not be written to the output file." << std::endl;
		}
	}

	return returnVal;

}


void FlashAmrMesh::readNFileVarsList() {
	//Reads in the list of unknown names from the HDF5 file-> These coorespond to the
	//four character variables output by FLASH like gamc, pres, dens, etc.

	//Also reads in the total amount of variables in the file. 

	if (!datasetAvailable("unknown names","required")) {
		return;
	}

		//Get coordinate scalars list.
		// H5std_string datasetName( "unknown names" );
		// DataSet varData = file->openDataSet( datasetName );
	
	//Get scalars list.
	hid_t dset = H5Dopen(inFile, "unknown names", H5P_DEFAULT);

	//Get dimensions of dataset.
	hsize_t dims_out[2];
	hid_t space = H5Dget_space( dset );

	H5Sget_simple_extent_dims( space, dims_out, NULL );
	
	hid_t memspace = H5Screate_simple(1, dims_out, NULL);

	hid_t strtype = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype, MESH_VAR_STRING_SIZE+1);

	//Create temp 1d array.
	char (*tempVarOut)[MESH_VAR_STRING_SIZE+1] = new char[dims_out[0]][MESH_VAR_STRING_SIZE+1];

	H5Dread( dset, strtype, memspace, space, H5P_DEFAULT, tempVarOut);

	
	//Copy variable names over.
	for (int i = 0; i < (int)dims_out[0]; i++) {
		// varNamesMap_iToS.insert(pair<std::string,int>(string(tempVarOut[i]),i));
		// varNamesMap_sToI.insert(pair<int,std::string>(i,string(tempVarOut[i])));

		//strcpy(variableNames[i],tempVarOut[i]);
		//cout << variableNames[i] << endl;

		fileVarsMap_sToI.insert(pair<std::string,int>(string(tempVarOut[i]),i));
		fileVarsMap_iToS.insert(pair<int,std::string>(i,string(tempVarOut[i])));

	}
	nFileVars = (int)dims_out[0];

	//Fill the rest of the char array with "" so we can count it easily later.
	for (int i = (int)dims_out[0]; i < FR_MAXVARS; i++) {
		strcpy(variableNames[i],"");
	}

	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);

	delete[] tempVarOut;

}

int FlashAmrMesh::getNFileVars() {
	return nFileVars;
}

int FlashAmrMesh::getNUserVars() {
	return nUserVars;
}

void FlashAmrMesh::addVariableNames(std::vector<std::string> fileVars, std::vector<std::string> auxVars) {
	// Reads in the names of both the requested variables found in the file (fileVars)
	// as well as the names of the new auxiliary (derived) variables. Space allocation
	// is done in allocateVarArray() and the reading of the file variables is done in both
	// readVariables() and readVarData().
	// 
	// Note that either ALL of the variables or a list of specific variables can be read in, but not
	// both. If all mesh variables are to be read in, fileVars should only contain one element: 'allVars'.
	// Otherwise, populate it with the specific mesh variables. The auxVars function does not need to be 
	// populated.

	int numInputMeshVars = (int)fileVars.size();

	// if (numInputMeshVars == 0) {
	// 	cout << "[flash reader] No variables specified. Exiting." << endl;
	// 	exit(EXIT_FAILURE);
	// }

	for (string variable : fileVars) {
		if (strlen(variable.c_str()) > 4 && variable != "allVars") {
			cout << "[flash_reader] Size of variable '" << variable << "' must be <= 4 characters." << endl;
			exit(EXIT_FAILURE);		
		}

		if (numInputMeshVars == 1 && variable == "allVars") {
			// Add all variables
			for (int i = 0; i < nFileVars; i++) {
				addUserVariable(fileVarsMap_iToS[i].c_str(), "file");
			}
		}
		else {
			// Add the specified variables
			addUserVariable(variable.c_str(), "file");	
		}
	}
	for (string variable : auxVars) {
		if (strlen(variable.c_str()) > 4) {
			cout << "[flash_reader] Size of variable '" << variable << "' must be <= 4 characters." << endl;
			exit(EXIT_FAILURE);		
		}

		// Pad to 4 characters like other FLASH variables
		std::string paddedVar = string(variable);
		paddedVar.insert(paddedVar.end(), 4 - paddedVar.size(), ' ');

		addUserVariable(paddedVar.c_str(), "derived");
	}
}

void FlashAmrMesh::allocateVarArray() 
{
	// Used to make sure we don't try to deallocate something that isn't there if the most 
	// basic constructor 'FlashAmrMesh()' is used and no variables are later allocated in main.
	variablesAllocated = true;

	// Reserve 16K of memory that can be deleted just in case we run out of memory
	char* _emergencyMemory = new char[16384];

	double sizeEstimate = (double)(sizeof(double)*getNUserVars()*nBlocksAssigned*nCellsVec[IAXIS]*nCellsVec[JAXIS]*nCellsVec[KAXIS]/(1024));
	if (!SUPPRESSOUTPUT) cout << "Attempting to allocate variable array." << endl;
	if (!SUPPRESSOUTPUT) cout << "Need approx. " << sizeEstimate << "kb of space on processor " << comm.Get_rank() << endl;
	for (int blk = 0; blk < nBlocksAssigned; blk++) {
		try
		{
			blocks[blk].data = new double***[getNUserVars()];
		}
		catch (bad_alloc ex)
		{
			// Delete the reserved memory so we can print an error message before exiting
    		delete[] _emergencyMemory;

			cout << "[FlashAmrMesh] Out of memory while allocating variable array. Exiting..." << endl;
			exit(EXIT_FAILURE);    		
		}

		for (int i = 0; i < getNUserVars(); i++) {
			try
			{
				blocks[blk].data[i] = new double**[nCellsVec[IAXIS]];
			}
			catch (bad_alloc ex)
			{
				// Delete the reserved memory so we can print an error message before exiting
	    		delete[] _emergencyMemory;

				cout << "[FlashAmrMesh] Out of memory while allocating variable array. Exiting..." << endl;
				exit(EXIT_FAILURE);    		
			}

			for (int j = 0; j < nCellsVec[IAXIS]; j++) {
				try
				{
					blocks[blk].data[i][j] = new double*[nCellsVec[JAXIS]];
				}
				catch (bad_alloc ex)
				{
					// Delete the reserved memory so we can print an error message before exiting
		    		delete[] _emergencyMemory;

					cout << "[FlashAmrMesh] Out of memory while allocating variable array. Exiting..." << endl;
					exit(EXIT_FAILURE);    		
				}

				for (int k = 0; k < nCellsVec[JAXIS]; k++) {
					try
					{
						blocks[blk].data[i][j][k] = new double[nCellsVec[KAXIS]];
					}
					catch (bad_alloc ex)
					{
						// Delete the reserved memory so we can print an error message before exiting
			    		delete[] _emergencyMemory;

						cout << "[FlashAmrMesh] Out of memory while allocating variable array. Exiting..." << endl;
						exit(EXIT_FAILURE);    		
					}
				}
			}
		}
	}

	// delete the reserve memory
	delete[] _emergencyMemory;
}

void FlashAmrMesh::deallocateVarArray() {

	for (int blk = 0; blk < nBlocksAssigned; blk++) {
		if (blocks[blk].data != NULL) {
			
			for (int i = 0; i < getNUserVars(); i++) {
				for (int j = 0; j < nCellsVec[IAXIS]; j++) {
					for (int k = 0; k < nCellsVec[JAXIS]; k++) {
						delete[] blocks[blk].data[i][j][k];
					}
					delete[] blocks[blk].data[i][j];
				}
				delete[] blocks[blk].data[i];
			}
			delete[] blocks[blk].data;
		}
		else {
			if (FORCEDEALLOCATECRASH) {
				cout << "[flash_reader] Failure in deallocateVarArray function." << endl;
				exit(EXIT_FAILURE);	
			}
		}
		
	}

}

int FlashAmrMesh::findVarIndex(const char requestedVar[]) {
	//Returns the index of the requested variable from the
	// data array.

	if (strlen(requestedVar) > 4) {
		cout << "[flash_reader] Size of variable '" << requestedVar << "' must be <= 4 characters." << endl;
		exit(EXIT_FAILURE);		
	}
	// Pad to 4 characters like other FLASH variables
	std::string paddedVar = string(requestedVar);
	paddedVar.insert(paddedVar.end(), 4 - paddedVar.size(), ' ');

	int returnIndex;

	// Make sure the variable is there
	std::map<std::string,int>::iterator iter = uVarsMap_sToI.find(paddedVar);
    if (iter == uVarsMap_sToI.end()) {
        cout << "[flash_reader] Could not find requested variable '" << paddedVar << "' in findVarIndex function. " << endl;
		cout << "Note that the requested variable must exist in the searchable array, must be four characters long, " << endl;
		cout << "and must be in lowercase." << endl;
		returnIndex = -1;
		//exit(EXIT_FAILURE); //Leaving this commented out to allow the caller to handle and explain error
    }
    else {
    	returnIndex = uVarsMap_sToI[paddedVar];
    }

	return returnIndex;

}

std::string FlashAmrMesh::findVarName(const int index) {
	std::string outStr;

	std::map<int,std::string>::iterator iter = uVarsMap_iToS.find(index);
	if (iter == uVarsMap_iToS.end()) {
		cout << "[flash_reader] Could not find requested variable at index" << index << " in findVarName function. " << endl;
		outStr = "-1";
		//exit(EXIT_FAILURE); //Leaving this commented out to allow the caller to handle and explain error
	}
	else {
		outStr = uVarsMap_iToS[index];
	}
	return outStr;
}


void FlashAmrMesh::readVarData(const int varIndex, const char requestedVar[] ) {
	//Internal function used by the variable read functions. This function
	// takes in the requested variable, searches for the corresponding index for
	// that variable in the array variableNames, and fills each block's data array
	// at that index.

	//std::cout << "Processor in readVarData: " << comm.Get_rank() << std::endl;
	//Get coordinate scalars list.
	
	if (!datasetAvailable(requestedVar,"required")){ 
		return;
	}

	hid_t dset = H5Dopen(inFile, requestedVar, H5P_DEFAULT);
	
	hid_t space = H5Dget_space( dset );

	hsize_t dims_out[4];
	H5Sget_simple_extent_dims( space, dims_out, NULL );

	//Define memory space
	hid_t memspace = H5Screate_simple(4, dims_out, NULL);

	//Create temp 1d array.
	double *tempVarOut;
	tempVarOut = new double[1*(int)dims_out[1]*(int)dims_out[2]*(int)dims_out[3]];

	hsize_t offset[4]; //hyperslab offset in file
	hsize_t count[4];  //size of hyperslab
	hsize_t offset_out[4];   // hyperslab offset in memory
	hsize_t count_out[4];    // size of the hyperslab in memory

	for (int blk = 0; blk < nBlocksAssigned; blk++) {
				
		offset[0] = mpiGetGlobalBlockID(blk);
		offset[1] = 0;
		offset[2] = 0;
		offset[3] = 0;

		count[0]  = 1;
		count[1]  = dims_out[1];
		count[2]  = dims_out[2];
		count[3]  = dims_out[3];

		H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL);

		/* http://www.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html */

	    offset_out[0] = 0;
	    offset_out[1] = 0;
	    offset_out[2] = 0;
	    offset_out[3] = 0;
	    count_out[0]  = 1;
	    count_out[1]  = dims_out[1];
	    count_out[2]  = dims_out[2];
	    count_out[3]  = dims_out[3];

	    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);

		H5Dread( dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, tempVarOut);

		//Read temp data into corresponding block coord arrays.
		int index;
		for (int i = 0; i < (int)dims_out[3]; i++) {
			for (int j = 0; j < (int)dims_out[2]; j++) {
				for (int k = 0; k < (int)dims_out[1]; k++) {
					index = k * (int)dims_out[2]*(int)dims_out[3] + j * (int)dims_out[3] + i;
					// std::cout << "i: " << i << " j: " << j << " k: " << k << std::endl;
					blocks[blk].data[varIndex][i][j][k] = tempVarOut[index];
					// std::cout << "Block: " << mpiGetGlobalBlockID(blk) << " Index: " << k   * dims_out[2]*dims_out[3] + 
					// 			    j   * dims_out[3] +
					// 			    i   << " Value: " << tempVarOut[k   * dims_out[2]*dims_out[3] + 
					// 			    j   * dims_out[3] +
					// 			    i   ] << std::endl;
				}
			}
		}

		// std::cout << blk << " " << mpiGetGlobalBlockID(blk) << "/" << nBlocks << " ";
		// std::cout << blocks[blk].data[index][0][0][0] << std::endl;

		
		
	}

	if (!SUPPRESSOUTPUT) cout << "Read variable: " << requestedVar << " at index: " << varIndex << endl;

	delete[] tempVarOut;

	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}

void FlashAmrMesh::initAuxVarData(const int varIndex, const char requestedVar[]) {
	// Initialize all of the data for the given auxiliary variable to zero.
	// This prevents the write function from outputting garbage to certain blocks when, 
	// for example, only the leaf blocks are written to.
	
	hsize_t dims_out[4];
	dims_out[0] = getNBlocks();
	dims_out[1] = getNcells(KAXIS);
	dims_out[2] = getNcells(JAXIS);
	dims_out[3] = getNcells(IAXIS);

	for (int blk = 0; blk < nBlocksAssigned; blk++) {
		int index;
		for (int i = 0; i < (int)dims_out[3]; i++) {
			for (int j = 0; j < (int)dims_out[2]; j++) {
				for (int k = 0; k < (int)dims_out[1]; k++) {
					index = k * (int)dims_out[2]*(int)dims_out[3] + j * (int)dims_out[3] + i;
					// std::cout << "i: " << i << " j: " << j << " k: " << k << std::endl;
					blocks[blk].data[varIndex][i][j][k] = 0;
					// std::cout << "Block: " << mpiGetGlobalBlockID(blk) << " Index: " << k   * dims_out[2]*dims_out[3] + 
					// 			    j   * dims_out[3] +
					// 			    i   << " Value: " << tempVarOut[k   * dims_out[2]*dims_out[3] + 
					// 			    j   * dims_out[3] +
					// 			    i   ] << std::endl;
				}
			}
		}
	}
	if (!SUPPRESSOUTPUT) cout << "Initialized aux variable: " << requestedVar << " at index: " << varIndex << endl;
}

void FlashAmrMesh::noWrite(const char inputVariable[]) {
	// All variables read in are written to a file (if we are writing to a file at all) by default.
	// This function allows the user to not output a specified variable.

	if (strlen(inputVariable) > 4) {
		cout << "[flash_reader] Size of variable '" << inputVariable << "' must be <= 4 characters." << endl;
		exit(EXIT_FAILURE);		
	}
	// Pad to 4 characters like other FLASH variables
	std::string paddedVar = string(inputVariable);
	paddedVar.insert(paddedVar.end(), 4 - paddedVar.size(), ' ');

	// Make sure the variable is there
	std::map<std::string,bool>::iterator iter = writeVarMap.find(paddedVar);
    if (iter == writeVarMap.end()) {
		cout << "[flash_reader] Could not find requested variable " << paddedVar << " in the noWrite function. " << endl;
		cout << "Note that the requested variable must exist in the searchable array, must be four characters long, " << endl;
		cout << "and must be in lowercase." << endl;
		exit(EXIT_FAILURE);
    }
    
    // for (iter = writeVarMap.begin(); iter != writeVarMap.end(); iter++) {
    // 	cout << iter->first << " " << iter->second << " " << endl;
    // }

    writeVarMap[paddedVar] = false;

    // for (iter = writeVarMap.begin(); iter != writeVarMap.end(); iter++) {
    // 	cout << iter->first << " " << iter->second << " " << endl;
    // }
    
}

void FlashAmrMesh::addUserVariable(const char inputVariable[], const char varType[]) {
	//Adds a variable to the userVariableNames array and
	// updates the number of user requested variables.
	// First the variables are added (in main) using this function,
	// and once all of the desired variables have been added, 
	// allocateVariables is called, then readVariables is called to read 
	// them all at once. Added variables must be four characters long. 
	//    Example call: addUserVariable("pres");

	//This function will automatically check to make sure the desired variable
	// exists in the file via the user variables map.

	if (nUserVars == 0) {
		//Initialize the list.
		for (int i = 0; i < FR_MAXVARS; i++) {
			strcpy(variableNames[i],"");
		}
	}

	// See if we've exceeded the limit for the number of variables
	if (nUserVars >= FR_MAXVARS) {
		cout << "[flash_reader] Could not add " << inputVariable << "." << endl;
		cout << "Too many variables. Can only have a maximum of " << nUserVars << " variables. " << endl;
		exit(EXIT_FAILURE);
	}

	// Different behavior depending on whether the variable is already in the HDF5 file.
	if (string(varType) == "file"){
		//Check to make sure the variable exists in the map.
		std::map<std::string,int>::iterator it;
		it = fileVarsMap_sToI.find(string(inputVariable));
		if (it == fileVarsMap_sToI.end()) {
			cout << "[flash_reader] Variable '" << inputVariable << "' does not exist in the file variables map." << endl;
			exit(EXIT_FAILURE);
		}
		else {
			//Add variable to list and update nUserVars. Note that this list can include derived variables too.
			strcpy(variableNames[nUserVars],inputVariable);
			nUserVars++;
			//Add variable to the list of present variables so we know it should be read.
			strcpy(variablesInFile[nReqFileVars],inputVariable);
			nReqFileVars++;
		}
	}
	else if (string(varType) == "derived") {
		//Add variable to list and update nUserVars. Note that this list can include derived variables too.
		strcpy(variableNames[nUserVars],inputVariable);
		nUserVars++;
		//Add variable to the list of derived variables so we know it should not be read.
		strcpy(variablesAux[nReqAuxVars],inputVariable);
		nReqAuxVars++;
	}
	else {
		cout << "[flash_reader] User variable input type not recognized. Acceptable inputs are 'file' for ";
	    cout << "variables present in the HDF5 file or 'derived' for new derived variables in LAVAflow." << endl;
	    exit(EXIT_FAILURE);
	}

	// Assume that we want to write this variable if we're creating an output file.
	// This can be changed by the user with the noWrite(variable name) function
	writeVarMap.insert(pair<std::string,bool>(std::string(inputVariable),true));
}

void FlashAmrMesh::readVariables(){
	//Reads the user desired variables into the data array. See the
	// addUserVariable function.

	/* Overall structure of variable reading:
		NOTE: The code does not have the functionality to use both addAllVariables and addUserVariables at
			  the same time. 

		First: Read in variables from the file. Put them inside the fileVarsMap(s).
					-Two maps: int --> string and string --> int.

		Second: If specific variables are being read in, addUserVariable(var) will be used.
				This will check to make sure the variable exists within fileVarsMap and will
				add the variable name to the array variableNames[]. 

				If all variables are being read in, addAllVariables() will be used.
				This will take all variables from the fileVarsMap and put them into variableNames[].
		
		Third:  Allocate the data array by the nUserVars value. This value is updated in both the 
				addAllVariables() function and the addUserVariable() function.

		Fourth: Take all the variables found in variableNames[] and read them in using readVarData.
				This function is responsible for that. It will read each variable one by one and add
				them to the uVarsMap map.
					-Two maps: int --> string and string --> int.
				These maps allow us to easily access values in the variableNames array by string or by int.
				However, they must be called using functions:
					findVarIndex(requestedVar[]) returns an index in variableNames[] corresponding to the requested variable.
					findVarName(index) returns a variable name for some index. Basically an array access at that index.

	*/

	// Read variables in the file
	int varIndex = 0;
	while (varIndex < nReqFileVars) {
		readVarData(varIndex,variablesInFile[varIndex]);

		//Insert variable into maps for later reference.
		uVarsMap_sToI.insert(pair<std::string,int>(string(variablesInFile[varIndex]),varIndex));
		uVarsMap_iToS.insert(pair<int,std::string>(varIndex,string(variablesInFile[varIndex])));
		varIndex++;
	}

	// Add aux variables to the maps
	int varAuxIndex = 0;
	while (varIndex < nUserVars) {
		initAuxVarData(varIndex,variablesAux[varAuxIndex]);

		//Insert variable into maps for later reference.
		uVarsMap_sToI.insert(pair<std::string,int>(string(variablesAux[varAuxIndex]),varIndex));
		uVarsMap_iToS.insert(pair<int,std::string>(varIndex,string(variablesAux[varAuxIndex])));
		varAuxIndex++;
		varIndex++;
	}
}

double FlashAmrMesh::getVariable(int blockId, int varIndex, int i, int j, int k) {
	return blocks[blockId].data[varIndex][i][j][k];
}

void FlashAmrMesh::deallocateVariables() {

	deallocateGid();
	if (variablesAllocated) {
		deallocateVarArray();
	}
	deallocateBlocks();

}

void FlashAmrMesh::generalRead() {
	//Reads in everything except for the variables. Each function
	// written is called here to fill in as much information as possible
	// for the grid and blocks. To read in the variables, either
	// readAllVariables can be used to read all of the variables, or
	// readPartialvariables can be used to read only some of the variables.
	// The list of user desired variables is generated using multiple calls
	// to addUserVariable. See that function for more details.

	if (!SUPPRESSOUTPUT) cout << "Reading file type, if it exists. " << endl;
	this->readFileType();
	//fileType = 7;
	if (!SUPPRESSOUTPUT) cout << "File type: " << fileType << endl;

	if (!SUPPRESSOUTPUT) cout << "Reading scalar lists. " << endl;
	if (fileType == 7) {
		this->readSimParameters();
		//Dim, k2d, and k3d are read in using readCoordinates.
		//See readCoordinates function for more details.
	}
	else {
		this->readInts("scalars");
		this->readStrings("scalars");
		this->readReals("scalars");
		this->readLogicals("scalars");
	}

	//Read in parameters after scalars so that any values with the same name
	//take the parameter values rather than the scalar values.
	if (!SUPPRESSOUTPUT) cout << "Reading parameter lists. " << endl;
	this->readInts("parameters");
	this->readStrings("parameters");
	this->readReals("parameters");
	this->readLogicals("parameters");

	this->readNBlocks();
	if (!SUPPRESSOUTPUT) cout << "Number of blocks: " << this->getNBlocks() << endl;

	//MPI stuff
	this->mpiAssignBlocks();

	this->allocateBlocks();
	if (!SUPPRESSOUTPUT) cout << "Blocks allocated. " << endl;

	this->readBlockNumber();
	if (!SUPPRESSOUTPUT) cout << "Block numbers assigned. " << endl;

	if (fileType != 7) {
		this->readDim();
	}

	this->readGeometryType();
	if (!SUPPRESSOUTPUT) cout << "Geometry Type: " << geometryName << endl;

	this->readCoordinates();
	if (!SUPPRESSOUTPUT) cout << "Coordinates read. " << endl;

	this->readCPUNumbers();
	if (!SUPPRESSOUTPUT) cout << "CPU numbers read. " << endl;

	if (!SUPPRESSOUTPUT) cout << "Dimensionality: " << this->getDim() << endl;

	this->readSizes();
	if (!SUPPRESSOUTPUT) cout << "Physical sizes of blocks read. " << endl;

	this->readBounds();
	if (!SUPPRESSOUTPUT) cout << "Block bounds read. " << endl;

	this->readNodeType();
	if (!SUPPRESSOUTPUT) cout << "Node types read. " << endl;

	this->readRefineLevel();
	if (!SUPPRESSOUTPUT) cout << "Refine levels read. " << endl;

	this->readTime();
	if (!SUPPRESSOUTPUT) cout << "Simulation time read. " << "Time: " << this->getTime() << endl;

	this->allocateGid();
	if (!SUPPRESSOUTPUT) cout << "Gid array allocated. " << endl;

	this->readGid();
	if (!SUPPRESSOUTPUT) cout << "Gid array read. " << endl;

	this->assignGidArrays();
	if (!SUPPRESSOUTPUT) cout << "Parent, neighbors, and children assigned. " << endl;

	this->readNFileVarsList();
	if (!SUPPRESSOUTPUT) cout << "List of variable names read. " << endl;
	if (!SUPPRESSOUTPUT) cout << "Number of variables in file: " << this->getNFileVars() << endl;

	this->readNcells();
	if (!SUPPRESSOUTPUT) cout << "Cells per block read in. Total cells per block: " << this->getTotalCells() << endl; 

	/*Allocation and reading of variables is done in main.*/

}

void FlashAmrMesh::openInFile() {
	MPI_Barrier(MPI_COMM_WORLD);
	if (!SUPPRESSOUTPUT) cout << "Opening input file on processor " << comm.Get_rank() << "...";
	inFile = H5Fopen (this->getGridInFileName().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (!SUPPRESSOUTPUT) cout << "Input file opened on processor" << comm.Get_rank() << "." << endl;
}

void FlashAmrMesh::closeInFile() {
	if (!SUPPRESSOUTPUT) cout << "Trying to close input file on processor " << comm.Get_rank() << "...";
	H5Fclose(inFile);
	if (!SUPPRESSOUTPUT) cout << "Input file closed on processor " << comm.Get_rank() << "." << endl;
}

void FlashAmrMesh::readInFile(vector<string> fileVars, vector<string> auxVars) {
	if (mpiGetID() == 0) {
		if (!SUPPRESSOUTPUT) std::cout << "Reading file" << std::endl;
	}
	openInFile();
	generalRead();
	addVariableNames(fileVars, auxVars);
	allocateVarArray();
	readVariables();
	// The file can't be closed until the end of execution since the 
	// write functions partially read the input datasets for metadata.
	if (mpiGetID() == 0) {
		if (!SUPPRESSOUTPUT) std::cout << "Reading complete" << std::endl;
	}
}

void FlashAmrMesh::createOutFile() {

	if (!SUPPRESSOUTPUT) cout << "Creating output file...";

    // Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);

    //Create a new file collectively and release property list identifier.
	outFile = H5Fcreate(this->getGridOutFileName().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	if (!SUPPRESSOUTPUT) cout << "Output file created." << endl;

}

void FlashAmrMesh::writeFileType() {
	// Open the file type from the input file 'file' and
	// output it in the output file 'outFile'. This function
	// is specific to file type '7' FLASH files.

	if (!datasetAvailable("file format version","required")) {
		return;
	}

    /* Read old dataset */
	hid_t dset_in = H5Dopen(inFile, "file format version", H5P_DEFAULT);
	hid_t space_in = H5Dget_space( dset_in );
	hsize_t dims_in[1];
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );
	hid_t memspace_in = H5Screate_simple(1, dims_in, NULL);
	int *fileType_in = new int[dims_in[0]];
	H5Dread( dset_in, H5T_NATIVE_INT, memspace_in, space_in, H5P_DEFAULT, fileType_in);

    /* Write new dataset */
    // There's redundant variables here, but it helps clarify what's being written.
	hid_t space_out = H5Dget_space( dset_in );
	hsize_t dims_out[1];
	dims_out[0] = dims_in[0];
	hid_t memspace_out = H5Screate_simple(space_out, dims_out, NULL);

    hid_t dset_out = H5Dcreate2(outFile, "file format version", H5T_NATIVE_INT, memspace_out,
    	                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, fileType_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);

	delete[] fileType_in;
	
	H5Sclose(memspace_in);
	H5Sclose(space_in);
	H5Dclose(dset_in);

	H5Sclose(memspace_out);
	H5Sclose(space_out);
	H5Dclose(dset_out);
}

void FlashAmrMesh::writeSimParameters() {
	// This function is for outputting file type 7 information

	typedef struct simParamStruct {
		//Instead of this file being formatted like the rest where there is a name
		// column and a value column, file format 7 has 8 columns with no names.
		// The values have to be read directly into the variables here.
		int totalBlocks;
		float time;
		float timeStep;
		float redShift;
		int numberOfSteps;
		int nxb;
		int nyb;
		int nzb;
	} simParamStruct;

	if (!datasetAvailable("simulation parameters","required")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "simulation parameters", H5P_DEFAULT);
	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );

	//"simulation parameters" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(simParamStruct));

	//NOTE: IF THERE ARE ANY EXTRA PARAMETERS HERE, THEY WILL NOT BE READ.
	//      THIS IS DUE TO THE WEIRD WAY THAT THE PARAMETERS WERE WRITTEN TO THE HDF5 FILE.
	H5Tinsert( memtype_in, "total blocks", 	HOFFSET(simParamStruct, totalBlocks), 	H5T_NATIVE_INT);
	H5Tinsert( memtype_in, "time", 			HOFFSET(simParamStruct, time), 			H5T_NATIVE_FLOAT);
	H5Tinsert( memtype_in, "timestep", 		HOFFSET(simParamStruct, timeStep), 		H5T_NATIVE_FLOAT);
	H5Tinsert( memtype_in, "redshift", 		HOFFSET(simParamStruct, redShift), 		H5T_NATIVE_FLOAT);
	H5Tinsert( memtype_in, "number of steps", 	HOFFSET(simParamStruct, numberOfSteps), H5T_NATIVE_INT);
	H5Tinsert( memtype_in, "nxb", 				HOFFSET(simParamStruct, nxb), 			H5T_NATIVE_INT);
	H5Tinsert( memtype_in, "nyb", 				HOFFSET(simParamStruct, nyb), 			H5T_NATIVE_INT);
	H5Tinsert( memtype_in, "nzb", 				HOFFSET(simParamStruct, nzb), 			H5T_NATIVE_INT);
	
	simParamStruct* simParamVals_in = new simParamStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, simParamVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(simParamStruct) );
	H5Tinsert( memtype_out, "total blocks", 	HOFFSET(simParamStruct, totalBlocks), 	H5T_NATIVE_INT);
	H5Tinsert( memtype_out, "time", 			HOFFSET(simParamStruct, time), 			H5T_NATIVE_FLOAT);
	H5Tinsert( memtype_out, "timestep", 		HOFFSET(simParamStruct, timeStep), 		H5T_NATIVE_FLOAT);
	H5Tinsert( memtype_out, "redshift", 		HOFFSET(simParamStruct, redShift), 		H5T_NATIVE_FLOAT);
	H5Tinsert( memtype_out, "number of steps", 	HOFFSET(simParamStruct, numberOfSteps), H5T_NATIVE_INT);
	H5Tinsert( memtype_out, "nxb", 				HOFFSET(simParamStruct, nxb), 			H5T_NATIVE_INT);
	H5Tinsert( memtype_out, "nyb", 				HOFFSET(simParamStruct, nyb), 			H5T_NATIVE_INT);
	H5Tinsert( memtype_out, "nzb", 				HOFFSET(simParamStruct, nzb), 			H5T_NATIVE_INT);

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(simParamStruct) );
	H5Tinsert( filetype_out, "total blocks", 	HOFFSET(simParamStruct, totalBlocks), 	H5T_NATIVE_INT);
	H5Tinsert( filetype_out, "time", 			HOFFSET(simParamStruct, time), 			H5T_NATIVE_FLOAT);
	H5Tinsert( filetype_out, "timestep", 		HOFFSET(simParamStruct, timeStep), 		H5T_NATIVE_FLOAT);
	H5Tinsert( filetype_out, "redshift", 		HOFFSET(simParamStruct, redShift), 		H5T_NATIVE_FLOAT);
	H5Tinsert( filetype_out, "number of steps", 	HOFFSET(simParamStruct, numberOfSteps), H5T_NATIVE_INT);
	H5Tinsert( filetype_out, "nxb", 				HOFFSET(simParamStruct, nxb), 			H5T_NATIVE_INT);
	H5Tinsert( filetype_out, "nyb", 				HOFFSET(simParamStruct, nyb), 			H5T_NATIVE_INT);
	H5Tinsert( filetype_out, "nzb", 				HOFFSET(simParamStruct, nzb), 			H5T_NATIVE_INT);

    hid_t dset_out = H5Dcreate(outFile, "simulation parameters", filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, simParamVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, simParamVals_in);
	H5Dclose(dset_in);
	H5Sclose(space_in);
	H5Tclose(memtype_in);

	delete[] simParamVals_in;
}

void FlashAmrMesh::writeSimInfo() {
	// This function is only called if it is file type 9 
	// (which is the default behavior).

	typedef struct simInfoStruct {
		// File format 9 sim info
		int fileFormatVersion;
		char setupCall[400];
		char fileCreationTime[MAX_STRING_LENGTH];
		char flashVersion[MAX_STRING_LENGTH];
		char buildDate[MAX_STRING_LENGTH];
		char buildDir[MAX_STRING_LENGTH];
		char buildMachine[MAX_STRING_LENGTH];
		char cFlags[400];
		char fFlags[400];
		char setupTimeStamp[MAX_STRING_LENGTH];
		char buildTimeStamp[MAX_STRING_LENGTH];
	} simInfoStruct;

	if (!datasetAvailable("sim info","required")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "sim info", H5P_DEFAULT);
	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );


	hid_t strtype80_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype80_in, MAX_STRING_LENGTH);
	hid_t strtype400_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype400_in, 400);

	//"simulation parameters" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(simInfoStruct));

	H5Tinsert( memtype_in, "file format version", HOFFSET(simInfoStruct, fileFormatVersion), H5T_NATIVE_INT  );
	H5Tinsert( memtype_in, "setup call",          HOFFSET(simInfoStruct, setupCall),         strtype400_in );
	H5Tinsert( memtype_in, "file creation time",  HOFFSET(simInfoStruct, fileCreationTime),  strtype80_in  );
	H5Tinsert( memtype_in, "flash version",       HOFFSET(simInfoStruct, flashVersion),      strtype80_in  );
	H5Tinsert( memtype_in, "build date",          HOFFSET(simInfoStruct, buildDate),         strtype80_in  );
	H5Tinsert( memtype_in, "build dir",           HOFFSET(simInfoStruct, buildDir),          strtype80_in  );
	H5Tinsert( memtype_in, "build machine",       HOFFSET(simInfoStruct, buildMachine),      strtype80_in  );
	H5Tinsert( memtype_in, "cflags",              HOFFSET(simInfoStruct, cFlags),            strtype400_in );
	H5Tinsert( memtype_in, "fflags",              HOFFSET(simInfoStruct, fFlags),            strtype400_in );
	H5Tinsert( memtype_in, "setup time stamp",    HOFFSET(simInfoStruct, setupTimeStamp),    strtype80_in  );
	H5Tinsert( memtype_in, "build time stamp",    HOFFSET(simInfoStruct, buildTimeStamp),    strtype80_in  );
	
	simInfoStruct* simInfoVals_in = new simInfoStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, simInfoVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(simInfoStruct) );
	H5Tinsert( memtype_out, "file format version", HOFFSET(simInfoStruct, fileFormatVersion), H5T_NATIVE_INT  );
	H5Tinsert( memtype_out, "setup call",          HOFFSET(simInfoStruct, setupCall),         strtype400_in );
	H5Tinsert( memtype_out, "file creation time",  HOFFSET(simInfoStruct, fileCreationTime),  strtype80_in  );
	H5Tinsert( memtype_out, "flash version",       HOFFSET(simInfoStruct, flashVersion),      strtype80_in  );
	H5Tinsert( memtype_out, "build date",          HOFFSET(simInfoStruct, buildDate),         strtype80_in  );
	H5Tinsert( memtype_out, "build dir",           HOFFSET(simInfoStruct, buildDir),          strtype80_in  );
	H5Tinsert( memtype_out, "build machine",       HOFFSET(simInfoStruct, buildMachine),      strtype80_in  );
	H5Tinsert( memtype_out, "cflags",              HOFFSET(simInfoStruct, cFlags),            strtype400_in );
	H5Tinsert( memtype_out, "fflags",              HOFFSET(simInfoStruct, fFlags),            strtype400_in );
	H5Tinsert( memtype_out, "setup time stamp",    HOFFSET(simInfoStruct, setupTimeStamp),    strtype80_in  );
	H5Tinsert( memtype_out, "build time stamp",    HOFFSET(simInfoStruct, buildTimeStamp),    strtype80_in  );

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(simInfoStruct) );
	H5Tinsert( filetype_out, "file format version", HOFFSET(simInfoStruct, fileFormatVersion), H5T_NATIVE_INT  );
	H5Tinsert( filetype_out, "setup call",          HOFFSET(simInfoStruct, setupCall),         strtype400_in );
	H5Tinsert( filetype_out, "file creation time",  HOFFSET(simInfoStruct, fileCreationTime),  strtype80_in  );
	H5Tinsert( filetype_out, "flash version",       HOFFSET(simInfoStruct, flashVersion),      strtype80_in  );
	H5Tinsert( filetype_out, "build date",          HOFFSET(simInfoStruct, buildDate),         strtype80_in  );
	H5Tinsert( filetype_out, "build dir",           HOFFSET(simInfoStruct, buildDir),          strtype80_in  );
	H5Tinsert( filetype_out, "build machine",       HOFFSET(simInfoStruct, buildMachine),      strtype80_in  );
	H5Tinsert( filetype_out, "cflags",              HOFFSET(simInfoStruct, cFlags),            strtype400_in );
	H5Tinsert( filetype_out, "fflags",              HOFFSET(simInfoStruct, fFlags),            strtype400_in );
	H5Tinsert( filetype_out, "setup time stamp",    HOFFSET(simInfoStruct, setupTimeStamp),    strtype80_in  );
	H5Tinsert( filetype_out, "build time stamp",    HOFFSET(simInfoStruct, buildTimeStamp),    strtype80_in  );

    hid_t dset_out = H5Dcreate(outFile, "sim info", filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, simInfoVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, simInfoVals_in);
	H5Dclose(dset_in);
	H5Sclose(space_in);
	H5Tclose(memtype_in);

	delete[] simInfoVals_in;
}

void FlashAmrMesh::writeInts(const char *which){
	// Open the file type from the input file 'file' and
	// output it in the output file 'outFile'. Change the type of 
	// read based on 'scalars' or 'parameters' input.

	typedef struct intScalarStruct {
		char name[MAX_STRING_LENGTH];
		int  value;
	} intScalarStruct;

	char datasetName[50] = {};
	if (string(which) == "scalars") {
		strcpy(datasetName,"integer scalars");
	}
	else if (string(which) == "parameters") {
		strcpy(datasetName,"integer runtime parameters");
	}

	if (!datasetAvailable(datasetName,"required")) {
		return;
	}

	//Get scalars list.
	hid_t dset_in = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );

	//Make specific type for the 80 character string.
	hid_t strtype_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype_in, MAX_STRING_LENGTH);

	//"integer scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(intScalarStruct) );
	H5Tinsert( memtype_in, "name", HOFFSET(intScalarStruct, name), strtype_in);
	H5Tinsert( memtype_in, "value", HOFFSET(intScalarStruct, value), H5T_NATIVE_INT);

	intScalarStruct* intVals_in = new intScalarStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, intVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(intScalarStruct) );
	H5Tinsert( memtype_out, "name", HOFFSET(intScalarStruct, name), strtype_in);
	H5Tinsert( memtype_out, "value", HOFFSET(intScalarStruct, value), H5T_NATIVE_INT);

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(intScalarStruct) );
	H5Tinsert( filetype_out, "name", HOFFSET(intScalarStruct, name), strtype_in);
	H5Tinsert( filetype_out, "value", HOFFSET(intScalarStruct, value), H5T_NATIVE_INT);

    hid_t dset_out = H5Dcreate(outFile, datasetName, filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, intVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, intVals_in);
	H5Tclose(memtype_in);
	H5Sclose(space_in);
	H5Dclose(dset_in);	
	H5Tclose(strtype_in);

	delete[] intVals_in;

}

void FlashAmrMesh::writeLogicals(const char *which){
	// Open the file type from the input file 'file' and
	// output it in the output file 'outFile'. Change the type of 
	// read based on 'scalars' or 'parameters' input.

	typedef struct logicalScalarStruct {
		char name[MAX_STRING_LENGTH];
		int  value;
	} logicalScalarStruct;

	char datasetName[50] = {};
	if (string(which) == "scalars") {
		strcpy(datasetName,"logical scalars");
	}
	else if (string(which) == "parameters") {
		strcpy(datasetName,"logical runtime parameters");
	}

	if (!datasetAvailable(datasetName,"optional")) {
		return;
	}

    /* Check to make sure that the logicals list is in the file. */

	/* Save old error handler */
	herr_t (*old_func)(void*);
	void* old_client_data;
	H5Eget_auto1(&old_func, &old_client_data);

	/* Turn off error handling */
	H5Eset_auto1(NULL, NULL);

	hid_t dset_in = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	if (dset_in < 0) { // Dataset doesn't exist.
		if (!SUPPRESSOUTPUT) cout << "No " << datasetName << " found." << endl;
		return;
	}

	/* Restore previous error handler */
	H5Eset_auto1(old_func, old_client_data);

	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );

	//Make specific type for the 80 character string.
	hid_t strtype_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype_in, MAX_STRING_LENGTH);

	//"integer scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(logicalScalarStruct) );
	H5Tinsert( memtype_in, "name", HOFFSET(logicalScalarStruct, name), strtype_in);
	H5Tinsert( memtype_in, "value", HOFFSET(logicalScalarStruct, value), H5T_NATIVE_INT);

	logicalScalarStruct* logicalVals_in = new logicalScalarStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, logicalVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(logicalScalarStruct) );
	H5Tinsert( memtype_out, "name", HOFFSET(logicalScalarStruct, name), strtype_in);
	H5Tinsert( memtype_out, "value", HOFFSET(logicalScalarStruct, value), H5T_NATIVE_INT);

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(logicalScalarStruct) );
	H5Tinsert( filetype_out, "name", HOFFSET(logicalScalarStruct, name), strtype_in);
	H5Tinsert( filetype_out, "value", HOFFSET(logicalScalarStruct, value), H5T_NATIVE_INT);

    hid_t dset_out = H5Dcreate(outFile, datasetName, filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, logicalVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, logicalVals_in);
	H5Tclose(memtype_in);
	H5Sclose(space_in);
	H5Dclose(dset_in);
	H5Tclose(strtype_in);

	delete[] logicalVals_in;
}

void FlashAmrMesh::writeReals(const char *which){
	//This routine returns a double corresponding to the 
	// requested scalar char array argument from either the 
	// real scalar list or the real runtime parameters list. 
	// These values go into the map object realScalarsMap.

	typedef struct realScalarStruct {
		char 	name[MAX_STRING_LENGTH];
		double  value;
	} realScalarStruct;


	char datasetName[50] = {};
	if (string(which) == "scalars") {
		strcpy(datasetName,"real scalars");
	}
	else if (string(which) == "parameters") {
		strcpy(datasetName,"real runtime parameters");
	}

	if (!datasetAvailable(datasetName,"required")) {
		return;
	}

	//Get scalars list.
	hid_t dset_in = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );

	//Make specific type for the 80 character string.
	hid_t strtype_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype_in, MAX_STRING_LENGTH);

	//"real scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(realScalarStruct) );
	H5Tinsert( memtype_in, "name", HOFFSET(realScalarStruct, name), strtype_in);
	H5Tinsert( memtype_in, "value", HOFFSET(realScalarStruct, value), H5T_NATIVE_DOUBLE);

	realScalarStruct* realVals_in = new realScalarStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, realVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(realScalarStruct) );
	H5Tinsert( memtype_out, "name", HOFFSET(realScalarStruct, name), strtype_in);
	H5Tinsert( memtype_out, "value", HOFFSET(realScalarStruct, value), H5T_NATIVE_DOUBLE);

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(realScalarStruct) );
	H5Tinsert( filetype_out, "name", HOFFSET(realScalarStruct, name), strtype_in);
	H5Tinsert( filetype_out, "value", HOFFSET(realScalarStruct, value), H5T_NATIVE_DOUBLE);

    hid_t dset_out = H5Dcreate(outFile, datasetName, filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, realVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, realVals_in);
	H5Dclose(dset_in);
	H5Tclose(memtype_in);
	H5Sclose(space_in);
	H5Tclose(strtype_in);

	delete[] realVals_in;
}

void FlashAmrMesh::writeStrings(const char *which){
	//This routine returns a string corresponding to the 
	// requested scalar char array argument from either the
	// string scalar list or the string runtime parameters list. 
	// These values go into the map object stringScalarsMap.

	typedef struct stringScalarStruct {
		char name[MAX_STRING_LENGTH];
		char value[MAX_STRING_LENGTH];
	} stringScalarStruct;

	char datasetName[50] = {};
	if (string(which) == "scalars") {
		strcpy(datasetName,"string scalars");
	}
	else if (string(which) == "parameters") {
		strcpy(datasetName,"string runtime parameters");
	}

	if (!datasetAvailable(datasetName,"required")) {
		return;
	}

	//Get scalars list.
	hid_t dset_in = H5Dopen(inFile, datasetName, H5P_DEFAULT);
	hsize_t dims_in[1];
	hid_t space_in = H5Dget_space( dset_in );
	H5Sget_simple_extent_dims( space_in, dims_in, NULL );

	//Make specific type for the 80 character string.
	hid_t strtype_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strtype_in, MAX_STRING_LENGTH);
	hid_t strval_in = H5Tcopy(H5T_C_S1);
	H5Tset_size (strval_in, MAX_STRING_LENGTH);

	//"string scalars" is a compound list. Have to deal with this differently:
	//Have to create a CompType object and read from that. 
	hid_t memtype_in = H5Tcreate(H5T_COMPOUND, sizeof(stringScalarStruct) );
	H5Tinsert( memtype_in, "name", HOFFSET(stringScalarStruct, name), strtype_in);
	H5Tinsert( memtype_in, "value", HOFFSET(stringScalarStruct, value), strval_in);

	stringScalarStruct* strVals_in = new stringScalarStruct[dims_in[0]];
	H5Dread( dset_in, memtype_in, H5S_ALL, H5S_ALL, H5P_DEFAULT, strVals_in);

    hid_t space_out = H5Dget_space( dset_in );

    // Refer to tutorial h5ex_t_cmpd.c for the reasoning behind this split. Basically
    // it's a placeholder if the sizes of what's going into the file needs to be slightly
    // different from what's in memory.
    hid_t memtype_out = H5Tcreate(H5T_COMPOUND, sizeof(stringScalarStruct) );
	H5Tinsert( memtype_out, "name", HOFFSET(stringScalarStruct, name), strtype_in);
	H5Tinsert( memtype_out, "value", HOFFSET(stringScalarStruct, value), strval_in);

    hid_t filetype_out = H5Tcreate(H5T_COMPOUND, sizeof(stringScalarStruct) );
	H5Tinsert( filetype_out, "name", HOFFSET(stringScalarStruct, name), strtype_in);
	H5Tinsert( filetype_out, "value", HOFFSET(stringScalarStruct, value), strval_in);

    hid_t dset_out = H5Dcreate(outFile, datasetName, filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, strVals_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    H5Dclose(dset_out);
    H5Tclose(filetype_out);
    H5Tclose(memtype_out);
    H5Sclose(space_out);

	H5Dvlen_reclaim (memtype_in, space_in, H5P_DEFAULT, strVals_in);
	H5Dclose(dset_in);
	H5Tclose(memtype_in);
	H5Sclose(space_in);
	H5Tclose(strtype_in);
	H5Tclose(strval_in);

	delete[] strVals_in;
}

void FlashAmrMesh::writeCoordinates() {

	if (!datasetAvailable("coordinates","required")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "coordinates", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[2];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	// cout << "dims_in: " << dims_in[0] << " " << dims_in[1] << endl;

	//Create temp 2d array.
	double *tempCoordsOut = new double[nBlocksAssigned*dims_in[1]];

	for (int i = 0; i < nBlocksAssigned; i++) {
		for (int j = 0; j < (int)dims_in[1]; j++) {
			int ind = i * (int)dims_in[1] + j;
			tempCoordsOut[ind] = blocks[i].coord[j];
			// cout << tempCoordsOut[ind] << " ";
		}
		// cout << endl;
	}

	/* Output */
	hsize_t dims_full[2]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];
	dims_full[1] = dims_in[1];

	// cout << "dims_full: " << dims_full[0] << " " << dims_full[1] << endl;
	

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(2, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, "coordinates", H5T_NATIVE_DOUBLE, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[2]; //hyperslab offset in new file
	hsize_t count_out[2];  //size of hyperslab
	offset_out[0] = mpiGetGlobalBlockID(0);
	offset_out[1] = 0;
	count_out[0] = nBlocksAssigned;
	count_out[1] = dims_in[1];

	// cout << "offset_out: " << offset_out[0] << " " << offset_out[1] << endl;
	// cout << "count_out: " << count_out[0]  << " " << count_out[1]  << endl;

	hid_t memspace_out = H5Screate_simple(2, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	herr_t status = H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

	// std::cout << "Status: " << status << std::endl;

	// MPI_Barrier(MPI_COMM_WORLD);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // std::cout << "Status 2: " << status << std::endl;

    // Write hyperslab to dataset
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_out, filespace_out,
	  	    	      plist_id, tempCoordsOut);

    // cout << "status 3: " << status << endl;

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempCoordsOut;

}

void FlashAmrMesh::writeCPUNumbers() {

	if (!datasetAvailable("processor number","optional")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "processor number", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[1];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	//Create temp 2d array.
	int *tempProcessorOut = new int[nBlocksAssigned];

	for (int i = 0; i < nBlocksAssigned; i++) {
		tempProcessorOut[i] = blocks[i].cpuNumber;
		// cout << tempProcessorOut[i] << " ";
	}
	// cout << endl;

	/* Output */
	hsize_t dims_full[1]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(1, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, "processor number", H5T_NATIVE_INT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[1]; //hyperslab offset in new file
	hsize_t count_out[1];  //size of hyperslab
	offset_out[0] = mpiGetGlobalBlockID(0);
	count_out[0] = nBlocksAssigned;

	hid_t memspace_out = H5Screate_simple(1, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write hyperslab to dataset
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_out, filespace_out,
	  	    	      plist_id, tempProcessorOut);

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempProcessorOut;

}

void FlashAmrMesh::writeSizes() {

	if (!datasetAvailable("block size","optional")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "block size", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[2];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	//Create temp 2d array.
	double *tempSizesOutD = new double[nBlocksAssigned*dims_in[1]];
	float *tempSizesOutF = new float[nBlocksAssigned*dims_in[1]];

	if (fileType == 7) {
		for (int i = 0; i < nBlocksAssigned; i++) {
			for (int j = 0; j < (int)dims_in[1]; j++) {
				tempSizesOutF[(int)dims_in[1]*i + j] = (float)blocks[i].size[j];
				// cout << tempSizesOutF[(int)dims_in[1]*i + j] << " ";
			}
			// cout << endl;
		}
	}
	else {
		for (int i = 0; i < nBlocksAssigned; i++) {
			for (int j = 0; j < (int)dims_in[1]; j++) {
				tempSizesOutD[(int)dims_in[1]*i + j] = blocks[i].size[j];
				// cout << tempSizesOutD[(int)dims_in[1]*i + j] << " ";
			}
			// cout << endl;
		}
	}

	/* Output */
	hsize_t dims_full[2]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];
	dims_full[1] = dims_in[1];

	// Create a data space that expects the total amount of data across all processors.
	hid_t dset_id;
	hid_t filespace;
	if (fileType == 7) {
		filespace = H5Screate_simple(2, dims_full, NULL);
    	dset_id = H5Dcreate(outFile, "block size", H5T_NATIVE_FLOAT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	else {
		filespace = H5Screate_simple(2, dims_full, NULL);
   		dset_id = H5Dcreate(outFile, "block size", H5T_NATIVE_DOUBLE, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	H5Sclose(filespace);
    
	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[2]; //hyperslab offset in new file
	hsize_t count_out[2];  //size of hyperslab
	offset_out[0] = mpiGetGlobalBlockID(0);
	offset_out[1] = 0;
	count_out[0] = nBlocksAssigned;
	count_out[1] = dims_in[1];

	hid_t memspace_out = H5Screate_simple(2, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write hyperslab to dataset
    if (fileType == 7) {
    	H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace_out, filespace_out,
	  	        	      plist_id, tempSizesOutF);    	
    }
    else {
    	H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_out, filespace_out,
	  		        	  plist_id, tempSizesOutD);
    }

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempSizesOutD;
	delete[] tempSizesOutF;

}

void FlashAmrMesh::writeBounds() {

	if (!datasetAvailable("bounding box","required")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "bounding box", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[3];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	//Create temp 2d array.
	double *tempBoundsOut = new double[nBlocksAssigned*dims_in[1]*dims_in[2]];

	for (int i = 0; i < nBlocksAssigned; i++) {
		for (int j = 0; j < (int)dims_in[1]; j++) {
			int index = (i*(int)dims_in[1]+j)*(int)dims_in[2];
			tempBoundsOut[index] = blocks[i].lowerBound[j];
			tempBoundsOut[index + 1] = blocks[i].upperBound[j];
			// cout << tempBoundsOut[index] << " " << tempBoundsOut[index + 1] << " ";
		}
		// cout << endl;
	}

	/* Output */
	hsize_t dims_full[3]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];
	dims_full[1] = dims_in[1];
	dims_full[2] = dims_in[2];

	// cout << "dims_full: " << dims_full[0] << " " << dims_full[1] << endl;
	

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(3, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, "bounding box", H5T_NATIVE_DOUBLE, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[3]; //hyperslab offset in new file
	hsize_t count_out[3];  //size of hyperslab
	offset_out[0] = mpiGetGlobalBlockID(0);
	offset_out[1] = 0;
	offset_out[2] = 0;
	count_out[0] = nBlocksAssigned;
	count_out[1] = dims_in[1];
	count_out[2] = dims_in[2];

	hid_t memspace_out = H5Screate_simple(3, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

	// std::cout << "Status: " << status << std::endl;

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // std::cout << "Status 2: " << status << std::endl;

    // Write hyperslab to dataset
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_out, filespace_out,
	  	    	      plist_id, tempBoundsOut);

    // cout << "status 3: " << status << endl;

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempBoundsOut;

}

void FlashAmrMesh::writeNodeType() {

	if (!datasetAvailable("node type","required")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "node type", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[1];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	//Create temp 2d array.
	int *tempNodeOut = new int[nBlocksAssigned];

	for (int i = 0; i < nBlocksAssigned; i++) {
		tempNodeOut[i] = blocks[i].nodeType;
		// cout << tempNodeOut[i] << " ";
	}
	// cout << endl;

	/* Output */
	hsize_t dims_full[1]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(1, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, "node type", H5T_NATIVE_INT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[1]; //hyperslab offset in new file
	hsize_t count_out[1];  //size of hyperslab
	offset_out[0] = mpiGetGlobalBlockID(0);
	count_out[0] = nBlocksAssigned;

	hid_t memspace_out = H5Screate_simple(1, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write hyperslab to dataset
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_out, filespace_out,
	  	    	      plist_id, tempNodeOut);

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempNodeOut;

}

void FlashAmrMesh::writeRefineLevel() {

	if (!datasetAvailable("refine level","required")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "refine level", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[1];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	//Create temp 2d array.
	int *tempRefineOut = new int[nBlocksAssigned];

	for (int i = 0; i < nBlocksAssigned; i++) {
		tempRefineOut[i] = blocks[i].refineLevel;
		// cout << tempRefineOut[i] << " ";
	}
	// cout << endl;

	/* Output */
	hsize_t dims_full[1]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(1, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, "refine level", H5T_NATIVE_INT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[1]; //hyperslab offset in new file
	hsize_t count_out[1];  //size of hyperslab
	offset_out[0] = mpiGetGlobalBlockID(0);
	count_out[0] = nBlocksAssigned;

	hid_t memspace_out = H5Screate_simple(1, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write hyperslab to dataset
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_out, filespace_out,
	  	    	      plist_id, tempRefineOut);

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempRefineOut;

}

void FlashAmrMesh::writeGid() {

	if (!datasetAvailable("gid","required")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "gid", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[2];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	//Create temp 2d array.
	int *tempGidOut = new int[nBlocks*dims_in[1]];

	for (int i = 0; i < nBlocks; i++) {
		for (int j = 0; j < (int)dims_in[1]; j++) {
			int ind = (int)dims_in[1]*i + j;
			tempGidOut[ind] = gid[i][j];
			// cout << tempGidOut[ind] << " ";
		}
		// cout << endl;
	}

	/* Output */
	hsize_t dims_full[2]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];
	dims_full[1] = dims_in[1];

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(2, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, "gid", H5T_NATIVE_INT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[2]; //hyperslab offset in new file
	hsize_t count_out[2];  //size of hyperslab
	offset_out[0] = 0;
	offset_out[1] = 0;
	count_out[0] = nBlocks;
	count_out[1] = dims_in[1];

	hid_t memspace_out = H5Screate_simple(2, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    // Write hyperslab to dataset
    if (mpiGetID() == 0) {
	    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_out, filespace_out,
		  	    	      plist_id, tempGidOut);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempGidOut;

}

void FlashAmrMesh::writeWhichChild() {
	// Reads in the 'Which Child' information and directly writes it to
	// the output file. This information is never read into the AMR Mesh
	// object itself.

	if (!datasetAvailable("which child","optional")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "which child", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[1];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	//MPI stuff
	hsize_t offset_in[1]; //hyperslab offset in file
	hsize_t count_in[1];  //size of hyperslab

	offset_in[0] = mpiGetGlobalBlockID(0);
	count_in[0]  = nBlocksAssigned;

	H5Sselect_hyperslab(memspace_in, H5S_SELECT_SET, offset_in, NULL, count_in, NULL);

	hid_t memspace_mem = H5Screate_simple(1, dims_in, NULL);

    hsize_t      offset_mem[1];   // hyperslab offset in memory
    hsize_t      count_mem[1];    // size of the hyperslab in memory
    offset_mem[0] = 0;
    count_mem[0]  = nBlocksAssigned;

    H5Sselect_hyperslab(memspace_mem, H5S_SELECT_SET, offset_mem, NULL, count_mem, NULL);

	//Create temp 2d array.
	int *tempWhichChildOut = new int[nBlocksAssigned];

	H5Dread( dset_in, H5T_NATIVE_INT, memspace_mem, memspace_in, H5P_DEFAULT, tempWhichChildOut);

	/* Output */
	hsize_t dims_full[1]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(1, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, "which child", H5T_NATIVE_INT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[1]; //hyperslab offset in new file
	hsize_t count_out[1];  //size of hyperslab
	offset_out[0] = mpiGetGlobalBlockID(0);
	count_out[0] = nBlocksAssigned;

	hid_t memspace_out = H5Screate_simple(1, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write hyperslab to dataset
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_out, filespace_out,
	  	    	      plist_id, tempWhichChildOut);

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_mem);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempWhichChildOut;	
}

void FlashAmrMesh::writeBFlags() {
	// Reads in the 'bflags' information and directly writes it to
	// the output file. This information is never read into the AMR Mesh
	// object itself.

	if (!datasetAvailable("bflags","optional")) {
		return;
	}

	hid_t dset_in = H5Dopen(inFile, "bflags", H5P_DEFAULT);
	hid_t memspace_in = H5Dget_space( dset_in );
	hsize_t dims_in[2];

	H5Sget_simple_extent_dims( memspace_in, dims_in, NULL );

	//MPI stuff
	hsize_t offset_in[2]; //hyperslab offset in file
	hsize_t count_in[2];  //size of hyperslab

	offset_in[0] = mpiGetGlobalBlockID(0);
	offset_in[1] = 0;
	count_in[0]  = nBlocksAssigned;
	count_in[1]  = dims_in[1];

	H5Sselect_hyperslab(memspace_in, H5S_SELECT_SET, offset_in, NULL, count_in, NULL);

	hid_t memspace_mem = H5Screate_simple(2, dims_in, NULL);

    hsize_t      offset_mem[2];   // hyperslab offset in memory
    hsize_t      count_mem[2];    // size of the hyperslab in memory
    offset_mem[0] = 0;
    offset_mem[1] = 0;
    count_mem[0]  = nBlocksAssigned;
    count_mem[1]  = dims_in[1];

    H5Sselect_hyperslab(memspace_mem, H5S_SELECT_SET, offset_mem, NULL, count_mem, NULL);

	//Create temp 2d array.
	int *tempBFlagsOut = new int[nBlocksAssigned*1];

	H5Dread( dset_in, H5T_NATIVE_INT, memspace_mem, memspace_in, H5P_DEFAULT, tempBFlagsOut);

	/* Output */
	hsize_t dims_full[2]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];
	dims_full[1] = dims_in[1];

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(2, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, "bflags", H5T_NATIVE_INT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	// Prepare the hyperslab containing this processor's data
	hsize_t offset_out[2]; //hyperslab offset in new file
	hsize_t count_out[2];  //size of hyperslab
	offset_out[0] = mpiGetGlobalBlockID(0);
	offset_out[1] = 0;
	count_out[0] = nBlocksAssigned;
	count_out[1] = dims_in[1];

	hid_t memspace_out = H5Screate_simple(2, count_out, NULL);
    
	hid_t filespace_out = H5Dget_space(dset_id);
	H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
		                                NULL, count_out, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write hyperslab to dataset
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace_out, filespace_out,
	  	    	      plist_id, tempBFlagsOut);

    H5Pclose(plist_id);
	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Dclose(dset_id);
    H5Sclose(memspace_mem);
    H5Sclose(memspace_in);	
	H5Dclose(dset_in);	

	delete[] tempBFlagsOut;	
}

void FlashAmrMesh::writeNFileVarsList() {
	// Writes the 'unknown names' list in the HDF5 file.

	// No dataset checking: No need for checking if there's a corresponding variable in the input file.

	// Count the number of variables to be written from the writeVarMap map.
	int numVarsToWrite = 0;
	for (int i = 0; i < nUserVars; i++) {
		if (writeVarMap[uVarsMap_iToS[i]]) {
			numVarsToWrite++;
		}
	}

	//Create temp 1d array.
	char (*tempVarOut)[MESH_VAR_STRING_SIZE] = new char[numVarsToWrite][MESH_VAR_STRING_SIZE];
	
	int count = 0;
	for (int i = 0; i < nUserVars; i++) {
		if (writeVarMap[uVarsMap_iToS[i]]) {
			strcpy(tempVarOut[count],uVarsMap_iToS[i].c_str());
			count++;
		}
	}

	//Get dimensions of dataset.
	hsize_t dims_out[2];
	dims_out[0] = numVarsToWrite;
	dims_out[1] = 1;

	hid_t filetype_out = H5Tcopy(H5T_C_S1);
	H5Tset_size (filetype_out, MESH_VAR_STRING_SIZE);
	hid_t memtype_out = H5Tcopy(H5T_C_S1);
	H5Tset_size (memtype_out, MESH_VAR_STRING_SIZE);

	hid_t space_out = H5Screate_simple (2, dims_out, NULL);

    hid_t dset_out = H5Dcreate(outFile, "unknown names", filetype_out, space_out, H5P_DEFAULT,
    	                       H5P_DEFAULT, H5P_DEFAULT);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mpiGetID() == 0) {
    	H5Dwrite(dset_out, memtype_out, H5S_ALL, H5S_ALL, plist_id, tempVarOut[0]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

	H5Dclose(dset_out);
	H5Sclose(space_out);
	H5Tclose(memtype_out);
	H5Tclose(filetype_out);

	delete[] tempVarOut;

}

void FlashAmrMesh::writeVarData(const char *outputVariable) {

	// No dataset checking: No need for checking if there's a corresponding variable in the input file.

	hsize_t dims_in[4];
	dims_in[0] = getNBlocks();
	dims_in[1] = getNcells(KAXIS);
	dims_in[2] = getNcells(JAXIS);
	dims_in[3] = getNcells(IAXIS);

	/* Output */
	hsize_t dims_full[4]; //Sum total size of the dataset across all processors.
	dims_full[0] = dims_in[0];
	dims_full[1] = dims_in[1];
	dims_full[2] = dims_in[2];
	dims_full[3] = dims_in[3];

	// Create a data space that expects the total amount of data across all processors.
	hid_t filespace = H5Screate_simple(4, dims_full, NULL);
    hid_t dset_id = H5Dcreate(outFile, outputVariable, H5T_NATIVE_DOUBLE, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

	//Create temp 2d array.
	double *tempVarOut = new double[1*dims_in[1]*dims_in[2]*dims_in[3]];

	hid_t memspace_out;
	hid_t filespace_out;
	hid_t plist_id;

	hsize_t count_out[4];    // size of the hyperslab in memory
	hsize_t offset_out[4];   // hyperslab offset in memory
	hsize_t dims_out[4];

	count_out[0] = 1; //nBlocksAssigned;
	count_out[1] = dims_in[1];
	count_out[2] = dims_in[2];
	count_out[3] = dims_in[3];

	dims_out[0] = 1;
	dims_out[1] = dims_in[1];
	dims_out[2] = dims_in[2];
	dims_out[3] = dims_in[3];

	double maxVarVal = -1e308;
	double minVarVal = 1e308;
    for (int blk = 0; blk < nBlocksAssigned; blk++) {
    	// std::cout << "Blocks assigned: " << nBlocksAssigned << std::endl;
		// std::cout << "Current block: " << blk << std::endl;

		// Prepare the hyperslab containing this processor's data
		offset_out[0] = mpiGetGlobalBlockID(blk);
		offset_out[1] = 0;
		offset_out[2] = 0;
		offset_out[3] = 0;

		//Read temp data into corresponding block coord arrays.
		int index;
		int var = uVarsMap_sToI[string(outputVariable)];
		for (int i = 0; i < (int)dims_out[3]; i++) {
			for (int j = 0; j < (int)dims_out[2]; j++) {
				for (int k = 0; k < (int)dims_out[1]; k++) {
					index = k * dims_out[2]*dims_out[3] + j * dims_out[3] + i;
					// std::cout << "i: " << i << " j: " << j << " k: " << k << std::endl;
					tempVarOut[index] = blocks[blk].data[var][i][j][k];
					// std::cout << "Block: " << mpiGetGlobalBlockID(blk) 
							  // << " Index: " << index 
							  // << " Value: " << tempVarOut[index] 
							  // << std::endl;
					maxVarVal = max(maxVarVal,tempVarOut[index]);
					minVarVal = min(minVarVal,tempVarOut[index]);
				}
			}
		}
		
		memspace_out = H5Screate_simple(4, count_out, NULL);
	    
		filespace_out = H5Dget_space(dset_id);
		H5Sselect_hyperslab(filespace_out, H5S_SELECT_SET, offset_out, 
			                                NULL, count_out, NULL);

	    plist_id = H5Pcreate(H5P_DATASET_XFER);
	    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

	    // Write hyperslab to dataset
	    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_out, filespace_out,
		  	    	      plist_id, tempVarOut);

	}
	double maxVarValAll;
	double minVarValAll;

	/* Attributes must be written collectively, so this is done on all
	   processors.
	*/
	MPI_Allreduce(&maxVarVal, &maxVarValAll, 1, MPI_DOUBLE, MPI_MAX, comm);
	MPI_Allreduce(&minVarVal, &minVarValAll, 1, MPI_DOUBLE, MPI_MIN, comm);

	/* Write the min/max attributes 
	   Note: This will take the max and min of ALL data, so if
	   only leaf blocks are written to in main, the min or max
	   may end up being 0 (the initialized value of a new variable). 
	*/	
	hsize_t dims_attr = 1;

	hid_t attrspace_id = H5Screate_simple(1, &dims_attr, NULL);
	hid_t attribute_id = H5Acreate(dset_id, "maximum", H5T_NATIVE_DOUBLE, attrspace_id,
		                      H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &maxVarValAll);
	

	H5Aclose(attribute_id);
	H5Sclose(attrspace_id);

	attrspace_id = H5Screate_simple(1, &dims_attr, NULL);
	attribute_id = H5Acreate(dset_id, "minimum", H5T_NATIVE_DOUBLE, attrspace_id,
		                      H5P_DEFAULT, H5P_DEFAULT);
	
	H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &minVarValAll);
	

	H5Aclose(attribute_id);
	H5Sclose(attrspace_id);


	H5Sclose(filespace_out);
	H5Sclose(memspace_out);
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    


	delete[] tempVarOut;

}

void FlashAmrMesh::writeVariables() {
	// Only write the variables the user wants to have written.
	// Default behavior is to write all variables, but the user
	// can use the noWrite() function in main to switch from writing
	// a variable to not writing a variable.

	for (int i = 0; i < nUserVars; i++) {
		if (writeVarMap[variableNames[i]]) {
			if (!SUPPRESSOUTPUT) cout << "Writing variable " << variableNames[i] << "...";
			writeVarData(variableNames[i]);
			if (!SUPPRESSOUTPUT) cout << "Variable " << variableNames[i] << "written." << endl;
		}
	}
}


void FlashAmrMesh::generalWrite() {
	// Writes variables similar to how they are read in the generalRead function.

	if (!SUPPRESSOUTPUT) cout << "Writing scalar lists. " << endl;
	if (fileType == 7) {
		if (!SUPPRESSOUTPUT) cout << "File type 7: Writing file type." << endl;
		this->writeFileType();
		if (!SUPPRESSOUTPUT) cout << "File type 7: Writing Simulation Parameters." << endl;
		this->writeSimParameters();
	}
	else {
		this->writeInts("scalars");
		this->writeStrings("scalars");
		this->writeReals("scalars");
		this->writeLogicals("scalars");
	}

	if (!SUPPRESSOUTPUT) cout << "Reading parameter lists. " << endl;
	this->writeInts("parameters");
	this->writeStrings("parameters");
	this->writeReals("parameters");
	this->writeLogicals("parameters");		

	this->writeNFileVarsList();
	if (!SUPPRESSOUTPUT) cout << "List of variable names written. " << endl;
	if (!SUPPRESSOUTPUT) cout << "Number of variables in file: " << this->getNUserVars() << endl;

	//TODO: This might need to only be in file type 9, or at least not in file type 7.
	//      See if this is present in all files that aren't type 7. If so, write it to the fileType variable.
	//      Make a function to check for the presence of datasets in the file.
	this->writeSimInfo();
	if (!SUPPRESSOUTPUT) cout << "Simulation info written. " << endl;

	/* These dataset are required in order for VisIt to be able to
	   read the HDF5 file. If they cannot be written, the code WILL crash
	   (see noWrite function, which is called at the top of each
	   of the below functions).
	*/

	this->writeCoordinates();
	if (!SUPPRESSOUTPUT) cout << "Coordinates written. " << endl;

	this->writeBounds();
	if (!SUPPRESSOUTPUT) cout << "Block bounds written. " << endl;

	this->writeNodeType();
	if (!SUPPRESSOUTPUT) cout << "Node types written. " << endl;

	this->writeRefineLevel();
	if (!SUPPRESSOUTPUT) cout << "Refine levels written. " << endl;

	this->writeGid();
	if (!SUPPRESSOUTPUT) cout << "Gid array written. " << endl;

	/* These datasets are NOT required for VisIt, but are written anyways.
	   If they cannot be written, the code does not crash (see noWrite function).
	*/

	this->writeCPUNumbers();
	if (!SUPPRESSOUTPUT) cout << "CPU numbers written. " << endl;

	this->writeSizes();
	if (!SUPPRESSOUTPUT) cout << "Physical sizes of blocks written. " << endl;

	this->writeWhichChild();
	if (!SUPPRESSOUTPUT) cout << "Which Child array written." << endl;

	this->writeBFlags();
	if (!SUPPRESSOUTPUT) cout << "BFlags array written." << endl;

	// TODO: Find a file type 7 flash file to test those functions.
}

void FlashAmrMesh::closeOutFile() {
	if (!SUPPRESSOUTPUT) cout << "Trying to close output file...";
	H5Fclose(outFile);
	if (!SUPPRESSOUTPUT) cout << "Closed output file." << endl;
}

void FlashAmrMesh::writeOutFile() {
	if (mpiGetID() == 0) {
		std::cout << "Writing file" << std::endl;
	}	
	createOutFile();
	generalWrite();
	writeVariables();
	closeOutFile();
	if (mpiGetID() == 0) {
		std::cout << "Writing complete" << std::endl;
	}	
}

//Tim's functions.

void FlashAmrMesh::fillGuardCells() {
	if (USEGUARDCELLS) {
		cout << "[flash_reader] Reader does not yet support fillGuardCells function." << endl;
		exit(EXIT_FAILURE);
	}
	else {
		return;
	}
}

void FlashAmrMesh::getBlkBoundBox(int blockId, double boundBox[2][MESH_MDIM]) {
	//  Returns bounding box for a certain block similar to how FLASH4 returns the bounding box.
	//
	//  Inputs: 
	//    blockId: id of the block (usually an array index)
	//
	//  Outputs:
	//    boundBox: 2D array returning the lower and upper bounds of the box.
	//
	// 
	//  Example:
	//    For 2D block with bounding coordinates (-0.5,0.0) for the bottom left corner
	//    and (0.5, 1.0) for the top right corner, the function will return 
    //      boundBox(LOWER, IAXIS) = -0.5
    //      boundBox(UPPER, IAXIS) = 0.5
    //      boundBox(LOWER, JAXIS) = 0.0
    //      boundBox(UPPER, JAXIS) = 1.0
    //      boundBox(LOWER, KAXIS) = 1 !returned as 1 because only 2 dims
    //      boundBox(UPPER, KAXIS) = 1 !returned as 1 because only 2 dims
    //    exactly as is specified in the FLASH API for this function, only with the added
    //    enums LOWER and UPPER.
	//
	//  For quick reference, enums are LOWER=0, UPPER=1, IAXIS=0, JAXIS=1, KAXIS=2
	//
	//  Example code in Main, using the same variables as the example at the top of this file:
	//
	/*	
		double blkBoundBox[2][MESH_MDIM];
		
		//Get physical size of the first block's bounding box and put it into blkBoundBox array.
		newGrid.getBlkBoundBox(0, blkBoundBox);
		
		std::cout << "Block bound box [0][0]: " << blkBoundBox[0][0] << std::endl;
		std::cout << "Block bound box [1][0]: " << blkBoundBox[1][0] << std::endl;
		std::cout << "Block bound box [0][1]: " << blkBoundBox[0][1] << std::endl;
		std::cout << "Block bound box [1][1]: " << blkBoundBox[1][1] << std::endl;
		std::cout << "Block bound box [0][2]: " << blkBoundBox[0][2] << std::endl;
		std::cout << "Block bound box [1][2]: " << blkBoundBox[1][2] << std::endl;
	*/

	
	if (MESH_MDIM >= 1) {
		boundBox[LOWER][IAXIS] = getBound(blockId, LOWER, IAXIS);
		boundBox[UPPER][IAXIS] = getBound(blockId, UPPER, IAXIS);
		for (int i = JAXIS; i <= KAXIS; i++) { //Fill the rest with 1's if they are not present.
			boundBox[LOWER][i] = 1;
			boundBox[UPPER][i] = 1;
		}
	}
	if (MESH_MDIM >= 2) {
		boundBox[LOWER][JAXIS] = getBound(blockId, LOWER, JAXIS);
		boundBox[UPPER][JAXIS] = getBound(blockId, UPPER, JAXIS);
		for (int i = KAXIS; i <= KAXIS; i++) {
			boundBox[LOWER][i] = 1;
			boundBox[UPPER][i] = 1;
		}
	}
	if (MESH_MDIM == 3) {
		boundBox[LOWER][KAXIS] = getBound(blockId, LOWER, KAXIS);
		boundBox[UPPER][KAXIS] = getBound(blockId, UPPER, KAXIS);
	}
	else{
		cout << "[flash_reader] Dimension used in getBlkBoundBox (MESH_MDIM) is not 1, 2, or 3. Exiting routine." << endl;
		exit(EXIT_FAILURE);
	}
}

void FlashAmrMesh::getBlkData(int blockID, GRIDDATASTRUCT gridDataStruct, const char structIndex[MESH_VAR_STRING_SIZE+1], 
	                     BEGINCOUNT beginCount, int startingPos[MESH_MDIM], double****& dataBlock 
	                     /* , int *dataSize */) {
	//Returns block data similar to how FLASH 4 returns block data. (Note that dataSize is not 
	// used, since this function passes a referenced block struct back to the caller, not an 
	// array of doubles like in FLASH 4. It will always read all of the data available for the 
	// desired block (except for guard cells), making dataSize unnecessary.)

	// Inputs:
	//		blockID:  id of the block to have variable data read into it.
	//
	//      gridDataStruct: Location of the data on the cells. 
	//                      -->Not yet implemented (defaults to CENTER).
	//
	//      structIndex: The variable being read in to dataBlock. A character array is used to
	//                   specify the variable. Format must be lowercase with only four characters. 
	//                   Example: "pres" or "dens".
	//
	//      beginCount:  Tells the reader whether or not guard cells should be included in the data read in.
	//                   -->Not yet implemented (defaults to INTERIOR since guard cells are not yet implemented). 
	//
	//      startingPos: MESH_MDIM size array containing the starting indexes of where to read in data.
	//                   This is basically implemented, but for now it defaults to (0,0,0).
	//
	// Outputs:
	//      dataBlock:   Data array that the data will be read into, passed by reference.
	//

	// Example: 
	/*
		const H5std_string FILE_NAME( "snb_hdf5_chk_0000" );

		int main() {

			FlashAmrMesh newGrid;


			newGrid.setGridFileName( FILE_NAME );
			
			newGrid.openFile();

			newGrid.generalRead();

			//Read all the variables
			//newGrid.readAllVariables();

			//Read only pressure and density, for example.
			newGrid.addUserVariable("pres");
			newGrid.addUserVariable("dens");
			newGrid.readPartialVariables();

			int numberOfBlocks = newGrid.getNBlocks();

			//Allocate data arrays for each block.
			double***** test;
			test = new double****[numberOfBlocks];

			for (int i = 0; i < numberOfBlocks; i++) {
				allocateVarArray(test[i], newGrid.getNFileVars(), newGrid.getNcells() );	
			}

			int startingPos[MESH_MDIM];
			startingPos[0] = 0;
			startingPos[1] = 0;
			startingPos[2] = 0;

			//Retrieve the data from each block.
			for (int i = 0; i < numberOfBlocks; i++) {
				newGrid.getBlkData(i, CENTER, "pres", INTERIOR, startingPos, test[i]);
			}

			//Confirm that something was read into the 0,0,0 cell for the corresponding 
			//pressue index
			int varIndex = newGrid.findVarIndex("pres");
			for (int i = 0; i < numberOfBlocks; i++) {
				std::cout << test[i][varIndex][0][0][0] << std::endl;
			}

			for (int i = 0; i < numberOfBlocks; i++) {
				deallocateVarArray(test[i], newGrid.getNFileVars(), newGrid.getNcells() );
			}

			newGrid.closeFile();

			return 0;

		}
	*/


	if (gridDataStruct != CENTER) {
		cout << "[flash_reader] Function getBlkData can only use the grid data struct CENTER currently. Defaulting to CENTER." << endl;
		gridDataStruct = CENTER;
	}

	if (startingPos[0] != 0 || startingPos[1] != 0 || startingPos[2] != 0) {
		cout << "[flash_reader] Function getBlkData can only use startingPos[] = 0 currently. Defaulting to 0. " << endl;
		startingPos[0] = 0;
		startingPos[1] = 0;
		startingPos[2] = 0;
	}

	if (beginCount != INTERIOR) {
		cout << "[flash_reader] Function getBlkData does not yet support counting guard cells, so only INTERIOR can be used. Defaulting to INTERIOR. " << endl;
		beginCount = INTERIOR;
	}


	int variableIndex = findVarIndex(structIndex);
	if (variableIndex == -1) {
		cout << "[flash_reader] Requested variable in getBlkData function not found. Exiting routine. " << endl;
		exit(EXIT_FAILURE);
	}

	//This conditional will be updated if there is a need to include functionality for 
	// different gridDataStructs and guard cells. Right now the conditional is meaningless;
	// it is more of a placeholder. These conditions would be required for this specific call, though.
	if (gridDataStruct == CENTER && beginCount == INTERIOR) {
		for (int k = startingPos[2]; k < getNcells(KAXIS); k++) {
			for (int j = startingPos[1]; j < getNcells(JAXIS); j++) {
				for (int i = startingPos[0]; i < getNcells(IAXIS); i++) {
					dataBlock[variableIndex][i][j][k] = getVariable(blockID, variableIndex, i, j, k);
				}
			}
		}
	}
}

void FlashAmrMesh::getBlkData(int blockID, GRIDDATASTRUCT gridDataStruct, int structIndex, 
	                     BEGINCOUNT beginCount, int startingPos[MESH_MDIM], double****& dataBlock 
	                     /* , int *dataSize */) {
	//Returns block data similar to how FLASH 4 returns block data. (Note that dataSize is not 
	// used, since this function passes a referenced block struct back to the caller, not an 
	// array of doubles like in FLASH 4. It will always read all of the data available for the 
	// desired block (except for guard cells), making dataSize unnecessary.)

	// Inputs:
	//		blockID:  id of the block to have variable data read into it.
	//
	//      gridDataStruct: Location of the data on the cells. 
	//                      -->Not yet implemented (defaults to CENTER).
	//
	//      structIndex: Integer index specifyig the variable being read in to dataBlock. 
	//
	//      beginCount:  Tells the reader whether or not guard cells should be included in the data read in.
	//                   -->Not yet implemented (defaults to INTERIOR since guard cells are not yet implemented). 
	//
	//      startingPos: MESH_MDIM size array containing the starting indexes of where to read in data.
	//                   This is basically implemented, but for now it defaults to (0,0,0).
	//
	// Outputs:
	//      dataBlock:   Data array that the data will be read into, passed by reference.
	//


	if (gridDataStruct != CENTER) {
		cout << "[flash_reader] Function getBlkData can only use the grid data struct CENTER currently. Defaulting to CENTER." << endl;
		gridDataStruct = CENTER;
	}

	if (startingPos[0] != 0 || startingPos[1] != 0 || startingPos[2] != 0) {
		cout << "[flash_reader] Function getBlkData can only use startingPos[] = 0 currently. Defaulting to 0. " << endl;
		startingPos[0] = 0;
		startingPos[1] = 0;
		startingPos[2] = 0;
	}

	if (beginCount != INTERIOR) {
		cout << "[flash_reader] Function getBlkData does not yet support counting guard cells, so only INTERIOR can be used. Defaulting to INTERIOR. " << endl;
		beginCount = INTERIOR;
	}


	if (structIndex == -1) {
		cout << "[flash_reader] Requested variable in getBlkData function not found. Exiting routine. " << endl;
		exit(EXIT_FAILURE);
	}

	//This conditional will be updated if there is a need to include functionality for 
	// different gridDataStructs and guard cells. Right now the conditional is meaningless;
	// it is more of a placeholder. These conditions would be required for this specific call, though.
	if (gridDataStruct == CENTER && beginCount == INTERIOR) {
		for (int k = startingPos[2]; k < getNcells(KAXIS); k++) {
			for (int j = startingPos[1]; j < getNcells(JAXIS); j++) {
				for (int i = startingPos[0]; i < getNcells(IAXIS); i++) {
					dataBlock[structIndex][i][j][k] = getVariable(blockID, structIndex, i, j, k);
				}
			}
		}
	}
}

void FlashAmrMesh::getBlkIndexLimits(int blockID, int blkLimits[2][MESH_MDIM], int blkLimitsGC[2][MESH_MDIM], GRIDDATASTRUCT gridDataStruct) {
	// Returns the block index limits for a given block. Returns both internal block limits as well as 
	// external, which include guard cells if they are used. 
	// Set whether or not guard cells are used by changing USEGUARDCELLS as well as NGUARD 
	// in the header file FlashAmrMesh.h.

	// Inputs:
	//     blockID: index of the block in the blocks array.
	//
	//     gridDataStruct: Location of the data on the cells. 
	//                  -->Not yet implemented (defaults to CENTER).
	//
	// Outputs:
	//     blkLimits[2][MESH_MDIM]: Internal block index limits 
	//
	//     blkLimitsGC[2][MESH_MDIM]: External block index limits (The same as internal if no guard cells are used)
	//

	/*	Example code in Main, using the same variables as the example at the top of this file:

		int blkBoundBoxInt[2][MESH_MDIM];
		int blkBoundBoxGCint[2][MESH_MDIM];

		newGrid.getBlkIndexLimits(0, blkBoundBoxInt, blkBoundBoxGCint, CENTER);
		std::cout << "Block index bound box [0][0]: " << blkBoundBoxInt[0][0] << std::endl;
		std::cout << "Block index bound box [1][0]: " << blkBoundBoxInt[1][0] << std::endl;
		std::cout << "Block index bound box [0][1]: " << blkBoundBoxInt[0][1] << std::endl;
		std::cout << "Block index bound box [1][1]: " << blkBoundBoxInt[1][1] << std::endl;
		std::cout << "Block index bound box [0][2]: " << blkBoundBoxInt[0][2] << std::endl;
		std::cout << "Block index bound box [1][2]: " << blkBoundBoxInt[1][2] << std::endl;

		std::cout << "GC Block index bound box [0][0]: " << blkBoundBoxInt[0][0] << std::endl;
		std::cout << "GC Block index bound box [1][0]: " << blkBoundBoxInt[1][0] << std::endl;
		std::cout << "GC Block index bound box [0][1]: " << blkBoundBoxInt[0][1] << std::endl;
		std::cout << "GC Block index bound box [1][1]: " << blkBoundBoxInt[1][1] << std::endl;
		std::cout << "GC Block index bound box [0][2]: " << blkBoundBoxInt[0][2] << std::endl;
		std::cout << "GC Block index bound box [1][2]: " << blkBoundBoxInt[1][2] << std::endl;

	*/

	int GX, GY, GZ;

	if (gridDataStruct != CENTER) {
		cout << "[flash_reader] Function getBlkIndexLimits can only use the grid data struct CENTER currently. Defaulting to CENTER." << endl;
		gridDataStruct = CENTER;
	}

	if (!USEGUARDCELLS) {
		GX = 0;
		GY = 0;
		GZ = 0;
	}
	else {
		GX = NGUARD;
		GY = NGUARD;
		GZ = NGUARD;
	}

	if (MESH_MDIM >= 1) {
  		blkLimits[LOWER][IAXIS] = GX; //lower bound index of the first interior cell of a block in xdir
  		blkLimits[UPPER][IAXIS] = GX + getNcells(IAXIS) - 1; //the upper bound index of the last interior cell of a block in xdir
  	
  		blkLimitsGC[LOWER][IAXIS] = 0; //lower bound index of the first cell in the entire block in xdir
  		blkLimitsGC[UPPER][IAXIS] = getNcells(IAXIS) + 2*GX; //upper bound index of the last cell in the entire block in xdir
		
		for (int i = JAXIS; i <= KAXIS; i++) { //Fill the rest with 1's if they are not present.
			blkLimits[LOWER][i] = 1;
			blkLimits[UPPER][i] = 1;

			blkLimitsGC[LOWER][i] = 1;
			blkLimitsGC[UPPER][i] = 1;
		}
	}

	if (MESH_MDIM >= 2) {
  		blkLimits[LOWER][JAXIS] = GY; 
  		blkLimits[UPPER][JAXIS] = GY + getNcells(JAXIS) - 1;
  		
  		blkLimitsGC[LOWER][JAXIS] = 0;
  		blkLimitsGC[UPPER][JAXIS] = getNcells(JAXIS) + 2*GY;
		
		for (int i = KAXIS; i <= KAXIS; i++) { //Fill the rest with 1's if they are not present.
			blkLimits[LOWER][i] = 1;
			blkLimits[UPPER][i] = 1;

			blkLimitsGC[LOWER][i] = 1;
			blkLimitsGC[UPPER][i] = 1;
		}

	}

	if (MESH_MDIM == 3) {
  		blkLimits[LOWER][KAXIS] = GZ;
  		blkLimits[UPPER][KAXIS] = GZ + getNcells(KAXIS) - 1; 
  	
  		blkLimitsGC[LOWER][KAXIS] = 0;
		blkLimitsGC[UPPER][KAXIS] = getNcells(KAXIS) + 2*GZ;
	}

}


void FlashAmrMesh::getBlkPhysicalSize(int blockID, double blockSize[MESH_MDIM]) {
	//	Returns the physical size of the block corresponding to the block ID.
	
	// Inputs:
	//     blockID: index of the block in the blocks array.
	//
	//
	// Outputs:
	//     blockSize[MESH_MDIM]: Physical sizes of the block in order IAXIS, JAXIS, KAXIS. 
	//
	
	/*	Example code in Main, using the same variables as the example at the top of this file:
		
		double physicalSize[MESH_MDIM];

		//Get physical size of the first block and put it into physicalSize array.
		newGrid.getBlkPhysicalSize(0, physicalSize);

		std::cout << "Physical size [0]: " << physicalSize[0] << std::endl;
		std::cout << "Physical size [1]: " << physicalSize[1] << std::endl;
		std::cout << "Physical size [2]: " << physicalSize[2] << std::endl;
	*/

	blockSize[IAXIS] = getSize(blockID,IAXIS);
	blockSize[JAXIS] = getSize(blockID,JAXIS);
	blockSize[KAXIS] = getSize(blockID,KAXIS);
}

void FlashAmrMesh::getBlkPtr(int blockID, double****& dataPtr, GRIDDATASTRUCT gridDataStruct) {
	//Gives dataPtr the address of the given blockID's data array in the blocks array.
	//Note that this will only return the data array for a corresponding blockID. NOT the block struct. 
	
	// Inputs:
	//     blockID: index of the block in the blocks array.
	//
	//     gridDataStruct: Location of the data on the cells. 
	//                  -->Not yet implemented (defaults to CENTER).
	//
	//
	// Outputs:
	//     dataPtr: 4-d array pointer used to directly read variable data (and only variable data) 
	//     from the block.
	//
	
	/* Example code in Main, using the same variables as the example at the top of this file:
		
		double ****dataPtr;

		//Get pointer to the first block's data
		newGrid.getBlkPtr(0, dataPtr, CENTER);

		//Print out the (0,0,0) pressure value of the first block's data
		int varIndex = newGrid.findVarIndex("pres");
		std::cout << dataPtr[newGrid.findVarIndex("pres")][0][0][0] << std::endl;
	*/
	
	if (gridDataStruct != CENTER) {
		cout << "[flash_reader] Function getBlkPtr can only use the grid data struct CENTER currently. Defaulting to CENTER." << endl;
		gridDataStruct = CENTER;
	}

	dataPtr = blocks[blockID].data;
}



void FlashAmrMesh::getDomainBoundBox(double boundBox[2][MESH_MDIM]) {
	//Get block IDs of minimum and maximum block locations for each direction.

	// Inputs:
	//     None
	//
	// Outputs:
	//     boundBox[2][MESH_MDIM]: Array contiaining the minimum and maximum physical locations
	//     on the grid for each direction.
	//
	
	/*  Example code in Main, using the same variables as the example at the top of this file:
	
		double domainBound[2][MESH_MDIM];

		newGrid.getDomainBoundBox(domainBound);

		std::cout << "Domain bound box [0][0]: " << domainBound[0][0] << std::endl;
		std::cout << "Domain bound box [1][0]: " << domainBound[1][0] << std::endl;
		std::cout << "Domain bound box [0][1]: " << domainBound[0][1] << std::endl;
		std::cout << "Domain bound box [1][1]: " << domainBound[1][1] << std::endl;
		std::cout << "Domain bound box [0][2]: " << domainBound[0][2] << std::endl;
		std::cout << "Domain bound box [1][2]: " << domainBound[1][2] << std::endl;
	*/

	double minX = realScalarsMap["xmin"];
	double minY = realScalarsMap["ymin"];
	double minZ = realScalarsMap["zmin"];
	
	double maxX = realScalarsMap["xmax"];
	double maxY = realScalarsMap["ymax"];
	double maxZ = realScalarsMap["zmax"];

	if (MESH_MDIM >= 1) {
		boundBox[LOWER][IAXIS] = minX;
		boundBox[UPPER][IAXIS] = maxX;
	}
	if (MESH_MDIM >= 2) {
		boundBox[LOWER][JAXIS] = minY;
		boundBox[UPPER][JAXIS] = maxY;
	}
	if (MESH_MDIM == 3) {
		boundBox[LOWER][KAXIS] = minZ;
		boundBox[UPPER][KAXIS] = maxZ;
	}
	else{
		cout << "[flash_reader] Dimension used in getBlkBoundBox (MESH_MDIM) is not 1, 2, or 3. Exiting routine." << endl;
		exit(EXIT_FAILURE);
	}
}

void FlashAmrMesh::releaseBlkPtr(int blockID, double****& blkPtr, GRIDDATASTRUCT gridDataStruct) {
	// This function doesn't actually do anything. The way that FLASH works isn't the same
	// as the way this code works. It is a 5-d array, and we use a blocks array and 4-d, so the 
	// release isn't necessary. Note that this is regarding the AMR mesh currently.

	// Inputs:
	//     blockID: index of the block in the blocks array.
	//
	//     gridDataStruct: Location of the data on the cells. 
	//                     -->Not yet implemented (defaults to CENTER).
	//
	//
	// Outputs:
	//     blkPtr: Releases the block pointer to the data array (not yet, but maybe when the code
	//             is updated later).
	//

	/* Example code in Main, using the same variables as the example at the top of this file:
		
		double ****dataPtr;

		//Get pointer to the first block's data
		newGrid.getBlkPtr(0, dataPtr, CENTER);

		//Print out the (0,0,0) pressure value of the first block's data
		int varIndex = newGrid.findVarIndex("pres");
		std::cout << dataPtr[newGrid.findVarIndex("pres")][0][0][0] << std::endl;
		
		//"Release" the block pointer
		releaseBlkPtr(0, dataPtr, CENTER);
	*/
	
	if (gridDataStruct != CENTER) {
		cout << "[flash_reader] Function releaseBlkPtr can only use the grid data struct CENTER currently. Defaulting to CENTER." << endl;
		gridDataStruct = CENTER;
	}

}


void FlashAmrMesh::getListOfBlocks(BLOCKTYPE blockType, int *listOfBlocks, int &count, int refinementLevel, 
							  double region_bndBox[2][MESH_MDIM], bool includePartialBlocks) {
	
	// TODO: Haven't tested this function yet. Should be returning the list of blocks assigned to
	//       one processor.

	//Does not deal with different processors
	double minX = realScalarsMap["xmin"];
	double minY = realScalarsMap["ymin"];
	double minZ = realScalarsMap["zmin"];
	
	double maxX = realScalarsMap["xmax"];
	double maxY = realScalarsMap["ymax"];
	double maxZ = realScalarsMap["zmax"];

	double epsilon = 1e-5;
	
	//want epsilon to be a fraction of the domain, not
	//just some fixed value for all possible domains.
	double epsilonX = epsilon * (maxX - minX); 
	double epsilonY = epsilon * (maxY - minY);
	double epsilonZ = epsilon * (maxZ - minZ);

	switch(blockType) {
		case ALL_BLKS:
			count = 0;
			for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) {
				listOfBlocks[i] = i; //Block index is no different than regular index.
				count++;
			}
			break;

		case IBDRY_BLKS: 
			count = 0;
			for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) {
				double coordX = getCoordinate(i,IAXIS);
				if (abs(coordX - minX) < epsilonX || abs(coordX - maxX) < epsilonX) {
					listOfBlocks[count] = i;
					count++;
				}
			}
			break;

		case JBDRY_BLKS:
			count = 0;
			for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) {
				double coordY = getCoordinate(i,JAXIS);
				if (abs(coordY - minY) < epsilonY || abs(coordY - maxY) < epsilonY) {
					listOfBlocks[count] = i;
					count++;
				}
			}
			break;

		case KBDRY_BLKS:
			count = 0;
			for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) {
				double coordZ = getCoordinate(i,KAXIS);
				if (abs(coordZ - minZ) < epsilonZ || abs(coordZ - maxZ) < epsilonZ) {
					listOfBlocks[count] = i;
					count++;
				}
			}
			break;

		case ANY_BDRY_BLKS:
			count = 0;
			for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) {
				double coordX = getCoordinate(i,IAXIS);
				double coordY = getCoordinate(i,JAXIS);
				double coordZ = getCoordinate(i,KAXIS);
				if (abs(coordX - minX) < epsilonX || abs(coordX - maxX) < epsilonX ||
					abs(coordY - minY) < epsilonY || abs(coordY - maxY) < epsilonY ||
					abs(coordZ - minZ) < epsilonZ || abs(coordZ - maxZ) < epsilonZ) 
				{
						listOfBlocks[count] = i;
						count++;
				}
			}

			break;


		case ACTIVE_BLKS:
			cout << "[flash_reader] ACTIVE_BLKS flag not yet implemented in function getListOfBlocks. Exiting program" << endl;
			exit(EXIT_FAILURE);
			break;

		case LEAF:
			count = 0;
			for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) 
            {   if(this->getNodeType(i) == 1)
                {
				    listOfBlocks[count] = i; //Block index is no different than regular index.
				    count++;
                }
                
			}
			break;

		case PARENT_BLK:
			cout << "[flash_reader] PARENT_BLK flag not yet implemented in function getListOfBlocks. Exiting program" << endl;
			exit(EXIT_FAILURE);
			break;

		case ANCESTOR:
			cout << "[flash_reader] ANCESTOR flag not yet implemented in function getListOfBlocks. Exiting program" << endl;
			exit(EXIT_FAILURE);
			break;

		case REFINEMENT:
			cout << "[flash_reader] REFINEMENT flag not yet implemented in function getListOfBlocks. Exiting program" << endl;
			exit(EXIT_FAILURE);
			break;

		case INREGION:
			//Deal with the blocks in region_bndBox
			if (region_bndBox != NULL) {
				if (includePartialBlocks) {
					count = 0;
					for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) {
						double coordX = getCoordinate(i,IAXIS);
						double coordY = getCoordinate(i,JAXIS);
						double coordZ = getCoordinate(i,KAXIS);
						double physicalSizeX = getSize(i,IAXIS);
						double physicalSizeY = getSize(i,JAXIS);
						double physicalSizeZ = getSize(i,KAXIS);

						if (coordX > (region_bndBox[LOWER][IAXIS] - physicalSizeX) && 
							coordX < (region_bndBox[UPPER][IAXIS] + physicalSizeX) &&
							coordY > (region_bndBox[LOWER][JAXIS] - physicalSizeY) && 
							coordY < (region_bndBox[UPPER][JAXIS] + physicalSizeY) &&
							coordZ > (region_bndBox[LOWER][KAXIS] - physicalSizeZ) && 
							coordZ < (region_bndBox[UPPER][KAXIS] + physicalSizeZ) )
						{
							listOfBlocks[count] = i;
							count++;
						}
					}
				}
				else {
					count = 0;
					for (int i = 0; i < getLocalNumBlocks() /*getNBlocks()*/; i++) {
						double coordX = getCoordinate(i,IAXIS);
						double coordY = getCoordinate(i,JAXIS);
						double coordZ = getCoordinate(i,KAXIS);
						if (coordX >= region_bndBox[LOWER][IAXIS] && 
							coordX <= region_bndBox[UPPER][IAXIS] &&
							coordY >= region_bndBox[LOWER][JAXIS] && 
							coordY <= region_bndBox[UPPER][JAXIS] &&
							coordZ >= region_bndBox[LOWER][KAXIS] && 
							coordZ <= region_bndBox[UPPER][KAXIS]) 
						{
							listOfBlocks[count] = i;
							count++;
						}
					}
				}
			}
			else {
				cout << "[flash_reader] region_bndBox array not used as argument to function getListOfBlocks. Exiting program" << endl;
				exit(EXIT_FAILURE);
			}
			break;
	}

	//Not currently implemented.
	// if (region_bndBox != NULL && blockType == REFINEMENT) {
	// 	//Provide refinement values for each block

	// }


}

void FlashAmrMesh::putBlkData(int blockID, GRIDDATASTRUCT gridDataStruct, const char structIndex[MESH_VAR_STRING_SIZE+1], 
	                     BEGINCOUNT beginCount, int startingPos[MESH_MDIM], double****& dataBlock 
	                     /* , int *dataSize */) {
	//Puts block data from main into the class. (Note that dataSize is not 
	// used, since the size of the array in the argument is arbitrary. 
	// It will always read all of the data available for the 
	// desired block (except for guard cells).

	// Inputs:
	//		blockID:  id of the block to have variable data read into it.
	//
	//      gridDataStruct: Location of the data on the cells. 
	//                      -->Not yet implemented (defaults to CENTER).
	//
	//      structIndex: The variable being read in to dataBlock. A character array is used to
	//                   specify the variable. Format must be lowercase with only four characters. 
	//                   Example: "pres" or "dens".
	//
	//      beginCount:  Tells the reader whether or not guard cells should be included in the data read in.
	//                   -->Not yet implemented (defaults to INTERIOR since guard cells are not yet implemented). 
	//
	//      startingPos: MESH_MDIM size array containing the starting indexes of where to read in data.
	//                   This is basically implemented, but for now it defaults to (0,0,0).
	//
	// Outputs:
	//      dataBlock:   Data array that the data will be read into, passed by reference.
	//

	// Example: 
	
	/*
	int main() {

		char inFileName[] = "snb_hdf5_chk_0000";
		char varName[] = "dens";
		//char outFileName[] = "test.txt";
		//char outFileType[] = "ascii";

		H5std_string FILE_NAME( inFileName );

		FlashAmrMesh newGrid;


		newGrid.setGridFileName( FILE_NAME );
		
		newGrid.openFile();

		newGrid.generalRead();

		//Read all the variables
		//newGrid.readAllVariables();

		//Read only some variables
		newGrid.addUserVariable(varName);
		newGrid.readPartialVariables();

		//int numberOfBlocks = newGrid.getNBlocks();

		double**** testBlock;

		//Allocate one block of data with values 1
		allocateVarArray(testBlock, newGrid.getNUserVars(), newGrid.getNcells());

		//Show some of the variables before assignment
		std::cout << "[Before] Block 0, var " << newGrid.findVarIndex(varName) << ", cell 0,0,0: " 
			 << newGrid.getVariable(0,newGrid.findVarIndex(varName),0,0,0) << std::endl;

		std::cout << "[Before] Block 0, var " << newGrid.findVarIndex(varName) << ", cell 1,1,0: " 
			 << newGrid.getVariable(0,newGrid.findVarIndex(varName),1,1,0) << std::endl;

		std::cout << "[Before] Block 0, var " << newGrid.findVarIndex(varName) << ", cell 14,14,0: " 
			 << newGrid.getVariable(0,newGrid.findVarIndex(varName),14,14,0) << std::endl;

		//Set all values in testBlock to 1.1.
		for (int i = 0; i < newGrid.getNcells(IAXIS); i++) {
			for (int j = 0; j < newGrid.getNcells(JAXIS); j++) {
				for (int k = 0; k < newGrid.getNcells(KAXIS); k++) {
					testBlock[0][i][j][k] = 1.1;
				}
			}
		}

		//Put block into grid.
		int startingPos[MESH_MDIM];
		startingPos[0] = 0;
		startingPos[1] = 0;
		startingPos[2] = 0;

		newGrid.putBlkData(0,CENTER,varName,INTERIOR,startingPos,testBlock);

		//Show some of the variables after assignment
		std::cout << "[After] Block 0, var " << newGrid.findVarIndex(varName) << ", cell 0,0,0: " 
			 << newGrid.getVariable(0,newGrid.findVarIndex(varName),0,0,0) << std::endl;

		std::cout << "[After] Block 0, var " << newGrid.findVarIndex(varName) << ", cell 1,1,0: " 
			 << newGrid.getVariable(0,newGrid.findVarIndex(varName),1,1,0) << std::endl;

		std::cout << "[After] Block 0, var " << newGrid.findVarIndex(varName) << ", cell 14,14,0: " 
			 << newGrid.getVariable(0,newGrid.findVarIndex(varName),14,14,0) << std::endl;




		newGrid.closeFile();

		return 0;

	}
	*/


	if (gridDataStruct != CENTER) {
		cout << "[flash_reader] Function getBlkData can only use the grid data struct CENTER currently. Defaulting to CENTER." << endl;
		gridDataStruct = CENTER;
	}

	if (startingPos[0] != 0 || startingPos[1] != 0 || startingPos[2] != 0) {
		cout << "[flash_reader] Function getBlkData can only use startingPos[] = 0 currently. Defaulting to 0. " << endl;
		startingPos[0] = 0;
		startingPos[1] = 0;
		startingPos[2] = 0;
	}

	if (beginCount != INTERIOR) {
		cout << "[flash_reader] Function getBlkData does not yet support counting guard cells, so only INTERIOR can be used. Defaulting to INTERIOR. " << endl;
		beginCount = INTERIOR;
	}

	//First assume that we're just searching through all the variables.
	int userArrayIndex = findVarIndex(structIndex);
//	if (!SUPPRESSOUTPUT) cout << "User Array Index " << userArrayIndex << endl;

	if (userArrayIndex == -1) {
		cout << "[flash_reader] Requested variable in getBlkData function not found. Exiting routine. " << endl;
		exit(EXIT_FAILURE);
	}


	//This conditional will be updated if there is a need to include functionality for 
	// different gridDataStructs and guard cells. Right now the conditional is meaningless;
	// it is more of a placeholder. These conditions would be required for this specific call, though.
	if (gridDataStruct == CENTER && beginCount == INTERIOR) {
		for (int k = startingPos[2]; k < getNcells(KAXIS); k++) {
			for (int j = startingPos[1]; j < getNcells(JAXIS); j++) {
				for (int i = startingPos[0]; i < getNcells(IAXIS); i++) {
					blocks[blockID].data[userArrayIndex][i][j][k] = dataBlock[userArrayIndex][i][j][k];
				}
			}
		}
	}
}


double FlashAmrMesh::getSingleCellVol(int blockID, BEGINCOUNT beginCount, int point[MESH_MDIM])
{
	// this function only works for Cartesian geometry
	if (getGeometryType() != CARTESIAN)
	{
		cout << "getSingleCellVol only works for Cartesian geometry...exiting" << endl;
		exit(EXIT_FAILURE);
	}

	double sideLengths[MESH_MDIM];
	getCellSideLengths(blockID, sideLengths);

	return sideLengths[IAXIS] * sideLengths[JAXIS] * sideLengths[KAXIS];
}


void FlashAmrMesh::getCellCoords(int axis, int blockID, int edge, bool guardcell, double* coordinates)
{

    
    int ncells = getNcells(axis);
    double ncellsD = ncells;
    int length, size;
    double boundBox[2][MESH_MDIM];
    double lowerBound, upperBound, dx;
    
    getBlkBoundBox(blockID, boundBox);
    lowerBound = boundBox[LOWER][axis];
    upperBound = boundBox[UPPER][axis];
    size = ncells; 
    dx = (upperBound-lowerBound) / ncellsD;

    
    if(guardcell)
    {
       lowerBound = lowerBound - double(NGUARD)*dx;
       size = size + NGUARD; 
    }

    
    switch(edge)
    {
        case CENTER:
        {
            coordinates[0] = lowerBound+(dx/2.0);
            break;
        }
        case LEFT_EDGE:
        {
            coordinates[0] = lowerBound;
            break;
        }
        case RIGHT_EDGE:
        {
            coordinates[0] = lowerBound + dx;
            break;
        }

        case FACE:
        {
            coordinates[0] = lowerBound;
            break;
        }
        default:
        {
            cout << "Get Coords : invalid edge specification, must be LEFT_EDGE  RIGHT_EDGE, or CENTER" << endl;
		    exit(EXIT_FAILURE);
        }
    }
        
    for( int i = 1; i <= size; i++)
    {
            coordinates[i] = coordinates[0] + double(i)*dx;
            
    }
    
    //cout << "Exiting Normally" <<endl;
    return;


}


void FlashAmrMesh::getCellBoundBox(int blockId, BEGINCOUNT beginCount, int i, int j, int k, double boundBox[2][MESH_MDIM])
{
	if (beginCount != INTERIOR)
	{
		std::cout << "There is currently no provision for guard cells in getCellBoundBox. beginCount must be INTERIOR..." << std::endl;
		exit(EXIT_FAILURE);
	}

	double boundBlock[2][MESH_MDIM];
	getBlkBoundBox(blockId, boundBlock);
	double cellLengthx = (boundBlock[UPPER][IAXIS] - boundBlock[LOWER][IAXIS]) / nCellsVec[IAXIS],
	cellLengthy = (boundBlock[UPPER][JAXIS] - boundBlock[LOWER][JAXIS]) / nCellsVec[JAXIS],
	cellLengthz = (boundBlock[UPPER][KAXIS] - boundBlock[LOWER][KAXIS]) / nCellsVec[KAXIS];

	boundBox[LOWER][IAXIS] = boundBlock[LOWER][IAXIS] + cellLengthx * i;
	boundBox[UPPER][IAXIS] = boundBlock[LOWER][IAXIS] + cellLengthx * (i+1);

	boundBox[LOWER][JAXIS] = boundBlock[LOWER][JAXIS] + cellLengthy * j;
	boundBox[UPPER][JAXIS] = boundBlock[LOWER][JAXIS] + cellLengthy * (j+1);

	boundBox[LOWER][KAXIS] = boundBlock[LOWER][KAXIS] + cellLengthx * k;
	boundBox[UPPER][KAXIS] = boundBlock[LOWER][KAXIS] + cellLengthx * (k+1);	
}

void FlashAmrMesh::getCellSideLengths(int blockID, double sideLengths[MESH_MDIM])
{
	// this function only works for Cartesian geometry
	if (getGeometryType() != CARTESIAN)
	{
		cout << "getCellSideLength only works for Cartesian geometry...exiting" << endl;
		exit(EXIT_FAILURE);
	}

	double blockLengths[MESH_MDIM];
	getBlkPhysicalSize(blockID, blockLengths);

	int nxb = getNcells(IAXIS);
	int nyb = getNcells(JAXIS);
	int nzb = getNcells(KAXIS);

	sideLengths[IAXIS] = blockLengths[IAXIS] / nxb;
	sideLengths[JAXIS] = blockLengths[JAXIS] / nyb;
	sideLengths[KAXIS] = blockLengths[KAXIS] / nzb;

}

void FlashAmrMesh::setVariable(int blockId, string structIndex, int i, int j, int k, double value) {
	int variableIndex = findVarIndex(structIndex.c_str());
	if (variableIndex == -1) {
		cout << "[flash_reader] Requested variable in getBlkData function not found. Exiting routine. " << endl;
		exit(EXIT_FAILURE);
	}
	blocks[blockId].data[variableIndex][i][j][k] = value;
}

int FlashAmrMesh::mpiGetID() {
	return comm.Get_rank();
}

int FlashAmrMesh::mpiGetProcs() {
	return comm.Get_size();
}

void FlashAmrMesh::mpiAssignBlocks() {
	if (nBlocks == 0) {
		cout << "[flash_reader] mpiAssignBlocks must come after the general read function." << std::endl;
		exit(EXIT_FAILURE);
	}

	int extra = nBlocks % comm.Get_size();

	//Gives the first few processors an extra block when
	//the number of blocks doesn't divide evenly.
	if (comm.Get_rank() < extra) {
		nBlocksAssigned = nBlocks / comm.Get_size() + 1;
	}
	else {
		nBlocksAssigned = nBlocks / comm.Get_size();
	}
}

int FlashAmrMesh::getLocalNumBlocks() {
	return nBlocksAssigned;
}

int FlashAmrMesh::mpiGetGlobalBlockID(int localBlockNumber) {
	if (nBlocksAssigned == 0) {
		cout << "[flash_reader] mpiGetGlobalBlockID must come after the mpiAssignBlocks function." << std::endl;
		exit(EXIT_FAILURE);
	}	

	int extra = nBlocks % comm.Get_size();

	if (comm.Get_rank() < extra) {
		return comm.Get_rank() * nBlocksAssigned + localBlockNumber;
	}
	else {
		return (extra * (nBlocksAssigned + 1)) // Account for the first few larger amounts
			   + ((comm.Get_rank() - extra) * nBlocksAssigned) //Add the rest of the smaller ones
			   + localBlockNumber; //Offset
	}
}

int FlashAmrMesh::mpiGetLocalBlockID(int globalBlockNumber) {
	if (nBlocksAssigned == 0) {
		cout << "[flash_reader] mpiGetLocalBlockID must come after the mpiAssignBlocks function." << std::endl;
		exit(EXIT_FAILURE);
	}

	int extra = nBlocks % comm.Get_size();

	if (comm.Get_rank() < extra) {
		return globalBlockNumber % nBlocksAssigned;
	}
	else {
		return (globalBlockNumber - (extra * (nBlocksAssigned + 1))) //Subtract larger amounts
			   % nBlocksAssigned; 
	}
}

void FlashAmrMesh::setLogical(const char *whichValue, const bool newValue)
{
	logicalScalarsMap[whichValue] = newValue;
	std::cout << "New value of usehydro: " << logicalScalarsMap[whichValue] << endl;
}

double FlashAmrMesh::getPointData(int blockID, GRIDDATASTRUCT gridDataStruct, int structIndex, BEGINCOUNT beginCount, int position[MESH_MDIM])
{
	if (beginCount != INTERIOR)
	{
		cout << "[flash_reader] beginCount must be set to INTERIOR in FlashAmrMesh::getPointData." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (gridDataStruct != CENTER)
	{
		cout << "[flash_reader] gridDataStruct must be set to CENTER in FlashAmrMesh::getPointData." << std::endl;
		exit(EXIT_FAILURE);
	}

	return blocks[blockID].data[structIndex][position[IAXIS]][position[JAXIS]][position[KAXIS]];
}

void FlashAmrMesh::getCellCoords(int blockID, int i, int j, int k, double cellCoords[MESH_MDIM])
{
	double cellBoundBox[2][MESH_MDIM];
	getCellBoundBox(blockID, INTERIOR, i, j, k, cellBoundBox);

	double sideLengths[MESH_MDIM];
    
    getCellSideLengths(blockID, sideLengths);

	for (int dim=0; dim<MESH_MDIM; dim++)
		cellCoords[dim] = cellBoundBox[LOWER][dim] + sideLengths[dim] / 2.;

}


// void amrGrid::mpiInit( int argc, char** argv ) {
// 	MPI::Init( argc, argv );
// 	mpiID = MPI::COMM_WORLD.Get_rank();
// 	procs = MPI::COMM_WORLD.Get_size();

// }

// void amrGrid::mpiFinalize() {
// 	MPI::Finalize();
// }

