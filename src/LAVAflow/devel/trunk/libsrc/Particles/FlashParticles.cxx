#include "FlashParticles.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "ArrayOperations.h"

using namespace std;

FlashParticles::FlashParticles()
{
	nUserVars = 0;
	nFileVars = 0;
	nParticlesLocal = 0;
    nParticlesGlobal = 0;
	particleData = NULL;
	comm = MPI::COMM_WORLD;
    simTime = 0;
}


FlashParticles::FlashParticles(const char *inFileName, const char **varNames, int nVars, int sortMode, std::string maskVar, int maskID, MPI::Intracomm recvComm)
{
	nUserVars = 0;
	nFileVars = 0;
	nParticlesLocal = 0;
    nParticlesGlobal = 0;
	comm = recvComm;
    this->sortMode = sortMode;
    particleData = NULL;

	filename = inFileName;
		
	// std::cout << "Initializing..." << std::endl;

    Initialize();

    // std::cout << "Initialized..." << std::endl;

    if (varNames != NULL)
    {
        std::vector<std::string> varNameList = CharToVecOfString(varNames, nVars);

    	for (auto variable : varNameList)
    		addUserVariable(variable.c_str());
    	
    	readVarData(maskVar, maskID);
    }
}

FlashParticles::~FlashParticles()
{
	closeFile();

	if (particleData != NULL)
        delete[] particleData[0];
	delete[] particleData;

}

void FlashParticles::openFile()
{
	std::ifstream infile(filename.c_str());

	if (infile.good())
	{

		h5File = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

		if (h5File < 0)
		{
			cout << "Error opening file...aborting" << endl;
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		cout << "Requested file does not exist...aborting" << endl;
		exit(EXIT_FAILURE);
	}
}

void FlashParticles::closeFile()
{
	H5Fclose(h5File);
}




void FlashParticles::readVarData(std::string maskVar, int maskID)
{
    if (maskID > -1 && sortMode != PART_SORT_TAG)
    {
        std::cerr << "Masking only supported when sorting by particle tag. Specified mask will have no effect." << std::endl;
    }
	
	string requestedVar = "tracer particles";

	// this represents the HDF5 dataset
	hid_t dset = H5Dopen(h5File, requestedVar.c_str(), H5P_DEFAULT);
	
	// get a copy of the file's dataspace
	// this represents the form of the raw data in the file
	hid_t space = H5Dget_space( dset );

	//  get the dimensionality of the space
	int nDims = H5Sget_simple_extent_ndims(space);

    // cout << "nDims: " << nDims << endl;

	// create a dataspace in memory
	// this is where the data from the file will be copied
	// its dimensions correspond to the file's dataspace
	hsize_t dsDims[nDims];
	H5Sget_simple_extent_dims( space, dsDims, NULL );

	// int nParticlesGlobal = dsDims[0];
    
    // allocate the data array
    particleData = new double*[nUserVars];
    particleData[0] = new double[nUserVars*nParticlesLocal];
    for (int i=1; i<nUserVars; i++)
        particleData[i] = particleData[i-1] + nParticlesLocal;

	
    // reading in the particles will depend on how they are sorted

    // this is only used if sorting by tag
    int *partIndThisProc = new int[nParticlesLocal];
    // int partIndThisProcTemp[nParticlesLocal];

    // if sorting by owner proc, we can use hyperslabs
    // figure out hyperslab bounds
    
    
    
    if (sortMode == PART_SORT_TAG)
    {
        // if sorting by tag, we need to construct an array of indices for this proc
        // these represent the rows in the hdf5 file that this proc should read in

        // read in just the tag IDs
        double *tagIDs = new double[nParticlesGlobal];

        // Define memory space
        // this describes the layout of the memory buffer
        hsize_t msDimsTags[] = {nParticlesGlobal, 1};
        hid_t memspaceTags = H5Screate_simple(nDims, msDimsTags, NULL);

        // point to the correct location to read from the file
        int tagColInd = fileVarsMap_sToI["tag"];
        hsize_t dsOffset[] = {0, tagColInd};
        hsize_t dsCount[] = {nParticlesGlobal, 1};

        H5Sselect_hyperslab(space, H5S_SELECT_SET, dsOffset, NULL, dsCount, NULL);
        H5Dread( dset, H5T_NATIVE_DOUBLE, memspaceTags, space, H5P_DEFAULT, tagIDs);

        // if masking is in effect, read in the masking data
        double *maskVals = new double[nParticlesGlobal];
        for (int i=0; i<nParticlesGlobal; i++)
            maskVals[i] = -1;


        if (maskID > -1)
        {
            // read in just the tag IDs
            
            // point to the correct location to read from the file
            int maskColInd = fileVarsMap_sToI[maskVar.c_str()];
            hsize_t dsOffsetMask[] = {0, maskColInd};
            
            H5Sselect_hyperslab(space, H5S_SELECT_SET, dsOffsetMask, NULL, dsCount, NULL);
            H5Dread( dset, H5T_NATIVE_DOUBLE, memspaceTags, space, H5P_DEFAULT, maskVals);
        }
        
        
        // this is holds the tag IDs for each particle found to belong to this proc
        // used for sorting
        int *tagID = new int[nParticlesLocal];
        int *partIndThisProcTemp = new int[nParticlesLocal];

        // initialize the tagID array
        // for (auto &i : tagID) i=-1;

        // std::cout << "Particle indices: " << beginParticleIndex << ", " << endParticleIndex << std::endl;

        // now we have all the tag IDs
        // record which particles have a tag ID belonging to this proc
        int partIndex = 0;
        for (int tagIndex=0; tagIndex < nParticlesGlobal; tagIndex++)
        {
            int tagIDInt = std::round(tagIDs[tagIndex]-1);

            if (tagIDInt >= beginParticleIndex && tagIDInt <= endParticleIndex)
            {
                
                if (maskID < 0 || (int)std::round(maskVals[tagIndex]) == maskID)
                {
                    partIndThisProcTemp[partIndex] = tagIndex;
                    tagID[partIndex] = tagIDInt;
                    partIndex++;
                }
            }
            // else
            //     std::cout << "Error in tag index: " << tagIDInt << std::endl;
        }

        nParticlesLocal = partIndex;

        // temporary check -- are there any uninitialized values?
        // int index = 0;
        // for (auto i : tagID)
        // {
        //     if (i == -1) std::cout << "Found uninitialized value in tagID at index " << index << std::endl;
        //     index++;
        // }

        // now we have all the row indices of the particles for this proc, along with 
        // associated tag IDs. now sort by tag ID
        
        // build the permutation vector
        std::vector<std::size_t> p(nParticlesLocal);
        
        // fill the permutation vector with sequential values from 0 to nParticlesLocal-1
        std::iota(p.begin(), p.end(), 0);


        std::sort(p.begin(), p.end(), 
            [&](std::size_t i, std::size_t j){ return tagID[i] < tagID[j]; });

        // apply the permutations to the particle IDs
        std::transform(p.begin(), p.end(), partIndThisProc,
        [&](std::size_t i){ return partIndThisProcTemp[i]; });

        delete[] tagID;
        delete[] partIndThisProcTemp;
        delete[] tagIDs;
        delete[] maskVals;

        // mask the values.
    }

	// cout << "Rank " << comm.Get_rank() << " indices: " << beginParticleIndex << "   " << endParticleIndex << endl;

    // cout << "Rank " << mpiRank << " rows: " << endl;

    // for (int i=0; i<10; i++)
    //     cout << partIndThisProc[i] << endl;


	// cout << "Total partices: " << nParticlesGlobal << endl;
	// cout << "Rank " << mpiRank << " particles this proc: " << nParticlesLocal << endl;


	// Define memory space
    // this describes the layout of the memory buffer
    hsize_t msDims[] = {nParticlesLocal, 1};
    hid_t memspace = H5Screate_simple(nDims, msDims, NULL);

    for (int uVarInd=0; uVarInd<nUserVars; uVarInd++)	
	{
		int varInd = fileVarsMap_sToI[variableNames[uVarInd]];
		
		// if sorting by owner proc, just select contiguous regions
        if (sortMode == PART_SORT_PROC)
        {
            // these hyperslabs represent regions of the dataspace
    		// so they allow us to read in only certain fields
    		hsize_t dsOffset[] = {beginParticleIndex, varInd};
    		hsize_t dsCount[] = {nParticlesLocal, 1};
    				
    		H5Sselect_hyperslab(space, H5S_SELECT_SET, dsOffset, NULL, dsCount, NULL);
    		// H5Sselect_hyperslab(memspace, H5S_SELECT_SET, msOffset, NULL, msCount, NULL);
        }
        // if sorting by tag, use the row indices calculated above
        else if (sortMode == PART_SORT_TAG)
        {
            hsize_t pointCoords[nParticlesLocal*2];

            // build the array describing the point coordinates
            for (int i=0; i<nParticlesLocal*2; i+=2)
            {
                pointCoords[i] = partIndThisProc[i/2];
                pointCoords[i+1] = varInd;
            }

            H5Sselect_elements( space, H5S_SELECT_SET, nParticlesLocal, pointCoords );
        }

		
        H5Dread( dset, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, particleData[uVarInd]);
		
		//Insert variable into maps for later reference.
		uVarsMap_sToI.insert(pair<std::string,int>(variableNames[uVarInd],uVarInd));
		uVarsMap_iToS.insert(pair<int,std::string>(uVarInd,variableNames[uVarInd]));
	}

    delete[] partIndThisProc;

	H5Dclose(dset);
	H5Sclose(space);
	H5Sclose(memspace);
}


void FlashParticles::readFieldNames()
{

	char datasetName[] = {"particle names"};

	//Get scalars list.
	hid_t dset = H5Dopen(h5File, datasetName, H5P_DEFAULT);
	
	//Get dimensions of dataset.
	hid_t space = H5Dget_space( dset );
	int nDims = H5Sget_simple_extent_ndims(space);

	hsize_t dims_out[nDims];
	H5Sget_simple_extent_dims( space, dims_out, NULL );

	nFileVars = dims_out[0];

	hid_t filetype = H5Dget_type (dset);
    size_t sdim = H5Tget_size (filetype);
    sdim++;

	char **fileVarNames = new char*[nFileVars];
	fileVarNames[0] = new char[nFileVars*sdim];
	for (int i=1; i<nFileVars; i++)
		fileVarNames[i] = fileVarNames[i-1] + sdim;
	
    hid_t memtype = H5Tcopy (H5T_C_S1);
    herr_t status = H5Tset_size (memtype, sdim);

	H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, fileVarNames[0]);

	//Copy variable names over.
	string whiteSpaces(" \t\f\v\n\r");
	string varNamesStrings[nFileVars];
	for (int i = 0; i < nFileVars; i++)
	{
		varNamesStrings[i] = string(fileVarNames[i]);
		size_t firstWhiteSpace = varNamesStrings[i].find_last_not_of(whiteSpaces);
		varNamesStrings[i].erase(firstWhiteSpace+1);

		fileVarsMap_sToI.insert(pair<std::string,int>(string(varNamesStrings[i]),i));
		fileVarsMap_iToS.insert(pair<int,std::string>(i,string(varNamesStrings[i])));
	}

	
	H5Dclose(dset);
	H5Sclose(space);	
}

int FlashParticles::findVarIndex(const char varName[])
{
	return uVarsMap_sToI[string(varName)];
}

std::string FlashParticles::findVarName(const int varIndex)
{
	return uVarsMap_iToS[varIndex];
}

void FlashParticles::addUserVariable(const char inputVariable[])
{
	string inputVar = inputVariable;

	if (nUserVars > nFileVars || nUserVars >= FR_MAXVARS) {
		cout << "[particle_reader] Could not add " << inputVariable << "." << endl;
		cout << "Too many variables. Can have a maximum of " << nUserVars << " variables. " << endl;
		exit(EXIT_FAILURE);
	}
	else {
		//Check to make sure the variable exists in the map.
		std::map<std::string,int>::iterator it;
		it = fileVarsMap_sToI.find(string(inputVariable));
		if (it == fileVarsMap_sToI.end()) {
			cout << "[particle_reader] Variable '" << inputVariable << "' does not exist in the file variables map." << endl;
			exit(EXIT_FAILURE);
		}
		else {
			//Add variable to list and update nUserVars.
			variableNames.push_back(inputVar);
			nUserVars++;
		}
	}
}

double FlashParticles::getParticleData(int partID, char* varName)
{
	int varIndex = findVarIndex(varName);
	return particleData[varIndex][partID];
}

double FlashParticles::getParticleData(int partID, int varIndex)
{
	return particleData[varIndex][partID];
}

double** FlashParticles::getDataPtr()
{
	return particleData;
}


void FlashParticles::readReals(const char *which)
{
    typedef struct realScalarStruct {
        char    name[MAX_STRING_LENGTH+1];
        double  value;
    } realScalarStruct;


    char datasetName[50] = {};
    if (string(which) == "scalars") {
        strcpy(datasetName,"real scalars");
    }
    else if (string(which) == "parameters") {
        strcpy(datasetName,"real runtime parameters");
    }


    //Get scalars list.
    hid_t dset = H5Dopen(h5File, datasetName, H5P_DEFAULT);
    
    //Get dimensions of dataset.
    hsize_t dims_out[2];
    hid_t space = H5Dget_space( dset );
    // hid_t ndims = H5Sget_simple_extent_dims( space, dims_out, NULL );
    H5Sget_simple_extent_dims( space, dims_out, NULL );

    int sizeOfScalarArray = dims_out[0];

    //Make specific type for the 80 character string.
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size (strtype, MAX_STRING_LENGTH+1);

    //"integer scalars" is a compound list. Have to deal with this differently:
    //Have to create a CompType object and read from that. 
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(realScalarStruct) );
    H5Tinsert( memtype, "name", HOFFSET(realScalarStruct, name), strtype);
    H5Tinsert( memtype, "value", HOFFSET(realScalarStruct, value), H5T_IEEE_F64LE);

    realScalarStruct out[sizeOfScalarArray];
    H5Dread( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);

    //Make array of strings for easy comparisons
    string strArray[sizeOfScalarArray];

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
}

void FlashParticles::readInts(const char *which)
{
    typedef struct intScalarStruct {
        char    name[MAX_STRING_LENGTH+1];
        int  value;
    } intScalarStruct;


    char datasetName[50] = {};
    if (string(which) == "scalars") {
        strcpy(datasetName,"integer scalars");
    }
    else if (string(which) == "parameters") {
        strcpy(datasetName,"integer runtime parameters");
    }


    //Get scalars list.
    hid_t dset = H5Dopen(h5File, datasetName, H5P_DEFAULT);
    
    //Get dimensions of dataset.
    hsize_t dims_out[2];
    hid_t space = H5Dget_space( dset );
    // hid_t ndims = H5Sget_simple_extent_dims( space, dims_out, NULL );
    H5Sget_simple_extent_dims( space, dims_out, NULL );

    int sizeOfScalarArray = dims_out[0];

    //Make specific type for the 80 character string.
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size (strtype, MAX_STRING_LENGTH+1);

    //"integer scalars" is a compound list. Have to deal with this differently:
    //Have to create a CompType object and read from that. 
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(intScalarStruct) );
    H5Tinsert( memtype, "name", HOFFSET(intScalarStruct, name), strtype);
    H5Tinsert( memtype, "value", HOFFSET(intScalarStruct, value), H5T_NATIVE_INT);

    intScalarStruct out[sizeOfScalarArray];
    H5Dread( dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);

    //Make array of strings for easy comparisons
    string strArray[sizeOfScalarArray];

    //Erase white space after name.
    string whiteSpaces(" \t\f\v\n\r");
    for (int i = 0; i < sizeOfScalarArray; i++) {
        strArray[i] = string(out[i].name);
        size_t firstWhiteSpace = strArray[i].find_last_not_of(whiteSpaces);
        strArray[i].erase(firstWhiteSpace+1);
    }

    for (int i = 0; i < sizeOfScalarArray; i++) {
        intScalarsMap.insert(pair<std::string,double>(strArray[i],out[i].value));
    }

    H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, out);
    H5Dclose(dset);
    H5Sclose(space);
    H5Tclose(memtype);
    H5Tclose(strtype);
}

double FlashParticles::getTime()
{
    return realScalarsMap["time"];
}


void FlashParticles::Initialize()
{
    // std::cout << "Filename in initialize: " << filename << std::endl;
    openFile();
    // std::cout << "Opened file" << std::endl;
    readFieldNames();
    // std::cout << "Read field namess" << std::endl;
    readReals("scalars");
    readInts("scalars");

    // figure out how many particles this processor should own
    int numProcs = comm.Get_size();
    int mpiRank = comm.Get_rank();

    nParticlesGlobal = intScalarsMap["globalnumparticles"];
    int leftOverParticles = nParticlesGlobal % numProcs;
    bool extraParticle = mpiRank < leftOverParticles;
    nParticlesLocal = int(nParticlesGlobal / double(numProcs)) + int(extraParticle);

    // figure out the offset in the file based on the rank
    // the particles are already sorted by owner processor in the file
    // so just divide the data into chunks
    beginParticleIndex = mpiRank*nParticlesLocal + int(!extraParticle) * leftOverParticles;
    endParticleIndex = beginParticleIndex + nParticlesLocal - 1;

    // cout << "In initialize: number of global parts: " << nParticlesGlobal << endl;
    // cout << "In initialize: number of local parts: " << nParticlesLocal << endl;
}
