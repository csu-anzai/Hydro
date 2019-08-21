
#include <cstdlib>
#include <iostream>
#include <string>
#include "FlashAmrMesh.h"
#include "Driver_main.h"
#include <string.h>

int LAVAFLOW_DRIVER_BASENM(int argc, char* argv[]) {
	/*

	ExportUniformData:
		Takes multiple command-line arguments to refine one variable in a given AMR mesh down to a uniform
		refinement level. Specifically intended for use with the BoxCounter3D or MultifracBoxCounter3D
		codes.

	Arguments:
		*File name including the path (e.g. /path/to/folder/cartesianShells_hdf5_chk_0000)
		*Variable name
		*Refinement level
		*Output file name
		*Output file type (ascii or binary)

	Output:
		*Ascii or Binary file with uniformly refined mesh data points in the current working directory.

	*/

	char arg_inFileName[MAX_FILENAME_LENGTH+1] = 
		"/home/pb12c/Research/LAVAflow/drivers/files/testData/cartesianShells/cartesianShells_hdf5_chk_0000";
	char arg_varName[MESH_VAR_STRING_SIZE+1] = "vari";
	char arg_refineLevel[MAX_STRING_LENGTH] = "4";
	char arg_outFileName[MAX_FILENAME_LENGTH] = "cartesianShells_0000.txt";
	char arg_outFileType[MAX_STRING_LENGTH] = "ascii";

	std::string defaultFlag("-default");

	std::cout << std::endl << std::endl;
	
	//Error message if the number of arguments is not enough or is too much.
	if (argc == 1) {
		std::cout << "ERROR: Command line arguments not present." << std::endl;
		std::cout << "Input command line arguments in the following format: flashFileName variableName refinementLevel outputFileName outputFileType" << std::endl;
		std::cout << "Example: ./exportUniformData /path/to/folder/cartesianShells_hdf5_chk_0000 vari 4 test.txt binary" << std::endl;
		std::cout << "Otherwise, use -default to use the current default arguments: " << std::endl;
		std::cout << "\tfileName: "         << arg_inFileName  << std::endl;
		std::cout << "\tvarName: "          << arg_varName     << std::endl;
		std::cout << "\tRefinement level: " << arg_refineLevel << std::endl;
		std::cout << "\toutputFileName: "   << arg_outFileName << std::endl;
		std::cout << "\toutputFileType: "   << arg_outFileType << std::endl;
		exit(EXIT_FAILURE);
	}
	if (argc > 1 && argc != 6 && std::string(argv[1]) != defaultFlag) {
		std::cout << "ERROR: Incorrect amount of arguments." << std::endl;
		std::cout << "Input command line arguments in the following format: flashFileName variableName outputFileName outputFileType" << std::endl;
		std::cout << "Example: ./exportUniformData /path/to/folder/cartesianShells_hdf5_chk_0000 vari 4 test.txt binary" << std::endl;
		std::cout << "Otherwise, use -default to use the current default arguments: " << std::endl;
		std::cout << "\tfileName: "         << arg_inFileName  << std::endl;
		std::cout << "\tvarName: "          << arg_varName     << std::endl;
		std::cout << "\tRefinement Level: " << arg_refineLevel << std::endl;
		std::cout << "\toutputFileName: "   << arg_outFileName << std::endl;
		std::cout << "\toutputFileType: "   << arg_outFileType << std::endl;	
		exit(EXIT_FAILURE);
	}
	//Default settings
	// Should really just be used for testing purposes.
	else if (std::string(argv[1]) == defaultFlag) {
		std::cout << "Using default arguments: " << std::endl;
		std::cout << "\tfileName: "         << arg_inFileName  << std::endl;
		std::cout << "\tvarName: "          << arg_varName     << std::endl;
		std::cout << "\tRefinement Level: " << arg_refineLevel << std::endl;
		std::cout << "\toutputFileName: "   << arg_outFileName << std::endl;
		std::cout << "\toutputFileType: "   << arg_outFileType << std::endl;	
	}

	else {
		char yesOrNo[2];

		strcpy(arg_inFileName,"");
		strcpy(arg_inFileName,argv[1]);

		strcpy(arg_varName,"");
		strcpy(arg_varName,argv[2]);

		strcpy(arg_refineLevel,"");
		strcpy(arg_refineLevel,argv[3]);

		strcpy(arg_outFileName,"");
		strcpy(arg_outFileName,argv[4]);

		strcpy(arg_outFileType,"");
		strcpy(arg_outFileType,argv[5]);
		
		std::cout << "Arguments read in are: " << std::endl;
		std::cout << "\tflashFileName: "    << arg_inFileName  << std::endl;
		std::cout << "\tvarName: "          << arg_varName     << std::endl;
		std::cout << "\tRefinement Level: " << arg_refineLevel << std::endl;
		std::cout << "\toutputFileName: "   << arg_outFileName << std::endl;
		std::cout << "\toutputFileType: "   << arg_outFileType << std::endl;
		std::cout << std::endl;
		std::cout << "Are the entered arguments correct? (y/n): ";
		std::cin >> yesOrNo;
		if (std::string(yesOrNo) != "y") {
			std::cout << "Exiting program to allow for entering correct arguments." << std::endl;
			exit(EXIT_FAILURE);
		}
	}


	std::cout << std::endl;
	std::cout << std::endl;

	//Create the specific variables needed for the AMR Mesh class in LAVAflow
	std::string inFileName(arg_inFileName);

	std::vector<std::string> meshVars;
	meshVars.push_back(arg_varName);

	int refineLevel = atoi(arg_refineLevel);

	char outFileName[MAX_FILENAME_LENGTH];
	strcpy(outFileName,arg_outFileName);

    char outFileHDF5Name[MAX_FILENAME_LENGTH];
	strcpy(outFileHDF5Name,arg_outFileName);

	char outFileType[MAX_STRING_LENGTH];
	strcpy(outFileType,arg_outFileType);

	std::cout << inFileName << " " << meshVars[0] << " " << refineLevel << " " << 
	             outFileName << " " << outFileType << std::endl;

	FlashAmrMesh mesh(inFileName, outFileHDF5Name, meshVars);

	double**** fineArray;
	
	double subdomainCoords[2][MESH_MDIM];

	subdomainCoords[LOWER][IAXIS] = 0;
	subdomainCoords[UPPER][IAXIS] = 0;
	subdomainCoords[LOWER][JAXIS] = 0;
	subdomainCoords[UPPER][JAXIS] = 0;
	subdomainCoords[LOWER][KAXIS] = 0;
	subdomainCoords[UPPER][KAXIS] = 0;

	mesh.refineToFinest(fineArray, subdomainCoords, refineLevel);
	std::cout << std::string(outFileType) << std::endl;

    if (std::string(outFileType) == "binary" || std::string(outFileType) == "ascii") {
	    mesh.printRefinedData(fineArray, outFileName, outFileType, meshVars[0].c_str());
    }
    else {
    	/* Updated HDF5 output. Does not used the specific format from 
    	   preintRefinedDataHDF5 that was originally written, but never used 
    	   for boxCounter3D or multifracBoxCounter3D. -pb
    	*/
		mesh.writeOutFileUniform(fineArray);
    }

	mesh.deallocateFinest(fineArray);

	return 0;

}