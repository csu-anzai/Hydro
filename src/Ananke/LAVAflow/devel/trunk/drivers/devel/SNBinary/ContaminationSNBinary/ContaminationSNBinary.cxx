
#include "CommonSNBinary.h"
using namespace std;

/* TODO: Look at resolution effects between the two models. 
	 	 Make table of each of the stripped masses
	 	 Make table of resolutions based on radii of our models, I think.*/

void printOut(char* fileName, double checkTime, double iron, double nickel) {

	std::ofstream outFile;
	outFile.open(fileName);

	outFile << checkTime << "," << iron << "," << nickel << std::endl;
	// std::cout << fileName << std::endl;
	// std::cout << checkTime << std::endl;
	// std::cout << mass << std::endl;

	outFile.close();
}

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv) {

	modelContainer models[NUMMODELS];
	readModels(models);

	/* Declare variables to be read */
	int nVars = 13;
	char** vars;
	int* varIndex;
	vars = new char*[nVars];
	for (int i = 0; i < nVars; i++) {
		vars[i] = new char[MESH_VAR_STRING_SIZE+1];
	}
	varIndex = new int[nVars];
	char inFileName[MAX_FILENAME_LENGTH];
	char outFileName[MAX_FILENAME_LENGTH];
	char inProgramFileName[MAX_FILENAME_LENGTH];
	char outProgramFileName[MAX_FILENAME_LENGTH];
	char outFileExtraName[MAX_FILENAME_LENGTH];

	/* --- Specify variables --- */

	MODEL currentModel = MS38_L200;
	// MODEL currentModel = MS54;
	// MODEL currentModel = MS63;
	// MODEL currentModel = MS7_L200;
	// MODEL currentModel = SG_L200;
	// MODEL currentModel = SY319_L200;
	// MODEL currentModel = SY428;

	int checkNum = MS38_L200_F;
	// int checkNum = MS54_F;
	// int checkNum = MS63_F;
	// int checkNum = MS7_L200_F;
	// int checkNum = SG_L200_F;
	// int checkNum = SY319_L200_F;
	// int checkNum = SY428_F;

	if (argc > 1) {
		currentModel = (MODEL)atoi(argv[1]);
		checkNum = atoi(argv[2]);
	}

	strcpy(vars[0],"dens");
	strcpy(vars[1],"gpot");
	strcpy(vars[2],"ener");
	strcpy(vars[3],"ms_2");
	strcpy(vars[4],"ms_3");
	strcpy(vars[5],"temp");
	strcpy(vars[6],"eint");
	strcpy(outFileExtraName, "_contamination.txt");


	strcpy(vars[7] ,"H1  ");
	strcpy(vars[8] ,"He4 ");
	strcpy(vars[9],"O16 ");
	strcpy(vars[10],"Si28");
	strcpy(vars[11],"Fe52");
	strcpy(vars[12],"Ni56");

	/* Useful Variable names:
		dens = Density
		gpot = Gravitational Potential
		ener = Kinetic Energy
		ms_1 = Ambient Medium Mass Scalar
		ms_2 = Ejecta Mass Scalar
		ms_3 = Companion Mass Scalar
		ms_4 = mass...Mass Scalar?
	*/	

	strcpy(inFileName, models[currentModel].checkpoint[checkNum]);
	strcpy(outFileName, models[currentModel].checkpoint[checkNum]);

	strcat(outFileName, outFileExtraName);

	strcpy(inProgramFileName, inFolderLocation);
	strcat(inProgramFileName, models[currentModel].modelFolder);
	
	strcpy(outProgramFileName, outFolderLocation);
	strcat(outProgramFileName, models[currentModel].modelFolder);
	strcat(outProgramFileName, "contamination/");

	strcat(inProgramFileName,inFileName);
	strcat(outProgramFileName,outFileName);
	
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "WARNING: THIS ANALYSIS PROGRAM NOT YET BUILT FOR MPI" << std::endl;
	std::cout << "Input file: " << inProgramFileName << std::endl;
	std::cout << "Output file: " << outProgramFileName << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	/* ------------------------- */

	string inFileNameStr(inProgramFileName);
	string outFileNameStr(outProgramFileName);

	/* Initialization */
	vector<string> meshVars;
	for (int i = 0; i < nVars; i++) {
		meshVars.push_back(vars[i]);
	}

    FlashAmrMesh mesh(inFileNameStr, outFileNameStr, meshVars);

	for (int i = 0; i < nVars; i++) {
		varIndex[i] = mesh.findVarIndex(vars[i]);
	}

	std::cout << std::endl << "--Entering Main Functions--" << std::endl << std::endl;
	int numberOfBlocks = mesh.getNBlocks();

	double***** blocksMain;

	//Allocate data arrays for each block.
	blocksMain = new double****[numberOfBlocks];

	for (int i = 0; i < numberOfBlocks; i++) {
		allocateVarArray(blocksMain[i], mesh.getNUserVars(), mesh.getNcells());
	}

	//Get leaf block IDs
	std::set<int> leafBlockIDs;

	for (int i = 0; i < numberOfBlocks; i++) {
		if (mesh.getNodeType(i) == FR_LEAF_NODE) {
			leafBlockIDs.insert(i);
		}
	}

	//Retrieve the data from each block.
	int startingPos[MESH_MDIM];
	startingPos[IAXIS] = 0;
	startingPos[JAXIS] = 0;
	startingPos[KAXIS] = 0;

	for (int i = 0; i < numberOfBlocks; i++) {
		for (int j = 0; j < nVars; j++) {
			mesh.getBlkData(i, CENTER, vars[j], INTERIOR, startingPos, blocksMain[i]);
		}
	}


	// ---------------
	// MAIN LOOP SETUP
	// ---------------

	int nCellsX = mesh.getNcells(IAXIS);
	int nCellsY = mesh.getNcells(JAXIS);
	int* nCellsArr = mesh.getNcells();

	// ---------
	// MAIN LOOP
	// ---------

	std::set<int>::iterator leaf;

	// double totalMass = 0;
	double totalSNEJ = 0;
	double totalOxygen = 0;
	double totalSilicon = 0;
	double totalIron = 0;
	double totalNickel = 0;

	/* Iterate over leaf blocks */
	for (leaf = leafBlockIDs.begin(); leaf != leafBlockIDs.end(); ++leaf) {
		
		/* Get block bounds */
		double blockBounds[2][MESH_MDIM];
		mesh.getBlkBoundBox(*leaf, blockBounds); //In terms of coordinates

		/* Iterate over cell in leaf blocks */
		for (int i = 0; i < nCellsX; i++) {
			for (int j = 0; j < nCellsY; j++) {
				double cellBounds[2][MESH_MDIM];
				getCellBoundBoxInMain(i, j, nCellsArr, blockBounds, cellBounds);

				 // Sum gravitational and kinetic energy inside cell to get total energy. 
				double totalEnergy = blocksMain[*leaf][varIndex[1]][i][j][0] + 
									 blocksMain[*leaf][varIndex[2]][i][j][0];
				
				// Make sure we're only summing up energies inside of the companion
				double compMassScalar = blocksMain[*leaf][varIndex[4]][i][j][0];

				if (totalEnergy * compMassScalar < 0) {

					double volume = PI * (cellBounds[UPPER][IAXIS]*cellBounds[UPPER][IAXIS]
										- cellBounds[LOWER][IAXIS]*cellBounds[LOWER][IAXIS]) 
									   * (cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);
		
					//double minX = cellBounds[LOWER][IAXIS];
					//double maxX	= cellBounds[UPPER][IAXIS];

					double density = blocksMain[*leaf][varIndex[0]][i][j][0];

					double snejMassScalar = blocksMain[*leaf][varIndex[3]][i][j][0];

					// double hydrogenScalar = blocksMain[*leaf][varIndex[7]][i][j][0];
					// double heliumScalar   = blocksMain[*leaf][varIndex[8]][i][j][0];
					double oxygenScalar   = blocksMain[*leaf][varIndex[9]][i][j][0];
					double siliconScalar  = blocksMain[*leaf][varIndex[10]][i][j][0];
					double ironScalar     = blocksMain[*leaf][varIndex[11]][i][j][0];
					double nickelScalar   = blocksMain[*leaf][varIndex[12]][i][j][0];
					
					// double temperature = blocksMain[*leaf][varIndex[5]][i][j][0];

					totalSNEJ    += density * volume * snejMassScalar;
					totalOxygen  += density * volume * oxygenScalar;
					totalSilicon += density * volume * siliconScalar;
					totalIron    += density * volume * ironScalar;
					totalNickel  += density * volume * nickelScalar;
					
				}

			}
		}
	}

	// std::cout << totalMass << std::endl;
	
	//Print out results.
	printOut(outProgramFileName, mesh.getTime(), totalIron, totalNickel);

	std::cout << "------------------------------" << std::endl;
	std::cout << "Checkpoint File: " << checkNum << std::endl;
	std::cout << "Physical Time: " << mesh.getTime() << std::endl;
	std::cout << "Total Captured Ejecta: " << totalSNEJ << std::endl;
	std::cout << "Total Oxygen: " << totalOxygen << std::endl;
	std::cout << "Total Silicon: " << totalSilicon << std::endl;
	std::cout << "Total Iron: " << totalIron << std::endl;
	std::cout << "Total Nickel: " << totalNickel << std::endl;
	std::cout << "Total Nickel (M_solar): " << totalNickel / 1.9891e33 << std::endl;
	std::cout << "Iron+Nickel (M_solar): " << (totalIron + totalNickel) / 1.9891e33 << std::endl;
	std::cout << "------------------------------" << std::endl;
	std::cout << std::endl;

	for (int i = 0; i < numberOfBlocks; i++) {
		deallocateVarArray(blocksMain[i], mesh.getNUserVars(), mesh.getNcells());
	}

	for (int i = 0; i < nVars; i++) {
		delete[] vars[i];
	}
	delete[] vars;
	delete[] varIndex;

	return 0;

}
