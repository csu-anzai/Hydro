
#include "CommonSNBinary.h"
using namespace std;

/* TODO: Look at resolution effects between the two models. 
	 	 Make table of each of the stripped masses
	 	 Make table of resolutions based on radii of our models, I think.*/

void printOut(char* fileName, double checkTime, double mass) {
	std::ofstream outFile;
	outFile.open(fileName);

	outFile << checkTime << "," << mass << std::endl;
	std::cout << fileName << std::endl;
	std::cout << checkTime << std::endl;
	std::cout << mass << std::endl;

	outFile.close();
}

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv) {

	modelContainer models[NUMMODELS];
	readModels(models);

	/* Declare variables to be read */
	int nVars = 5;
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
	
	double gToSolarMass = 5.02739933e-34;

	/* --- Specify variables --- */

	MODEL currentModel;
	int checkNum;
	if (argc == 1) {
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
	}
	else if (argc == 3) {
		currentModel = (MODEL)atoi(argv[1]);
		checkNum = atoi(argv[2]);
	}

	strcpy(vars[0],"dens");
	strcpy(vars[1],"gpot");
	strcpy(vars[2],"ener");
	strcpy(vars[3],"ms_2");
	strcpy(vars[4],"ms_3");
	strcpy(outFileExtraName, "_stripped_mass.txt");

	/* Useful Variable names:
		dens = Density
		gpot = Gravitational Potential
		ener = Total energy (Kinetic + Internal)
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
	strcat(outProgramFileName, "strippedMass/");

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

	double totalMass = 0;
	double snejMass = 0;

	/* Iterate over leaf blocks */
	for (leaf = leafBlockIDs.begin(); leaf != leafBlockIDs.end(); ++leaf) {
		
		/* Get block bounds */
		double blockBounds[2][MESH_MDIM];
		mesh.getBlkBoundBox(*leaf, blockBounds);


		/* Iterate over cell in leaf blocks */
		for (int i = 0; i < nCellsX; i++) {
			for (int j = 0; j < nCellsY; j++) {
				double cellBounds[2][MESH_MDIM];
				getCellBoundBoxInMain(i, j, nCellsArr, blockBounds, cellBounds);

				/* Sum gravitational and kinetic energy inside cell to get total energy. */
				double totalEnergy = blocksMain[*leaf][varIndex[1]][i][j][0] + 
									 blocksMain[*leaf][varIndex[2]][i][j][0];

				/* Currently checking unbound mass July 31st, 2014 */

				

				double volumeOuter = PI * cellBounds[UPPER][IAXIS] * cellBounds[UPPER][IAXIS] * 
								(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);
				double volumeInner = PI * cellBounds[LOWER][IAXIS] * cellBounds[LOWER][IAXIS] * 
								(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);

				double volume = volumeOuter - volumeInner;

				if (blocksMain[*leaf][varIndex[3]][i][j][0] > 0.9) {
					snejMass += blocksMain[*leaf][varIndex[0]][i][j][0] * blocksMain[*leaf][varIndex[3]][i][j][0] * volume;
				}

				if (totalEnergy > 0) {

					double volumeOuter = PI * cellBounds[UPPER][IAXIS] * cellBounds[UPPER][IAXIS] * 
									(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);
					double volumeInner = PI * cellBounds[LOWER][IAXIS] * cellBounds[LOWER][IAXIS] * 
									(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);

					double volume = volumeOuter - volumeInner;
					
					double density = blocksMain[*leaf][varIndex[0]][i][j][0];
					
					double compMassScalar = blocksMain[*leaf][varIndex[4]][i][j][0];

					// double area = (cellBounds[UPPER][IAXIS] - cellBounds[LOWER][IAXIS]) *
					// 			  (cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);
					
					// double r = cellBounds[LOWER][IAXIS];

					// Volume element = dTheta (which is 4pi) * density of cell * radius * dr * dz
					// which is equal to dTheta * density * radius * dx * dy
					// Analogous to area of torus via integral calculation.

					// totalMass += compMassScalar * 4 * PI * density * r * area;
					totalMass += density * compMassScalar * volume;

				}

			}
		}
		
	}

	// std::cout << "Initial mass: " << 
	//std::cout << "Total remaining mass: " << totalMass << " (g)" << std::endl;
	//std::cout << "Total remaining mass: " << totalMass * gToSolarMass << " (m_solar)" << std::endl;
	//std::cout << "Simulation time: " << mesh.getTime() << std::endl;
	//std::cout << "Unbound mass: " << totalMass * gToSolarMass << " (m_solar)" << std::endl;
	
	std::cout << "Output: " << mesh.getTime() << " " << totalMass * gToSolarMass << std::endl;
	std::cout << std::endl;

	std::cout << "Snej total mass: " << snejMass << std::endl;


	printOut(outProgramFileName, mesh.getTime(), totalMass * gToSolarMass);

	for (int i = 0; i < numberOfBlocks; i++) {
		deallocateVarArray(blocksMain[i], mesh.getNUserVars(), mesh.getNcells());
	}

	for (int i = 0; i < nVars; i++) {
		delete[] vars[i];
	}
	delete[] vars;
	delete[] varIndex;

	// mesh.closeFile();

	return 0;

}
