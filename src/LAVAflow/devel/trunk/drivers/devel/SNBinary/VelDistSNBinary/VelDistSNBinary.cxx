
#include "CommonSNBinary.h"
using namespace std;

/* TODO: Look at resolution effects between the two models. 
	 	 Make table of each of the stripped masses
	 	 Make table of resolutions based on radii of our models, I think.*/

void printOut(char* fileName, double* velBinValues, double* velValuesComp, double* velValuesSnej, 
			  double* velValuesH1, double* velValuesHe4, double* velValuesO16, double* velValuesSi28,
			  double* velValuesFe52, int bins) {
	std::ofstream outFile;
	outFile.open(fileName);

	for (int i = 0; i < bins; i++) {
		outFile << velBinValues[i] << "," << 
				   velValuesComp[i] << "," << 
				   velValuesSnej[i] << "," << 
				   velValuesH1[i] << "," << 
				   velValuesHe4[i] << "," << 
				   velValuesO16[i] << "," << 
				   velValuesSi28[i] << "," << 
				   velValuesFe52[i] << std::endl;
	}

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
	
	int bins = 1000;
	// For linear scale, use 1e6
	// double minVel = 1e6;

	// For log scale, use 1e0
	double minVel = 1e0;
	double maxVel = 5e9;
	
	
	double* velBinValues;
	velBinValues = new double[bins];

	// // Linear spacing
	// double velSpacing = (maxVel - minVel)/double(bins);
	// for (int i = 0; i < bins; i++) {
	// 	velBinValues[i] = minVel + i * velSpacing;
	// }

	// Log spacing
	double velSpacing = pow((maxVel/minVel),(1/double(bins)));
	for (int i = 0; i < bins; i++) {
		velBinValues[i] = minVel * pow(velSpacing,i);
	}

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

	strcpy(vars[0],"dens");
	strcpy(vars[1],"gpot");
	strcpy(vars[2],"ener");
	strcpy(vars[3],"ms_2");
	strcpy(vars[4],"ms_3");
	strcpy(vars[5],"velx");
	strcpy(vars[6],"vely");
	strcpy(vars[7],"velz");

	strcpy(vars[8] ,"H1  ");
	strcpy(vars[9] ,"He4 ");
	strcpy(vars[10],"O16 ");
	strcpy(vars[11],"Si28");
	strcpy(vars[12],"Fe52");




	strcpy(outFileExtraName, "_vel_dist_log.txt");

	/* Useful Variable names:
		dens = Density
		gpot = Gravitational Potential
		ener = total energy 
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
	strcat(outProgramFileName, "velDist/");

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
	double* velHistogramComp;
	double* velHistogramSnej;
	double* velHistogramH1;
	double* velHistogramHe4;
	double* velHistogramO16;
	double* velHistogramSi28;
	double* velHistogramFe52;
	velHistogramComp = new double[bins];
	velHistogramSnej = new double[bins];
	velHistogramH1 = new double[bins];
	velHistogramHe4 = new double[bins];
	velHistogramO16 = new double[bins];
	velHistogramSi28 = new double[bins];
	velHistogramFe52 = new double[bins];

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

				/* Sum gravitational and total energy inside cell to get the binding energy. */
				double totalEnergy = blocksMain[*leaf][varIndex[1]][i][j][0] + 
									 blocksMain[*leaf][varIndex[2]][i][j][0];

				if (totalEnergy > 0) {
					/* We're outside of the potential well */

					double volumeOuter = PI * cellBounds[UPPER][IAXIS] * cellBounds[UPPER][IAXIS] * 
									(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);
					double volumeInner = PI * cellBounds[LOWER][IAXIS] * cellBounds[LOWER][IAXIS] * 
									(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);

					double volume = volumeOuter - volumeInner;
					
					double density = blocksMain[*leaf][varIndex[0]][i][j][0];
					
					double snejMassScalar = blocksMain[*leaf][varIndex[3]][i][j][0];
					double compMassScalar = blocksMain[*leaf][varIndex[4]][i][j][0];
					double hydrogenScalar = blocksMain[*leaf][varIndex[8]][i][j][0];
					double heliumScalar   = blocksMain[*leaf][varIndex[9]][i][j][0];
					double oxygenScalar   = blocksMain[*leaf][varIndex[10]][i][j][0];
					double siliconScalar  = blocksMain[*leaf][varIndex[11]][i][j][0];
					double ironScalar     = blocksMain[*leaf][varIndex[12]][i][j][0];

					// double temperature = blocksMain[*leaf][varIndex[5]][i][j][0];

					double velx = blocksMain[*leaf][varIndex[5]][i][j][0];
					double vely = blocksMain[*leaf][varIndex[6]][i][j][0];
					double velz = blocksMain[*leaf][varIndex[7]][i][j][0];
					double velMag = sqrt(velx*velx + vely*vely + velz*velz);

					// totalMass += density * compMassScalar * volume;
					/* Add to the correct point on the mass histogram, which is a function of temperature. */
					int i = 0;
					while (velMag >= velBinValues[i]) {
						i++;
					}

					velHistogramComp[i] += density * compMassScalar * volume;
					velHistogramSnej[i] += density * snejMassScalar * volume;
					velHistogramH1[i]   += density * hydrogenScalar * volume;
					velHistogramHe4[i]  += density * heliumScalar * volume;
					velHistogramO16[i]  += density * oxygenScalar * volume;
					velHistogramSi28[i] += density * siliconScalar * volume;
					velHistogramFe52[i] += density * ironScalar * volume;
				
				}
			}
		}
	}

	// std::cout << totalMass << std::endl;
	
	//Print out results.
	printOut(outProgramFileName, velBinValues, velHistogramComp, velHistogramSnej, 
			 velHistogramH1, velHistogramHe4, velHistogramO16, velHistogramSi28,
			 velHistogramFe52, bins);

	std::cout << "Velocity, Companion, Snej, H1, He4, O16, Si28, Fe52" << std::endl;
	for (int i = 0; i < bins; i++) {
		std::cout << velBinValues[i] << "," << 
				     velHistogramComp[i] << "," << 
				     velHistogramSnej[i] << "," << 
				     velHistogramH1[i] << "," << 
				     velHistogramHe4[i] << "," << 
				     velHistogramO16[i] << "," << 
				     velHistogramSi28[i] << "," << 
				     velHistogramFe52[i] << std::endl;
	}

	for (int i = 0; i < numberOfBlocks; i++) {
		deallocateVarArray(blocksMain[i], mesh.getNUserVars(), mesh.getNcells());
	}

	for (int i = 0; i < nVars; i++) {
		delete[] vars[i];
	}
	delete[] vars;
	delete[] varIndex;
	delete[] velBinValues;
	delete[] velHistogramComp;
	delete[] velHistogramSnej;
	delete[] velHistogramH1;
	delete[] velHistogramHe4;
	delete[] velHistogramSi28;
	delete[] velHistogramO16;
	delete[] velHistogramFe52;

	return 0;

}
