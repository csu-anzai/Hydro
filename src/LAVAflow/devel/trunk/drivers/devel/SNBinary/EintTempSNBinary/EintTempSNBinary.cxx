
#include "CommonSNBinary.h"
using namespace std;

/* TODO: Look at resolution effects between the two models. 
	 	 Make table of each of the stripped masses
	 	 Make table of resolutions based on radii of our models, I think.*/

void printOut(char* fileName, double* temperatureValues, double* eintValuesCompanion, 
			  double* eintValuesEjecta, double* eintValuesEjectaUnaffected, int bins) {
	
	std::ofstream outFile;
	outFile.open(fileName);

	for (int i = 0; i < bins; i++) {
		outFile << temperatureValues[i] << "," << eintValuesCompanion[i] << ","
				<< eintValuesEjecta[i] << "," << eintValuesEjectaUnaffected[i] << std::endl;
	}

	outFile.close();
}

void printOutBoth(char* fileName, double* temperatureValues, double* eintValues, 
			 double** massValues, int bins){
	std::ofstream outFile;
	outFile.open(fileName);

	for (int i = 0; i < bins; i++) {
		for (int j = 0; j < bins; j++) {
			outFile << temperatureValues[i] << "," << eintValues[j] << "," << massValues[i][j] << std::endl;
		}
	}

	outFile.close();
}

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv) {

	modelContainer models[NUMMODELS];
	readModels(models);

	/* Declare variables to be read */
	int nVars = 7;
	char** vars;
	int* varIndex;
	vars = new char*[nVars];
	for (int i = 0; i < nVars; i++) {
		vars[i] = new char[MESH_VAR_STRING_SIZE+1];
	}
	varIndex = new int[nVars];
	char inFileName[MAX_FILENAME_LENGTH];
	char inProgramFileName[MAX_FILENAME_LENGTH];
	
	char outFileName[MAX_FILENAME_LENGTH];
	char outProgramFileName[MAX_FILENAME_LENGTH];
	char outFileExtraName[MAX_FILENAME_LENGTH];

	// char outFileNameEint[MAX_FILENAME_LENGTH];
	// char outProgramFileNameEint[MAX_FILENAME_LENGTH];
	// char outFileExtraNameEint[MAX_FILENAME_LENGTH];

	// char outFileNameBoth[MAX_FILENAME_LENGTH];
	// char outProgramFileNameBoth[MAX_FILENAME_LENGTH];
	// char outFileExtraNameBoth[MAX_FILENAME_LENGTH];
	
	/* Parameters */
	
	// MODEL currentModel = MS38_L50;
	// MODEL currentModel = MS38_L100;
	MODEL currentModel = MS38_L200;
	// MODEL currentModel = MS54;
	// MODEL currentModel = MS63;
	// MODEL currentModel = MS7_L200;
	// MODEL currentModel = SG_L50;
	// MODEL currentModel = SG_L100;
	// MODEL currentModel = SG_L200;
	// MODEL currentModel = SY319_L50;
	// MODEL currentModel = SY319_L100;
	// MODEL currentModel = SY319_L200;
	// MODEL currentModel = SY428;

	/* INITIAL TIMES */
	// int checkNum = MS38_L50_I;
	// int checkNum = MS38_L100_I;
	// int checkNum = MS38_L200_I;
	// int checkNum = MS54_I;
	// int checkNum = MS63_I;
	// int checkNum = MS7_L200_I;
	// int checkNum = SG_L50_I;
	// int checkNum = SG_L100_I;
	// int checkNum = SG_L200_I;
	// int checkNum = SY319_L50_I;
	// int checkNum = SY319_L100_I;
	// int checkNum = SY319_L200_I;
	// int checkNum = SY428_I;

	/* HOTTEST TIME (not exactly, but close to it for every model) 
	   The times for the non-L200 models are chosen based on the 
	   L200 60th checkpoint times because the simulation times are not 
	   tethered to the checkpoint times.
	*/
	// int checkNum = 49; //MS38_L50
	// int checkNum = 53; //MS38_L100
	int checkNum = 60; //MS38_L200

	// int checkNum = 42; //SG_L50
	// int checkNum = 43; //SG_L100
	// int checkNum = 60; //SG_L200

	// int checkNum = 49; //SY319_L50
	// int checkNum = 54; //SY319_L100
	// int checkNum = 60; //SY319_L200


	if (argc > 1) {
		currentModel = (MODEL)atoi(argv[1]);
		checkNum = atoi(argv[2]);
	}

	double minTemp = 1e3;
	double maxTemp = 1.2e10;
	int bins = 100;
	
	// // linear spacing
	// double tempSpacing = (maxTemp - minTemp)/double(bins);
	// double* temperatureValues;
	// temperatureValues = new double[bins];
	// for (int i = 0; i < bins; i++) {
	// 	temperatureValues[i] = minTemp + i * tempSpacing;
	// }

	// Log spacing
	double* temperatureValues;
	temperatureValues = new double[bins];
	double tempSpacing = pow((maxTemp/minTemp),(1/double(bins)));
	for (int i = 0; i < bins; i++) {
		temperatureValues[i] = minTemp * pow(tempSpacing,i);
	}

	/* --- Specify variables --- */


	strcpy(vars[0],"dens");
	strcpy(vars[1],"gpot");
	strcpy(vars[2],"ener");
	strcpy(vars[3],"ms_2");
	strcpy(vars[4],"ms_3");
	strcpy(vars[5],"temp");
	strcpy(vars[6],"eint");
	strcpy(outFileExtraName, "_eint_temp_L200.txt");
	// strcpy(outFileExtraNameEint, "_mass_eint_dens.txt");
	// strcpy(outFileExtraNameBoth, "_mass_both.txt");

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
	// strcpy(outFileNameEint, models[currentModel].checkpoint[checkNum]);
	// strcpy(outFileNameBoth, models[currentModel].checkpoint[checkNum]);
	
	strcpy(inProgramFileName, inFolderLocation);
	strcat(inProgramFileName, models[currentModel].modelFolder);
	strcat(inProgramFileName, inFileName);

	strcat(outFileName, outFileExtraName);
	strcpy(outProgramFileName, outFolderLocation);
	strcat(outProgramFileName, models[currentModel].modelFolder);
	strcat(outProgramFileName, "eintTemp/");
	strcat(outProgramFileName, outFileName);

	// strcat(outFileNameEint, outFileExtraNameEint);
	// strcpy(outProgramFileNameEint, outFolderLocation);
	// strcat(outProgramFileNameEint, models[currentModel].modelFolder);
	// strcat(outProgramFileNameEint, "massTemp/");
	// strcat(outProgramFileNameEint,outFileNameEint);

	// strcat(outFileNameBoth, outFileExtraNameBoth);
	// strcpy(outProgramFileNameBoth, outFolderLocation);
	// strcat(outProgramFileNameBoth, models[currentModel].modelFolder);
	// strcat(outProgramFileNameBoth, "massTemp/");
	
	// strcat(outProgramFileNameBoth,outFileNameBoth);
	
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "WARNING: THIS ANALYSIS PROGRAM NOT YET BUILT FOR MPI" << std::endl;
	std::cout << "Input file: " << inProgramFileName << std::endl;
	std::cout << "Output file Temp: " << outProgramFileName << std::endl;
	// std::cout << "Output file Eint-Dens: " << outProgramFileNameEint << std::endl;
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

	double* eintHistogramCompanion; 
	double* eintHistogramEjectaAll;  //All of the ejecta, including that which interacts with companion.
	double* eintHistogramEjectaUnaffected; //Bottom half of the domain, doubled.
	
	eintHistogramCompanion = new double[bins];
	eintHistogramEjectaAll = new double[bins];
	eintHistogramEjectaUnaffected = new double[bins];

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

				/* We're restricting the domain of interest to a subdomain surrounding
				   the most energetic part of the domain, defined with minX, maxX, minY, maxY. */

				// if (   cellBounds[LOWER][IAXIS] >= minX && cellBounds[UPPER][IAXIS] <= maxX
				// 	&& cellBounds[LOWER][JAXIS] >= minY && cellBounds[UPPER][JAXIS] <= maxY) {
					/* We're within the region of interest */

					double volumeOuter = PI * cellBounds[UPPER][IAXIS] * cellBounds[UPPER][IAXIS] * 
									(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);
					double volumeInner = PI * cellBounds[LOWER][IAXIS] * cellBounds[LOWER][IAXIS] * 
									(cellBounds[UPPER][JAXIS] - cellBounds[LOWER][JAXIS]);

					double volume = volumeOuter - volumeInner;
					
					double density = 		blocksMain[*leaf][varIndex[0]][i][j][0];
					double snejMassScalar = blocksMain[*leaf][varIndex[3]][i][j][0];
					double compMassScalar = blocksMain[*leaf][varIndex[4]][i][j][0];
					double temperature = 	blocksMain[*leaf][varIndex[5]][i][j][0];
					double eint = 			blocksMain[*leaf][varIndex[6]][i][j][0];

					// totalMass += density * compMassScalar * area;

					int ii = 0;
					while (temperature > temperatureValues[ii] && ii < bins) {
						ii++;
					}
					eintHistogramCompanion[ii] += density * volume * eint * compMassScalar;
					eintHistogramEjectaAll[ii] += density * volume * eint * snejMassScalar;
					if (cellBounds[UPPER][JAXIS] <= 0) { //Only add the unaffected bottom half of the ejecta
						eintHistogramEjectaUnaffected[ii] += 2 * density * volume * eint * snejMassScalar;
					}

					// specific internal energy, energy density, and internal energy
					// energy density is energy/gram
					// specific energy is energy/volume
					// energy is the energy density * density * volume

					// /* Add to the correct point on the mass histogram, which is a function of temperature. */
					// int ii = 0;
					// while (temperature > temperatureValues[ii] && ii < bins) {
					// 	ii++;
					// }
					// /* ONLY USING THE COMPANION'S MASS CURRENTLY */
					// //massHistogram[i] += density * compMassScalar * area;
					// // eintHistogramCompanion[ii] += density * snejMassScalar * area;
					
					// if (temperature <= temperatureValues[bins-1] && temperature >= temperatureValues[0]) {
					// 	eintHistogramCompanion[ii] += density * volume;
					// }

					// int jj = 0;
					// while (eint*density > eintValues[jj] && jj < bins) {
					// 	jj++;
					// }
					// if (eint*density <= eintValues[bins-1] && eint*density >= eintValues[0]) {
					// 	massHistogramEint[jj] += density * volume;
					// }

					// if (temperature <= temperatureValues[bins-1] && temperature >= temperatureValues[0]
					// 	&& eint*density <= eintValues[bins-1] && eint*density >= eintValues[0]) {
					// 	massHistogramBoth[ii][jj] += density * volume;
					// }
				// }
			}
		}
	}
	
	//Print out results.
	printOut(outProgramFileName, temperatureValues, eintHistogramCompanion, 
			 eintHistogramEjectaAll, eintHistogramEjectaUnaffected, bins);
	//printOut(outProgramFileNameEint, eintValues, massHistogramEint, bins);
	//printOutBoth(outProgramFileNameBoth, temperatureValues, eintValues, massHistogramBoth, bins);

    std::cout << "Temperature, Companion Eint, Ejecta Eint, Unaffected Ejecta Eint" << std::endl;
	for (int i = 0; i < bins; i++) {
		std::cout << temperatureValues[i] << "," << eintHistogramCompanion[i] << ","
				  << eintHistogramEjectaAll[i] << "," << eintHistogramEjectaUnaffected[i] << std::endl;
	}

	for (int i = 0; i < numberOfBlocks; i++) {
		deallocateVarArray(blocksMain[i], mesh.getNUserVars(), mesh.getNcells());
	}

	for (int i = 0; i < nVars; i++) {
		delete[] vars[i];
	}
	delete[] vars;
	delete[] varIndex;
	delete[] temperatureValues;

	delete[] eintHistogramCompanion;
	delete[] eintHistogramEjectaUnaffected;
	delete[] eintHistogramEjectaAll;

	return 0;

}
