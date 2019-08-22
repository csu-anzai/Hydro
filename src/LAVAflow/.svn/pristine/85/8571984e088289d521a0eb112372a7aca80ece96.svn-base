#include "CommonSNBinary.h"

void getCheckNames(modelContainer &outModel, int modelBegin, int modelEnd) {
	for (int i = modelBegin; i <= modelEnd; i++) {
		//Make string to build output file name.
		char checkpointFileName[MAX_FILENAME_LENGTH];
		strcpy(checkpointFileName,"snb_hdf5_chk_");
		//Fill in the preceding zeroes.
		if (i < 1000) {
			strcat(checkpointFileName,"0");
		}
		if (i < 100) {
			strcat(checkpointFileName,"0");
		}
		if (i < 10) {
			strcat(checkpointFileName,"0");
		}
		char intStr[5];
		sprintf(intStr, "%d", i);
		strcat(checkpointFileName,intStr);
		//std::cout << checkpointFileName << std::endl;

		//Copy the string into the array.
		strcpy(outModel.checkpoint[i], checkpointFileName);
	}
}

void readModels(modelContainer *outModels) {

	/* MS38 - l50 */

	//Folder information
	strcpy(outModels[MS38_L50].modelFolder,  "w7_ms38_l50-20070509/");
	//Checkpoint file information
	getCheckNames(outModels[MS38_L50], MS38_L50_I, MS38_L50_F);

	/* MS38 - l100 */
	strcpy(outModels[MS38_L100].modelFolder, "w7_ms38_l100-20070509/");
	getCheckNames(outModels[MS38_L100], MS38_L100_I, MS38_L100_F);

	/* MS38 - l200 */
	strcpy(outModels[MS38_L200].modelFolder, "w7_ms38_l200-20070509/");
	getCheckNames(outModels[MS38_L200], MS38_L200_I, MS38_L200_F);

	/* MS54 - l200 */
	strcpy(outModels[MS54].modelFolder,      "w7_ms54_l200-20070509/");
	getCheckNames(outModels[MS54], MS54_I, MS54_F);

	/* MS63 - l200 */
	strcpy(outModels[MS63].modelFolder,      "w7_ms63_l200-20070509/");
	getCheckNames(outModels[MS63], MS63_I, MS63_F);

	/* MS7 - l50 */
	strcpy(outModels[MS7_L50].modelFolder,   "w7_ms7_l50-20070509/");
	getCheckNames(outModels[MS7_L50], MS7_L50_I, MS7_L50_F);

	/* MS7 - l100 */
	strcpy(outModels[MS7_L100].modelFolder,  "w7_ms7_l100-20070509/");
	getCheckNames(outModels[MS7_L100], MS7_L100_I, MS7_L100_F);

	/* MS7 - l200 */
	strcpy(outModels[MS7_L200].modelFolder,  "w7_ms7_l200-20070509/");
	getCheckNames(outModels[MS7_L200], MS7_L200_I, MS7_L200_F);

	/* SG - L50 */
	strcpy(outModels[SG_L50].modelFolder,    "w7_sg_l50-20070509/");
	getCheckNames(outModels[SG_L50], SG_L50_I, SG_L50_F);

	/* SG - L100 */
	strcpy(outModels[SG_L100].modelFolder,   "w7_sg_l100-20070509/");
	getCheckNames(outModels[SG_L100], SG_L100_I, SG_L100_F);

	/* SG - L200 */
	strcpy(outModels[SG_L200].modelFolder,   "w7_sg_l200-20070509/");
	getCheckNames(outModels[SG_L200], SG_L200_I, SG_L200_F);

	/* SY319 - L50 */
	strcpy(outModels[SY319_L50].modelFolder, "w7_sy319_l50-20070509/");
	getCheckNames(outModels[SY319_L50], SY319_L50_I, SY319_L50_F);

	/* SY319 - L100 */
	strcpy(outModels[SY319_L100].modelFolder,"w7_sy319_l100-20070509/");
	getCheckNames(outModels[SY319_L100], SY319_L100_I, SY319_L100_F);

	/* SY319 - L200 */
	strcpy(outModels[SY319_L200].modelFolder,"w7_sy319_l200-20070509/");
	getCheckNames(outModels[SY319_L200], SY319_L200_I, SY319_L200_F);

	/* SY428 */
	strcpy(outModels[SY428].modelFolder,     "w7_sy428_l200-20070509/");
	getCheckNames(outModels[SY428], SY428_I, SY428_F);
}

void getCellBoundBoxInMain(int i,int j,/*int k,*/ int *nCellsArray, 
					       double blockBounds[2][MESH_MDIM], double cellBounds[2][MESH_MDIM]) {
	
	double boxLengthX = blockBounds[UPPER][IAXIS] - blockBounds[LOWER][IAXIS];
	double boxLengthY = blockBounds[UPPER][JAXIS] - blockBounds[LOWER][JAXIS];
	// double boxLengthZ = blockBounds[UPPER][KAXIS] - blockBounds[LOWER][KAXIS];

	double cellLengthX = boxLengthX / double(nCellsArray[IAXIS]);
	double cellLengthY = boxLengthY / double(nCellsArray[JAXIS]);
	// double cellLengthZ = boxLengthZ / nCellsArray[KAXIS];

	cellBounds[LOWER][IAXIS] = blockBounds[LOWER][IAXIS] + double(i) * cellLengthX;
	cellBounds[LOWER][JAXIS] = blockBounds[LOWER][JAXIS] + double(j) * cellLengthY;
	// cellBounds[LOWER][KAXIS] = blockBounds[LOWER][KAXIS] + (j*nCellsArray[KAXIS] + i) * cellLengthZ;

	cellBounds[UPPER][IAXIS] = blockBounds[LOWER][IAXIS] + double(i+1) * cellLengthX;
	cellBounds[UPPER][JAXIS] = blockBounds[LOWER][JAXIS] + double(j+1) * cellLengthY;
	// cellBounds[UPPER][KAXIS] = blockBounds[UPPER][KAXIS] + (j*nCellsArray[KAXIS] + i + 1) * cellLengthZ;

}

/* 
	Functions not used in the AMR class, but still needed for functions in main.
	Most notably allocating data arrays outside of the class.
*/

void allocateVarArray(double****& inData, int numVars, int *nCellsVec) {

	inData = new double***[numVars];

	for (int i = 0; i < numVars; i++) {
		inData[i] = new double**[nCellsVec[IAXIS]];

		for (int j = 0; j < nCellsVec[IAXIS]; j++) {
			inData[i][j] = new double*[nCellsVec[JAXIS]];

			for (int k = 0; k < nCellsVec[JAXIS]; k++) {
				inData[i][j][k] = new double[nCellsVec[KAXIS]];
			}
		}
	}
}

void deallocateVarArray(double****& inData, int numVars, int *nCellsVec) {

	if (inData != NULL) {
		for (int i = 0; i < numVars; i++) {
			for (int j = 0; j < nCellsVec[IAXIS]; j++) {
				for (int k = 0; k < nCellsVec[JAXIS]; k++) {
					delete[] inData[i][j][k];
				}
				delete[] inData[i][j];
			}
			delete[] inData[i];
		}
		delete[] inData;
	}
	else {
		if (FORCEDEALLOCATECRASH) {
			std::cout << "[flash_reader] Failure in deallocateVarArray function outside of class." << std::endl;
			exit(EXIT_FAILURE);	
		}
	}
}
