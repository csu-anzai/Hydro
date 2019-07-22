#include "FlashAmrMesh.h"
#include "Driver_main.h"

/* Constants */
#define PI 3.14159265
#define MAX_CHECKPOINT 300

enum MODEL {MS38_L50, MS38_L100, MS38_L200, MS54, MS63, 
			MS7_L50, MS7_L100, MS7_L200, SG_L50, SG_L100, SG_L200, 
			SY319_L50, SY319_L100, SY319_L200, SY428, NUMMODELS};
// ------------------------------------------------------
// When adding a model above, be sure to add it before 
// NUMMODELS. 
// ------------------------------------------------------

//Initial checkpoints
const int MS38_L50_I=0;
const int MS38_L100_I=0;
const int MS38_L200_I=0;
const int MS54_I=0;
const int MS63_I=0;
const int MS7_L50_I=0;
const int MS7_L100_I=0;
const int MS7_L200_I=0;
const int SG_L50_I=0;
const int SG_L100_I=0;
const int SG_L200_I=0;
const int SY319_L50_I=0;
const int SY319_L100_I=0;
const int SY319_L200_I=0;
const int SY428_I=0;

//Final checkpoints
const int MS38_L50_F=172;
const int MS38_L100_F=184;
const int MS38_L200_F=213;
const int MS54_F=245; 
const int MS63_F=203; 
const int MS7_L50_F=141;
const int MS7_L100_F=146;
const int MS7_L200_F=208;
const int SG_L50_F=149; 
const int SG_L100_F=155; 
const int SG_L200_F=234; 
const int SY319_L50_F=165;
const int SY319_L100_F=186;
const int SY319_L200_F=205;
const int SY428_F=213;

	/* Change this based on where the files are located. */
//Desktop
const char inFolderLocation[MAX_FILENAME_LENGTH] = {"/data1/Data/SNB_2007_Data/"};
const char outFolderLocation[MAX_FILENAME_LENGTH] = {"/data1/Data/SNB_2007_Data/"};

//NERSC
// const char inFolderLocation[MAX_FILENAME_LENGTH] = {"/global/project/projectdirs/incite2/snb/2007/"};
// const char outFolderLocation[MAX_FILENAME_LENGTH] = {"/global/homes/p/pboehner/analysisResults/"};

enum ORIGINTYPE {USERDEF, HIGHESTDENSITY};
/* Above: Enum for choice of where the origin lies. 
     -Userdef: Use user-defined values, usually near the origin of the grid
     -HighestDensity: Let the code find the most dense point and set the origin where it's projected
        on the symmetry axis (Y).
*/



struct modelContainer{
	char checkpoint[MAX_CHECKPOINT][MAX_FILENAME_LENGTH];
	char modelFolder[MAX_FILENAME_LENGTH];
};

void getCheckNames(modelContainer &outModel, int modelBegin, int modelEnd);
void readModels(modelContainer *outModels);
void getCellBoundBoxInMain(int i,int j,/*int k,*/ int *nCellsArray, 
					       double blockBounds[2][MESH_MDIM], double cellBounds[2][MESH_MDIM]);
void allocateVarArray(double****& inData, int numVars, int *nCellsVec);
void deallocateVarArray(double****& inData, int numVars, int *nCellsVec);