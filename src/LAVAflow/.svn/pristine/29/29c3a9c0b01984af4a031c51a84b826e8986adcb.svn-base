//TEST TEST TEST
#include <iostream>
#include <string>
#include <vector>

#include "FlashAmrMesh.h"
#include "Driver_main.h"
#include "LayerAverage.h"
using namespace std;

int LAVAFLOW_DRIVER_BASENM(int argc, char** argv)
{

	// LavaMPI mpiObj(argc, argv);

    string filename("/home/rjl09c/AnankeREPO/LAVAflow/devel/trunk/drivers/files/testData/cartesianShells/cartesianShells_hdf5_chk_0000");
	vector<string> meshVars = {"vari"};
	FlashAmrMesh mesh(filename, meshVars);
    int densIndex = mesh.findVarIndex("vari");
    LayerAverage trial1;
    trial1.calcLayerAverage(&mesh, IAXIS, 4, densIndex);





}
