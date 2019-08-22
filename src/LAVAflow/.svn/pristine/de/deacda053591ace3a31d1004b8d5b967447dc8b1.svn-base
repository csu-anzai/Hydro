//****h* LAVAflow/LayerAverage
//  NAME
//    Driver
//
//  DESCRIPTION
//    LayerAveraging driver that computes the average of different layers
//    for specified quantaties on the mesh. These different layers are 
//    seperated by bins that are either user-defined or automatically generated
//    based on user requests for a specific number or spacing. These
//    bin averages are then stored in a classData structure, along with
//    other information used in the computation, such as the volume of each
//    bin, the number of bins, and the total volume of the domain. 
//
//  METHODS
//     LayerAverage::LayerAverage
//     Driver::~LayerAverage
//     Driver::calcLayerAverage
//****



#include<vector>
#include "FlashAmrMesh.h"
using namespace std;

struct classData 
{
	char dir;
	int nbins;
	double binsize;
    vector <double> binFaces;
	vector <double> binAvg;
    vector <double> binVol;
	

};


class LayerAverage
{




	public:
		LayerAverage();
		void calcLayerAverage(FlashAmrMesh*,int, int, int);
		void calcLayerAverage(FlashAmrMesh*,int, double, int);
		void calcLayerAverage(FlashAmrMesh* , int, double*, int);
         
		~LayerAverage();
        classData averageData;
    private:


	


};

