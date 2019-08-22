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


class ShellAverage
{
	public:
		ShellAverage();
		void calcShellAverage(FlashAmrMesh*,int, int, int);
		void calcShellAverage(FlashAmrMesh*,int, double, int);
		void calcShellAverage(FlashAmrMesh* , int, vector<double>&, int);
         
		~ShellAverage();
        classData averageData;
    private:


	


};
