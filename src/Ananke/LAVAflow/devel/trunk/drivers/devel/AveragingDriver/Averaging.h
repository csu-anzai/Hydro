#include<vector>
#include "FlashAmrMesh.h"

struct averageData 
{
	char dir;
	int nbins;
	double binsize;
	vector<double> binAverages;
	

};


class ShellAverage
{
	public:
		ShellAverage();
		void calcShellAverage(int, int);
		void calcShellAverage(int, double);
		void calcShellAverage(int, vector<double>&);
		~ShellAverage();


	private:


};
