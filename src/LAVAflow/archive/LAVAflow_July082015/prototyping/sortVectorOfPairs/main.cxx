#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

bool myfunction( const pair<int, char>& i, const pair<int, char>& j ) {
    if( i.first < j.first ) return true;
    if( j.first < i.first ) return false;
    return j.second < i.second;
}

int main(void)
{
	vector<double> originalData;
	vector<pair<double,int> > toBeSorted;

	// Set original data
	originalData.push_back(1.0);
	originalData.push_back(0.0);
	originalData.push_back(2.0);


	// construct toBeSorted
	for(int ind=0; ind<originalData.size(); ind++)
	{
		toBeSorted.push_back(pair<double,int>(originalData[ind],ind));
	}


	// Print unsorted toBeSorted
	std::cout<<"Unsorted toBeSorted:"<<std::endl;
	for (std::vector<pair<double,int> >::iterator i = toBeSorted.begin(); i != toBeSorted.end(); ++i)
	{
		std::cout<<"\t"<<(*i).first<<"\t"<<(*i).second<<std::endl;
	}

	sort(toBeSorted.begin(), toBeSorted.end(), myfunction);


	std::cout<<"Sorted toBeSorted:"<<std::endl;
	for (std::vector<pair<double,int> >::iterator i = toBeSorted.begin(); i != toBeSorted.end(); ++i)
	{
		std::cout<<"\t"<<(*i).first<<"\t"<<(*i).second<<std::endl;
	}



}