#include <list>
#include <iostream>
using namespace std;

int main(void)
{

	list<int> data;
	list<int>::iterator fIt;
	list<int>::reverse_iterator rIt;
	list<int>::iterator finalIt;
	int count = 0;


	// Add some UNQUE elements to list
	data.push_back(0);
	data.push_back(1);
	data.push_back(2);
	data.push_back(3);
	data.push_back(4);

	// Loop from the front and the back simultaneously to find an element
	int valToFind = 4;

	for(fIt = data.begin(), rIt = data.rbegin(); fIt!=data.end(), rIt!=data.rend(); ++fIt, ++rIt)
	{
		count++;
		if((*fIt)==valToFind)
		{
			cout<<"Found value with forward iterator!"<<endl;
			finalIt = fIt;
			break;
		}

		if((*rIt)==valToFind)
		{
			cout<<"Found value with reverse iterator!"<<endl;
			finalIt = fIt;
			break;
		}
	}

	cout<<"Number of times through loop: "<<count<<endl;

	// Removing element




}