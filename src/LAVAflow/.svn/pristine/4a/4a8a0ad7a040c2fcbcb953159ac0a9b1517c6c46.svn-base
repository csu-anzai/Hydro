#include <list>
#include <iostream>
using namespace std;

int main(void)
{
	// Create list of elements
	list<int> L;
	L.push_back(1);
	L.push_back(2);
	L.push_back(-1);
	L.push_back(-1);
	L.push_back(3);
	L.push_back(-1);
	L.push_back(4);
	L.push_back(-1);


	// Iterate over the list and remove all elements <0
	list<int>::iterator it = L.begin();
	int whileCounter = 0;
	while(it != L.end())
	{
		if((*it) < 0)
		{
			it = L.erase(it);
		}
		else
		{
			++it;
		}
		whileCounter++;
	}


	// Print results
	std::cout<<"Number of times through while loop: "<<whileCounter<<std::endl;
	std::cout<<"List elements: "<<std::endl;
	for(it = L.begin(); it!=L.end(); ++it)
	{
		std::cout<<*it<<std::endl;
	}



}