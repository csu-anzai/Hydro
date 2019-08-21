#include <list>
#include <iostream>
using namespace std;

class Element
{
public:
	int value;
	Element(int v)
	{
		this->value = v;
	}
};

int main(void)
{

	list<Element*> data;
	list<Element*>::iterator fIt;
	list<Element*>::reverse_iterator rIt;


	// Add some UNQUE elements to list

	Element* e0 = new Element(0);
	Element* e1 = new Element(1);
	Element* e2 = new Element(2);

	data.push_back(e0);
	data.push_back(e1);
	data.push_back(e2);

	// Loop from the front and the back simultaneously to find an element
	Element* valToFind = e1;

	cout<<"Removing element using list::remove"<<endl;
	data.remove(valToFind);


	cout<<"List contains: "<<endl;
	for(fIt=data.begin(); fIt!=data.end(); ++fIt)
	{
		cout<<(*fIt)->value<<endl;
	}

	cout<<"Is removed element deallocated?"<<endl;
	cout<<valToFind->value<<endl;





}