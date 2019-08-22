#include "ScalarList.h"
#include <iostream>
#include <ostream>
#include <map>
#include <fstream>
#include <sstream>

namespace Util
{


bool ScalarList::addScalar(double value, const char* id)
{
	// this->identifiers.push_back(std::string(id));
	// this->values.push_back(value);
	// std::string idStr(id);
	this->dictionary.insert(std::pair<std::string,double>(id,value));
	this->numberOfEntries++;
}

double ScalarList::getScalar(const char* id)
{
	double value = 1e999;

	// Find the value in the map
	std::map<std::string,double>::iterator it;
	it = this->dictionary.find(id);

	if(it != this->dictionary.end())
	{
		value = (*it).second;
	}

	return value;

}

bool ScalarList::print(std::string filename)
{
	std::ofstream out;
	out.open(filename.c_str());

	if(!out.good())
	{
		std::cerr<<"[ScalarList::read] Unable to open file \""<<filename<<"\"!"<<std::endl;
		return false;
	}

	return this->print(out);
}


bool ScalarList::print(std::ostream& out)
{


	// Save old format flags
	std::ios_base::fmtflags ff;
	ff = out.flags();

	// Set output flags
	out.precision(15);
	out.setf(std::ios::scientific, std::ios::floatfield);

	// Loop over all elements and write them to the output
	// for(int i=0; i<numberOfEntries; i++)
	// {
	// 	out<<this->identifiers[i]<<this->delimiter<<this->values[i]<<std::endl;
	// }

	std::map<std::string,double>::iterator it;
	for(it = this->dictionary.begin(); it != this->dictionary.end(); ++it)
	{
		out<<(*it).first<<this->delimiter<<(*it).second<<std::endl;
	}

	// Reinstate old format flags
	out.setf(ff);


}


bool ScalarList::read(std::string filename)
{
	std::ifstream in;
	in.open(filename.c_str());

	if(!in.good())
	{
		std::cerr<<"[ScalarList::read] Unable to open file \""<<filename<<"\"!"<<std::endl;
		return false;
	}
		
	return this->read(in);

}


bool ScalarList::read(std::istream& in)
{

	// std::cout<<"[ScalarList::read] Entering function"<<std::endl;
	while(in.good())
	{	
		std::stringstream ss;
		std::string line;
		std::getline(in,line);

		if(line.empty())
		{
			continue;
		}

		std::string id;
		double value;

		ss<<line;
		ss>>id>>value;

		// Store data
		this->addScalar(value,id.c_str());

		// std::cout<<id<<"\t"<<value<<std::endl;

	}

	// std::cout<<"[ScalarList::read] Leaving function"<<std::endl;
	return true;

}


};