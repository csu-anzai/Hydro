#ifndef LAVA_UTIL_SCALARLIST_H
#define LAVA_UTIL_SCALARLIST_H

#include <ostream>
#include <map>
#include <string>


namespace Util
{

class ScalarList
{
public:
	ScalarList()
	{
		this->numberOfEntries = 0;
		this->delimiter = "\t";
	}

	/* data */

	int numberOfEntries;
	std::map<std::string,double> dictionary;

	std::string delimiter;

	int getSize() { return this->numberOfEntries; }
	double getScalar(const char* id);

	bool addScalar(double value, const char* id);
	bool print(std::ostream& out);
	bool print(std::string filename);
	bool read(std::istream& in);
	bool read(std::string filename);

};

};

#endif