#include <iostream>
#include <fstream>
#include "../libsrc/Utilities/LAVAUtil.h"


int main(void)
{

	Util::ScalarList outputList, inputList;

	// Add data
	outputList.addScalar(0.123,"a");
	outputList.addScalar(1.0,"b");

	// Write to file
	outputList.print("scalarlist.test");


	// Read from file
	inputList.read("scalarlist.test");

	// Print to screen
	inputList.print(std::cout);




}