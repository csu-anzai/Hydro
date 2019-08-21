
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "FlashAmrMesh.h"
#include "Driver_main.h"
#include <string>

using namespace std;
// #include "H5Cpp.h"

int LAVAFLOW_DRIVER_BASENM(int argc, char* argv[]) {


	if (argc != 2) {
		std::cout << "Need file name for program to run." << std::endl;
		exit(EXIT_FAILURE);
	}

	char arg_inFileName[MAX_FILENAME_LENGTH];
	strcpy(arg_inFileName,argv[1]);
	string inFileName(arg_inFileName);

	// Not using one of the other constructors because variables are
	// required for those.
	FlashAmrMesh mesh;
	mesh.setGridFileName(inFileName);
	mesh.openInFile();
	mesh.generalRead();

	//Print out time, timestep, and stepnumber to file.

	std::cout << std::setprecision(6) << std::scientific;

	double timeVal;
	double timeStep;
	double stepNumber;

	timeVal = mesh.getReal("time");

	if (mesh.getFileType() == 7) {
		stepNumber = mesh.getInt("number of steps");
		timeStep = mesh.getReal("timestep");
	}
	else {
		stepNumber = mesh.getInt("nstep");
		timeStep = mesh.getReal("dt");
	}
	
	std::string timeValStr;
	std::string stepNumberStr;
	std::string timeStepStr;

	std::stringstream convert;
	convert.precision(6);

	convert << std::scientific << timeVal;
	timeValStr = convert.str();
	convert.str("");

	convert << std::scientific << timeStep;
	timeStepStr = convert.str();
	convert.str("");

	convert << std::scientific << stepNumber;
	stepNumberStr = convert.str();
	convert.str("");

	std::string output;
	output.reserve(100);
	output += "JSON={";
		output += "";
			output += "\"time\"";
		output += ":";
			output += timeValStr;
		output += "";
	output += ",";
		output += "";
			output += "\"timeStep\"";
		output += ":";
			output += timeStepStr;
		output += "";
	output += ",";
		output += "";
			output += "\"stepNumber\"";
		output += ":";
			output += stepNumberStr;
		output += "";
	output += "}";

	std::cout << std::setprecision(6) << std::scientific;
	std::cout << output << std::endl;

	return 0;

}