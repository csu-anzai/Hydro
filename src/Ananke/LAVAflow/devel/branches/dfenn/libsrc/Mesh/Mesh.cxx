#include "Mesh.h"
#include <string>

Mesh::Mesh() {
	
}

Mesh::~Mesh() {
}

void Mesh::setGridFileName(std::string inFile){
	fileName = inFile;
}

std::string Mesh::getGridFileName() {
	return fileName;
}