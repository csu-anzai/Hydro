#include <hdf5.h>
#include <mpi.h>
#include <map>
#include <string>
#include <cstring>
#include <vector>
#include "../Mesh/FlashIO.h"


#ifndef LAVA_PARTICLES_H
#define LAVA_PARTICLES_H

#define MESH_VAR_STRING_SIZE 4
#define MAX_STRING_LENGTH 80
#define FR_MAXVARS 100

#define PART_SORT_PROC 1
#define PART_SORT_TAG  2


class FlashParticles : public FlashIO
{
	public:
	int nParticlesLocal;
	int nParticlesGlobal;

	std::vector<std::string> variableNames; /* Names of stored variables */
	double **particleData;
	int nUserVars;
	int nFileVars;
	int sortMode;
	double simTime;

	std::string   filename;
	
	std::map<int,std::string>     		fileVarsMap_iToS;
	std::map<std::string,int>			fileVarsMap_sToI;
	std::map<int,std::string>     		uVarsMap_iToS;
	std::map<std::string,int>			uVarsMap_sToI;
	std::map<std::string,int>     		intScalarsMap;
	std::map<std::string,double>   		realScalarsMap;

	hid_t h5File;

	MPI::Intracomm comm;

private:
	int beginParticleIndex;
	int endParticleIndex;

	public:

	FlashParticles();
	FlashParticles(const char* , const char**, int, int = PART_SORT_PROC, std::string="", int=-1, MPI::Intracomm = MPI::COMM_WORLD);
	~FlashParticles();


	void Initialize();
	void openFile();
	void closeFile();
	void readVarData(std::string ="", int=-1);

	void readFieldNames();
	std::string findVarName(const int );
	void addUserVariable(const char inputVariable[]);

	double getParticleData(int, char*);
	double getParticleData(int , int );
	double** getDataPtr();

	void readReals(const char*);
	void readInts(const char *);

	virtual double getTime();
	virtual int findVarIndex(const char []);
	
};

#endif
