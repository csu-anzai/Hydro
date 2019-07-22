std::string basePartFilename = "tburn_hdf5_part_";
vector<string> partVars = {"tag"};
std:vector<std::string> partFilenames(numTimeSeps+1);

std::vector<double> remainingPartFrac(numTimeSeps);


// load particles. 
std::vector<const char*> varsChar = VecStringToVecChar(partVars);
FlashParticles particles(partFilenames[index].c_str(), &varsChar[0], partVars.size(), PART_SORT_PROC);
int tagInd = particles.findVarIndex("tag");

// count how many belonged to each original cluster
double** partData = particles.getDataPtr();
int nPartsGlobal = particles.nParticlesGlobal;

vector<int> numPartsByTag(2, 0);
maxTagEncountered = 1;

for (int i=0; i<nPartsGlobal; i++)
{
    int tag = partData[tagInd][i];
    if (tag >= maxTagEncountered)
        numPartsByTag.resize(tag+1, 0);

    numPartsByTag[tag]++;
}


