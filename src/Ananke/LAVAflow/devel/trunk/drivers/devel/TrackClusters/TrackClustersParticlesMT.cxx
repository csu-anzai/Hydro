#include <string>
#include <vector>
#include "TrackClusters.h"
#include "FlashParticles.h"

using namespace std;

const double massTracerThreshold = 0.2;

void TrackClustersParticlesMT(std::vector<string> &partFilenames, std::vector<string> &meshFilenames, std::string originalClusterFilename, ConfigData &cfgData)
{
    vector<string> partVars = {"tag", "posx", "posy", "posz"};
    std::vector<std::string> clusterMeshVar = {"clst"};

    vector<int> numOrigPartsByTag(2, 0);
    int nClusters = 0;

    std::vector<double> timeSeps;
    std::vector<std::vector<double>> remainingParticleFracPerFile;

    int nOldClusters;
    double origSimTime;
    std::vector<double> origClusterMasses;
    GetOrigMeshInfo(originalClusterFilename, nOldClusters, origSimTime, origClusterMasses);

    int fileIndex = 0;
    for (auto thisFile : partFilenames)
    {

        // cout << "Loading file " << thisFile << endl;
        // load particles. 
        std::vector<const char*> varsChar = VecStringToVecChar(partVars);
        FlashParticles particles(thisFile.c_str(), &varsChar[0], partVars.size(), PART_SORT_PROC);
        int tagInd = particles.findVarIndex("tag");
        int posXInd = particles.findVarIndex("posx");
        int posYInd = particles.findVarIndex("posy");
        int posZInd = particles.findVarIndex("posz");

        // count how many belonged to each original cluster
        double** partData = particles.getDataPtr();
        int nPartsGlobal = particles.nParticlesGlobal;

        if (fileIndex == 0)
        {
            int maxTagEncountered = 1;
            for (int i=0; i<nPartsGlobal; i++)
            {
                int tag = partData[tagInd][i];
                if (tag >= maxTagEncountered)
                {
                    numOrigPartsByTag.resize(tag+1, 0);
                    maxTagEncountered = tag;
                }

                numOrigPartsByTag[tag]++;
            }

            nClusters = numOrigPartsByTag.size() - 2;
        }

        // load the mesh
        FlashAmrMesh mesh(meshFilenames[fileIndex], clusterMeshVar);
        int clusterIDMeshIndex = mesh.findVarIndex("clst");

        double deltaT = mesh.getTime() - origSimTime;
        timeSeps.push_back(deltaT);

        int nCellsBlockX = mesh.getInt("nxb"),
            nCellsBlockY = mesh.getInt("nyb"),
            nCellsBlockZ = mesh.getInt("nzb");

        int nBlocksI = mesh.getInt("nblockx");

        auto blockLinearIndexMap = BuildBlockLinearIndexMap(mesh);

        // get domain bounds and cell side lengths
        double domainBoundBox[2][MESH_MDIM];
        mesh.getDomainBoundBox(domainBoundBox);

        // only works on a uniform mesh!
        double sideLengths[MESH_MDIM];
        mesh.getCellSideLengths(0, sideLengths);

        // this will hold the number of particles still belonging to a cluster, indexed by particle tag ID
        std::vector<int> currentParticleCounts(nClusters+2, 0);

        // for each particle
        for (int partIndex=0; partIndex < nPartsGlobal; partIndex++)
        {
            // find corresponding mesh cell
            int partX = partData[posXInd][partIndex];
            int partY = partData[posYInd][partIndex];
            int partZ = partData[posZInd][partIndex];

            int domainI = int((partX - domainBoundBox[LOWER][IAXIS]) / sideLengths[IAXIS]),
                domainJ = int((partY - domainBoundBox[LOWER][JAXIS]) / sideLengths[JAXIS]),
                domainK = int((partZ - domainBoundBox[LOWER][KAXIS]) / sideLengths[KAXIS]);

            // get the row-major linear block index, assuming an equal number of cells per dimension in the block, and a cubic domain
            int linBlockIndex = GetLinearBlockIndex( domainI,  domainJ,  domainK, nCellsBlockX, nBlocksI);
            int blockID = blockLinearIndexMap[linBlockIndex];

            int blockI = domainI % nCellsBlockX,
                blockJ = domainJ % nCellsBlockY,
                blockK = domainK % nCellsBlockZ;

            // get tracer ID
            double ****solnData;
            //Get pointer to the block's data
            mesh.getBlkPtr(blockID, solnData, CENTER);

            double clstMassTracer = solnData[clusterIDMeshIndex][blockI][blockJ][blockK];

            int particleTag = partData[tagInd][partIndex];

            // if the value of the mass tracer is less than half that of the particle, assume it no longer belongs to the cluster
            if (abs(clstMassTracer - particleTag) < particleTag * massTracerThreshold)
            {
                currentParticleCounts[particleTag]++;
            }
        } // end loop over particles

        std::vector<double> remainingParticleFrac;

        // cout << "Particle fractions remaining: " << endl;
        for (int i=2; i < currentParticleCounts.size(); i++)
        {
            double particleFrac = currentParticleCounts[i] / double(numOrigPartsByTag[i]);
            remainingParticleFrac.push_back(particleFrac);
            // cout << "Cluster " << i << ": " << particleFrac << "\n";
        }
        // cout << endl;

        remainingParticleFracPerFile.push_back(remainingParticleFrac);

        fileIndex++;
        
    } // end loop over files

    ofstream clusterPartFrac("clusterPartFracs.dat", ios::out);

    const int colWidth = 22;

    // column headers
    const char *vinit[] = {"deltaT", "particleID", "origclusterPartFrac"};
    std::vector<std::string> colHeaders(vinit, end(vinit)); // definition

    clusterPartFrac << std::left;
    for (auto thisColHeader : colHeaders)
        clusterPartFrac << std::setw(colWidth) << thisColHeader.c_str();

    clusterPartFrac << '\n';

    for (int fileIndex=0; fileIndex < timeSeps.size(); fileIndex++)
    {
        for (int i=0; i<nClusters; i++)
            clusterPartFrac <<   std::left << 
                                std::setw(colWidth) << timeSeps[fileIndex] <<   
                                std::setw(colWidth) << i+2 <<
                                std::setw(colWidth) << remainingParticleFracPerFile[fileIndex][i] << 
                                endl;
    }

    clusterPartFrac.close();
}
