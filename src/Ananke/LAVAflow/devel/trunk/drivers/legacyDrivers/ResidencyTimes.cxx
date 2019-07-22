#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Operators/SurfaceAreaOperator.h"
#include "../libsrc/Operators/shellAveraging/ShellAverageOperator.h"
#include "../libsrc/Operators/shellAveraging/planeShells/ShellAveragePlaneOperator.h"
#include "../libsrc/Operators/surfaceAveraging/planarSurface/SurfaceAveragePlaneOperator.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/includes/LAVAconstants.h"
#include "../libsrc/Utilities/LAVAUtil.h"
#include "../libsrc/Clustering/FloodFill/FloodFill.h"
#include "../libsrc/Particles/Particles.h"
#include "../libsrc/Statistics/Histogram/Histogram.h"

#include <string>
#include <vtkSphereSource.h>
#include <vtkRegularPolygonSource.h>
#include <vtkSphericalTransform.h>
#include <vtkTransformFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkRectilinearGridGeometryFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLRectilinearGridWriter.h>

#include <fstream>
#include <algorithm>
#include <sstream>

using namespace Statistics;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


void determineShockRadius(int nDim, int currentCS, std::string filenameIn, double**& shockRadius, double*& phiEdges, double*& thetaEdges, int& nPhi, int& nTheta)
{


	Mesh data(nDim, currentCS, filenameIn);

	// Get domain and mesh information
	int meshDimensions[3]; // Number of cells in each dimension
	data.getDataDimension(meshDimensions);

	nPhi = meshDimensions[1];
	nTheta = nDim>2 ? meshDimensions[2] : 1;

	if(shockRadius == NULL)
	{
		std::cout<<"Allocating shockRadius"<<std::endl;
		shockRadius = new double*[nPhi];

		for(int i=0; i<nPhi; i++)
		{
			shockRadius[i] = new double[nTheta];
			for(int j=0; j<nTheta; j++)
			{
				shockRadius[i][j] = 0.0;
			} 
		}
	}

	vtkSmartPointer<vtkDataArray> shock = data["shock"];

	for(int p = 0; p<nPhi; ++p)
	{
		for(int t = 0; t<nTheta; ++t)
		{

			// std::cout<<"p = "<<p<<"\tt= "<<t<<std::endl;
			
			int startingIndex[3];
			startingIndex[XAXIS]	= 0;
			startingIndex[YAXIS]	= p;
			startingIndex[ZAXIS]	= t;
			std::vector<int> cellIndices = data.getCellsAlongRay(XAXIS, startingIndex);

			bool foundFirst = false;

			for (std::vector<int>::reverse_iterator i = cellIndices.rbegin(); i != cellIndices.rend(); ++i)
			{
				
				bool isShocked = shock->GetTuple1(*i) > 0.5; // The data is either 0 or 1, so you need to pick something reasonable

				// This section of if statements determines if we've made it inside the shock
				if(isShocked)
				{
					if(!foundFirst)
					{
						foundFirst = true;
					}
				}
				else
				{
					// If we've already found the first element, then not being shocked means we left the shock
					// Therefore, we should start considering that we're in the GR
					if(foundFirst)
					{
						double cellBounds[6];
						data.getDataSet()->GetCellBounds(*i,cellBounds);
						shockRadius[p][t] = cellBounds[0];
						break;
					}

				}
			}
		}
	}


	// Copy rectilinear coordinates into phiEdges and thetaEdges
	vtkSmartPointer<vtkDataArray> phiCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetYCoordinates();
	vtkSmartPointer<vtkDataArray> thetaCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetZCoordinates();

	// Allocate memory for edges
	if(phiEdges == NULL)
	{
		std::cout<<"Allocating phiEdges"<<std::endl;
		phiEdges = new double[nPhi];
	}
	if(thetaEdges == NULL)
	{
		std::cout<<"Allocating thetaEdges"<<std::endl;
		thetaEdges = new double[nTheta+1];
	}

	// Do the copies
	for(int i=0; i<nPhi+1; i++)
	{
		phiEdges[i] = phiCoords->GetTuple1(i);
	}

	if(nDim>2)
	{
		for(int i=0; i<nTheta+1; i++)
		{
			thetaEdges[i] = thetaCoords->GetTuple1(i);
		}
	}
	else
	{
		thetaEdges[0] = -vtkMath::Pi();
		thetaEdges[1] = vtkMath::Pi();
	}



	return;
}



int main(int argc, char** argv)
{

	std::string filenameInBase, fileNameOutBase;
	int currentCS = CS_SPHERE;
	int currAxis   = XAXIS;
	int nDim = 3;
	double MSOLAR = 1.9891e33; // solar mass in grams
	double GRAVCONST = 6.67259e-8; // cm^3/g/s
	double unitFOE = 1.e51;
	double explosionCriterion = 1.e48;
	int fileNumStart = 0;
	int fileNumFinish = 0;
	std::string filenameFmtString = "%s.%s.%04d.vtk";
	std::string filenameFmtStringNoVTK = "%s.%s.%04d";

	double densMax = 1.e10;

	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		// filenameIn = "../Tests/files/averaging/b155d2a3s123SG1MHwr2LB.plt.0010.vtk";
		// filenameIn = "/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/b163d2a3s530SG1MHWr2LBp.plt.0200.vtk";
		filenameInBase = "/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/b163d2a3s530SG1MHWr2LBp";
		fileNameOutBase = "/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/output/b163d2a3s530SG1MHWr2LBp";
		fileNumStart = 0;
		fileNumFinish = 500;
		nDim = 2;

		
		// filenameIn = "/data1/sne/HOTB/3d/b157d3a3s401SG1MHWr3LB/b157d3a3s401SG1MHWr3LB.plt.0100.vtk";
		// nDim = 3;

		// fileNameOutBase = filenameIn;
	}
	else if(argc==5)
	{
		/* Ex:
		filenameInBase = "/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/b163d2a3s530SG1MHWr2LBp";
		fileNameOutBase = "/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp/output/b163d2a3s530SG1MHWr2LBp";
		*/
		filenameInBase = argv[1];
		fileNameOutBase = argv[2];

		// Read start and end numbers
		std::stringstream startNumSS, endNumSS;
		startNumSS<<argv[3];
		endNumSS<<argv[4];

		startNumSS>>fileNumStart;
		endNumSS>>fileNumFinish;

		std::cout<<"Input base:  "<<filenameInBase<<std::endl;
		std::cout<<"Output base: "<<fileNameOutBase<<std::endl;
		std::cout<<"Start number:  "<<fileNumStart<<std::endl;
		std::cout<<"Finish number: "<<fileNumFinish<<std::endl;

		// fileNumStart = 0;
		// fileNumFinish = 500;
	}

	// Histogram variables
	int 	numTimeBins = 250;
	int 	numEintBins = 40;
	double	timeRange[] = {0.0, 0.25};
	double	eintAlpha = 1.0e-3;
	double	eintBeta  = 10.0;

	double  ebindRange[] = {1.e16,1.e20};
	int 	numEbindBins = 40;

	double	radRange[] = {1.e6,1.e10};
	int 	numRadBins = 40;


	double*	avgEintPerTimeAll;
	double*	avgEintPerTimeInst;
	int* 	timeCounterAll;
	int* 	timeCounterInst;

	avgEintPerTimeAll 	= new double[numTimeBins];
	avgEintPerTimeInst 	= new double[numTimeBins];
	timeCounterAll 	= new int[numTimeBins];
	timeCounterInst 	= new int[numTimeBins];
	for(int i=0; i<numTimeBins; i++)
	{
		avgEintPerTimeAll[i] = 0.0;
		avgEintPerTimeInst[i] = 0.0;
		timeCounterAll[i] = 0;
		timeCounterInst[i] = 0;
	}

	// Particle arrays
	double** shockRadius = NULL;
	double*  phiEdges = NULL;
	double*  thetaEdges = NULL;
	double*  residencyTime = NULL;
	double*  eintInitial = NULL;
	double*  eintFinal = NULL;
	bool*    isInGR = NULL;
	bool*    beenInGR = NULL;
	int nPhi = 0;
	int nTheta = 0;

	// Additional data required for PV-work diagram investigation
	double** pvData;
	int nPVDataCats = 10;
	enum{MAXE, INITE, DE, ETIME, ETRES, ERAD, MINRAD, RTIME, RTRES, REINT};


	double simTimePrev = 0.0;
	double simTimeStart = 0.0;
	int nPartInGR = 0;
	int nPartInGRold = 0;
	bool firstPass = true;

	for(int fileNum = fileNumStart; fileNum <= fileNumFinish; fileNum++)
	{
		char filenameMesh[512];
		sprintf(filenameMesh,filenameFmtString.c_str(),	filenameInBase.c_str(), "plt", fileNum);

		char filenamePart[512];
		sprintf(filenamePart,filenameFmtString.c_str(), filenameInBase.c_str(), "prt", fileNum);

		char filenameOut[512];
		sprintf(filenameOut,(filenameFmtStringNoVTK+".residency").c_str(), fileNameOutBase.c_str(), "prt", fileNum);

		char filenameOutPV[512];
		sprintf(filenameOutPV,(filenameFmtStringNoVTK+".pvdata").c_str(), fileNameOutBase.c_str(), "prt", fileNum);


		std::cout<<"\n\n"<<std::endl;
		std::cout<<"Processing file \""<<filenameMesh<<"\""<<std::endl;

		// Determine the shock radius as a function of (phi,theta)
		determineShockRadius(nDim, currentCS, filenameMesh, shockRadius, phiEdges, thetaEdges, nPhi, nTheta);

		// Determine the angular spacing. The following work assumes a uniform spacing in angle
		double dPhi = phiEdges[1]-phiEdges[0];
		double dTheta = thetaEdges[1]-thetaEdges[0];

		// Load in particles
		Particles parts(currentCS, filenamePart);
		int nParts = parts.getNumberOfParticles();
		double simTime = parts.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);
		double simDt;// = parts.getDataSet()->GetFieldData()->GetArray("DT")->GetTuple1(0);

		// parts.getDataSet()->PrintSelf(std::cout,vtkIndent(0));

		if(firstPass)
		{
			simTimePrev = simTime;
			simTimeStart = simTime;
		}

		simDt = simTime-simTimePrev;


		std::cout<<"Sim time: "<<simTime<<std::endl;
		std::cout<<"Sim dt:   "<<simDt<<std::endl;


		// Allocate arrays if necessary
		if(firstPass)
		{
			isInGR = new bool[nParts];
			beenInGR = new bool[nParts];
			residencyTime = new double[nParts];
			eintInitial = new double[nParts];
			eintFinal   = new double[nParts];
			pvData = new double*[nParts];
			for(int i=0; i<nParts; i++)
			{
				isInGR[i] = false;
				beenInGR[i] = false;
				residencyTime[i] = 0.0;
				eintInitial[i] = 0.0;
				eintFinal[i] = 0.0;

				pvData[i] = new double[nPVDataCats];
				for(int j=0; j<nPVDataCats; j++)
				{
					pvData[i][j] = 0.0;
				}
			}
		}

		// Create histogram container
		Histogram hist;
		Histogram::HistObj* tresHistAll = hist.createStaticHistogram(timeRange, numTimeBins, "residencyTimeAll", Histogram::SPACING_LIN);
		Histogram::HistObj* tresHistInstant = hist.createStaticHistogram(timeRange, numTimeBins, "residencyTimeInstant", Histogram::SPACING_LIN);

		// Reinitialize some arrays

		for(int i=0; i<numTimeBins; i++)
		{
			avgEintPerTimeAll[i] = 0.0;
			avgEintPerTimeInst[i] = 0.0;
			timeCounterAll[i] = 0;
			timeCounterInst[i] = 0;
		}

		std::cout<<"Number of particles: "<<nParts<<std::endl;

		// Determine which particles are in the GR
		// The gain region for this method is defined as being below the shock (r<r_Shock), and having a density less than densMax (dens<densMax)
		// densMax should be ~1e10
		vtkSmartPointer<vtkDataArray> dens = parts["dem"];
		vtkSmartPointer<vtkDataArray> eint = parts["eim"];
		vtkSmartPointer<vtkDataArray> ebind = parts["eim"] + 0.5*(parts["pvx"]*parts["pvx"] + parts["pvy"]*parts["pvy"] + parts["pvz"]*parts["pvz"]) + parts["gpm"];
		int cntLostParticles = 0;
		nPartInGR = 0;
		for(int i=0; i<nParts; i++)
		{
			// Get particle location and density
			double ptPos[3];
			double ptDens;

			parts.getParticlePosition(i,ptPos);
			ptDens = dens->GetTuple1(i);

			// Determine the index into the phi-coordinate that the particle lies in
			int phiInd = 0;
			phiInd = int((ptPos[1]-phiEdges[0])/dPhi);

			// Determine the index into the theta-coordinate that the particle lies in
			int thetaInd = 0;
			if(nDim>2)
			{
				thetaInd = int((ptPos[2]-thetaEdges[0])/dTheta);
			}
			else
			{
				thetaInd = 0;
			}

			// NOTE: Sometimes, the integration of particles can push them outside of the physical domain
			// These guys are lost, so we'll just ignore them.
			if(phiInd<0 || phiInd>=nPhi || thetaInd<0 || thetaInd>=nTheta)
			{
				isInGR[i] = false;
				cntLostParticles++;
				continue;
			}


			// Get the shock position for this (phi,theta) pair
			double rShock = shockRadius[phiInd][thetaInd];

			// Determine if the particle is in the gain region
			bool isBelowShock = ptPos[0] < rShock;
			bool isBelowDensMax = ptDens < densMax;
			bool isInGainRegion = isBelowShock && isBelowDensMax;

			// If the particle used to be in the gain region (isInGR==true), and is currently in the gain region (isInGainRegion==true)
			// add the timestep to the particle's residency time
			if(isInGR[i] && isInGainRegion)
			{
				// Determine if this is the first step that the particle is in the gain region
				// If so, record it's initial internal energy
				if(!beenInGR[i])
				{
					beenInGR[i] = true;
					eintInitial[i] = eint->GetTuple1(i);

					// set initial pvData values
					pvData[i][MAXE] = eintInitial[i];
					pvData[i][INITE] = eintInitial[i];
					pvData[i][DE] = 0.0;
					pvData[i][ETIME] = simTime;
					pvData[i][ETRES] = simDt;
					pvData[i][ERAD]  = ptPos[0];
					pvData[i][MINRAD] = ptPos[0];
					pvData[i][RTIME] = simTime;
					pvData[i][RTRES] = simDt;
					pvData[i][REINT] = eintInitial[i];

				}

				// Increment residency time
				residencyTime[i] += simDt;

				// Compute most recent eint
				eintFinal[i] = eint->GetTuple1(i);

				// Update the instantaneous residency time histogram
				int tresInd = hist.addDatum(tresHistInstant, residencyTime[i]);

				// update average array
				avgEintPerTimeInst[tresInd] += eintFinal[i] - eintInitial[i];

				// update counter 
				timeCounterInst[tresInd]++;

				// Update pvData
				if(eintFinal[i]>pvData[i][MAXE])
				{
					pvData[i][MAXE] = eintFinal[i];
					pvData[i][DE] = pvData[i][MAXE]-pvData[i][INITE];
					pvData[i][ETIME] = simTime;
					pvData[i][ETRES] = residencyTime[i];
					pvData[i][ERAD]  = ptPos[0];
				}

				if(ptPos[0]<pvData[i][MINRAD])
				{
					pvData[i][MINRAD] = ptPos[0];
					pvData[i][RTIME] = simTime;
					pvData[i][RTRES] = residencyTime[i];
					pvData[i][REINT] = eintFinal[i];
				}
			}

			// Update the total residency time histogram
			if(beenInGR[i])
			{
				int tresInd = hist.addDatum(tresHistAll, residencyTime[i]);

				// update average array
				avgEintPerTimeAll[tresInd] += eintFinal[i] - eintInitial[i];

				// update counter 
				timeCounterAll[tresInd]++;
			}

			// Set particle's past GR state (from the next time slice's perspective)
			isInGR[i] = isInGainRegion;

			// Increase counter
			if(isInGainRegion) nPartInGR++;

		}

		// print residency time histogram to screen
		// hist.printObject(std::cout);

		// Finish computing the average eint at each residency time and determine the maximum value
		double avgEintPerTimeAllMax = 0.0;
		double avgEintPerTimeInstMax = 0.0;
		for (int i=0; i<numTimeBins; i++)
		{
			if(timeCounterAll[i]!=0)
			{
				avgEintPerTimeAll[i] /= timeCounterAll[i];
				// std::cout<<"avgEintPerTimeAll["<<i<<"]: "<<avgEintPerTimeAll[i]<<std::endl;
			}
			if(timeCounterInst[i]!=0)
			{
				avgEintPerTimeInst[i] /= timeCounterInst[i];
				// std::cout<<"avgEintPerTimeInst["<<i<<"]: "<<avgEintPerTimeInst[i]<<std::endl;
			}
			if(avgEintPerTimeAll[i]>avgEintPerTimeAllMax)
			{
				avgEintPerTimeAllMax = avgEintPerTimeAll[i];
			}
			if(avgEintPerTimeInst[i]>avgEintPerTimeInstMax)
			{
				avgEintPerTimeInstMax = avgEintPerTimeInst[i];
			}
		}


		if(avgEintPerTimeAllMax<= 0.0)
		{
			avgEintPerTimeAllMax = 1.e-6;
		}
		if(avgEintPerTimeInstMax<= 0.0)
		{
			avgEintPerTimeInstMax = 1.e-6;
		}


		std::cout<<"avgEintPerTimeAllMax: "<<avgEintPerTimeAllMax<<std::endl;
		std::cout<<"avgEintPerTimeInstMax: "<<avgEintPerTimeInstMax<<std::endl;

		// Create deltaEint bins
		double eintRangeAll[]  = {eintAlpha*avgEintPerTimeAllMax, eintBeta*avgEintPerTimeAllMax};
		double eintRangeInst[] = {eintAlpha*avgEintPerTimeInstMax, eintBeta*avgEintPerTimeInstMax};
		std::vector<Histogram::HistObj*> eintHistsAll;
		std::vector<Histogram::HistObj*> eintHistsInst;
		for(int i=0; i<numTimeBins; i++)
		{
			char histName[256];
			sprintf(histName,"deltaEint_%04d_All",i);
			eintHistsAll.push_back(hist.createStaticHistogram(eintRangeAll, numEintBins, histName, Histogram::SPACING_LOG));
		}		
		for(int i=0; i<numTimeBins; i++)
		{
			char histName[256];
			sprintf(histName,"deltaEint_%04d_Inst",i);
			eintHistsInst.push_back(hist.createStaticHistogram(eintRangeInst, numEintBins, histName, Histogram::SPACING_LOG));
		}

		// Create radius histograms
		std::vector<Histogram::HistObj*> radHistsAll;
		std::vector<Histogram::HistObj*> radHistsInst;
		for(int i=0; i<numTimeBins; i++)
		{
			char histName[256];
			sprintf(histName,"radialPos_%04d_All",i);
			radHistsAll.push_back(hist.createStaticHistogram(radRange, numRadBins, histName, Histogram::SPACING_LOG));
		}		
		for(int i=0; i<numTimeBins; i++)
		{
			char histName[256];
			sprintf(histName,"radialPos_%04d_Inst",i);
			radHistsInst.push_back(hist.createStaticHistogram(radRange, numRadBins, histName, Histogram::SPACING_LOG));
		}

		// Create ebind histograms
		std::vector<Histogram::HistObj*> ebindHistsAll;
		std::vector<Histogram::HistObj*> ebindHistsInst;
		for(int i=0; i<numTimeBins; i++)
		{
			char histName[256];
			sprintf(histName,"ebind_%04d_All",i);
			ebindHistsAll.push_back(hist.createStaticHistogram(ebindRange, numEbindBins, histName, Histogram::SPACING_LOG));
		}		
		for(int i=0; i<numTimeBins; i++)
		{
			char histName[256];
			sprintf(histName,"ebind_%04d_Inst",i);
			ebindHistsInst.push_back(hist.createStaticHistogram(ebindRange, numEbindBins, histName, Histogram::SPACING_LOG));
		}

		// Compute histograms
		for(int i=0; i<nParts; i++)
		{
			int tresInd = hist.getBinIndex(tresHistAll, residencyTime[i]);
			double eintDeltaTmp = eintFinal[i] - eintInitial[i];
			double eintDeltaAllTmp = eintDeltaTmp;
			double eintDeltaInstTmp = eintDeltaTmp;
			double ebindTmp = ebind->GetTuple1(i);

			double partPos[3];
			parts.getParticlePosition(i,partPos);
			double partRad = partPos[0];

			// If the eint value for this particle is <alpha*eintAvgMax, set the value to alpha*eintAvgMax
			if(eintDeltaTmp<eintAlpha*avgEintPerTimeAllMax)
			{
				eintDeltaAllTmp = eintAlpha*avgEintPerTimeAllMax;
			}
			if(eintDeltaTmp<eintAlpha*avgEintPerTimeInstMax)
			{
				eintDeltaInstTmp = eintAlpha*avgEintPerTimeInstMax;
			}

			// Total histograms
			if(beenInGR[i])
			{
				hist.addDatum(eintHistsAll[tresInd],  eintDeltaAllTmp);
				hist.addDatum(radHistsAll[tresInd],   partRad);
				hist.addDatum(ebindHistsAll[tresInd], ebindTmp);
			}

			// Instantaneous histograms
			if(isInGR[i])
			{
				hist.addDatum(eintHistsInst[tresInd],  eintDeltaAllTmp);
				hist.addDatum(radHistsInst[tresInd],   partRad);
				hist.addDatum(ebindHistsInst[tresInd], ebindTmp);
			}
		}

		// Write to file
		// Open file for output
		std::cout<<"Writing output file: "<<filenameOut<<std::endl;
		std::ofstream output;
		output.open(filenameOut);
		if(output.is_open())
		{
			output<<"simTimeStart"<<"\t"<<simTimeStart<<std::endl;
			output<<"simTime"<<"\t"<<simTime<<std::endl;
			output<<"simDt"<<"\t"<<simDt<<std::endl;
			output<<"averageEintAllMax"<<"\t"<<avgEintPerTimeAllMax<<std::endl;
			output<<"averageEintInstMax"<<"\t"<<avgEintPerTimeInstMax<<std::endl;
			hist.printObject(output);
			output.close();
		}
		else
		{
			std::cerr<<"[ResidencyTimes.cxx] ERROR: Unable to open output file \""<<filenameOut<<"\". Exiting now."<<std::endl;
			exit(-1);
		}

		// Write pvData to file
		std::cout<<"Writing output file: "<<filenameOutPV<<std::endl;
		output.open(filenameOutPV);
		if(output.is_open())
		{
			output<<"simTimeStart"<<"\t"<<simTimeStart<<std::endl;
			output<<"simTime"<<"\t"<<simTime<<std::endl;
			output<<"simDt"<<"\t"<<simDt<<std::endl;
			output	<<"partNum"<<"\t"
					<<"inGRnow"<<"\t"
					<<"maxEint"<<"\t"
					<<"initEint"<<"\t"
					<<"deltaEint"<<"\t"
					<<"timeEint"<<"\t"
					<<"tresEint"<<"\t"
					<<"radiusEint"<<"\t"
					<<"minRadius"<<"\t"
					<<"timeRadius"<<"\t"
					<<"tresRadius"<<"\t"
					<<"eintRadius";
			output<<std::endl;
			for(int i=0; i<nParts; i++)
			{
				if(beenInGR[i])
				{
					output	<<i<<"\t"
							<<isInGR[i]<<"\t"
							<<pvData[i][MAXE]<<"\t"
							<<pvData[i][INITE]<<"\t"
							<<pvData[i][DE]<<"\t"
							<<pvData[i][ETIME]<<"\t"
							<<pvData[i][ETRES]<<"\t"
							<<pvData[i][ERAD]<<"\t"
							<<pvData[i][MINRAD]<<"\t"
							<<pvData[i][RTIME]<<"\t"
							<<pvData[i][RTRES]<<"\t"
							<<pvData[i][REINT];
					output<<std::endl;
				}	
			}

			output.close();
		}
		else
		{
			std::cerr<<"[ResidencyTimes.cxx] ERROR: Unable to open output file \""<<filenameOutPV<<"\". Exiting now."<<std::endl;
			exit(-1);
		}

#if 0
		// Create histogram
		if(!firstPass)
		{
			Statistics::Histogram hist;
			// hist.processFixedWidth(residencyTime, nParts, 1e-3, "test");
			double binRange[2]; binRange[0] = 0.0; binRange[1] = 0.25;
			hist.processStaticHistogram(residencyTime, nParts, binRange, 250, "test");
			hist.printObject(std::cout);

			// Open file for output
			std::cout<<"Writing output file: "<<filenameOut<<std::endl;
			std::ofstream output;
			output.open(filenameOut);
			if(output.is_open())
			{
				output<<"simTimeStart"<<"\t"<<simTimeStart<<std::endl;
				output<<"simTime"<<"\t"<<simTime<<std::endl;
				output<<"simDt"<<"\t"<<simDt<<std::endl;
				hist.printObject(output);
				output.close();
			}
			else
			{
				std::cerr<<"[ResidencyTimes.cxx] ERROR: Unable to open output file \""<<filenameOut<<"\". Exiting now."<<std::endl;
				exit(-1);
			}
		}
#endif

		// Print some stats
		std::cout<<"Current number of particles in GR:  "<<nPartInGR<<std::endl;
		std::cout<<"Previous number of particles in GR: "<<nPartInGRold<<std::endl;
		std::cout<<"Percent change: "<<100.0*(nPartInGR-nPartInGRold)/nPartInGR<<std::endl;
		std::cout<<"Lost particles: "<<cntLostParticles<<std::endl;

		nPartInGRold = nPartInGR;
		simTimePrev = simTime;

		if(firstPass) firstPass = false;
	}



	/* -------------------------------------
					Cleanup
	   ------------------------------------- */
	   	
	delete[] shockRadius;
	delete[] phiEdges;
	delete[] thetaEdges;
	delete[] residencyTime;
	delete[] isInGR;
	delete[] beenInGR;
	delete[] eintInitial;
	delete[] eintFinal;
	delete[] avgEintPerTimeAll;
	delete[] avgEintPerTimeInst;
	delete[] timeCounterAll;
	delete[] timeCounterInst;
	delete[] pvData;


}
