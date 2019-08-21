#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <list>

#include "vtkDataArray.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkMath.h"

#include "../libsrc/Particles/Particles.h"
#include "../libsrc/includes/LAVAconstants.h"
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/Utilities/LAVAUtil.h"


void createDissipationGridVar(Mesh &data, int nDim);

int main(int argc, char** argv)
{

	std::cout<<"Hello world!"<<std::endl;
	

	/* ==================================================

		Read particles
	
	================================================== */
	int nDim = 2;

	int fileNum = 300;
	std::string dirInFmt = "/data1/sne/HOTB/2d/%s";
	std::string dirOutFmt = "/data1/sne/HOTB/2d/%s/output";
	std::string modelName = "b163d2a3s530SG1MHWr2LBp";
	// std::string directory = "/data1/sne/HOTB/2d/b163d2a3s132SG1MHWr2LBp";
	// std::string modelName = "b163d2a3s132SG1MHWr2LBp";
	char dirIn[256], dirOut[256];
	std::string pltFileName, prtFileName, expFileName, fileNameOutBase;
	
	if(argc==1)
	{
		std::cout<<"No arguments passed!"<<std::endl;
		
		char pltFileNameChar[256], prtFileNameChar[256], expFileNameChar[256];
		// Create file names
		std::sprintf(dirIn,dirInFmt.c_str(),modelName.c_str());
		std::sprintf(dirOut,dirOutFmt.c_str(),modelName.c_str());
		std::sprintf(prtFileNameChar,"%s/%s.prt.%04d.vtk",		dirIn,	modelName.c_str(),	fileNum);
		std::sprintf(expFileNameChar,"%s/%s.plt.%04d.expstat",	dirOut,	modelName.c_str(),	fileNum);
		std::sprintf(pltFileNameChar,"%s/%s.plt.%04d.vtk",		dirIn,	modelName.c_str(),	fileNum);

		pltFileName = pltFileNameChar;
		prtFileName = prtFileNameChar;
		expFileName = expFileNameChar;
		fileNameOutBase = "evsf";
	}
	else if(argc==5)
	{

		pltFileName = argv[1];
		prtFileName = argv[2];
		fileNameOutBase = argv[3];
		expFileName = argv[4];
	}
	else
	{
		std::cerr<<"Wrong number of arguments!"<<std::endl;
	}


	// Read in scalar list from explosion statistics file
	// and determine the gain radius and the minimum shock radius
	Util::ScalarList expStat;
	expStat.read(expFileName);
	double gainRadius = expStat.getScalar("grMax");
	double shockRadiusMin = expStat.getScalar("shockMin");

	// Read the particle file
	Particles parts(CS_SPHERE,prtFileName);
	parts.createCoordinateDataArrays("r","phi","theta");
	vtkSmartPointer<vtkDataArray> radii = parts["r"];
	int nParticlesTotal = parts.getNumberOfParticles();
	double simTime = parts.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);



	// Create cartesian coordinates for each particle
	parts.addVariable("x",0.0);
	parts.addVariable("y",0.0);
	parts.addVariable("z",0.0);
	parts.addVariable("pvxCart",0.0);
	parts.addVariable("pvyCart",0.0);
	parts.addVariable("pvzCart",0.0);

	for(int p=0; p<nParticlesTotal; p++)
	{
		double r = parts["r"]->GetTuple1(p);
		double phi = parts["phi"]->GetTuple1(p);
		double theta = parts["theta"]->GetTuple1(p);
		double vr = parts["pvx"]->GetTuple1(p);
		double vp = parts["pvy"]->GetTuple1(p);
		double vt = parts["pvz"]->GetTuple1(p);

		double x = r*cos(theta)*sin(phi);
		double y = r*sin(theta)*sin(phi);
		double z = r*cos(phi);

		double vx = vr*sin(phi)*cos(theta) + vp*cos(phi)*cos(theta) - vt*sin(theta);
		double vy = vr*sin(phi)*sin(theta) + vp*cos(phi)*sin(theta) + vt*cos(theta);
		double vz = vr*cos(theta) + vp*sin(phi);

		parts["x"]->SetTuple1(p,x);
		parts["y"]->SetTuple1(p,y);
		parts["z"]->SetTuple1(p,z);
		parts["pvxCart"]->SetTuple1(p,vx);
		parts["pvyCart"]->SetTuple1(p,vy);
		parts["pvzCart"]->SetTuple1(p,vz);
	}


	/* ==================================================

		Read grid data
	
	================================================== */


	std::cout<<"Reading grid data from file \""<<pltFileName<<"\"...";
	std::cout.flush();

	Mesh gridData(nDim, CS_SPHERE, pltFileName);
	gridData.createCoordinateDataArrays("r","phi","theta");

	std::cout<<"done!"<<std::endl;

	std::cout<<"Computing strain double dot product...";
	std::cout.flush();

	createDissipationGridVar(gridData, nDim);
	
	std::cout<<"done!"<<std::endl;
	

	/* ==================================================

		Determine viable particles
		(1)	This includes particles that are inside the GR
	
	================================================== */
	std::list<int> particleIndices;
	parts.addVariable("mask",0.0);


	for(int p=0; p<nParticlesTotal; p++)
	{
		double r = parts["r"]->GetTuple1(p);

		if(r>=gainRadius && r<=shockRadiusMin)
		{
			// Add to list of viable indices
			particleIndices.push_back(p);

			// Set mask to 1.0
			parts["mask"]->SetTuple1(p,1.0);
		}
	}

	int nParticlesViable = particleIndices.size();
	std::cout<<"Number of viable particles: "<<nParticlesViable<<std::endl;

	/* ==================================================

		Remove the mean component from the velocities
	
	================================================== */
	double velxMean 	= parts.mean("pvx","mask");
	double velyMean 	= parts.mean("pvy","mask");
	double velzMean 	= parts.mean("pvz","mask");		
	double velxCartMean = parts.mean("pvxCart","mask");
	double velyCartMean = parts.mean("pvyCart","mask");
	double velzCartMean = parts.mean("pvzCart","mask");


	for(std::list<int>::iterator it = particleIndices.begin(); it != particleIndices.end(); ++it)
	{
		parts["pvx"]->SetTuple1(*it,parts["pvx"]->GetTuple1(*it)-velxMean);
		parts["pvy"]->SetTuple1(*it,parts["pvy"]->GetTuple1(*it)-velyMean);
		parts["pvz"]->SetTuple1(*it,parts["pvz"]->GetTuple1(*it)-velzMean);
		
		parts["pvxCart"]->SetTuple1(*it,parts["pvxCart"]->GetTuple1(*it)-velxCartMean);
		parts["pvyCart"]->SetTuple1(*it,parts["pvyCart"]->GetTuple1(*it)-velyCartMean);
		parts["pvzCart"]->SetTuple1(*it,parts["pvzCart"]->GetTuple1(*it)-velzCartMean);
	}

	/* ==================================================

		Prepare the bins for the cumulative averaging
		of pair-wise data as it streams in. 

		We take the bin range minimum to be 0 (no separation)
		and the bin range maximum to be twice the shock radius

		We estimate the number of (linearly spaced) bins
		as nBins = sqrt(nParticlesViable*nParticlesViable) = nParticlesViable
		via Square-root choice. This may want to be changed later
	
		This gives a bin width, h, via h=(max-min)/nParticlesViable

		We'll allow the actual bin center to vary, by also
		taking the cumulative average of the particle separation
		for those found to lie in the bin

		Therefore, for each bin we have to store:
			(1)	Bin edges (2 doubles, static)
			(2)	Bin center (1 double, dynamic)
			(3)	Bin values (N doubles, dynamic)

		We also need to record the number of times each bin has been accessed
			(1) Number of Updates (nBin ints, dynamic)

	
	================================================== */
	int nBins = nParticlesViable;
	int nQuantities = 6; // autocorrelation r/phi/theta/dotprod, p=2, p=3
	int nDoubles = 2+1+nQuantities; // edges + center + data

	// Allocate storage
	int* binAccesses = new int[nBins];
	for(int i=0; i<nBins; i++)
	{
		binAccesses[i] = 0;
	}

	double** bins;
	bins = new double*[nBins];
	for(int i=0; i<nBins; i++)
	{
		bins[i] = new double[nDoubles];
		for(int j=0; j<nQuantities; j++)
		{
			bins[i][j] = 0.0;
		}
	}

	// Set the bin edges
	double binEdgeMin = 0.0;
	double binEdgeMax = 2.0*shockRadiusMin;

	double binWidth = (binEdgeMax-binEdgeMin)/double(nBins);

	bins[0][0] = binEdgeMin;
	bins[nBins-1][1] = binEdgeMax;
	for(int i=0; i<nBins-1; i++)
	{
		bins[i][1] = binEdgeMin + binWidth*double(i+1);
		bins[i+1][0] = bins[i][1];
	}


	/* ==================================================

		Compute the p=2 Eulerian structure functions
		(1)	For each particle I
			(a)	For each particle J
				(i)		Determine if the particle J to the "right"
						of particle I. If so, compute the velocity 
						difference and separation
				(ii)	Compute the local constribution to the structure
						function S(r)
				(ii) 	Store this result as a pair<double,double>

	================================================== */

	double strainAvg = 0.0;

	std::vector<std::pair<double,double> > SxData;
	std::vector<std::pair<double,double> > SyData;
	std::vector<std::pair<double,double> > SzData;
	std::vector<std::pair<double,double> > SrData;
	std::list<int>::iterator partIt, otherIt;

	vtkSmartPointer<vtkDataArray> x = parts["x"];
	vtkSmartPointer<vtkDataArray> y = parts["y"];
	vtkSmartPointer<vtkDataArray> z = parts["z"];
	vtkSmartPointer<vtkDataArray> r = parts["r"];
	vtkSmartPointer<vtkDataArray> phi = parts["phi"];
	vtkSmartPointer<vtkDataArray> theta = parts["theta"];
	vtkSmartPointer<vtkDataArray> velx = parts["pvx"];
	vtkSmartPointer<vtkDataArray> vely = parts["pvy"];
	vtkSmartPointer<vtkDataArray> velz = parts["pvz"];
	vtkSmartPointer<vtkDataArray> velxCart = parts["pvxCart"];
	vtkSmartPointer<vtkDataArray> velyCart = parts["pvyCart"];
	vtkSmartPointer<vtkDataArray> velzCart = parts["pvzCart"];


	int counter = 0;

	for(partIt = particleIndices.begin(); partIt != particleIndices.end(); ++partIt)
	{

		// std::cout<<int(double(nParticlesViable*nParticlesViable)/100.0)<<std::endl;
		if(counter%1000==0)
		{
			std::cout<<"Percent complete: "<<double(counter)/double(nParticlesViable)*100<<std::endl;
		}


		int part = *partIt;
		double partX = x->GetTuple1(part);
		double partY = y->GetTuple1(part);
		double partZ = z->GetTuple1(part);

		double partVX = velxCart->GetTuple1(part);
		double partVY = velyCart->GetTuple1(part);
		double partVZ = velzCart->GetTuple1(part);

		// Sample the strain at the current particle
		double partPt[3];
		partPt[0] = r->GetTuple1(part);
		partPt[1] = phi->GetTuple1(part);
		partPt[2] = theta->GetTuple1(part);
		strainAvg += gridData.interpolate(partPt,"energyDissipation",0);

		for(otherIt = particleIndices.begin(); otherIt != particleIndices.end(); ++otherIt)
		{
			int other = *otherIt;

			// std::cout<<part<<"\t"<<other<<std::endl;
			
			// Don't look at self
			if(part==other)
				continue;


			double otherX = x->GetTuple1(other);
			double otherY = y->GetTuple1(other);
			double otherZ = z->GetTuple1(other);
			double otherVX = velxCart->GetTuple1(other);
			double otherVY = velyCart->GetTuple1(other);
			double otherVZ = velzCart->GetTuple1(other);

			// Compute the vector between the two particles
			double coordDiff[3];
			coordDiff[0] = otherX - partX;
			coordDiff[1] = otherY - partY;
			coordDiff[2] = otherZ - partZ;

			// Compute the separation length
			double separation = coordDiff[0]*coordDiff[0] +
								coordDiff[1]*coordDiff[1] +
								coordDiff[2]*coordDiff[2];

			separation = sqrt(separation);

			// Normalize the vector between the particles
			coordDiff[0] /= separation;
			coordDiff[1] /= separation;
			coordDiff[2] /= separation;

			// Compute the velocity difference
			double velDiff[3];
			velDiff[0] = otherVX-partVX;
			velDiff[1] = otherVY-partVY;
			velDiff[2] = otherVZ-partVZ;

			// Project the velocity onto the separation vector
			double diffProj = 	coordDiff[0]*velDiff[0] + 
								coordDiff[1]*velDiff[1] + 
								coordDiff[2]*velDiff[2];

			// Store the results
			// SrData.push_back(std::pair<double,double>(separation,diffProj));

			// Determine the bin number
			int binNumber = int(floor(separation/binWidth));

			if(binNumber>=nBins || binNumber<0)
			{
				std::cerr<<"Bin number is outside of range! binNumber = "<<binNumber<<std::endl;
			}

			// Increase the bin accessor count
			binAccesses[binNumber]++;
			int N = binAccesses[binNumber];

			// Update cumulative averages
			// Bin center
			bins[binNumber][2] = (separation+(N-1)*bins[binNumber][2])/double(N);
			// Autocorrelation
			double aX = velx->GetTuple1(part)*velx->GetTuple1(other);
			double aY = vely->GetTuple1(part)*vely->GetTuple1(other);
			double aZ = velz->GetTuple1(part)*velz->GetTuple1(other);
			bins[binNumber][3] = (aX+(N-1)*bins[binNumber][3])/double(N);
			bins[binNumber][4] = (aY+(N-1)*bins[binNumber][4])/double(N);
			bins[binNumber][5] = (aZ+(N-1)*bins[binNumber][5])/double(N);
			bins[binNumber][6] = ((aX+aY+aZ)+(N-1)*bins[binNumber][6])/double(N);
			// p=2
			bins[binNumber][7] = (diffProj*diffProj+(N-1)*bins[binNumber][7])/double(N);

			// p=3
			bins[binNumber][8] = (diffProj*diffProj*diffProj+(N-1)*bins[binNumber][8])/double(N);



		}

		counter++;

	}
	// Finish computing the average strain
	strainAvg /= double(nParticlesViable);

	std::cout<<"Average strain: "<<strainAvg<<std::endl;



	/* ==================================================

		Determine the kinematic viscosity as a function
		of separation distance.

		nu(r) = dissipation(r)/(2*average strain)

		The dissipation quantity is obtained from 
		Kolmogorov's structure function relations.

		For the second order longitudinal structure function (p=2)
		the dissipation as a function of separation (r) is given by

		edis(r) = <u_r(r)^2>^(3/2)/(C_2^(3/2)*r)

		We assume C_2 is ~O(1)
	
	================================================== */
	double* kinVisc = new double[nBins];

	for(int i=0; i<nBins; i++)
	{
		kinVisc[i] = (pow(bins[i][7], 1.5)/bins[i][2])/(2.0*strainAvg);
	}


	/* ==================================================

		Setup output stream
	
	================================================== */
	// std::string outputFileNameFmt = "eulerian_statistics.%04d.x.dat";
	// char outputFileName[256];
	// std::sprintf(outputFileName,outputFileNameFmt.c_str(),fileNum);

	// std::ofstream output;
	// output.open(outputFileName);
	// output<<"separation"<<"\t"<<"difference"<<std::endl;

	// output.precision(15);
	// output.setf(std::ios::scientific,std::ios::floatfield);

	// std::vector<std::pair<double,double> >::iterator it;
	// for(it = SrData.begin(); it != SrData.end(); it++)
	// {
	// 	output<<(*it).first<<"\t"<<(*it).second<<std::endl;
	// }
	// output.close();

	std::string outputFileNameFmt = "%s.evsf";
	char outputFileName[256];
	std::sprintf(outputFileName,outputFileNameFmt.c_str(),fileNameOutBase.c_str());

	std::ofstream output;
	output.open(outputFileName);
	output<<"simTime"<<"\t"<<simTime<<std::endl;
	output<<"strainAvg"<<"\t"<<strainAvg<<std::endl;
	output<<"nParticlesViable"<<"\t"<<nParticlesViable<<std::endl;
	output 	<<"binLeft"<<"\t"
			<<"binRight"<<"\t"
			<<"binCenter"<<"\t"
			<<"binAccesses"<<"\t"
			<<"autocorR"<<"\t"
			<<"autocorPhi"<<"\t"
			<<"autocorTheta"<<"\t"
			<<"autocorDotProd"<<"\t"
			<<"p=2"<<"\t"
			<<"p=3"<<"\t"
			<<"kinVisc"<<std::endl;


	output.precision(15);
	output.setf(std::ios::scientific,std::ios::floatfield);

	for(int i=0; i<nBins; i++)
	{
		output 	<<bins[i][0]<<"\t"
				<<bins[i][1]<<"\t"
				<<bins[i][2]<<"\t"
				<<binAccesses[i]<<"\t"
				<<bins[i][3]<<"\t"
				<<bins[i][4]<<"\t"
				<<bins[i][5]<<"\t"
				<<bins[i][6]<<"\t"
				<<bins[i][7]<<"\t"
				<<bins[i][8]<<"\t"
				<<kinVisc[i]<<std::endl;
	}
	output.close();

	/* ==================================================

		Cleanup
	
	================================================== */

	delete[] bins;
	delete[] binAccesses;


	// vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	// writer->SetFileName("particles.vtk");
	// writer->SetInputData(partsInitial.getDataSet());
	// writer->Write();

	std::cout<<"Finished"<<std::endl;
}


void createDissipationGridVar(Mesh &data, int nDim)
{
	/* ==================================================

		Compute the strain rate tensor
	
	================================================== */
	data.addVariable("sinPhi");
	data.addVariable("cosPhi");
	data.addVariable("cotPhi");

	if(nDim==2)
	{
		data.addVariable("velz",0.0);
	}

	vtkSmartPointer<vtkDataArray> uR = data["velx"];
	vtkSmartPointer<vtkDataArray> uP = data["vely"];
	vtkSmartPointer<vtkDataArray> uT = data["velz"];
	vtkSmartPointer<vtkDataArray> r = data["r"];
	vtkSmartPointer<vtkDataArray> phi = data["phi"];


	for(int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{
		data["sinPhi"]->SetTuple1(i, sin(phi->GetTuple1(i)));
		data["cosPhi"]->SetTuple1(i, cos(phi->GetTuple1(i)));
		data["cotPhi"]->SetTuple1(i, tan(vtkMath::Pi()/2.0 - phi->GetTuple1(i)));
	}


	// Compute Srr
	data.firstDerivative(XAXIS, "velx", "Srr");

	// Compute Spp
	data.firstDerivative(YAXIS, "vely", "DUpDp");
	data.storeData((data["DUpDp"]+data["velx"])/data["r"],"Spp");


	// Compute Stt
	if(nDim==3)
	{
		data.firstDerivative(ZAXIS, "velz", "DUtDt");
	}
	else
	{
		data.addVariable("DUtDt",0.0);
	}

	data.storeData(	(data["DUtDt"] + data["sinPhi"]*data["velx"] + data["cosPhi"]*data["vely"])/(data["r"]*data["sinPhi"]), "Stt");


	// Compute Srp
	data.firstDerivative(YAXIS,"velx","DUrDp");
	data.firstDerivative(XAXIS,"vely","DUpDr");
	data.storeData(	0.5*(data["DUrDp"]/data["r"] + data["DUpDr"] - data["vely"]/data["r"]), "Srp");


	// Compute Spt
	if(nDim==3)
	{
		data.firstDerivative(ZAXIS, "vely", "DUpDt");
	}
	else
	{
		data.addVariable("DUpDt",0.0);
	}
	data.firstDerivative(YAXIS,"velz","DUtDp");
	data.storeData(	0.5*(( (data["DUpDt"]/data["sinPhi"]) + data["DUtDp"] - data["velz"]*data["cotPhi"])/data["r"]), "Spt");

	std::cout<<"Finished computing Spt"<<std::endl;


	// Compute Srt
	if(nDim==3)
	{
		data.firstDerivative(ZAXIS, "velx", "DUrDt");
	}
	else
	{
		data.addVariable("DUrDt",0.0);
	}
	data.firstDerivative(XAXIS,"velz","DUtDr");
	data.storeData(	0.5*(data["DUrDt"]/(data["r"]*data["sinPhi"]) + data["DUtDr"] - data["velz"]/data["r"]), "Srt");

	std::cout<<"Finished computing Srt"<<std::endl;


	/* ==================================================

		Compute the energy dissipation rate
	
	================================================== */
	// data.storeData(data["Srr"] + data["Spp"] + data["Stt"] + (2.0*(data["Srp"] + data["Srt"] + data["Spt"])), "energyDissipation");
	data.addVariable("energyDissipation",0.0);

	for(int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{

		double edis = 0.0;

		edis += data["Srr"]->GetTuple1(i)*data["Srr"]->GetTuple1(i);
		edis += data["Spp"]->GetTuple1(i)*data["Spp"]->GetTuple1(i);
		edis += data["Stt"]->GetTuple1(i)*data["Stt"]->GetTuple1(i);


		edis += 4.0*data["Srp"]->GetTuple1(i)*data["Srp"]->GetTuple1(i);
		edis += 4.0*data["Srt"]->GetTuple1(i)*data["Srt"]->GetTuple1(i);
		edis += 4.0*data["Spt"]->GetTuple1(i)*data["Spt"]->GetTuple1(i);

		data["energyDissipation"]->SetTuple1(i,edis);

	}
}