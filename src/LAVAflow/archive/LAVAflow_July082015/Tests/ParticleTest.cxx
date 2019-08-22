#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <list>

#include "vtkDataArray.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkUnstructuredGridWriter.h"

#include "../libsrc/Particles/Particles.h"
#include "../libsrc/includes/LAVAconstants.h"
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/Utilities/LAVAUtil.h"

int determineViableParticles(std::string directory, std::string modelName, int startNum, int endNum, std::list<int>& particleIndices)
{

	char prtFileName[256], expFileName[256];

	bool firstPass =  true;

	for(int fileNum = startNum; fileNum <= endNum; fileNum++)
	{

		std::cout<<double(fileNum-startNum)/double(endNum-startNum)<<"... ";
		std::cout.flush();

		// Create file names
		std::sprintf(prtFileName,"%s/%s.prt.%04d.vtk",directory.c_str(),modelName.c_str(),fileNum);
		std::sprintf(expFileName,"%s/output/%s.plt.%04d.expstat",directory.c_str(),modelName.c_str(),fileNum);


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

		// If it's the first file, add the particles to the list that lie in the gain region
		if(firstPass)
		{

			for(int p=0; p<parts.getNumberOfParticles(); p++)
			{
				double partRadius = radii->GetTuple1(p);

				// If the particle lies in the gain region
				if(partRadius >= gainRadius && partRadius <= shockRadiusMin)
				{
					// Add the particle index to the list
					particleIndices.push_back(p);
				}

			}
			firstPass = false;
		}
		// If this is not the first file, we need to remove particles that do not fall inside the gain region
		else
		{
			std::list<int>::iterator it = particleIndices.begin();
			while(it != particleIndices.end())
			{
				int p = *it;

				double partRadius = radii->GetTuple1(p);


				if(partRadius < gainRadius || partRadius > shockRadiusMin)
				{
					it = particleIndices.erase(it);
				}
				else
				{
					++it;
				}

			}

		}

	} // end loop over files

	std::cout<<std::endl;
	// Return the number of viable particles over the entire temporal region of interest
	return particleIndices.size();
}


int main(void)
{

	std::cout<<"Hello world!"<<std::endl;
	
	// Setup information for time sequence analysis
	int startingFileNumber = 250;
	int currentFileNumber = startingFileNumber+1;
	int finalFileNumber = 500;
	int numberOfFiles = finalFileNumber - currentFileNumber + 1;

	std::cout<<"Number of files to use: "<<numberOfFiles<<std::endl;



	// Setup information for output
	std::string outputFileNameFmt = "lagrangian_autocorrelation.%04d.dat";
	char outputFileName[256];
	std::sprintf(outputFileName,outputFileNameFmt.c_str(),startingFileNumber);

	std::ofstream output;
	output.open(outputFileName);
	output<<"t"<<"\t"<<"tau"<<"\t"<<"Rx"<<"\t"<<"Ry"<<"\t"<<"Rz"<<"\t"
		  <<"D1x"<<"\t"<<"D1y"<<"\t"<<"D1z"<<"\t"
		  <<"D2x"<<"\t"<<"D2y"<<"\t"<<"D2z"<<"\t"
		  <<"D3x"<<"\t"<<"D3y"<<"\t"<<"D3z"
		  <<std::endl;

	output.precision(15);
	output.setf(std::ios::scientific,std::ios::floatfield);

	// Initialize some variables
	double timeInitial = 0.0;
	double timeCurrent = 0.0;

	/* ==================================================

		Determine viable particles for the entire run
	
	================================================== */
	

	std::string directory = "/data1/sne/HOTB/2d/b163d2a3s530SG1MHWr2LBp";
	std::string modelID = "b163d2a3s530SG1MHWr2LBp";
	std::list<int> particleIndices;

	std::cout<<"Finding viable particles for temporal extent..."; std::cout.flush();
	int nViableParticles = determineViableParticles(directory, modelID, startingFileNumber, finalFileNumber, particleIndices);
	std::cout<<"done!"<<std::endl;
	std::cout<<"Number of viable particles: "<<nViableParticles<<std::endl;

	// std::list<int>::iterator pIt;
	// for(pIt = particleIndices.begin(); pIt != particleIndices.end(); pIt++)
	// {
	// 	std::cout<<*pIt<<std::endl;
	// }


	/* ==================================================

		Get initial particles
	
	================================================== */
	
	std::string formatStr = "/data1/sne/HOTB/2d/%s/%s.prt.%04d.vtk";
	char filename[256];

	// Get appropriate file for this step
	std::sprintf(filename,formatStr.c_str(),modelID.c_str(),modelID.c_str(),startingFileNumber);

	// Load in the particles
	Particles partsInitial(CS_SPHERE,filename);
	timeInitial = partsInitial.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);

	// Mask the particles based on the initial particles in the gain region that stay the entire time
	partsInitial.addVariable("mask",0.0);
	vtkSmartPointer<vtkDataArray> mask = partsInitial["mask"];

	std::list<int>::iterator pIt;
	for(pIt = particleIndices.begin(); pIt != particleIndices.end(); pIt++)
	{
		mask->SetTuple1(*pIt,1.0);
	}


	std::cout<<"Particle information: "<<std::endl;
	std::cout<<"\tNumber of particles: "<<partsInitial.getDataSet()->GetNumberOfPoints()<<std::endl;
	std::cout<<"\tNumber after masking:"<<partsInitial.sum("mask")<<std::endl;


	/* ==================================================

		Begin the time sequence analysis
	
	================================================== */

	// Print some information
	std::cout<<"Initial file: "<<filename<<std::endl;
	std::cout<<"Initial time: "<<timeInitial<<std::endl;

	// Compute the primed particle velocities
	double varOrig[3];
	varOrig[0] = partsInitial.variance("pvx","mask");
	varOrig[1] = partsInitial.variance("pvy","mask");
	varOrig[2] = partsInitial.variance("pvz","mask");


	std::cout<<"Initial Variances: "<<std::endl;
	std::cout<<"\tx: "<<varOrig[0]<<std::endl;
	std::cout<<"\ty: "<<varOrig[1]<<std::endl;
	std::cout<<"\tz: "<<varOrig[2]<<std::endl;

	partsInitial.storeData(partsInitial["pvx"]-partsInitial.mean("pvx","mask"),"pvxPrime");
	partsInitial.storeData(partsInitial["pvy"]-partsInitial.mean("pvy","mask"),"pvyPrime");
	partsInitial.storeData(partsInitial["pvz"]-partsInitial.mean("pvz","mask"),"pvzPrime");

	// Get pointers to the data sets
	vtkSmartPointer<vtkDataArray> pvxPrimeOrig = partsInitial["pvxPrime"];
	vtkSmartPointer<vtkDataArray> pvyPrimeOrig = partsInitial["pvyPrime"];
	vtkSmartPointer<vtkDataArray> pvzPrimeOrig = partsInitial["pvzPrime"];
	vtkSmartPointer<vtkDataArray> pvxOrig = partsInitial["pvx"];
	vtkSmartPointer<vtkDataArray> pvyOrig = partsInitial["pvy"];
	vtkSmartPointer<vtkDataArray> pvzOrig = partsInitial["pvz"];


	/* ==================================================

		Loop over the remaining particle files
		and compute the time difference and
		needed quantities
	
	================================================== */
	int counter = 0;
	for(int fileNum = startingFileNumber; fileNum <= finalFileNumber; fileNum++)
	{
		// Get the current file name
		std::sprintf(filename,formatStr.c_str(),modelID.c_str(),modelID.c_str(),fileNum);
		std::cout<<"Processing file: "<<filename<<std::endl;

		// Get the current files
		Particles partsCurrent(CS_SPHERE,filename);
		timeCurrent = partsCurrent.getDataSet()->GetFieldData()->GetArray("TIME")->GetTuple1(0);

		// Store the particle mask that was calculated at the start
		partsCurrent.storeData(partsInitial["mask"],"mask");

		// std::cout<<"Number of particles in current mask: "<<partsCurrent.sum("mask")<<std::endl;

		// get dt and store it in tau
		double tau = timeCurrent - timeInitial;


		/* ==================================================

			Compute the Lagrangian autocorrelation
			function 
			(1)	Compute the primed velocities
			(2)	Compute the product of the primed velocities for
					the original data and the current data
			(3)	Compute the mean of that product
			(4)	Normalize by the variance
		
		================================================== */
		// Compute the primed particle velocities
		double varCurr[3];
		varCurr[0] = partsCurrent.variance("pvx","mask");
		varCurr[1] = partsCurrent.variance("pvy","mask");
		varCurr[2] = partsCurrent.variance("pvz","mask");

		// std::cout<<"Variances: "<<std::endl;
		// std::cout<<"\tx: "<<varCurr[0]<<std::endl;
		// std::cout<<"\ty: "<<varCurr[1]<<std::endl;
		// std::cout<<"\tz: "<<varCurr[2]<<std::endl;

		// Compute primed velocities of current particles
		partsCurrent.storeData(partsCurrent["pvx"]-partsCurrent.mean("pvx","mask"),"pvxPrime");
		partsCurrent.storeData(partsCurrent["pvy"]-partsCurrent.mean("pvy","mask"),"pvyPrime");
		partsCurrent.storeData(partsCurrent["pvz"]-partsCurrent.mean("pvz","mask"),"pvzPrime");

		// Get pointers to the data sets
		vtkSmartPointer<vtkDataArray> pvxPrimeCurr = partsCurrent["pvxPrime"];
		vtkSmartPointer<vtkDataArray> pvyPrimeCurr = partsCurrent["pvyPrime"];
		vtkSmartPointer<vtkDataArray> pvzPrimeCurr = partsCurrent["pvzPrime"];

		vtkSmartPointer<vtkDataArray> pvxCurr = partsCurrent["pvx"];
		vtkSmartPointer<vtkDataArray> pvyCurr = partsCurrent["pvy"];
		vtkSmartPointer<vtkDataArray> pvzCurr = partsCurrent["pvz"];

		// First order structure function
		double D1[3];
		D1[0] = 0.0;
		D1[1] = 0.0;
		D1[2] = 0.0;

		// Second order structure function
		double D2[3];
		D2[0] = 0.0;
		D2[1] = 0.0;
		D2[2] = 0.0;

		// Third order structure function
		double D3[3];
		D3[0] = 0.0;
		D3[1] = 0.0;
		D3[2] = 0.0;

		// Autocorrelation
		double R[3];
		R[0] = 0.0;
		R[1] = 0.0;
		R[2] = 0.0;

		for(int i=0; i<partsInitial.getNumberOfParticles(); i++)
		{
			if(mask->GetTuple1(i)>0.5)
			{
				R[0] += pvxPrimeOrig->GetTuple1(i)*pvxPrimeCurr->GetTuple1(i);
				R[1] += pvyPrimeOrig->GetTuple1(i)*pvyPrimeCurr->GetTuple1(i);
				R[2] += pvzPrimeOrig->GetTuple1(i)*pvzPrimeCurr->GetTuple1(i);

				double diff[3];
				diff[0] = pvxCurr->GetTuple1(i) - pvxOrig->GetTuple1(i);
				diff[1] = pvyCurr->GetTuple1(i) - pvyOrig->GetTuple1(i);
				diff[2] = pvzCurr->GetTuple1(i) - pvzOrig->GetTuple1(i); 

				D1[0] += diff[0];
				D1[1] += diff[1];
				D1[2] += diff[2];

				D2[0] += diff[0]*diff[0];
				D2[1] += diff[1]*diff[1];
				D2[2] += diff[2]*diff[2];

				D3[0] += diff[0]*diff[0]*diff[0];
				D3[1] += diff[1]*diff[1]*diff[1];
				D3[2] += diff[2]*diff[2]*diff[2];

			}
		}

		// Finish computing the autocorrelation
		for(int i=0; i<3; i++)
		{
			// Divide by the number of particles
			R[i] /= partsInitial.sum("mask");
			D1[i] /= partsInitial.sum("mask");
			D2[i] /= partsInitial.sum("mask");
			D3[i] /= partsInitial.sum("mask");

			// normalize by the product of the sqrt of variances
			R[i] /= sqrt(varOrig[i]*varCurr[i]);
		}


		// Write to output file
		

		output 	<<timeCurrent<<"\t"
				<<tau<<"\t"
				<<R[0]<<"\t"
				<<R[1]<<"\t"
				<<R[2]<<"\t"
				<<D1[0]<<"\t"
				<<D1[1]<<"\t"
				<<D1[2]<<"\t"
				<<D2[0]<<"\t"
				<<D2[1]<<"\t"
				<<D2[2]<<"\t"
				<<D3[0]<<"\t"
				<<D3[1]<<"\t"
				<<D3[2]
				<<std::endl;

		// Increase counter
		counter++;

	}


	output.close();




	// vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	// writer->SetFileName("particles.vtk");
	// writer->SetInputData(partsInitial.getDataSet());
	// writer->Write();

	std::cout<<"Finished"<<std::endl;
}