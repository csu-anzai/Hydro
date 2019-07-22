#include <iostream>
#include <vector>
#include <complex>
#include "../libsrc/Math/LAVAMath.h"
#include "../libsrc/includes/LAVAconstants.h"

int main(void)
{

	int LMin = 0;
	int LMax = 30;
	double* Espec = new double[LMax+1];

	for(int i=0; i<=LMax; i++){ Espec[i] = 0.0; }

	double thetaMin = 0.0;
	double thetaMax = 2.0*PI;
	double phiMin = 0.0;
	double phiMax = PI;
	double theta = 0.0;
	double phi = 0.0;

	int nPtsTheta = 50;
	int nPtsPhi = 50;

	double dPhi = (phiMax-phiMin)/double(nPtsPhi-1);
	double dTheta = (thetaMax-thetaMin)/double(nPtsTheta-1);


	// Test the factorial function
	for(int i=0; i<=20; i++)
	{
		std::cout<<i<<"\t"<<Math::factorial(i)<<std::endl;
	}

	// Move over multipole order
	for(int l=LMin; l<=LMax; l++)
	{
		std::cout<<l<<std::endl;
		// Loop over degree
		for(int m=-l; m<=l; m++)
		{
			// std::cout<<"\t"<<m<<std::endl;

			std::complex<double> sum(0,0);
			int counter = 0;
			// Loop over cells to compute the sum
			for(int t = 0; t < nPtsTheta; t++)
			{
			
				theta = thetaMin + dTheta*t;

				for(int p = 0; p < nPtsPhi; p++)
				{
			

					phi = phiMin + dPhi*p;

					// Get the spherical harmonic
					std::complex<double> sh = Math::sphericalHarmonic(l,m,theta,phi);	
					
					// std::cout<<"sh: "<<sh<<std::endl;

					// Take the conjugate
					sh = std::conj(sh);

					// Add to multipole sum
					sum += sh*sin(phi)*dPhi*dTheta;

					counter++;
				}
			}

			// std::cout<<"Visited "<<counter<<" cells"<<std::endl;
			// Perform the norm and square operations
			double energyValue = sqrt(std::real(sum)*std::real(sum) + std::imag(sum)*std::imag(sum));
			energyValue *= energyValue;

			// Store in Espec
			Espec[l] += energyValue;
		}
	}



	// Print Espec
	std::cout<<"L\tE"<<std::endl;
	for(int l=LMin; l<=LMax; l++)
	{
		std::cout<<l<<"\t"<<Espec[l]<<"\t"<<sqrt(Espec[l])<<std::endl;
	}



	// Cleanup
	delete [] Espec;

}