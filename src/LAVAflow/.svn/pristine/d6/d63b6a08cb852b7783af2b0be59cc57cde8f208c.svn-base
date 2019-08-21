#include "sphericalHarmonic.h"
#include "LAVAMath.h"
#include "LAVAconstants.h"
#include <iostream>
#include <complex>
#include <cmath>



std::complex<double> Math::sphericalHarmonic(int L, int M, double theta, double phi)
{

	int Mabs = abs(M);

	std::complex<double> Ylm(0.0,0.0);

	// Get the Legendre polynomial
	double x = cos(phi);
	double Plm = Math::legendre(L,Mabs,x);

	// Get the exponential contribution
	std::complex<double> expContrib = exp(std::complex<double>(0.0,double(Mabs)*theta));

	// std::cout<<"Plm = "<<Plm<<std::endl;

	// Get the coefficient
	// double coeff = sqrt((2*L+1)*double(Math::factorial(L-Mabs))/(4.0*PI*double(Math::factorial(L+Mabs))));


	// double coeffDenom = 1.0;
	// double coeffNumer = 1.0;


	// coeffDenom = 4.0*PI;
	// for(int i=1; i<=(L+Mabs); i++)
	// {
	// 	coeffDenom *= double(i);
	// }
	// coeffNumer = 2.0*L+1.0;
	// for(int i=1; i<=(L-Mabs); i++)
	// {
	// 	coeffNumer *= double(i);
	// }

	// // for(int i = (L-Mabs+1); i<= (L+Mabs); i++)
	// // {
	// // 	coeffDenom *= double(i);
	// // }

	// if(fabs(coeffDenom)>1.0e300)
	// {
	// 	std::cout<<"[Math::sphericalHarmonic] coeffDenom is getting dangerously large: "<<coeffDenom<<std::endl;
	// }


	// double coeff = sqrt( ((2*L+1)/PI) * (1.0/coeffDenom) );


	// double coeff = sqrt(coeffNumer/coeffDenom);

	// Compute the coefficient the log way
	double coeffLog = 0.0;

	coeffLog += log((2.0*L+1.0)/(4.0*PI));

	for(int x = (L-Mabs+1); x <= (L+Mabs); x++)
	{
		coeffLog -= log(double(x));
	} 

	// Handle the square root
	coeffLog *= 0.5;

	double coeff = exp(coeffLog);


// std::cout<<"l = "<<L<<"\tcoef = "<<coeff<<"\t coeffLog = "<<coeffLog<<std::endl;

	// Combine the terms
	Ylm = coeff*Plm*expContrib;

	// If M is negative, we need to modify Ylm
	if(M<0)
	{
		if(Mabs%2==0)
			Ylm = std::conj(Ylm);
		else
			Ylm = -1.0*std::conj(Ylm);
	}

	// Return the spherical harmonic value
	return Ylm;

}
