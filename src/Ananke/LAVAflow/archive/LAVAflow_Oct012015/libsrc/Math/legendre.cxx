#include "legendre.h"
#include <iostream>

#include <cmath>

// Evaluate the associated, non-normalized Legendre polynomial of order L and degree M at x

double Math::legendre(int L, int M, double x)
{
	// Check that we have acceptable inputs
	double Plm = 0.0;
	if(M<0 || M>L || fabs(x)>1.0)
	{
		Plm = -1e99;
		std::cerr<<"[Math::legendre] ERROR: Incorrect input parameters!"<<std::endl;
		return Plm;
	}

	// Compute Pmm
	double Pmm = 1.0;

	if(M>0)
	{
		double somx2 = sqrt((1.0-x)*(1.0+x));
		double fact = 1.0;
		for(int i=1; i<=M; i++)
		{
			Pmm = Pmm*(-fact*somx2);
			fact = fact + 2.0;
		}
	}


	if(L==M)
	{
		Plm = Pmm;
		return Plm;
	}
	else
	{        
		double Pmmp1 = x*(2.0*M+1.0)*Pmm;

		if(L == (M+1))
		{
			Plm = Pmmp1;
			return Plm;
		}
		else
		{
			double Pll = 0.0;
			for(int ll=M+2; ll<=L; ll++)
			{
				Pll = (x*(2.0*ll-1.0)*Pmmp1 - (ll+M-1.0)*Pmm)/double(ll-M);
				Pmm = Pmmp1;
				Pmmp1 = Pll;
			}
			Plm = Pll;
			return Plm;
		}
	}

}
