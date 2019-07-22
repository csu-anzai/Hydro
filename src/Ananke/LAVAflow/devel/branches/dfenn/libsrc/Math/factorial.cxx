#include "factorial.h"

long Math::factorial(int value)
{

	long i, product;

	product = 1;
	for(i=1; i<=value; i++)
	{
		product *= i;
	}
	return product;
}