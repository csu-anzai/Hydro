#include "Geometry.h"
#include "../Math/LAVAMath.h"
#include "../includes/LAVAconstants.h"
#include <iostream>

#include "Primitives/Edge3.h"
#include "Primitives/Vec3.h"



bool Geometry::intersectionLinePlane(Edge3* line, Vec3* pPlane, Vec3* normalPlane,  Vec3** pIntersection)
{
	double errorMax;
	errorMax= ERRMAX;
	bool isIntersecting = false;

	// std::cout<<"[Geometry::intersectionLinePlane] Entering function"<<std::endl;

	// make sure the points arent coincident
	if(line->getLength()<errorMax)
	{
	    std::cerr<<"[Geometry::intersectionLinePlane] (Warning) Points are indistinguishable"<<std::endl;
	    return false;
	}

	double denom = Math::dot(*(line->getNormalizedRay()), *normalPlane);

	// std::cout<<line<<"\t"<<pPlane<<"\t"<<normalPlane<<"\t"<<pIntersection<<std::endl;

	if((abs(denom)/line->getLength())<ERRMAX)
	{
	    // fprintf(2,'[Warning] Line lies in the plane\n');
	    return false;
	}
	else
	{

		double tParam = Math::dot(*pPlane-*(line->first()), *normalPlane)/denom;
		
		// std::cout<<line->first()<<std::endl;
		// line->first()->print(std::cout);
		// line->second()->print(std::cout);
		// std::cout<<pIntersection<<std::endl;

		**pIntersection = *line->first() + tParam*(*line->getNormalizedRay());

	    if(tParam>-errorMax && tParam < 1+errorMax)
	    {
			isIntersecting = true;
	    }
	    else
	    {    
			isIntersecting = false;
	    }
	}

	// std::cout<<"[Geometry::intersectionLinePlane] Leaving function"<<std::endl;

	return isIntersecting;

}