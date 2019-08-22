#include "Geometry.h"
#include "../Math/LAVAMath.h"
#include "../includes/LAVAconstants.h"
#include <iostream>

#include "Primitives/Vec3.h"

bool Geometry::isPointInfrontOfPlane(Vec3 point, Vec3 planePoint, Vec3 planeNormal)
{

	bool isInFront = false;

	Vec3 diffPt(point-planePoint);
	diffPt.normalize();

	double distFromPlane = Math::dot(diffPt, planeNormal);
	if(distFromPlane>=ERRMAX)
	{
		isInFront = true;
	}

	return isInFront;
}
