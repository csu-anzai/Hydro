#ifndef LAVA_EDGE3_H
#define LAVA_EDGE3_H

#include "Vec3.h"
#include <ostream>

class Edge3
{

public:

	Edge3(Vec3 vStart, Vec3 vEnd);

	Edge3(Vec3 vStart, Vec3 normalizedRayIn, double lengthIn);

	~Edge3();

	/* data */

protected:

	Vec3* vertices[2];
	Vec3* ray;
	Vec3* normalizedRay;
	double length;


	void computeRay();
	void computeLength();

public:

	Vec3*  first();
	Vec3*  second();
	Vec3*  operator[](const int i);
	Vec3*  getRay();
	Vec3*  getNormalizedRay();
	double getLength();
	void   print(std::ostream& out, int indent);

};



#endif