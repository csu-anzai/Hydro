#ifndef LAVA_MATH_CROSS_H
#define LAVA_MATH_CROSS_H

// #include "../Geometry/Primitives/VecTN.h"
#include "../Geometry/Primitives/Vec3.h"

namespace Math
{

inline Vec3 crossProduct(Vec3 A, Vec3 B);

};

// 3D cross product
Vec3 Math::crossProduct(Vec3 A, Vec3 B)
{
	Vec3 res;

	res[0] = A[1]*B[2] - A[2]*B[1];
	res[1] = A[2]*B[0] - A[0]*B[2];
	res[2] = A[0]*B[1] - A[1]*B[0];

	return res;
}

#endif