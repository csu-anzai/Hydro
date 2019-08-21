#ifndef LAVA_MATH_DOT_H
#define LAVA_MATH_DOT_H

// #include "../Geometry/Primitives/VecTN.h"
#include "../Geometry/Primitives/VecTN.h"

namespace Math
{

template <class T, int N>
T dot(VecTN<T,N> a, VecTN<T,N> b);
template <class T, int N>
T dot(VecTN<T,N>* a, VecTN<T,N>* b);


};



template <class T, int N>
T Math::dot(VecTN<T,N> a, VecTN<T,N> b)
{
	return dot(&a, &b);
}


template <class T, int N>
T Math::dot(VecTN<T,N>* a, VecTN<T,N>* b)
{
	T sum = 0.;
	for (int i = 0; i < N; i++)
	{
		sum += (*a)[i]*(*b)[i];
	}
	return sum;
}

#endif