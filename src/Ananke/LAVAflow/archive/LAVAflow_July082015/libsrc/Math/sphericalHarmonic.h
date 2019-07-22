#ifndef LAVA_MATH_SPHEREHARM_H
#define LAVA_MATH_SPHEREHARM_H

#include <complex>



namespace Math
{

// template<class T> class complex;
// 0<=L<Inf
// -L<=M<=L
// Angles in radians
std::complex<double> sphericalHarmonic(int L, int M, double theta, double phi);

};

#endif