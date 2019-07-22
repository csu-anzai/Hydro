#include "VecTN.h"
#include "Vec3.h"
#include "Edge3.h"
#include <iostream>
using namespace std;

// typedef VectTN<double,3> VectD3;
// typedef VectTN<int,3> VectI3;


int main(void)
{

	cout<<"Hello world!"<<endl;
	// VectTN<double,5> v3(3,1.0,2.0,3.0);
	// VectD3 vd3(3,6.0,7.0,8.0);
	// VectD3 vi3(1,2,3);
	Vec3 v1;
	Vec3 v2(2.0,4.0,-10.0);
	Vec3 sub;
	sub = 2*v2;

	// v2.normalize();
	v2.print(cout);
	sub.print(cout);

	cout<<"dot product: "<<dot(v2,sub)<<endl;

	Edge3 e(v2,sub);
	e.print(cout);

	cout<<"dot product: "<<dot(e.first(),e.second())<<endl;

}