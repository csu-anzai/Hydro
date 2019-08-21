#include "Edge3.h"
#include <iostream>

Edge3::Edge3(Vec3 vStart, Vec3 vEnd)
{
	vertices[0] = new Vec3(vStart);
	vertices[1] = new Vec3(vEnd);
	this->computeRay();
	this->computeLength();

}

Edge3::Edge3(Vec3 vStart, Vec3 normalizedRayIn, double lengthIn)
{
	vertices[0] = new Vec3(vStart);
	this->normalizedRay = new Vec3(normalizedRayIn);
	this->length = lengthIn;
	this->ray = new Vec3((*this->normalizedRay)*this->length);
	vertices[1] = new Vec3(*vertices[1] + *this->ray);
}

Edge3::~Edge3()
{
	delete vertices[0];
	delete vertices[1];
	delete ray;
	delete normalizedRay;
}

void Edge3::computeRay()
{
	this->ray = new Vec3(*vertices[1]-*vertices[0]);
	this->normalizedRay = new Vec3(*this->ray);
	this->normalizedRay->normalize();
}

void Edge3::computeLength()
{
	this->length = this->ray->norm();
}

Vec3* Edge3::first()
{
	return this->vertices[0];
}


Vec3* Edge3::second()
{
	return this->vertices[1];
}

Vec3* Edge3::operator[](const int i)
{
	switch(i)
	{
		case 1:
			return this->first();
		case 2:
			return this->second();
		default:
    			std::cerr<<"[Edge3->operator[]] Accessor value of "<<i<<" is greater than max value "<<1<<" (2 elements total)"<<std::endl;
    			return NULL;

	}
}


Vec3*  Edge3::getRay()
{
	return this->ray;
}
Vec3*  Edge3::getNormalizedRay()
{
	return this->normalizedRay;
}
double Edge3::getLength()
{
	return this->length;
}

void Edge3::print(std::ostream& out, int indent=0)
{

	out<<std::string(indent,'\t')<<"Edge information:\n";
	out<<std::string(indent+1,'\t')<<"Starting Vertex: \n";
	out<<std::string(indent+2,'\t');
	this->first()->print(out);
	out<<std::string(indent+1,'\t')<<"Ending Vertex:\n";
	out<<std::string(indent+2,'\t');
	this->second()->print(out);
	out<<std::string(indent+1,'\t')<<"Normalized Ray:\n";
	out<<std::string(indent+2,'\t');
	this->getNormalizedRay()->print(out);
	out<<std::string(indent+1,'\t')<<"Length:\n";
	out<<std::string(indent+2,'\t')<<this->getLength();
	out<<"\n\n";

	out<<std::endl;
}