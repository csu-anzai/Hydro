#ifndef LAVA_HEFACE_H
#define LAVA_HEFACE_H

// #include "HEEdge.h"
#include <cstdlib>
#include <ostream>
#include <string>
#include <vector>
#include "../Primitives/Vec3.h"

class HEEdge;

class HEFace
{
public:
	HEFace()
	{
		this->edge = NULL;
	}
	~HEFace()
	{
		// this->deleteFace();
	}

	/* data */
	HEEdge* edge;
	Vec3 normal;

	void print(std::ostream& out, int indent=0);
	void getAllEdges(std::vector<HEEdge*>& edgeList);
	bool deleteFace();
	void computeNormal();
	Vec3 getNormal();


};

#endif