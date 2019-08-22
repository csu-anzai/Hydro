#ifndef LAVA_GEOMETRY_H
#define LAVA_GEOMETRY_H

#include "Primitives/Edge3.h"
#include "Primitives/Vec3.h"

namespace Geometry
{
	bool intersectionLinePlane(Edge3* line, Vec3* pPlane, Vec3* normalPlane,  Vec3** pIntersection);
	int intersectionPolyhedronPlane(double** vEdgeStart, double** vEdgeEnd, int nEdges, double* pPlane, double* normalPlane, double** pIntersection);
	// void intersectionPolygonPlane(void);
	void clipPolygonUsingPlane(void);
	// void isPointInfrontOfPlane(void);
	bool isPointInfrontOfPlane(Vec3 point, Vec3 planePoint, Vec3 planeNormal);
};

#endif