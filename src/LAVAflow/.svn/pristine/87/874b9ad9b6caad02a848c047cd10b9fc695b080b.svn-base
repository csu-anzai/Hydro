#include <math.h>
#include <stdio.h>
#include "iso.h" 

int main()
{
	int i,j,k;
	int xdim = 110;
	int ydim = 110;
	int zdim = 110;
	double xmin = -1;
	double xmax = 1;
	double ymin = -1;
	double ymax = 1;
	double zmin = -1;
	double zmax = 1;
	double xval,yval,zval,dist;
	double scaling = 1;
	double xcent = (xmax+xmin)/2;
	double ycent = (ymax + ymin)/2; 
	double zcent	= (zmax+zmin)/2;	
	double scalfield [xdim*ydim*zdim];
	double xint = (xmax-xmin)/xdim;
	double yint = (ymax - ymin)/ydim;
	double zint = (zmax - zmin)/zdim;
	double area = 0;
	double threshold = 0.4;
	

	zval = zmin;
	for (k = 0; k <= zdim - 1; k ++)
	{		
		yval = ymin;
		for (j = 0; j <= ydim -1; j ++)
		{
			xval = xmin;
			for(i = 0; i <= xdim -1; i ++)
			{
				dist = sqrt(pow((xval-xcent),2)+pow((yval-ycent),2)+pow((zval-zcent),2));				
				scalfield[i+(j*xdim)+(k*xdim*ydim)] = dist;
				xval = xval + xint;
			}
			yval = yval + yint;
		}
		zval = zval + zint;
	}
		
	iso_surface(scalfield, &xdim, &ydim, &zdim, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax, &threshold, &area);
	
}


