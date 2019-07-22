#include <vtkType.h>
#include <vtkIdList.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkCurvatures.h>

#include "BaseOperator.h"
#include "CurvatureOperator.h"

CurvatureOperator::CurvatureOperator()
{
	curvatures=vtkSmartPointer<vtkCurvatures>::New();
	curvatures->SetCurvatureTypeToMean();
}

double distance2D(double* a, double* b)
{
	double delX = b[0]-a[0];
	double delY = b[1]-a[1];
	return sqrt(delX*delX+delY*delY); 
}

vtkSmartPointer<vtkIdList> getConnectedPoints(vtkSmartPointer<vtkPolyData> mesh, int point)
{
	vtkSmartPointer<vtkIdList> connectedVertices =
		vtkSmartPointer<vtkIdList>::New();	
	vtkSmartPointer<vtkIdList> cellIdList =
	      vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(point, cellIdList);
	//cout<<"Number of cell Ids "<<cellIdList->GetNumberOfIds()<<endl;

	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{

		vtkSmartPointer<vtkIdList> pointIdList =
		vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);
	//	cout<<"Number of point Ids "<<pointIdList->GetNumberOfIds()<<endl;

		if(pointIdList->GetId(0) != point)
		{
			connectedVertices->InsertNextId(pointIdList->GetId(0));
		}
		else
		{
			connectedVertices->InsertNextId(pointIdList->GetId(1));
		}
	}
	return connectedVertices;
}

vtkSmartPointer<vtkDataSet> CurvatureOperator::process(const vtkSmartPointer<vtkDataSet> ds)
{
	vtkSmartPointer<vtkPolyData> pd = vtkPolyData::SafeDownCast(ds);
	//If pd is null then we it's not a polydata object
	//and therfore we can't compute the Curvature.
	if(!pd)
	{
		cout<<"input to CurvatureOperator is not polydata!"<<endl;
		return pd;
	}
	
	vtkSmartPointer<vtkPolyData> retval = vtkSmartPointer<vtkPolyData>::New();
	retval->DeepCopy(pd);
	
	int numPts = pd->GetNumberOfPoints();
	// Empty array check
	if (numPts==0||(pd->GetNumberOfPolys()==0 && pd->GetNumberOfLines()==0))
	{
		cout<<"Input to CurvatureOperator has no polygons or points!"<<endl;
		return retval;
	}

	// vtkData
	vtkSmartPointer<vtkIdList> vertices = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> neighborVertices = vtkSmartPointer<vtkIdList>::New();
	//Holds current verteces index
	int v[3];
	//Stores each vertex xyz values.
	double points[3][3];
	double lengthBackward=0.0;
	double lengthForward=0.0;
	double thetaBackward=0.0;
	double thetaForward=0.0;
	double theta=0.0;
	double delBackward=0.0;
	double delForward=0.0;

	//retval->BuildLinks();
	//data init
	//int face = 0;
	int pointID=0;
	int numPoints = retval->GetNumberOfPoints();
	int nv = 0;
	//get the first face to determine the dimensionality
	retval->GetCellPoints(0,vertices);
	nv = vertices->GetNumberOfIds();
	
	//we have a 2d mesh
	if(nv==2||nv>3)
	{
		vtkDoubleArray* meanCurvature = vtkDoubleArray::New();
		meanCurvature->SetName("Mean_Curvature");
		meanCurvature->SetNumberOfComponents(1);
		meanCurvature->SetNumberOfTuples(numPoints);
		double *meanCurvatureData = meanCurvature->GetPointer(0);
		int curvatureOne = 0;
		//loop through each point 
		int totalPointsWithNeighbors=0;
		for (pointID = 0; pointID < numPoints; pointID++)
		{
			neighborVertices=getConnectedPoints(retval,pointID);
	
			if(neighborVertices->GetNumberOfIds()==0 || neighborVertices->GetNumberOfIds()==1)
			{
				cout<<"2 neighbors not found curvature is 0"<<endl;
				//no neighbors = no curvature
				meanCurvatureData[pointID]=0.0;
				continue;
			}
		
			retval->GetPoint(neighborVertices->GetId(0),points[0]);
			retval->GetPoint(pointID,points[1]);
			retval->GetPoint(neighborVertices->GetId(1),points[2]);
			//cout<<"x1 = "<<points[0][0]<<" x2 = "<<points[1][0]<<" x3 = "<<points[2][0]<<endl;
			//cout<<"y1 = "<<points[0][1]<<" y2 = "<<points[1][1]<<" y3 = "<<points[2][1]<<endl;
#if 0
			lengthBackward=distance2D(points[0],points[1]);
			lengthForward=distance2D(points[1],points[2]);
			thetaBackward=std::atan(std::abs((points[0][0]-points[1][0])/(points[0][1]-points[1][1])));
			thetaForward=std::atan(std::abs((points[2][0]-points[1][0])/(points[2][1]-points[1][1])));
			theta=0.5 * (thetaBackward+thetaForward);
			delBackward =std::abs(thetaBackward-theta);
			delForward =std::abs(thetaForward-theta);
			meanCurvatureData[pointID]=0.5*((delBackward/(lengthBackward))+(delForward/(lengthForward)));
			cout<<"length b="<<lengthBackward<<" length f ="<<lengthForward<<endl;
			cout<<"theta b="<<thetaBackward<<" theta f ="<<thetaForward<<" theta = "<<theta<<endl;
			cout<<"delta b="<<delBackward<<" delta f ="<<delForward<<" meanCurvature = "<<meanCurvatureData[pointID]<<endl<<endl;
			//cout<<"two neighbors"<<endl;
#else

			double tan1[2];
			double tan2[2];
			double delXBackward =(points[1][0]-points[0][0]);
			double delYBackward =(points[1][1]-points[0][1]);
			double delXForward =(points[2][0]-points[1][0]);
			double delYForward =(points[2][1]-points[1][1]);
			lengthBackward=sqrt(delXBackward*delXBackward+delYBackward*delYBackward);
			lengthForward=sqrt(delXForward*delXForward+delYForward*delYForward);
			tan1[0]=delXBackward/lengthBackward;
			tan1[1]=delYBackward/lengthBackward;
			tan2[0]=delXForward/lengthForward;
			tan2[1]=delYForward/lengthForward;
			theta = std::acos(tan1[0]*tan2[0]+tan1[1]*tan2[1]);
			double tmp1[2];
			double tmp2[2];
			//tmp1[0]=delXBackward/2.;
			//tmp1[1]=delYBackward/2.;
			//tmp2[0]=delXForward/2.;
			//tmp2[1]=delYForward/2.;
			//tmp1[0]=(point[0][0]+point[1][0])/2.;
			//tmp1[1]=(point[0][1]+point[1][1])/2.;
			//tmp2[0]=(point[2][0]+point[1][0])/2.;
			//tmp2[1]=(point[2][1]+point[1][1])/2.;
			//double dist = sqrt(tmp1[0]*tmp1[0]+tmp1[1]*tmp1[1])+sqrt(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]);
			double dist =(lengthBackward+lengthForward)/2.;

			//double dist = sqrt(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1])-sqrt(tmp1[0]*tmp1[0]+tmp1[1]*tmp1[1]);
			meanCurvatureData[pointID] = theta/dist;
			cout<<"theta = "<<theta<<endl;
			cout<<"length forward = "<<lengthForward<<" length backward = "<<lengthBackward<<endl;
			cout<<"Mean Curvature = "<<meanCurvatureData[pointID]<<endl;
			if(meanCurvatureData[pointID]>0.9999 &&meanCurvatureData[pointID]<1.001)
				curvatureOne++;
#endif
			totalPointsWithNeighbors++;
		}
		cout<<"Points with 2 neighbors "<<totalPointsWithNeighbors<<" Total number of points = "<<numPoints<<endl;
		cout<<"Points with approximately 1.0 curvature "<<curvatureOne<<" Total number of points = "<<numPoints<<endl;
		retval->GetPointData()->AddArray(meanCurvature);
		retval->GetPointData()->SetActiveScalars("Mean_Curvature");

		return retval;
	}
	curvatures->SetInputData(pd);
	curvatures->Update();
	return curvatures->GetOutput();
}
void CurvatureOperator::useMeanCurvature()
{
curvatures->SetCurvatureTypeToMean();
}
void CurvatureOperator::useMaximumCurvature()
{
curvatures->SetCurvatureTypeToMaximum();
}
void CurvatureOperator::useMinimumCurvature()
{
curvatures->SetCurvatureTypeToMinimum();
}
