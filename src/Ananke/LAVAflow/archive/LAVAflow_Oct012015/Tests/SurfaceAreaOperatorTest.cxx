#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Operators/SurfaceAreaOperator.h"
#include "../libsrc/Operators/ContourOperator.h"

#include <vtkSphereSource.h>
#include <vtkRegularPolygonSource.h>
#include <vtkSphericalTransform.h>
#include <vtkTransformFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkVersion.h>
//#include <math.h>

//double center[] = {2000.0,0.0,0.0};
double center[] = {0.0,0.0,0.0};
int circleRes = 10;
int spherePhiRes = 10;
int sphereThetaRes = 5;

/*vtkSmartPointer<vtkPolyData> generateCircle(int circleRes, double radius, double* center)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	//points->SetNumberOfPoints(circleRes-1);
	double angle = 0.0;
	double delAngle = (360./circleRes)*(vtkMath::Pi()/180.);
	for(int i = 0; i<(circleRes-1); i++)
	{
		double point[3];
		point[0]=cos(angle)*radius;
		point[1]=sin(angle)*radius;
		point[2]=0.0;
		points->InsertNextPoint(point);
		angle+=delAngle;
	}
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	for(int i = 0; i<circleRes; i++)
	{
		vtkSmartPointer<vtkPolygon> line = vtkSmartPointer<vtkPolygon>::New(); 
		line->GetPointIds()->SetNumberOfIds(2);
		
		cells->InsertNextCel(line);
	}
}*/

int main(int argc, char** argv)
{
	double radius=1.0;
	Renderer renderer(800,600);
	if(argc>1)
		radius=atof(argv[1]);
	printf("vtkVersion  = %s\n",vtkVersion::GetVTKVersion());
	cout<<"radius = "<<radius<<endl;
	VTKReader reader;
	std::string filename = "/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0200.vtk";
	vtkDataSet* ds = reader.readFile(filename.c_str());
	vtkRectilinearGrid* rg = vtkRectilinearGrid::SafeDownCast(ds);
	//ds->PrintSelf(cout,vtkIndent(0));
	//It seems contouring algorithms require point data.
	vtkSmartPointer<vtkCellDataToPointData> c2p = vtkSmartPointer<vtkCellDataToPointData>::New();
	//c2p->SetInputData(ds);
	c2p->SetInputData(rg);
	c2p->Update();
	vtkSmartPointer<vtkDataSet> pointDS=c2p->GetOutput();
	vtkSmartPointer<vtkSphericalTransform> transformer = vtkSmartPointer<vtkSphericalTransform>::New();
	//vtkSmartPointer<vtkTransformPolyDataFilter> tpd = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	vtkSmartPointer<vtkTransformFilter> tpd = vtkSmartPointer<vtkTransformFilter>::New();
	tpd->SetTransform(transformer);
	tpd->SetInputData(pointDS);
	tpd->Update();
	vtkSmartPointer<vtkDataSet> cartisianDS=tpd->GetOutput();

	ContourOperator contourOp;

        double isoValues[] = {50000.};
	contourOp.setIsoLevels(1, isoValues);
	//contourOp.setArrayName(vtkDataObject::FIELD_ASSOCIATION_POINTS,"dens");
	SurfaceAreaOperator surfaceAreaOp;
	double sphereExact=4.*vtkMath::Pi()*radius*radius;
	for(int i = 1; i < 20; i++)
	{
		vtkSmartPointer<vtkSphereSource> sphereGenerator = vtkSmartPointer<vtkSphereSource>::New();
		sphereGenerator->SetRadius(radius);
		sphereGenerator->SetCenter(center);
		sphereGenerator->SetPhiResolution(spherePhiRes);
		sphereGenerator->SetThetaResolution(sphereThetaRes);
		sphereGenerator->Update();
		vtkSmartPointer<vtkPolyData> sphereDS = sphereGenerator->GetOutput();
		surfaceAreaOp.process(sphereDS);
		cout<<"resolution = "<<spherePhiRes*sphereThetaRes<<endl;
		cout<<"Estimated Sphere Surface Area = "<<surfaceAreaOp.getSurfaceArea()<<" Exact Surface Area = "<<sphereExact<<endl;
		cout<<"Error ="<<std::abs(surfaceAreaOp.getSurfaceArea()-sphereExact)/sphereExact<<endl;
		spherePhiRes+=5;
		sphereThetaRes+=10;
	}

	double circleExact=2.*vtkMath::Pi()*radius;
	for(int i = 1; i < 20; i++)
	{
		vtkSmartPointer<vtkRegularPolygonSource> circleGenerator = vtkSmartPointer<vtkRegularPolygonSource>::New();
		circleGenerator->SetGeneratePolyline(0);
		circleGenerator->SetNumberOfSides(circleRes);
		circleGenerator->SetRadius(radius);
		circleGenerator->SetCenter(center);
		circleGenerator->Update();
		vtkSmartPointer<vtkPolyData> circleDS = circleGenerator->GetOutput();
		surfaceAreaOp.process(circleDS);
		cout<<"resolution = "<<circleRes<<endl;
		cout<<"Estimated circle circumference = "<<surfaceAreaOp.getSurfaceArea()<<" Exact circumference = "<<circleExact<<endl;
		cout<<"Error ="<<std::abs(surfaceAreaOp.getSurfaceArea()-(circleExact))/circleExact<<endl;
		circleRes+=10;
	}
	
	vtkSmartPointer<vtkDataSet> contourDS = contourOp.process(cartisianDS);
	surfaceAreaOp.process(contourDS);
	cout<<"Estimated surface area of hotb = "<<surfaceAreaOp.getSurfaceArea()<<endl;

	return 0;
}
