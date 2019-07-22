#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Readers/FlashReader.h"
#include "../libsrc/Operators/ContourOperator.h"
#include "../libsrc/Operators/CurvatureOperator.h"

#include <vtkSphericalTransform.h>
#include <vtkSphereSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransformFilter.h>
#include <vtkFillHolesFilter.h>
#include <vtkUniformGrid.h>
#include <vtkOverlappingAMR.h>
#include <vtkRegularPolygonSource.h>
#include <vtkRectilinearGrid.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>

#include <vtkCellDataToPointData.h>
#include <vtkVersion.h>
#include <vtkMath.h>

vtkSmartPointer<vtkPolyData> generateCircle(int circleRes, double radius, double* center, int sweepRes=0)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	//points->SetNumberOfPoints(circleRes-sweepRes);
	double angle = 0.0;
	double delAngle = ((vtkMath::Pi()*2.)/circleRes);//*(vtkMath::Pi()/180.);
	for(int i = 0; i<(circleRes-sweepRes); i++)
	{
		double point[3];
		point[0]=cos(angle)*radius+center[0];
		point[1]=sin(angle)*radius+center[1];
		point[2]=center[2];
		points->InsertNextPoint(point);
		angle+=delAngle;
	}
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	for(int i = 0; i<circleRes-sweepRes; i++)
	{
		vtkSmartPointer<vtkPolyLine> line = vtkSmartPointer<vtkPolyLine>::New(); 
		line->GetPointIds()->SetNumberOfIds(2);
		line->GetPointIds()->SetId(0,i);
		if(i==circleRes-sweepRes-1&&!sweepRes)
		{
			cout<<"at final point !"<<endl;
			line->GetPointIds()->SetId(1,0);
		}
		else if(i==circleRes-sweepRes-1&&sweepRes)
		{
			cout<<"at final point with sweepRes = "<<sweepRes<<"!"<<endl;
			break;
		}
		else
		{
			line->GetPointIds()->SetId(1,i+1);
		}
		cells->InsertNextCell(line);
	}
	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();	
	pd->SetPoints(points);
	pd->SetLines(cells);
	return pd;
}
int main(int argc, char** argv)
{
	printf("vtkVersion  = %s\n",vtkVersion::GetVTKVersion());
	//std::string filename = "/home/asy10/VTKData/Data/RectGrid2.vtk";
	//std::string filename = "/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0050.vtk";
	std::string filename = "/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0200.vtk";
	//std::string filename = "/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0253.vtk";
	std::string flashFilename = "/home/asy10/circularDataset.hdf5";

	VTKReader reader;
	FlashReader flashReader;
	NullSelector selector;
	Renderer renderer(800,600);
	ContourOperator contourOp;
        CurvatureOperator curveOp;
        //double isoValues[] = {0.5,500.,5000., 50000.};
        //double isoValues[] = {0.5,500.,5000., 50000., 100000.};
        double isoValues[] = {50000.};
        //double isoValues[] = {8.03907e+9};
	//contourOp.setIsoLevels(4, isoValues);
	//contourOp.setIsoLevels(5, isoValues);
	contourOp.setIsoLevels(1, isoValues);
	//contourOp.setArrayName(vtkDataObject::FIELD_ASSOCIATION_POINTS,"dens");
	contourOp.setArrayName(vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,"dens");

	vtkDataSet* ds = reader.readFile(filename.c_str());
	vtkRectilinearGrid* rg = vtkRectilinearGrid::SafeDownCast(ds);
	//ds->PrintSelf(cout,vtkIndent(0));
	//It seems contouring algorithms require point data.
	vtkSmartPointer<vtkCellDataToPointData> c2p = vtkSmartPointer<vtkCellDataToPointData>::New();
	//c2p->SetInputData(ds);
	c2p->SetInputData(rg);
	c2p->Update();
	vtkSmartPointer<vtkDataSet> pointDS=c2p->GetOutput();
	//cout<<"Point DS-----"<<endl;
	//pointDS->PrintSelf(cout,vtkIndent(5));
	//cout<<"End Point DS-----"<<endl;

	double range[2];
	pointDS->GetScalarRange(range);
	cout<<"PointDS min = "<<range[0]<<" max = "<<range[1]<<endl;
	#if 0 
	vtkSmartPointer<vtkRectilinearGridGeometryFilter> rggf = vtkSmartPointer<vtkRectilinearGridGeometryFilter>::New();
	rggf->SetExtent(rg->GetExtent());
	rggf->SetInputData(pointDS);
	rggf->Update();
	vtkSmartPointer<vtkDataSet> plyDS = rggf->GetOutput();
	cout<<"PLY DS-----"<<endl;
	plyDS->PrintSelf(cout,vtkIndent(5));
	cout<<"End PLY DS-----"<<endl;
	plyDS->GetScalarRange(range);
	cout<<"PLYDS min = "<<range[0]<<" max = "<<range[1]<<endl;

	vtkSmartPointer<vtkRectilinearGridToTetrahedra> rggf = vtkSmartPointer<vtkRectilinearGridToTetrahedra>::New();
	rggf->SetInputData(pointDS);
	rggf->Update();
	vtkSmartPointer<vtkDataSet> tetraDS = reinterpret_cast<vtkDataSet*>(rggf->GetOutput());
	cout<<"tetra DS-----"<<endl;
	tetraDS->PrintSelf(cout,vtkIndent(5));
	cout<<"End tetra DS-----"<<endl;
	tetraDS->GetScalarRange(range);
	cout<<"tetraDS min = "<<range[0]<<" max = "<<range[1]<<endl;
	#endif 

	vtkSmartPointer<vtkSphericalTransform> transformer = vtkSmartPointer<vtkSphericalTransform>::New();
	//vtkSmartPointer<vtkTransformPolyDataFilter> tpd = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	vtkSmartPointer<vtkTransformFilter> tpd = vtkSmartPointer<vtkTransformFilter>::New();
	tpd->SetTransform(transformer);
	tpd->SetInputData(pointDS);
	//tpd->SetInputData(opDS);
	//tpd->SetInputData(plyDS);
	//tpd->SetInputData(tetraDS);
	tpd->Update();
	vtkSmartPointer<vtkDataSet> cartisianDS=tpd->GetOutput();
	//cout<<"Cartisian DS-----"<<endl;
	//cartisianDS->PrintSelf(cout,vtkIndent(5));
	//cout<<"End Cartisian DS-----"<<endl;
	cartisianDS->GetScalarRange(range);
	cout<<"CartisianDS min = "<<range[0]<<" max = "<<range[1]<<endl;

	//vtkSmartPointer<vtkPolyData> polydataDS = vtkPolyData::SafeDownCast(cartisianDS);
	//polydataDS->BuildLinks(10000);
	//vtkSmartPointer<vtkDataSet> selectDS = selector.select(pointDS);
	//selectDS->PrintSelf(cout,vtkIndent(0));
	//vtkSmartPointer<vtkDataSet> opDS = contourOp.process(selectDS);	
	

	double center[] = {0.0,0.0,0.0};
	int circleRes = 20;
	int sweepRes = 3;
	double radius = 0.25;
#if 0
	vtkSmartPointer<vtkRegularPolygonSource> circleGenerator = vtkSmartPointer<vtkRegularPolygonSource>::New();
	circleGenerator->SetGeneratePolyline(1);
	circleGenerator->SetGeneratePolygon(0);
	circleGenerator->SetNumberOfSides(circleRes);
	circleGenerator->SetRadius(radius);
	circleGenerator->SetCenter(center);
	circleGenerator->Update();
	
	vtkSmartPointer<vtkPolyData> circleDS = circleGenerator->GetOutput();
#elif 1
	vtkSmartPointer<vtkPolyData> circleDS = generateCircle(circleRes,radius,center,sweepRes);
#else

	//The following is used just to test the 2d curvature. We have a 2d flash dataset with 1 block of cells.
	vtkDataSet* flashDS = vtkOverlappingAMR::SafeDownCast(flashReader.readFile(flashFilename.c_str()))->GetDataSet(0,0);
	flashDS->GetCellData()->SetScalars(flashDS->GetCellData()->GetArray(0));
	cout<<"Flash DS-----"<<endl;
	flashDS->PrintSelf(cout,vtkIndent(5));
	cout<<"End Flash DS-----"<<endl;
	flashDS->GetScalarRange(range);
	cout<<"Flash DS min = "<<range[0]<<" max = "<<range[1]<<endl;
	isoValues[0]=0.2;
	contourOp.setIsoLevels(1, isoValues);
	contourOp.setArrayName(vtkDataObject::FIELD_ASSOCIATION_POINTS,"vari");
	//contourOp.setArrayName(vtkDataObject::FIELD_ASSOCIATION_CELLS,"vari");

	c2p->SetInputData(flashDS);
	c2p->Update();
	vtkSmartPointer<vtkDataSet> flashPointsDS =  c2p->GetOutput();

	vtkSmartPointer<vtkPolyData> circleDS = vtkPolyData::SafeDownCast(contourOp.process(flashPointsDS));
#endif
	//vtkSmartPointer<vtkPolyData> circleDS = vtkPolyData::SafeDownCast(contourOp.process(flashPointsDS));

	cout<<"Circle DS-----"<<endl;
	circleDS->PrintSelf(cout,vtkIndent(5));
	cout<<"End Circle DS-----"<<endl;
	circleDS->GetScalarRange(range);
	cout<<"Circle DS min = "<<range[0]<<" max = "<<range[1]<<endl;

	vtkSmartPointer<vtkDataSet> circleCurveDS = curveOp.process(circleDS);
	//vtkSmartPointer<vtkDataSet> circleCurveDS = circleDS;
	cout<<"Circle Curve DS-----"<<endl;
	circleCurveDS->PrintSelf(cout,vtkIndent(5));
	cout<<"End Circle Curve DS-----"<<endl;
	//tpd->SetInputData(curveDS);
	//tpd->Update();
	//vtkSmartPointer<vtkDataSet> cartisianCurvatureDS=tpd->GetOutput();
	double circCurveRange[2];
	circleCurveDS->GetScalarRange(circCurveRange);
	cout<<"Circle Curvature min = "<<circCurveRange[0]<<" max = "<<circCurveRange[1]<<endl;

	vtkSmartPointer<vtkDataSet> opDS = contourOp.process(cartisianDS);	
	//vtkSmartPointer<vtkDataSet> opDS = contourOp.process(polydataDS);	
	//cout<<"Contour DS-----"<<endl;
	//opDS->PrintSelf(cout,vtkIndent(5));
	//cout<<"End Contour DS-----"<<endl;
	opDS->GetScalarRange(range);
	cout<<"ContourDS min = "<<range[0]<<" max = "<<range[1]<<endl;
	//opDS->PrintSelf(cout,vtkIndent(10));

	#if 0
	double center[] = {0.0,0.0,0.0};
	//int spherePhiRes = 100;//30;
	//int sphereThetaRes = 50;//15;
	int spherePhiRes = 500;//30;
	int sphereThetaRes = 250;//15;
	//double radius = 7.46567e+08;
	double radius = 1.0;
	cout<<"radius = "<<radius<<endl;
	vtkSmartPointer<vtkSphereSource> sphereGenerator = vtkSmartPointer<vtkSphereSource>::New();
	sphereGenerator->SetRadius(radius);
	sphereGenerator->SetCenter(center);
	sphereGenerator->SetPhiResolution(spherePhiRes);
	sphereGenerator->SetThetaResolution(sphereThetaRes);
	sphereGenerator->Update();
	vtkSmartPointer<vtkPolyData> sphereDS = sphereGenerator->GetOutput();
#endif

	//double range[2];
	//opDS->GetScalarRange(range);
	//cout<<"min = "<<range[0]<<" max = "<<range[1]<<endl;
	//renderer.renderPseudoColor(opDS,min,max);
	//vtkSmartPointer<vtkCurvatures> curveFilter = vtkSmartPointer<vtkCurvatures>::New();
	//curveFilter->SetInputData(opDS);
	//curveFilter->SetCurvatureTypeToMean();
	//curveFilter->Update();
	//vtkSmartPointer<vtkDataSet> curveDS = curveFilter->GetOutput();
	//curveDS->PrintSelf(cout,vtkIndent(0));
	//vtkSmartPointer<vtkFillHolesFilter> fillHoles = vtkSmartPointer<vtkFillHolesFilter>::New();
	//fillHoles->SetHoleSize(7.1e8);
	//fillHoles->SetInputData(opDS);
	//fillHoles->Update();
	//vtkSmartPointer<vtkDataSet> fillHolesDS = fillHoles->GetOutput();
	//vtkSmartPointer<vtkCurvatures> curveFilter = vtkSmartPointer<vtkCurvatures>::New();
	//curveFilter->SetInputData(circleDS);
	//curveFilter->SetInputData(opDS);
	//curveFilter->SetInputData(sphereDS);
	//curveFilter->SetInputData(fillHolesDS);
	//curveFilter->SetCurvatureTypeToMean();
	//curveFilter->SetCurvatureTypeToMaximum();
	//curveFilter->Update();
	vtkSmartPointer<vtkDataSet> curveDS = curveOp.process(opDS);
	cout<<"Curve DS-----"<<endl;
	curveDS->PrintSelf(cout,vtkIndent(5));
	cout<<"End Curve DS-----"<<endl;
	//tpd->SetInputData(curveDS);
	//tpd->Update();
	//vtkSmartPointer<vtkDataSet> cartisianCurvatureDS=tpd->GetOutput();
	curveDS->GetScalarRange(range);
	cout<<"Curvature min = "<<range[0]<<" max = "<<range[1]<<endl;
	//range[0]=mean-1.e-9;
	//range[1]=mean+1.e-9;
	//range[0]=-1.5e-6;
	//range[1]=1.5e-6;
	//range[0]=0.0;
	//range[1]=1.5e-9;
	range[0]=1.e-9;
	range[1]=1.5e-9;
	//range[1]=1.7e-9;
	cout<<"viewing min = "<<range[0]<<" max = "<<range[1]<<endl;
	//surfaceAreaOp.process(cartisianDS);
	//surfaceAreaOp.process(opDS);
	//cout<<"Estimated Surface Area cartisian coordinates is "<<surfaceAreaOp.getSurfaceArea()<<" Assuming radius = 7.46567e+08 then exact should be "<<4.*vtkMath::Pi()*7.46567e+08*7.46567e+08<<endl;
	//renderer.renderPseudoColor(opDS,range[0]-1.,range[1]+1.);
	//renderer.renderPseudoColor(cartisianCurvatureDS,range[0],range[1]);
	//renderer.renderPseudoColor(curveDS,range[0]+1e-9,range[1]-1e-9);
	//renderer.renderPseudoColor(curveDS,range[1]-2.e-9,range[1]);
	//renderer.renderPseudoColor(curveDS,range[0],range[1]);
	//renderer.renderPseudoColor(curveDS,range[0],range[1]);
	renderer.renderPseudoColor(circleCurveDS,circCurveRange[0]-.1,circCurveRange[1]+.1);
	//renderer.renderPseudoColor(circleCurveDS,.99,1.01);
	//renderer.renderPoints(circleCurveDS,circCurveRange[0],circCurveRange[1]);
	//renderer.renderPseudoColor(selectDS);
	//renderer2.renderPseudoColor(ds);
	renderer.interact();
	//renderer.renderPseudoColor(ds);
	return 0;
}
