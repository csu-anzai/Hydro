#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Operators/MultiplyOperator.h"
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridGeometryFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkSphericalTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVersion.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <float.h>

int main(int argc, char** argv)
{
	
	double scale = 2.0;
	if(argc>1)
		scale= atof(argv[1]);
	printf("vtkVersion  = %s\n",vtkVersion::GetVTKVersion());
	//std::string filename = "/home/asy10/VTKData/Data/RectGrid2.vtk";
	//std::string filename = "/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0050.vtk";
	std::string filename = "/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0200.vtk";
	
	VTKReader reader;
	NullSelector selector;
	Renderer renderer(800,600);
	Renderer renderer2(800,600);
	MultiplyOperator multiOp;
	multiOp.setScale(scale);
	vtkDataSet* ds = reader.readFile(filename.c_str());
	vtkSmartPointer<vtkSphericalTransform> transformer = vtkSmartPointer<vtkSphericalTransform>::New();
	vtkSmartPointer<vtkTransformPolyDataFilter> tpd = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	tpd->SetTransform(transformer);
	ds->PrintSelf(cout,vtkIndent(0));
	vtkRectilinearGrid* rg = vtkRectilinearGrid::SafeDownCast(ds);
	vtkSmartPointer<vtkRectilinearGridGeometryFilter> rggf = vtkSmartPointer<vtkRectilinearGridGeometryFilter>::New();
	vtkSmartPointer<vtkCellDataToPointData> c2p = vtkSmartPointer<vtkCellDataToPointData>::New();
	c2p->SetInputData(rg);
	c2p->Update();
	
	rggf->SetExtent(rg->GetExtent());
	rggf->SetInputData(c2p->GetOutput());
	rggf->Update();

	tpd->SetInputData(rggf->GetOutput());
	tpd->Update();
	vtkDataSet* newDS = reinterpret_cast<vtkDataSet*>(tpd->GetOutput());
	newDS->PrintSelf(cout,vtkIndent(0));
	//Debugging
	vtkDataArray* da = rg->GetCellData()->GetScalars();
	double total=0.0;
	unsigned long num = 0;
	if(da)
	{
		for(vtkIdType j=0;j<da->GetNumberOfTuples(); j++)
		{
			for(vtkIdType k = 0; k<da->GetNumberOfComponents(); k++)
			{
				double val = da->GetComponent(j,k);
				total+=val;
				num+=1;
			//	if(val>1000000000000000.0||val<0.0)
		//			printf("cell data ---- val = %f\n",val);
			}
		}
	}
	printf("cell data ---- Total = %f, num = %ld, mean = %f\n",total,num,total/num);
	//Debugging
	da = c2p->GetOutput()->GetPointData()->GetScalars();
	total=0.0;
	num = 0;
	if(da)
	{
		for(vtkIdType j=0;j<da->GetNumberOfTuples(); j++)
		{
			for(vtkIdType k = 0; k<da->GetNumberOfComponents(); k++)
			{
				double val = da->GetComponent(j,k);
				total+=val;
				num+=1;
				//if(val>1000000000000000.0||val<0.0)
				//	printf("val = %f\n",val);
			}
		}
	}
	printf("c2p ---- Total = %f, num = %ld, mean = %f\n",total,num,total/num);

	//Debugging
	da = rggf->GetOutput()->GetPointData()->GetScalars();
	total=0.0;
	num = 0;
	if(da)
	{
		for(vtkIdType j=0;j<da->GetNumberOfTuples(); j++)
		{
			for(vtkIdType k = 0; k<da->GetNumberOfComponents(); k++)
			{
				double val = da->GetComponent(j,k);
				total+=val;
				num+=1;
				//if(val>1000000000000000.0||val<0.0)
				//	printf("val = %f\n",val);
			}
		}
	}
	printf("PolyData ---- Total = %f, num = %ld, mean = %f\n",total,num,total/num);
	//Debugging
	da = newDS->GetPointData()->GetScalars();
	total=0.0;
	num = 0;
	double min=DBL_MAX;
	double max=-DBL_MAX;
	if(da)
	{
		for(vtkIdType j=0;j<da->GetNumberOfTuples(); j++)
		{
			for(vtkIdType k = 0; k<da->GetNumberOfComponents(); k++)
			{
				double val = da->GetComponent(j,k);
				total+=val;
				num+=1;
				if(val<min)
					min=val;
				if(val>max)
					max=val;
			}
		}
	}
	printf("final data ---- Total = %f, num = %ld, mean = %f\n",total,num,total/num);
	//Debugging

	
	//vtkPoints* pts=NULL;
	//ds->GetPoints(pts);
	/*double* x=NULL;
	double y[3];
	for(vtkIdType i = 0; i<ds->GetNumberOfPoints(); i++)
	{
		x=ds->GetPoint(i);
		transformer->TransformPoint(x,y);
		//cout<<"x("<<x[0]<<","<<x[1]<<","<<x[2]<<")"<<endl;
		cout<<"y("<<y[0]<<","<<y[1]<<","<<y[2]<<")"<<endl;
		*(x)=y[0];
		*(x+1)=y[1];
		*(x+2)=y[2];
		x=ds->GetPoint(i);
		cout<<"x("<<x[0]<<","<<x[1]<<","<<x[2]<<")"<<endl;
	}
	ds->PrintSelf(cout,vtkIndent(0));*/
	//tpd->SetInputData(ds);
	//tpd->SetTransform(transformer);
	//vtkDataSet* newDS = tpd->GetOutput();
	//newDS->PrintSelf(cout,vtkIndent(0));
	//vtkSmartPointer<vtkDataSet> selectDS = selector.select(ds);
	vtkSmartPointer<vtkDataSet> selectDS = selector.select(newDS);
	vtkSmartPointer<vtkDataSet> opDS = multiOp.process(selectDS);	
	renderer.renderPseudoColor(opDS,min,max);
	//renderer.renderPseudoColor(ds);
	//renderer.renderPseudoColor(selectDS);
	//renderer2.renderPseudoColor(ds);
	renderer2.renderPseudoColor(newDS,min,max);
	renderer.interact();
	renderer2.interact();
	//renderer.renderPseudoColor(ds);
	/*vtkSmartPointer<vtkRectilinearGridReader> reader =
		vtkSmartPointer<vtkRectilinearGridReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkSmartPointer<vtkDataSetMapper> mapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	//mapper->SetInputData(reinterpret_cast<vtkDataSet*>(reader->GetOutput()));
	mapper->SetInputConnection(reader->GetOutputPort());
	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	//actor->GetProperty()->SetRepresentationToWireframe();
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderer->AddActor(actor);
	renderer->SetBackground(.3,.6,.3);
	renderWindow->Render();
	renderWindowInteractor->Start();*/
	return 0;
}
