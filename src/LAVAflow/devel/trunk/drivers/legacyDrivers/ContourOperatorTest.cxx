#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Operators/ContourOperator.h"

#include <vtkSphericalTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkVersion.h>

int main(int argc, char** argv)
{
	printf("vtkVersion  = %s\n",vtkVersion::GetVTKVersion());
	//std::string filename = "/home/asy10/VTKData/Data/RectGrid2.vtk";
	//std::string filename = "/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0050.vtk";
	std::string filename = "/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0200.vtk";
	
	VTKReader reader;
	NullSelector selector;
	Renderer renderer(800,600);
	ContourOperator contourOp;
        double isoValues[] = {0.5,500.,5000., 50000.};
	contourOp.setIsoLevels(4, isoValues);
	contourOp.setArrayName(vtkDataObject::FIELD_ASSOCIATION_POINTS,"dens");

	vtkDataSet* ds = reader.readFile(filename.c_str());
	//ds->PrintSelf(cout,vtkIndent(0));
	//It seems contouring algorithms require point data.
	vtkSmartPointer<vtkCellDataToPointData> c2p = vtkSmartPointer<vtkCellDataToPointData>::New();
	c2p->SetInputData(ds);
	c2p->Update();
	vtkSmartPointer<vtkDataSet> pointDS=c2p->GetOutput();
	
	vtkSmartPointer<vtkDataSet> selectDS = selector.select(pointDS);
	selectDS->PrintSelf(cout,vtkIndent(0));
	vtkSmartPointer<vtkDataSet> opDS = contourOp.process(selectDS);	
	opDS->PrintSelf(cout,vtkIndent(10));
	double range[2];
	opDS->GetScalarRange(range);
	cout<<"min = "<<range[0]<<" max = "<<range[1]<<endl;
	//renderer.renderPseudoColor(opDS,min,max);
	vtkSmartPointer<vtkSphericalTransform> transformer = vtkSmartPointer<vtkSphericalTransform>::New();
	vtkSmartPointer<vtkTransformPolyDataFilter> tpd = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	tpd->SetTransform(transformer);
	tpd->SetInputData(opDS);
	tpd->Update();
	vtkSmartPointer<vtkDataSet> cartisianDS=tpd->GetOutput();
	//renderer.renderPseudoColor(opDS,range[0]-1.,range[1]+1.);
	//renderer.renderPseudoColor(cartisianDS,range[0]-1.,range[1]+1.);
	renderer.renderPseudoColor(cartisianDS);
	//renderer.renderPseudoColor(selectDS);
	//renderer2.renderPseudoColor(ds);
	renderer.interact();
	//renderer.renderPseudoColor(ds);
	return 0;
}
