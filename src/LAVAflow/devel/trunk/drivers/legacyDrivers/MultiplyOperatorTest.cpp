#include "../Selectors/NullSelector.h"
#include "../Rendering/Renderer.h"
#include "../Readers/VTKReader.h"
#include "../Operators/MultiplyOperator.h"

int main(int argc, char** argv)
{
	
	double scale = 2.0;
	if(argc>1)
		scale= atof(argv[1])
	
	VTKReader reader;
	NullSelector selector;
	Renderer renderer(800,600);
	MultiplyOperator multiOp;
	vtkDataSet* ds = reader.read("/data3/nsf/database/models/astrophysics/supernovae/core-collapse/explosion/hydrodynamics/3/hotb/b18d3a3s123SG2r3/b18d3a3s123SG2r3.plt.0150.vtk")
	
	vtkSmartPointer<vtkDataSet> selectDS = selector(ds);
	vtkSmartPointer<vtkDataSet> opDS = multiOp(selectDS);	
	renderer->renderPsuedoColor(opDS);
	renderer->interact();
	return 0;
}
