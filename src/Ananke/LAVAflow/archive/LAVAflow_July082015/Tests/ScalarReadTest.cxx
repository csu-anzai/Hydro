#include "../libsrc/Selectors/NullSelector.h"
#include "../libsrc/Operators/SurfaceAreaOperator.h"
#include "../libsrc/Operators/shellAveraging/ShellAverageOperator.h"
#include "../libsrc/Operators/shellAveraging/planeShells/ShellAveragePlaneOperator.h"
#include "../libsrc/Operators/surfaceAveraging/planarSurface/SurfaceAveragePlaneOperator.h"
#include "../libsrc/Rendering/Renderer.h"
#include "../libsrc/Readers/VTKReader.h"
#include "../libsrc/Mesh/Mesh.h"
#include "../libsrc/includes/LAVAconstants.h"

#include "../libsrc/Math/LAVAMath.h"

#include <string>
#include <vtkRectilinearGrid.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLRectilinearGridWriter.h>

#include <fstream>
#include <algorithm>
#include <complex>


int main(void)
{
	int nDim = 2;
	int currentCS = CS_CART;
	Mesh data(nDim, currentCS, "../Tests/files/averaging/structpts_constant_2d.vtk");

	data.getDataSet()->PrintSelf(std::cout,vtkIndent(0));


	std::vector<std::string> scalarNames;
	scalarNames.push_back("vari2");
	scalarNames.push_back("vari3");
	Mesh dataSelect(nDim, currentCS, "../Tests/files/averaging/structpts_constant_2d.vtk",scalarNames);

	dataSelect.getDataSet()->PrintSelf(std::cout,vtkIndent(0));
}