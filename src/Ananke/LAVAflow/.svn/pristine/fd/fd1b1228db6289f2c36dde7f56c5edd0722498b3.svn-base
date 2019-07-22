#ifndef RENDERER_H
#define RENDERER_H
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkLogLookupTable.h>


class Renderer
{
vtkSmartPointer<vtkActor> actor;
vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> rWin;
vtkSmartPointer<vtkRenderWindowInteractor> rInteract;
vtkSmartPointer<vtkDataSetMapper> dsMapper;
vtkSmartPointer<vtkPolyDataMapper> polyMapper;
vtkSmartPointer<vtkLogLookupTable> llTable;

public:
	Renderer(unsigned int width, unsigned int height);
	void interact();
	void renderPseudoColor(vtkSmartPointer<vtkDataSet> ds, double min = 0.0, double max=1.0);
	void renderPoints(vtkSmartPointer<vtkDataSet> ds, double min = 0.0, double max=1.0, double pointSize = 10.0);
	void renderContours(vtkSmartPointer<vtkDataSet> ds);
};

#endif
