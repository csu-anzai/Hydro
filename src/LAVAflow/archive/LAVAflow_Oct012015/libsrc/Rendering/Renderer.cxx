#include "Renderer.h"
#include <vtkProperty.h>

Renderer::Renderer(unsigned int width, unsigned int height)
{
	actor = vtkSmartPointer<vtkActor>::New();
	renderer = vtkSmartPointer<vtkRenderer>::New();
	rWin =vtkSmartPointer<vtkRenderWindow>::New();
	rInteract = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	rWin->SetSize(width,height);
	rWin->AddRenderer(renderer);
	rInteract->SetRenderWindow(rWin);
	dsMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	polyMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	llTable = vtkSmartPointer<vtkLogLookupTable>::New();
	//polyMapper->SetLookupTable(llTable);
	//polyMapper->GetLookupTable()->PrintSelf(cout,vtkIndent(0));
	renderer->AddActor(actor);
}
void Renderer::renderPseudoColor(vtkSmartPointer<vtkDataSet> ds, double min, double max)
{
	//dsMapper->SetInputData(ds);
	//dsMapper->SetInput(ds);
	//dsMapper->Update();
	polyMapper->SetScalarRange(min,max);
	polyMapper->SetInputData(vtkPolyData::SafeDownCast(ds));
	polyMapper->Update();
	actor->SetMapper(polyMapper);
	rWin->Render();
}
void Renderer::renderPoints(vtkSmartPointer<vtkDataSet> ds, double min, double max, double pointSize)
{
	//dsMapper->SetInputData(ds);
	//dsMapper->SetInput(ds);
	//dsMapper->Update();
	polyMapper->SetScalarRange(min,max);
	polyMapper->SetInputData(vtkPolyData::SafeDownCast(ds));
	polyMapper->Update();
	actor->SetMapper(polyMapper);
	actor->GetProperty()->SetPointSize(pointSize);
	int rep = actor->GetProperty()->GetRepresentation();
	cout<<"Representation type was "<<actor->GetProperty()->GetRepresentationAsString()<<endl;
	actor->GetProperty()->SetRepresentationToPoints();
	cout<<"Representation type is "<<actor->GetProperty()->GetRepresentationAsString()<<endl;
	rWin->Render();
	actor->GetProperty()->SetRepresentation(rep);
}
void Renderer::renderContours(vtkSmartPointer<vtkDataSet> ds)
{
	polyMapper->SetInputData(vtkPolyData::SafeDownCast(ds));
	//dsMapper->SetInput(ds);
	polyMapper->Update();
	actor->SetMapper(polyMapper);
	rWin->Render();
}
void Renderer::interact()
{
	rInteract->Start();
}
