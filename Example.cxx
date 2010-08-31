#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataNormals.h>

#include "vtkSuperPoints.h"

int main(int, char*[])
{
  vtkSmartPointer<vtkSphereSource> reader =
      vtkSmartPointer<vtkSphereSource>::New();
  reader->SetThetaResolution(20);
  reader->SetPhiResolution(20);
  reader->Update();
  
  vtkSmartPointer<vtkPolyDataNormals> normalGenerator =
    vtkSmartPointer<vtkPolyDataNormals>::New();
  normalGenerator->SetInputConnection(reader->GetOutputPort());
  normalGenerator->Update();

  vtkSmartPointer<vtkSuperPoints> superPoints =
    vtkSmartPointer<vtkSuperPoints>::New();
  superPoints->SetInputConnection(normalGenerator->GetOutputPort());
  superPoints->SetAutomaticRadiusRatio(15);
  superPoints->Update();

  {
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(argv[2]);
  writer->SetInputConnection(superPoints->GetOutputPort(0));
  writer->Write();
  }

  {
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(argv[3]);
  writer->SetInputConnection(superPoints->GetOutputPort(1));
  writer->Write();
  }
  return EXIT_SUCCESS;
}
