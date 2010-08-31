#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include "vtkSuperPoints.h"


int main(int argc, char* argv[])
{
  if(argc != 4)
    {
    std::cerr << "Required arguments: InputFile.vtp OutputSuperPointLabels.vtp OuputSuperPointCenters.vtp" << std::endl;
    return EXIT_FAILURE;
    }
    
  vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  vtkSmartPointer<vtkSuperPoints> superPoints =
    vtkSmartPointer<vtkSuperPoints>::New();
  superPoints->SetInputConnection(reader->GetOutputPort());
  superPoints->SetAutomaticRadiusRatio(5);
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
