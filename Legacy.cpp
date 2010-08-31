#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkKdTree.h>
#include <vtkIdList.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>

bool IsDone(std::vector<unsigned int> &used);
int GetUnusedIndex(std::vector<unsigned int> &used);
bool IsUsed(std::vector<unsigned int> &used, unsigned int index);

int main(int argc, char* argv[])
{
  vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  //decide on a reasonable e-sphere size
  double bounds[6];
  polyData->GetBounds(bounds);

  double delx = bounds[1] - bounds[0];
  double dely = bounds[3] - bounds[2];
  double delz = bounds[5] - bounds[4];
  double minDim = std::min(delx, std::min(dely,delz));

  double eRadius = minDim / 100.;

  vtkPoints* points = reader->GetOutput()->GetPoints();

  //create a vector to keep track of the points that are already assigned to a superpoint
  std::vector<unsigned int> used(points->GetNumberOfPoints(), 0);

  //Create a kd tree
  vtkSmartPointer<vtkKdTree> kdTree =
    vtkSmartPointer<vtkKdTree>::New();
  kdTree->BuildLocatorFromPoints(points);

  while(!IsDone(used))
    {
    int index = GetUnusedIndex(used);

    vtkIdType k = 1; // Start by finding the 2 closest points to the randomly selected point (the first closest is the point itself!)
    vtkSmartPointer<vtkIdList> result =
      vtkSmartPointer<vtkIdList>::New();
    do
      {
      k++;
      double testPoint[3];
      points->GetPoint(index, testPoint);



      kdTree->FindClosestNPoints(k, testPoint, result);

      }while(!IsUsed(used, result->GetId(k-1)));
    }

  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName("output.vtp");
//  writer->SetInputConnection(glyphFilter->GetOutputPort());
  writer->Write();

  return EXIT_SUCCESS;
}

bool IsDone(std::vector<unsigned int> &used)
{
  int vecSum = std::accumulate(used.begin(), used.begin() + used.size(), 0);
  if(vecSum == used.size())
    {
    return true;
    }

  return false;
}

int GetUnusedIndex(std::vector<unsigned int> &used)
{
  if(IsDone(used))
    {
    //all the indices are already used!
    return -1;
    }

  std::vector<unsigned int> indices;

  for(unsigned int i = 0; i < used.size(); i++)
    {
    indices.push_back(i);
    }

  std::random_shuffle(indices.begin(), indices.end());

  for(unsigned int i = 0; i < used.size(); i++)
    {
    if(!IsUsed(used, indices[i]))
      {
      return indices[i];
      }
    }

  return -1; //should never get here
}

bool IsUsed(std::vector<unsigned int> &used, unsigned int index)
{
  if(used[index] != 0)
    {
    return true;
    }

  return false;
}
