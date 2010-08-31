#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkIdList.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>

#include <vtkFullyConnectedGraphFilter.h>
#include <vtkGraphBFSIterator.h>

#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>

bool IsDone(std::vector<unsigned int> &used);
int GetUnusedIndex(std::vector<unsigned int> &used);
bool IsUsed(std::vector<unsigned int> &used, unsigned int index);

int main(int argc, char* argv[])
{
  /*
  vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  */

  vtkSmartPointer<vtkSphereSource> reader =
      vtkSmartPointer<vtkSphereSource>::New();
  reader->Update();

  //decide on a reasonable e-sphere size
  double bounds[6];
  reader->GetOutput()->GetBounds(bounds);

  double delx = bounds[1] - bounds[0];
  double dely = bounds[3] - bounds[2];
  double delz = bounds[5] - bounds[4];
  double minDim = std::min(delx, std::min(dely,delz));

  double eRadius = minDim / 100.;

  vtkPoints* points = reader->GetOutput()->GetPoints();

  //create a vector to keep track of the points that are already assigned to a superpoint
  std::vector<unsigned int> used(points->GetNumberOfPoints(), 0);

  //Create a kd tree
  vtkSmartPointer<vtkKdTreePointLocator> kdTree =
    vtkSmartPointer<vtkKdTreePointLocator>::New();
  kdTree->SetDataSet(reader->GetOutput());
  kdTree->BuildLocator();

  int superPointID = 0;

  while(!IsDone(used))
    {
    //find a point that we have not yet assigned to a region
    int index = GetUnusedIndex(used);
    double queryPoint[3];
    points->GetPoint(index, queryPoint);

    //find all the points around the query point
    vtkSmartPointer<vtkIdList> result =
      vtkSmartPointer<vtkIdList>::New();
    kdTree->FindPointsWithinRadius(eRadius, queryPoint, result);

    //extract the points around the query point
    vtkSmartPointer<vtkIdTypeArray> ids =
      vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);

    for(unsigned int i = 0; i < result->GetNumberOfIds(); i++)
      {
      ids->InsertNextValue(result->GetId(i));
      }

    vtkSmartPointer<vtkSelectionNode> selectionNode =
        vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::POINT);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);

    vtkSmartPointer<vtkSelection> selection =
        vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);

    vtkSmartPointer<vtkExtractSelection> extractSelection =
        vtkSmartPointer<vtkExtractSelection>::New();

    extractSelection->SetInput(0, reader->GetOutput());
    extractSelection->SetInput(1, selection);
    extractSelection->Update();

    //create a fully connected graph on the points in a e-neighborhood of the query point
    vtkSmartPointer<vtkFullyConnectedGraphFilter> fullyConnectedFilter =
      vtkSmartPointer<vtkFullyConnectedGraphFilter>::New();
    fullyConnectedFilter->SetInputConnection(extractSelection->GetOutputPort());
    fullyConnectedFilter->Update();

    vtkSmartPointer<vtkGraphBFSIterator> graphIterator =
      vtkSmartPointer<vtkGraphBFSIterator>::New();
    graphIterator->SetGraph(fullyConnectedFilter->GetOutput());
    graphIterator->SetStartVertex(index);
    graphIterator->Initialize();

    while(iterator->HasNext())
      {
      vtkIdType nextVertex = iterator->Next();
      }

    superPointID++;

    } //end while(!IsDone(used))

  /*
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName("output.vtp");
//  writer->SetInputConnection(glyphFilter->GetOutputPort());
  writer->Write();
  */

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
