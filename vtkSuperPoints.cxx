#include "vtkSuperPoints.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"

#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkIdList.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>

#include <vtkFullyConnectedGraphFilter.h>
#include <vtkGraphBFSIterator.h>
#include <vtkUnstructuredGridToGraph.h>
#include <vtkNearestNeighborGraphFilter.h>
#include <vtkConnectGraph.h>
#include <vtkGraphVertexDataConditionalIterator.h>

#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>

#define UNASSIGNED -1

vtkStandardNewMacro(vtkSuperPoints);

vtkSuperPoints::vtkSuperPoints()
{
  this->SetNumberOfOutputPorts(2);

  this->NNRadius = 1.0;
  this->AutomaticRadius = 1;
  this->AutomaticRadiusRatio = 10;
}

int vtkSuperPoints::RequestData(vtkInformation *vtkNotUsed(request),
                                             vtkInformationVector **inputVector,
                                             vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

  // get the input and ouptut
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

      
  vtkDataArray* normalsInput = input->GetPointData()->GetNormals();
  if(!normalsInput)
    {
    std::cerr << " no input normals" << std::endl;
    exit(-1);
    }
  //vtkDataArray* normalsInput = input->GetPointData()->GetArray("Normals");
  
  vtkPolyData *outputLabeled = vtkPolyData::SafeDownCast(
    outInfo0->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *outputCenters = vtkPolyData::SafeDownCast(
    outInfo1->Get(vtkDataObject::DATA_OBJECT()));

  if(this->AutomaticRadius)
    {
    //decide on a reasonable RBNN radius
    double bounds[6];
    input->GetBounds(bounds);

    double delx = bounds[1] - bounds[0];
    double dely = bounds[3] - bounds[2];
    double delz = bounds[5] - bounds[4];
    std::cout << "delx: " << delx << " dely: " << dely << " delz: " << delz << std::endl;
    double minDim = std::min(delx, std::min(dely,delz));

    this->NNRadius = minDim / static_cast<double>(this->AutomaticRadiusRatio);
    std::cout << "Automatic radius: " << this->NNRadius << std::endl;
    }

  //create a vector to keep track of the points that are already assigned to a superpoint
  std::vector<int> superPointLabels(input->GetNumberOfPoints(), UNASSIGNED);

 
  //Create a kd tree
  vtkSmartPointer<vtkKdTreePointLocator> kdTree =
    vtkSmartPointer<vtkKdTreePointLocator>::New();
  kdTree->SetDataSet(input);
  kdTree->BuildLocator();

  int superPointID = 0;
  vtkSmartPointer<vtkPoints> superPointCenters =
    vtkSmartPointer<vtkPoints>::New();

  //for(vtkIdType pointID = 0; pointID < input->GetNumberOfPoints(); ++pointID)
  for(vtkIdType pointID = 0; pointID < 4000; ++pointID)
    {
    if(pointID % 1000 == 0)
      {
      std::cout << "point " << pointID << " out of " << input->GetNumberOfPoints() << std::endl;
      }
    //std::cout << "point " << pointID << std::endl;
    if(pointID == 1313)
    {
      std::cout << "here";
    }
    
    if(superPointLabels[pointID] != UNASSIGNED)
      {
      continue;
      }

    double queryPoint[3];
    input->GetPoint(pointID, queryPoint);
    superPointCenters->InsertNextPoint(queryPoint);

    //find all the points around the query point
    vtkSmartPointer<vtkIdList> neighbors =
      vtkSmartPointer<vtkIdList>::New();
    kdTree->FindPointsWithinRadius(this->NNRadius, pointID, neighbors, true);

    //std::cout << "There are " << neighbors->GetNumberOfIds() << " neighbors." << std::endl;
    
    //extract the points around the query point
    vtkSmartPointer<vtkIdTypeArray> ids =
      vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);

    for(vtkIdType i = 0; i < neighbors->GetNumberOfIds(); i++)
      {
      if(superPointLabels[neighbors->GetId(i)] == UNASSIGNED)
        {
        ids->InsertNextValue(neighbors->GetId(i));
        }
      }

    if(ids->GetNumberOfTuples() == 0)
      {
      continue; //all points in this neighborhood already belong to a superpixel
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
    extractSelection->SetInput(0, input);
    extractSelection->SetInput(1, selection);
    extractSelection->Update();

    vtkDataArray* normalsAfter = vtkDataArray::SafeDownCast(vtkUnstructuredGrid::SafeDownCast(extractSelection->GetOutput())->GetPointData()->GetNormals());
    
    vtkIdTypeArray* originalIds =
      vtkIdTypeArray::SafeDownCast(vtkUnstructuredGrid::SafeDownCast(extractSelection->GetOutput())->GetPointData()->GetArray("vtkOriginalPointIds"));

    //std::cout << "There are " << vtkUnstructuredGrid::SafeDownCast(extractSelection->GetOutput())->GetNumberOfPoints() << " extracted points." << std::endl;
    
    
    // Create a fully connected graph on the points in a e-neighborhood of the query point
    /*
    vtkSmartPointer<vtkUnstructuredGridToGraph> unstructuredGridToGraph =
      vtkSmartPointer<vtkUnstructuredGridToGraph>::New();
    unstructuredGridToGraph->SetInputConnection(extractSelection->GetOutputPort());
    unstructuredGridToGraph->Update();

    vtkSmartPointer<vtkFullyConnectedGraphFilter> fullyConnectedFilter =
      vtkSmartPointer<vtkFullyConnectedGraphFilter>::New();
    fullyConnectedFilter->SetInputConnection(unstructuredGridToGraph->GetOutputPort());
    fullyConnectedFilter->Update();
    */
    
    
    // Create a NN graph on the extracted points
    vtkSmartPointer<vtkNearestNeighborGraphFilter> nnFilter =
      vtkSmartPointer<vtkNearestNeighborGraphFilter>::New();
    if(vtkUnstructuredGrid::SafeDownCast(extractSelection->GetOutput())->GetNumberOfPoints() < nnFilter->GetKNeighbors() + 2)
      {
      continue; //there will be some unlabeled points
      }
      
    nnFilter->SetInputConnection(extractSelection->GetOutputPort());
    nnFilter->Update();
    
    // Connect the NN graph
    vtkSmartPointer<vtkConnectGraph> connectFilter =
      vtkSmartPointer<vtkConnectGraph>::New();
    connectFilter->SetInputConnection(nnFilter->GetOutputPort());
    connectFilter->Update();
    
    
    vtkIdType query = -1;
    //map from original id to current id in selected graph
    for(vtkIdType j = 0; j < originalIds->GetNumberOfTuples(); j++)
      {
      if(originalIds->GetValue(j) == pointID)
        {
        query = j;
	//std::cout << "Original Id " << pointID << " is now ID " << j << " in extracted points" << std::endl;
        break;
        }
      }
    if(query == -1)
      {
      std::cerr << "Original Id " << pointID << " was not found in extracted points - this should never happen!" << std::endl;
      exit(-1);
      }

    vtkSmartPointer<vtkGraphVertexDataConditionalIterator> graphIterator =
      vtkSmartPointer<vtkGraphVertexDataConditionalIterator>::New();
    graphIterator->SetGraph(connectFilter->GetOutput());
    //graphIterator->SetGraph(fullyConnectedFilter->GetOutput());
    graphIterator->SetStartVertex(query);
    graphIterator->SetConditionTypeToNormalAngle();
    graphIterator->SetComparisonDirectionToLessThan();
    graphIterator->SetAcceptableValue(10);
    graphIterator->Initialize();
    
    superPointLabels[pointID] = superPointID;

    while(graphIterator->HasNext())
      {
      vtkIdType nextVertex = graphIterator->Next();
      superPointLabels[originalIds->GetValue(nextVertex)] = superPointID;
      }

    superPointID++;

    } //end for

  //color points
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(input->GetPoints());

  vtkSmartPointer<vtkIntArray> labelArray =
    vtkSmartPointer<vtkIntArray>::New();
  labelArray->SetNumberOfComponents(1);
  labelArray->SetName("SuperPointLabels");

  for(unsigned int i = 0; i < superPointLabels.size(); i++)
    {
    labelArray->InsertNextValue(superPointLabels[i]);
    }

  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputConnection(polydata->GetProducerPort());
  glyphFilter->Update();

  outputLabeled->ShallowCopy(glyphFilter->GetOutput());

  outputLabeled->GetPointData()->SetScalars(labelArray);

  vtkSmartPointer<vtkPolyData> centersPolyData =
    vtkSmartPointer<vtkPolyData>::New();
  centersPolyData->SetPoints(superPointCenters);

  vtkSmartPointer<vtkVertexGlyphFilter> centersGlyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  centersGlyphFilter->SetInputConnection(centersPolyData->GetProducerPort());
  centersGlyphFilter->Update();

  outputCenters->ShallowCopy(centersGlyphFilter->GetOutput());

  // Average normals
  vtkDataArray* inputNormals =
    vtkDataArray::SafeDownCast(input->GetPointData()->GetNormals());

  vtkSmartPointer<vtkDoubleArray> superNormals =
    vtkSmartPointer<vtkDoubleArray>::New();
  superNormals->SetNumberOfComponents(3);
  superNormals->SetName("Normals");

  if(inputNormals)
    {
    for(vtkIdType superID = 0; superID < centersPolyData->GetNumberOfPoints(); superID++)
      {
      double x=0; double y=0; double z=0;
      int count = 0;
      for(vtkIdType point = 0; point < outputLabeled->GetNumberOfPoints(); point++)
        {
        if(superPointLabels[point] == superID)
          {
          double n[3];
          inputNormals->GetTuple(point, n);
          x += n[0];
          y += n[1];
          z += n[2];
          count++;
          }
        }
      x /= static_cast<double>(count);
      y /= static_cast<double>(count);
      z /= static_cast<double>(count);

      double superNormal[3];
      superNormal[0] = x;
      superNormal[1] = y;
      superNormal[2] = z;
      superNormals->InsertNextTupleValue(superNormal);
      }

    outputCenters->GetPointData()->SetNormals(superNormals);
    }
  else
    {
    std::cerr << "No normals on input!" << std::endl;
    }


  // Average colors
  vtkUnsignedCharArray* inputColors =
    vtkUnsignedCharArray::SafeDownCast(input->GetPointData()->GetArray("Colors"));

  vtkSmartPointer<vtkUnsignedCharArray> superColors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  superColors->SetNumberOfComponents(3);
  superColors->SetName("Colors");

  if(inputColors)
    {
    for(vtkIdType superID = 0; superID < centersPolyData->GetNumberOfPoints(); superID++)
      {
      double r=0; double g=0; double b=0;
      int count = 0;
      for(vtkIdType point = 0; point < outputLabeled->GetNumberOfPoints(); point++)
        {
        if(superPointLabels[point] == superID)
          {
          unsigned char c[3];
          inputColors->GetTupleValue(point, c);
          r += c[0];
          g += c[1];
          b += c[2];
          count++;
          }
        }
      r /= static_cast<double>(count);
      g /= static_cast<double>(count);
      b /= static_cast<double>(count);

      unsigned char superColor[3];
      superColor[0] = static_cast<unsigned char>(r);
      superColor[1] = static_cast<unsigned char>(g);
      superColor[2] = static_cast<unsigned char>(b);
      superColors->InsertNextTupleValue(superColor);
      }
    outputCenters->GetPointData()->SetScalars(superColors);
    }
  else
    {
    std::cerr << "No colors on input!" << std::endl;
    }
    
  return 1;
}


//----------------------------------------------------------------------------
void vtkSuperPoints::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "NNRadius: (" << this->NNRadius << ")\n";
  os << indent << "AutomaticRadius: (" << this->AutomaticRadius << ")\n";
  os << indent << "AutomaticRadiusRatio: (" << this->AutomaticRadiusRatio << ")\n";
}