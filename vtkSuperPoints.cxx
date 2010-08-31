#include "vtkSuperPoints.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"

#include <vtkMath.h>
#include <vtkLookupTable.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkIdFilter.h>
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
#include <vtkGraphWriter.h>
#include <vtkExtractSelectedGraph.h>

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

// For testing only
#include <vtkGraphToPolyData.h>
#include <vtkXMLPolyDataWriter.h>

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

  vtkSmartPointer<vtkIdFilter> idFilter = 
    vtkSmartPointer<vtkIdFilter>::New();
  idFilter->SetInputConnection(input->GetProducerPort());
  idFilter->SetIdsArrayName("OriginalIds");
  idFilter->Update();
     
  vtkDataArray* normalsInput = input->GetPointData()->GetNormals();
  if(!normalsInput)
    {
    std::cerr << " no input normals" << std::endl;
    exit(-1);
    }
    
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

  // Build a graph on all of the input points

  vtkSmartPointer<vtkNearestNeighborGraphFilter> nnFilter =
    vtkSmartPointer<vtkNearestNeighborGraphFilter>::New();
  nnFilter->SetInputConnection(idFilter->GetOutputPort());
  nnFilter->Update();
  
  // Connect the NN graph
  vtkSmartPointer<vtkConnectGraph> connectFilter =
    vtkSmartPointer<vtkConnectGraph>::New();
  connectFilter->SetInputConnection(nnFilter->GetOutputPort());
  connectFilter->Update();
  
  {
  vtkSmartPointer<vtkGraphWriter> graphWriter = 
    vtkSmartPointer<vtkGraphWriter>::New();
  graphWriter->SetFileName("NNgraph.graph");
  graphWriter->SetInputConnection(connectFilter->GetOutputPort());
  graphWriter->Write();
  
  vtkSmartPointer<vtkGraphToPolyData> graphToPolyData = 
    vtkSmartPointer<vtkGraphToPolyData>::New();
  graphToPolyData->SetInputConnection(connectFilter->GetOutputPort());
  graphToPolyData->Update();
  
  vtkSmartPointer<vtkXMLPolyDataWriter> graphPolyDataWriter = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  graphPolyDataWriter->SetInputConnection(graphToPolyData->GetOutputPort());
  graphPolyDataWriter->SetFileName("NNgraph.vtp");
  graphPolyDataWriter->Write();
    
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

  for(vtkIdType pointID = 0; pointID < input->GetNumberOfPoints(); ++pointID)
  //for(vtkIdType pointID = 0; pointID < 4000; ++pointID)
    {
    if(pointID % 1000 == 0)
      {
      std::cout << "point " << pointID << " out of " << input->GetNumberOfPoints() << std::endl;
      }
    //std::cout << "point " << pointID << std::endl;
    
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
      continue; // All points in this neighborhood already belong to a superpixel
      }

    vtkSmartPointer<vtkSelectionNode> selectionNode =
      vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::VERTEX);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);

    vtkSmartPointer<vtkSelection> selection =
      vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);
    
    vtkSmartPointer<vtkExtractSelectedGraph> extractSelection =
      vtkSmartPointer<vtkExtractSelectedGraph>::New();
    extractSelection->SetInput(0, connectFilter->GetOutput());
    extractSelection->SetInput(1, selection);
    extractSelection->Update();
    
    //std::cout << "There are " << extractSelection->GetOutput()->GetNumberOfVertices() << " extracted points." << std::endl;
    
    vtkDataSetAttributes* vertexData = extractSelection->GetOutput()->GetVertexData();
        
    vtkDataArray* normalsAfter = vtkDataArray::SafeDownCast(vertexData->GetNormals());
    
    vtkIdTypeArray* originalIds =
      vtkIdTypeArray::SafeDownCast(vertexData->GetArray("OriginalIds"));

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
    graphIterator->SetGraph(extractSelection->GetOutput());
    //graphIterator->SetGraph(connectFilter->GetOutput());
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
    
    
  // Setup the colors array
  
  vtkSmartPointer<vtkLookupTable> lut = 
    vtkSmartPointer<vtkLookupTable>::New();
  lut->SetNumberOfTableValues(superPointCenters->GetNumberOfPoints());
  
  vtkSmartPointer<vtkUnsignedCharArray> superColors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  superColors->SetNumberOfComponents(3);
  superColors->SetName("SuperColors");
  
  vtkMath::RandomSeed(time(NULL));
  
  for(unsigned int i = 0; i < superPointCenters->GetNumberOfPoints(); i++)
    {
    lut->SetTableValue(i, vtkMath::Random(0.0,1.0), vtkMath::Random(0.0,1.0), vtkMath::Random(0.0,1.0));
    }
    
  for(unsigned int i = 0; i < superPointLabels.size(); i++)
    {
    double randomRGBA[4];
    lut->GetTableValue(superPointLabels[i], randomRGBA);
  
    unsigned char randomColor[3];
    randomColor[0] = static_cast<unsigned char>(255. * randomRGBA[0]);
    randomColor[1] = static_cast<unsigned char>(255. * randomRGBA[1]);
    randomColor[2] = static_cast<unsigned char>(255. * randomRGBA[2]);
  
    std::cout << "Converted " << (int)randomRGBA[0] << " " << (int)randomRGBA[1] << " " << (int)randomRGBA[2]
	      << " to " << (int)randomColor[0] << " " << (int)randomColor[1] << " " << (int)randomColor[2] << std::endl;
    superColors->InsertNextTupleValue(randomColor);
    }

  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputConnection(polydata->GetProducerPort());
  glyphFilter->Update();

  outputLabeled->ShallowCopy(glyphFilter->GetOutput());

  outputLabeled->GetPointData()->SetScalars(labelArray);
  outputLabeled->GetPointData()->AddArray(superColors);

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

  vtkSmartPointer<vtkUnsignedCharArray> averageColors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  averageColors->SetNumberOfComponents(3);
  averageColors->SetName("AverageColors");

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

      unsigned char averageColor[3];
      averageColor[0] = static_cast<unsigned char>(r);
      averageColor[1] = static_cast<unsigned char>(g);
      averageColor[2] = static_cast<unsigned char>(b);
      averageColors->InsertNextTupleValue(averageColor);
      }
    outputCenters->GetPointData()->SetScalars(averageColors);
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