#ifndef __vtkSuperPoints_h
#define __vtkSuperPoints_h

#include "vtkPolyDataAlgorithm.h"

class vtkSuperPoints : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSuperPoints,vtkAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSuperPoints *New();

  vtkSetMacro(NNRadius, double);
  vtkSetMacro(AutomaticRadius, int);
  vtkSetMacro(AutomaticRadiusRatio, int);

protected:
  vtkSuperPoints();
  ~vtkSuperPoints(){}

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkSuperPoints(const vtkSuperPoints&);  // Not implemented.
  void operator=(const vtkSuperPoints&);  // Not implemented.

  double NNRadius;
  int AutomaticRadius;
  int AutomaticRadiusRatio;
};

#endif