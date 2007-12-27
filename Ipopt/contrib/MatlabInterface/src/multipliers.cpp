#include "multipliers.h"

// Constants.
// -----------------------------------------------------------------
const char* lowerBoundMultipliersLabel = "zl";
const char* upperBoundMultipliersLabel = "zu";
const char* constraintMultipliersLabel = "lambda";

// Function definitions for clas Multipliers.
// -----------------------------------------------------------------
Multipliers::Multipliers (const mxArray*& ptr) {
  mxArray* p;

  // First check to see whether the MATLAB array is a structure
  // array. If not, throw an exception.
  if (!mxIsStruct(ptr))
    throw MatlabException("Matlab array must be a structure array");
    
  // Get the multipliers corresponding to the lower bounds on the
  // optimization variables.
  p = mxGetField(ptr,0,lowerBoundMultipliersLabel);
  if (p == 0)
    throw MatlabException("MATLAB multipliers input does not have the \
correct fields");
  zl = new Matrix(p);

  // Get the multipliers corresponding to the upper bounds on the
  // optimization variables.
  p = mxGetField(ptr,0,upperBoundMultipliersLabel);
  if (p == 0)
    throw MatlabException("MATLAB multipliers input does not have the \
correct fields");
  zu = new Matrix(p);  

  // Get the multipliers corresponding to the upper bounds on the
  // optimization variables.
  p = mxGetField(ptr,0,constraintMultipliersLabel);
  if (p == 0)
    throw MatlabException("MATLAB multipliers input does not have the \
correct fields");
  lambda = new Matrix(p);    
}

Multipliers::Multipliers (mxArray*& ptr, int n, int m) {
  const char* fieldnames[3];
  mxArray*    p;

  zl     = 0;
  zu     = 0;
  lambda = 0;

  // Create the Matlab struct array.
  fieldnames[0] = lowerBoundMultipliersLabel;
  fieldnames[1] = upperBoundMultipliersLabel;
  fieldnames[2] = constraintMultipliersLabel;
  ptr           = mxCreateStructMatrix(1,1,3,fieldnames);

  // Set up the multipliers corresponding to the lower bounds on the
  // optimization variables.
  zl = new Matrix(p,1,n);
  mxSetField(ptr,0,lowerBoundMultipliersLabel,p);

  // Set up the multipliers corresponding to the upper bounds on the
  // optimization variables.
  zu = new Matrix(p,1,n);
  mxSetField(ptr,0,upperBoundMultipliersLabel,p);
  
  // Set up the multipliers corresponding to the constraints.
  lambda = new Matrix(p,1,m);
  mxSetField(ptr,0,constraintMultipliersLabel,p);
} 

Multipliers::Multipliers (const Multipliers& source) {
  zl     = new Matrix(*source.zl);
  zu     = new Matrix(*source.zu);
  lambda = new Matrix(*source.lambda);
}

Multipliers::~Multipliers() {
  if (zl)     delete zl;
  if (zu)     delete zu;
  if (lambda) delete lambda;
}
