#include "iterate.h"
#include "matlabexception.h"

// Function definitions.
// -----------------------------------------------------------------
void copymemory (const double* source, double* dest, int n) {
  memcpy(dest,source,sizeof(double)*n);
}

// Function definitions for class Iterate.
// -----------------------------------------------------------------
Iterate::Iterate (const mxArray* source) 
  : nv(0), ptr(0) {
  const mxArray* p;  // Pointer to a MATLAB array.

  // Duplicate the MATLAB array.
  ptr = mxDuplicateArray(source);

  // Compute the number of optimization variables.
  if (mxIsCell(ptr)) {

    // The MATLAB array is a cell array. Repeat for each cell.
    int n = mxGetNumberOfElements(p);
    for (int i = 0; i < n; i++) {
      p = mxGetCell(ptr,i);  // Get the ith cell.
      if (!mxIsDouble(p))
	throw MatlabException("The initial iterate must be either an array \
in DOUBLE precision, or a cell array in which each cell is an array in \
DOUBLE precision");      
      nv += mxGetNumberOfElements(p);
    }
  } else {
    
    // The MATLAB array should be a numeric array.
    if (!mxIsDouble(ptr))
      throw MatlabException("The initial iterate must be either an array \
in DOUBLE precision, or a cell array in which each cell is an array in \
DOUBLE precision");
    nv = mxGetNumberOfElements(ptr);
  }
}

Iterate::Iterate (const Iterate& source) 
  : nv(source.nv), ptr(mxDuplicateArray(source.ptr) { }

Iterate::~Iterate() {
  if (ptr) mxDestroyArray(ptr);
}

void Iterate::inject (const double* x) {
  if (mxIsCell(ptr)) {

    // The MATLAB array is a cell array. Repeat for each cell.
    mxArray* p; 
    int      m;
    int      n = mxGetNumberOfElements(ptr);
    for (int i = 0; i < n; i++) {
      p = mxGetCell(ptr,i);  // Get the ith cell.
      m = mxGetNumberOfElements(p);
      copymemory(x,mxGetPr(p),m);
      x += m;
    }
  } else

    // The MATLAB array is a numeric array.
    copymemory(x,mxGetPr(ptr),nv);
}

void Iterate::copyto (double* x) const {
  if (mxIsCell(ptr)) {

    // The MATLAB array is a cell array. Repeat for each cell.
    const mxArray* p; 
    int            m;
    int            n = mxGetNumberOfElements(ptr);
    for (int i = 0; i < n; i++) {
      p = mxGetCell(ptr,i);  // Get the ith cell.
      m = mxGetNumberOfElements(p);
      copymemory(mxGetPr(p),x,m);
      x += m;
    }
  } else

    // The MATLAB array is a numeric array.
    copymemory(mxGetPr(ptr),x,nv);
}

