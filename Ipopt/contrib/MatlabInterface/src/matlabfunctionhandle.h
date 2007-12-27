#ifndef INCLUDE_MATLABFUNCTIONHANDLE
#define INCLUDE_MATLABFUNCTIONHANDLE

#include "mex.h"

// Class MatlabFunctionHandle.
// -----------------------------------------------------------------
class MatlabFunctionHandle {
public:

  // This is the default constructor. It stores a null function
  // handle.
  MatlabFunctionHandle();

  // This constructor accepts as input a pointer to a Matlab array,
  // which must be of the function handle class.
  explicit MatlabFunctionHandle (const mxArray* ptr);

  // The copy constructor makes a full copy of the source object.
  MatlabFunctionHandle (const MatlabFunctionHandle& source);
  
  // The destructor.
  ~MatlabFunctionHandle();

  // Conversion operator for pointer to MATLAB array.
  operator mxArray* () const { return f; };  

  // Returns false if and only if we have the null function handle.
  operator bool() const { return f != 0; };

protected:
  mxArray* f;  // The MATLAB array storing information concerning 
               // the function handle.

  // The copy assignment operator is not proper, thus remains protected.
  MatlabFunctionHandle& operator= (const MatlabFunctionHandle& source) 
  { return *this; };
};

#endif
