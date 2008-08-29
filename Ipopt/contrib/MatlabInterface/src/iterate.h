#ifndef INCLUDE_ITERATE
#define INCLUDE_ITERATE

#include "mex.h"
#include <string.h>

// Function declarations.
// -----------------------------------------------------------------
void copymemory (const double* source, double* dest, int n);

// Class Iterate.
// -----------------------------------------------------------------
// An object of this class stores all the information we would like 
// to keep track of regarding the variables in our optimization problem.
class Iterate {
public:

  // This constructor initializes the object basically by creating a
  // duplicate of the specified MATLAB array.
  Iterate (const mxArray* source);

  // The copy constructor makes a full, independent copy of the source
  // object.
  Iterate (const Iterate& source);

  // The destructor.
  ~Iterate();

  // Return the number of variables.
  int numvars() const { return nv; };

  // Copy the elements from the source array.
  void inject (const double* x);

  // Copy the elements to the location in memory pointed to by "dest".
  // It is assumed that sufficient memory is allocated for the
  // destination.
  void copyto (double* x) const;

  // Convert the Iterate object to a MATLAB array.
  mxArray* operator() { return ptr; };

protected:
  int      nv;   // The number of optimization variables.
  mxArray* ptr;  // All the information is stored in a MATLAB array.
};

#endif
