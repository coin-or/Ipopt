// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 15, 2008

#ifndef INCLUDE_ITERATE
#define INCLUDE_ITERATE

#include "mex.h"
#include <string.h>

// Function definitions.
// ---------------------------------------------------------------
template <class Type> void copymemory (const Type* source, Type* dest, int n) {
  memcpy(dest,source,sizeof(Type)*n);
}

// Class Iterate.
// -----------------------------------------------------------------
// An object of this class stores all the information we would like to
// keep track of regarding the variables in our optimization
// problem. An Iterate object is constructed from a MATLAB array,
// which is either an array in double precision, or a cell array with
// elements that are arrays in double precision. There are two methods
// of interest: a function that copies the elements from a source
// array to the Iterate object (inject), and a function that copies
// the Itereate to a destination array (copyto).
class Iterate {
public:

  // This constructor initializes the object basically by creating a
  // duplicate of the specified MATLAB array. It is important to note
  // that the Iterate object does *not* create an independent copy of
  // the MATLAB array. As such, it will not destroy the MATLAB arrays
  // when the object is destroyed.
  explicit Iterate (mxArray* ptr);

  // The destructor.
  ~Iterate() { };

  // Return the number of variables.
  friend int numvars (const Iterate& x) { return x.nv; };

  // Copy the elements from the source array.
  void inject (const double* x);

  // Copy the elements to the location in memory pointed to by "dest".
  // It is assumed that sufficient memory is allocated for the
  // destination.
  void copyto (double* x) const;

  // Convert the Iterate object to a MATLAB array.
  operator       mxArray*()       { return ptr; };
  operator const mxArray*() const { return ptr; };

  // This function dynamically creates a new array containing all the
  // information about the iterate as specified by the MATLAB
  // array. It is up to the user to deallocate the memory associated
  // with the newly allocated data array. The return value is the size
  // of the array.
  static int getMatlabData (const mxArray* ptr, double*& data);

protected:
  int      nv;   // The number of optimization variables.
  mxArray* ptr;  // All the information is stored in a MATLAB array.
};

#endif
