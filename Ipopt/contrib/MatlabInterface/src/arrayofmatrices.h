// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_ARRAYOFMATRICES
#define INCLUDE_ARRAYOFMATRICES

#include "array.h"
#include "matlabmatrix.h"
#include "mex.h"

// Class ArrayOfMatrices.
// -----------------------------------------------------------------
class ArrayOfMatrices : public Array<Matrix*> {
public:
    
  // This version of the constructor behaves just like its parent.
  explicit ArrayOfMatrices (int length) 
    : Array<Matrix*>(length) { };
    
  // This constructor creates an array of matrices from a Matlab
  // array. It accepts either a matrix in double precision, or a cell
  // array with entries that are matrices.
  explicit ArrayOfMatrices (const mxArray* ptr);

  // This constructor creates an array of matrices from a collection
  // of Matlab arrays. The Matlab arrays must be matrices.
  ArrayOfMatrices (const mxArray* ptrs[], int numptrs);

  // This constructor creates an array of matrices and the
  // associated Matlab structures. The Matlab structures are
  // matrices. The second input argument acts as a template for the
  // creation of the matrices, but the data from "model" is not
  // actually copied into the new ArrayOfMatrices object. It is up
  // to the user to make sure that the array of mxArray pointers has
  // enough room for the pointers to the Matlab arrays.
  ArrayOfMatrices (mxArray* ptrs[], const ArrayOfMatrices& model);

  // This constructor creates an array of matrices using the second
  // input argument as a model. The input argument "data" contains
  // the element data. Note that the information is NOT copied from
  // the model!
  ArrayOfMatrices (double* data, const ArrayOfMatrices& model);

  // The copy constructor makes a shallow copy of the source object.
  ArrayOfMatrices (const ArrayOfMatrices& source);
    
  // The destructor.
  ~ArrayOfMatrices();
  
  // Copy assignment operator that observes the same behaviour as
  // the Array copy assignment operator.
  ArrayOfMatrices& operator= (const ArrayOfMatrices& source);
  
  // Copy the data from the location in memory pointed to by "source"
  // to the matrices in the given order. It is assumed that the input
  // points to a valid source of data.
  void inject (const double* source);

  // Copy the elements from all the matrices to the location in
  // memory pointed to by "dest". It is assumed that sufficient
  // memory is allocated for the destination.
  void copyto (double* dest) const;

  // Returns true if the two objects have the exact same dimensions.
  bool operator== (const ArrayOfMatrices& a) const;
  bool operator!= (const ArrayOfMatrices& a) const 
  { return !(*this == a); };
  
  // Return the total number of elements in all the matrices.
  int numelems() const;
  
protected:
  static int getnummatlabmatrices (const mxArray* ptr);
};

#endif
