// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_SPARSEMATRIX
#define INCLUDE_SPARSEMATRIX

#include "array.h"
#include "mex.h"

// Type definitions.
// -----------------------------------------------------------------
#ifdef MWINDEXISINT
// This line is needed for versions of MATLAB prior to 7.3.
typedef int mwIndex;
#endif

// Function declarations.
// ---------------------------------------------------------------
int getSparseMatrixSize     (const mxArray* ptr);
int isSparseLowerTriangular (const mxArray* ptr);

// class SparseMatrixStructure
// ---------------------------------------------------------------
// An object of class SparseMatrixStructure stores information about
// the structure of a sparse matrix. It does not store the actual
// values of the matrix entries.
//
// WARNING: Starting with version 7.3, MATLAB can handle 64-bit
// addressing, and the authors of MATLAB have modified the
// implementation of sparse matrices to reflect this change. However,
// I convert all the row and column indices in the sparse matrix to
// signed integers, and this could potentially cause problems when
// dealing with large, sparse matrices on 64-bit platforms with MATLAB
// version 7.3 or greater.
class SparseMatrixStructure {
public:

  // This constructor takes as input a Matlab array. It it points to a
  // valid sparse matrix, it will store all the information pertaining
  // to the sparse matrix structure. If "makeCopy" is true, then the
  // object will obtain an independent copy of the sparse matrix
  // structure. If not, the object will be dependent on the data in
  // memory.
  explicit SparseMatrixStructure (const mxArray* ptr, 
				  bool makeCopy = false);
    
  // The copy constructor makes a shallow copy of the source object.
  SparseMatrixStructure (const SparseMatrixStructure& source);
    
  // The destructor.
  ~SparseMatrixStructure();
    
  // Get the height and width of the matrix.
  int height() const { return h; };
  int width () const { return w; };

  // Return the number of non-zero entries.
  int size() const { return nnz; };

  // Return the number of non-zero entries in the cth column.
  int size (int c) const;

  // Upon completion of this function, cols[i] contains the column
  // index for the ith element, and rows[i] contains the row index for
  // the ith element. It is assumed that "cols" and "rows" have
  // sufficient space to store this information. This routine is most
  // useful for converting the Matlab sparse matrix format into the
  // IPOPT format.
  void getColsAndRows (int* cols, int* rows) const;

  // Copy the matrix entries in a sensible manner while preserving the
  // structure of the destination. In order to preserve the structure
  // of the destination, it is required that its set of non-zero
  // entries be a (non-strict) superset of the non-zero entries of the
  // source.
  friend void copyElems (const SparseMatrixStructure& sourceStructure,
			 const SparseMatrixStructure& destStructure,
			 const double* sourceValues, double* destValues);

protected:
  mwIndex* jc;      // See mxSetJc in the MATLAB documentation.
  mwIndex* ir;      // See mxSetIr in the MATLAB documentation.
  int      nnz;     // The number of non-zero elements.
  int      h;       // The height of the matrix. 
  int      w;       // The width of the matrix.
  bool     owner;   // Whether or not the object has ownership of the
                    // "jc" and "ir" matrices.

  // The copy assignment operator is kept hidden because we don't
  // want it to be used.
  SparseMatrixStructure& operator= (const SparseMatrixStructure& source)
  { return *this; };
};

#endif
