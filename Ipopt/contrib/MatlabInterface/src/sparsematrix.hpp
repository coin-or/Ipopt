// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_SPARSEMATRIX
#define INCLUDE_SPARSEMATRIX

#include "mex.h"

// Type definitions.
// -----------------------------------------------------------------
// This line is needed for versions of MATLAB prior to 7.3.
#ifdef MWINDEXISINT
typedef int mwIndex;
#endif

// class SparseMatrix
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
class SparseMatrix {
public:

  // This constructor takes as input a Matlab array. It it points to a
  // valid sparse matrix in double precision, it will store all the
  // information pertaining to the sparse matrix structure. It is up
  // to the user to ensure that the MATLAB array is a sparse,
  // symmetric matrix with row indices in increasing order as the
  // nonzero elements appear in the matrix. Note that a SparseMatrix
  // object retains a completely independent copy of the sparse matrix
  // information by duplicating the data from the specified MATLAB
  // array.
  explicit SparseMatrix (const mxArray* ptr);

  // The destructor.
  ~SparseMatrix();
    
  // Get the height and width of the matrix.
  friend int height (const SparseMatrix& A) { return A.h; };
  friend int width  (const SparseMatrix& A) { return A.w; };

  // The first function returns the total number of non-zero entries.
  // The second function returns the number of non-zero entries in the
  // cth column.
  int numelems ()      const { return nnz; };
  int numelems (int c) const;

  // Upon completion of this function, cols[i] contains the column
  // index for the ith element, and rows[i] contains the row index for
  // the ith element. It is assumed that "cols" and "rows" have
  // sufficient space to store this information. This routine is most
  // useful for converting the Matlab sparse matrix format into the
  // IPOPT format.
  void getColsAndRows (int* cols, int* rows) const;

  // Copy the matrix entries in a sensible manner while preserving the
  // structure of the destination. In order to preserve the structure
  // of the destination, it is required that the source set of
  // non-zero entries be a subset of the destination non-zero
  // entries. On success, the value true is returned.
  bool copyto (SparseMatrix& dest) const;

  // Copy the values of the nonzero elements to the destination array
  // which of course must be of the proper length.
  void copyto (double* dest) const;

  // Returns the number of nonzeros in the sparse matrix.
  static int getSizeOfSparseMatrix (const mxArray* ptr);

  // Returns true if and only if the sparse matrix is symmetric and
  // lower triangular.
  static bool isLowerTri (const mxArray* ptr);

  // For the proper functioning of a sparse matrix object, it is
  // necessary that the row indices be in increasing order.
  static bool inIncOrder (const mxArray* ptr);

protected:
  int      h;    // The height of the matrix. 
  int      w;    // The width of the matrix.
  int      nnz;  // The number of non-zero elements.
  mwIndex* jc;   // See mxSetJc in the MATLAB documentation.
  mwIndex* ir;   // See mxSetIr in the MATLAB documentation.
  double*  x;    // The values of the non-zero entries.
};

#endif
