// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "sparsematrix.hpp"
#include "matlabexception.hpp"
#include "iterate.hpp"

// Function definitions for class SparseMatrix.
// ---------------------------------------------------------------
SparseMatrix::SparseMatrix (const mxArray* ptr) 
  : jc(0), ir(0), x(0) {

  // Get the height, width and number of non-zeros.
  h   = (int) mxGetM(ptr);
  w   = (int) mxGetN(ptr);
  nnz = getSizeOfSparseMatrix(ptr);

  // Copy the row and column indices, and the values of the nonzero
  // entries.
  jc = new mwIndex[w+1];
  ir = new mwIndex[nnz];
  x  = new double[nnz];
  copymemory(mxGetJc(ptr),jc,w+1);
  copymemory(mxGetIr(ptr),ir,nnz);
  copymemory(mxGetPr(ptr),x,nnz);
}

SparseMatrix::~SparseMatrix() {
  if (jc) delete[] jc;
  if (ir) delete[] ir;
  if (x)  delete[] x;
}

int SparseMatrix::numelems (int c) const {
  return (int) jc[c+1] - (int) jc[c];
}

void SparseMatrix::getColsAndRows (int* cols, int* rows) const {
       
  // Repeat for each column in the matrix, then repeat for each
  // non-zero entry in the current column.
  for (int c = 0, i = 0; c < w; c++)
    for ( ; i < (int) jc[c+1]; i++) {
      cols[i] = (int) c;
      rows[i] = (int) ir[i];
    }
}

bool SparseMatrix::copyto (SparseMatrix& dest) const {
  bool match, matchrow, matchcol;  // Loop variables.

  // Initialize the destination values to zero, because we might not
  // have corresponding nonzero entries from the source for some of
  // the destination nonzero entries.
  for (int i = 0; i < dest.nnz; i++)
    dest.x[i] = 0;

  // In order to properly copy the elements from one sparse matrix to
  // another, the non-zero elements of the destination matrix must be
  // a superset of the collection of non-zero elements in the source
  // matrix. If not, the copy operation is invalid, and the function
  // returns false. Repeat for each column.
  int i = 0;  // Index of element in source.
  int j = 0;  // Index of element in destination.
  for (int c = 0; c < dest.w; c++) {

    // Repeat for each non-zero element in the destination column.
    for ( ; j < (int) dest.jc[c+1]; j++) {

      // Let's check to see if the source column matches the
      // destination column, and the source row matches the
      // destination row. The first line checks to see if the row
      // indices match. The second line checks to see if the column
      // indices match. (If the source matrix is valid, there is no
      // need to check whether the column index of the source entry is
      // LESS THAN the destination column index since we assume that
      // there are less entries in the source, and under that
      // assumption we move faster through the source matrix than we
      // do through the destination matrix.)
      matchrow = (ir[i] == dest.ir[j]);
      matchcol = (i >= (int) jc[c]) && (i < (int) jc[c+1]);
      match    = matchrow && matchcol;

      // If the row & column indices match, then we copy the source
      // entry value and move on to the next non-zero entry in the
      // source.
      dest.x[j] = match * x[i];
      i        += match;
    }      
  }

  // If we've reached the end of the loop and we haven't visited all
  // the non-zero entries in the source matrix, it means that the
  // source matrix is invalid, and we throw an error.
  if (i < nnz)
    return false;
  else
    return true;
}

void SparseMatrix::copyto (double* dest) const {
  copymemory(x,dest,nnz);
}

// Function definitions for static members of class SparseMatrix.
// -----------------------------------------------------------------
int SparseMatrix::getSizeOfSparseMatrix (const mxArray* ptr) {
  
  // Get the width (the number of columns) of the matrix.
  int w = (int) mxGetN(ptr);
  
  // The precise number of non-zero elements is contained in the
  // last entry of the jc array. (There is one jc entry for each
  // column in the matrix, plus an extra one.)
  mwIndex* jc = mxGetJc(ptr);
  return (int) jc[w];    
}  

bool SparseMatrix::isLowerTri (const mxArray* ptr) {

  // Get the height and width of the matrix.
  int h = (int) mxGetM(ptr);
  int w = (int) mxGetN(ptr);
  
  // Check whether the sparse matrix is symmetric.
  if (h != w)
    return false;

  // A sparse lower triangular matrix has the property that the
  // column indices are always less than or equal to the row
  // indices.
  bool     b  = true;  // The return value.
  mwIndex* jc = mxGetJc(ptr);
  mwIndex* ir = mxGetIr(ptr);
  for (int c = 0, i = 0; c < w; c++)
    for ( ; i < (int) jc[c+1]; i++)
      b = b && (c <= (int) ir[i]);
  
  return b;
}

bool SparseMatrix::inIncOrder (const mxArray* ptr) {
  bool     b  = true;
  int      w  = (int) mxGetN(ptr);
  mwIndex* jc = mxGetJc(ptr);
  mwIndex* ir = mxGetIr(ptr);

  for (int c = 0, i = 0; c < w; c++)
    if (jc[c+1] > jc[c]) {
      i++;
      for ( ; i < (int) jc[c+1]; i++)
	b = b && (ir[i] > ir[i-1]);
    }

  return b;
}
