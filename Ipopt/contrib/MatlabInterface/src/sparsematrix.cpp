// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "sparsematrix.h"
#include "matlabexception.h"

// Function definitions.
// -----------------------------------------------------------------
int getSparseMatrixSize (const mxArray* ptr) {
  if (!mxIsSparse(ptr))
    throw MatlabException("Matlab array must be a sparse matrix");
  
  // Get the width (the number of columns) of the matrix.
  int w = (int) mxGetN(ptr);
  
  // The precise number of non-zero elements is contained in the
  // last entry of the jc array. (There is one jc entry for each
  // column in the matrix, plus an extra one.)
  mwIndex* jc = mxGetJc(ptr);
  return jc[w];    
}  

int isSparseLowerTriangular (const mxArray* ptr) {
  bool islt = true;  // The return value.

  // Get the dimensions of the matrix.
  int height = (int) mxGetM(ptr);
  int width  = (int) mxGetN(ptr);
  
  // Check whether the Matlab array is a proper sparse, N x N
  // matrix.
  if (!mxIsSparse(ptr))
    throw MatlabException("Matlab array must be a sparse matrix");
  if ((mxGetNumberOfDimensions(ptr) != 2) || (height != width))
    throw MatlabException("Matlab array must be an N x N matrix");

  // A sparse lower triangular matrix has the property that the
  // column indices are always less than or equal to the row
  // indices.
  mwIndex* jc = mxGetJc(ptr);
  mwIndex* ir = mxGetIr(ptr);
  for (int c = 0, i = 0; c < width; c++)
    for ( ; i < (int) jc[c+1]; i++)
      islt = islt && (c <= (int) ir[i]);
  
  return islt;
}

// Function definitions for class SparseMatrixStructure
// ---------------------------------------------------------------
SparseMatrixStructure::SparseMatrixStructure (const mxArray* ptr,
					      bool makeCopy) {

  if (mxGetNumberOfDimensions(ptr) != 2)
    throw MatlabException("Matlab array must be a matrix");
  if (!mxIsSparse(ptr))
    throw MatlabException("Matlab array must be a sparse matrix");
  if (!mxIsDouble(ptr))
    throw MatlabException("Matlab array must be in double precision");
  
  jc = 0;
  ir = 0;

  // Get the row and column indices.
  mwIndex* jcMatlab = mxGetJc(ptr);
  mwIndex* irMatlab = mxGetIr(ptr);

  // Get the height and width of the matrix.
  h = (int) mxGetM(ptr);
  w = (int) mxGetN(ptr);

  // Get the number of non-zero entries.
  nnz = getSparseMatrixSize(ptr);
  
  if (makeCopy) {
    owner = true;
    
    // Copy the row and column indices.
    jc = new mwIndex[w+1];
    ir = new mwIndex[nnz];
    copymemory(jcMatlab,jc,w+1);
    copymemory(irMatlab,ir,nnz);
  }
  else {
    owner = false;
    jc    = jcMatlab;
    ir    = irMatlab;
  }
  
  // For the proper functioning of a sparse matrix object, it is
  // necessary that the row indices be in increasing order.
  bool inIncOrder = true;
  for (int c = 0, i = 0; c < w; c++)
    if (size(c)) {
      i++;
      for ( ; i < (int) jc[c+1]; i++)
	inIncOrder = inIncOrder & (ir[i] > ir[i-1]);
    }    
  if (!inIncOrder)
    throw MatlabException("The rows in the sparse matrix are not in \
increasing order, as required");
}

SparseMatrixStructure::SparseMatrixStructure 
(const SparseMatrixStructure& source) {
  h     = source.h;
  w     = source.w;
  nnz   = source.nnz;
  jc    = source.jc;
  ir    = source.ir;
  owner = false;
}
  
SparseMatrixStructure::~SparseMatrixStructure() {
  if (owner) {
    if (jc) 
      delete[] jc;
    if (ir) 
      delete[] ir;
  }
}

int SparseMatrixStructure::size (int c) const {
  return jc[c+1] - jc[c];
}

void SparseMatrixStructure::getColsAndRows (int* cols, int* rows) 
const {
       
  // Repeat for each column in the matrix, then repeat for each
  // non-zero entry in the current column.
  for (int c = 0, i = 0; c < w; c++)
    for ( ; i < (int) jc[c+1]; i++) {
      cols[i] = (int) c;
      rows[i] = (int) ir[i];
    }
}

void copyElems (const SparseMatrixStructure& sourceStructure,
		const SparseMatrixStructure& destStructure,
		const double* sourceValues, double* destValues) {
  bool match, matchrow, matchcol;  // Loop variables.

  // First, make sure that the dimensions of the soure and destination match.
  if (sourceStructure.height() != destStructure.height())
    throw MatlabException("Unable to copy sparse matrix; source has \
a different height than destination");
  if (sourceStructure.width() != destStructure.width())
    throw MatlabException("Unable to copy sparse matrix; source has \
a different width than destination");
  if (sourceStructure.size() > destStructure.size())
    throw MatlabException("Unable to copy sparse matrix; source has \
more non-zero elements than destination");

  // Initialize the destination values to zero.
  for (int i = 0; i < destStructure.size(); i++)
    destValues[i] = 0;

  // In order to properly copy the elements from one sparse matrix to
  // another, the non-zero elements of the source matrix must be a
  // subset of the collection of non-zero elements in the
  // destination. If not, the copy operation is invalid.
  //
  // Repeat for each column. The index of the current element from
  // the source matrix is represented by i, and the index of the
  // current element in the destination is represented by j.
  int i = 0, j = 0;
  for (int c = 0; c < destStructure.w; c++) {

    // Repeat for each non-zero element in the destination column.
    for ( ; j < (int) destStructure.jc[c+1]; j++) {

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
      matchrow = (sourceStructure.ir[i] == destStructure.ir[j]);
      matchcol = (i >= (int) sourceStructure.jc[c]) && 
	         (i <  (int) sourceStructure.jc[c+1]);
      match    = matchrow && matchcol;

      // If the row & column indices match, then we copy the source
      // entry value and move on to the next non-zero entry in the
      // source.
      destValues[j] = match * sourceValues[i];
      i            += match;
    }      
  }

  // If we've reached the end of the loop and we haven't visited all
  // the non-zero entries in the source matrix, it means that the
  // source matrix is invalid, and we throw an error.
  if (i < sourceStructure.size())
    throw MatlabException("Error copying sparse matrix; the collection \
of non-zero entries in the source matrix is not a subset of the non-zero \
entries present in the destination");
}
