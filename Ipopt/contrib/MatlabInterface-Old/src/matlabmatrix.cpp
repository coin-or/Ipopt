// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabmatrix.h"
#include "matlabexception.h"

// Function definitions.
// -----------------------------------------------------------------
double* getMatlabMatrixDouble (const mxArray* ptr) {
  if (mxGetNumberOfDimensions(ptr) != 2)
    throw MatlabException("Matlab array must be a matrix");
  if (!mxIsDouble(ptr))
    throw MatlabException("Matlab array must be of type double");
  return mxGetPr(ptr);
}

double* createMatlabMatrixDouble (mxArray*& ptr, int height, int width) {
  ptr = mxCreateDoubleMatrix(height,width,mxREAL);
  return mxGetPr(ptr);
}
  
// Function definitions for class Matrix.
// -----------------------------------------------------------------
Matrix::Matrix (int height, int width)
  : Array<double>(height*width) { 
  h = height;
  w = width;
}

Matrix::Matrix (double* data, int height, int width) 
  : Array<double>(data,height*width) {
  h = height;
  w = width;
  }

Matrix::Matrix (const mxArray* ptr) 
  : Array<double>(getMatlabMatrixDouble(ptr),mxGetNumberOfElements(ptr)) {
  h = mxGetM(ptr);
  w = mxGetN(ptr);
}

Matrix::Matrix (mxArray*& ptr, int height, int width) 
  : Array<double>(createMatlabMatrixDouble(ptr,height,width),height*width) {
  h = height;
  w = width;
}

Matrix::Matrix (const Matrix& source)
  : Array<double>(source) {
  h = source.h;
  w = source.w;
}

Matrix& Matrix::operator= (const Matrix& source) {
  inject(source);
  return *this;
}

bool Matrix::operator== (const Matrix& X) const {
  return (h == X.h) && (w == X.w);
}

double& Matrix::entry (int r, int c) {
  return elems[h*c + r];
}

double Matrix::entry (int r, int c) const {
  return elems[h*c + r];
}
