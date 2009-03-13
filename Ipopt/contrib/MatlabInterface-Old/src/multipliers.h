// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 22, 2007

#ifndef INCLUDE_MULTIPLIERS
#define INCLUDE_MULTIPLIERS

#include "mex.h"
#include "matlabmatrix.h"
#include "arrayofmatrices.h"

// Class Multipliers.
// -----------------------------------------------------------------
// This class reserves storage for the Lagrange multipliers associated
// with a constrained, nonlinear program. There are three types of
// Lagrange multipliers: those associated with the upper bounds on the
// primal variables, those associated with the lower bounds, and those
// associated with the equality and inequality constraints. Of
// particular interest is the fact that one of the constructors
// accesses the information from a MATLAB structure. The structure
// must be created with the following fields: zl, zu, lambda. There is
// another constructor that creates a new MATLAB structure with those
// fields. See the descriptions of the constructors below for more
// information.
class Multipliers {
public:

  // Read the values of the multipliers from the specified MATLAB
  // structure. See the comments above for more information as to the
  // form the MATLAB structure is expected to take.
  explicit Multipliers (const mxArray*& ptr);

  // Create a set of multipliers for n variables and m constraints. It
  // creates a MATLAB struct array as a side effect.
  Multipliers (mxArray*& ptr, int n, int m);

  // The copy constructor makes a shallow copy of the data.
  Multipliers (const Multipliers& source);

  // The destructor.
  ~Multipliers();

  // Access the multipliers.
  const Matrix& lowerbounds() const { return *zl;     };
  const Matrix& upperbounds() const { return *zu;     };
  const Matrix& constraints() const { return *lambda; };
        Matrix& lowerbounds()       { return *zl;     };
        Matrix& upperbounds()       { return *zu;     };
        Matrix& constraints()       { return *lambda; };

protected:
  Matrix* zl;     // The Lagrange multipliers corresponding to the
	          // lower bounds on the optimization variables.
  Matrix* zu;     // The Lagrange multipliers corresponding to the
	          // upper bounds on the optimization variables.
  Matrix* lambda; // The Lagrange multipliers associated with the
		  // equality and inequality constraints.
};

#endif
