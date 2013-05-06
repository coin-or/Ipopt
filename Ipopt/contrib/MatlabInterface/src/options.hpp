// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

#ifndef INCLUDE_OPTIONS
#define INCLUDE_OPTIONS

#include "mex.h"
#include "iterate.hpp"
#include "ipoptoptions.hpp"

// Class Options.
// -----------------------------------------------------------------
// This class processes the options input from MATLAB.
class Options {
public:

  // The constructor expects as input a point to a MATLAB array, in
  // particular a structure array with the appropriate fields. Note
  // that the Options object does *not* possess an independent copy of
  // some of the MATLAB data (such as the auxiliary data).
  Options (const Iterate& x, Ipopt::IpoptApplication& app, 
	   const mxArray* ptr);
  
  // The destructor.
  ~Options();

  // Get the number of variables and the number of constraints.
  friend int numvars        (const Options& options) { return options.n; };
  friend int numconstraints (const Options& options) { return options.m; };

  // Access the lower and upper bounds on the variables and constraints. 
  const double* lowerbounds () const { return lb; };
  const double* upperbounds () const { return ub; };
  const double* constraintlb() const { return cl; };
  const double* constraintub() const { return cu; };

  // Access the IPOPT options object.
  const IpoptOptions ipoptOptions() const { return ipopt; };

  // Access the Lagrange multpliers.
  const double* multlb    () const { return zl;     };
  const double* multub    () const { return zu;     };
  const double* multconstr() const { return lambda; };

protected:
  int            n;       // The number of optimization variables.
  int            m;       // The number of constraints.
  double*        lb;      // Lower bounds on the variables.
  double*        ub;      // Upper bounds on the variables.
  double*        cl;      // Lower bounds on constraints.
  double*        cu;      // Upper bounds on constraints.
  double*        zl;      // Lagrange multipliers for lower bounds.
  double*        zu;      // Lagrange multipliers for upper bounds.
  double*        lambda;  // Lagrange multipliers for constraints.
  IpoptOptions   ipopt;   // The IPOPT options.

  // These are helper functions used by the class constructor.
  static double* loadLowerBounds      (int n, const mxArray* ptr, 
				       double neginfty);
  static double* loadUpperBounds      (int n, const mxArray* ptr, 
				       double posinfty);
  static int     loadConstraintBounds (const mxArray* ptr, double*& cl, 
				       double*& cu, double neginfty,
				       double posinfty);
  static void    loadMultipliers      (int n, int m, const mxArray* ptr, 
				       double*& zl, double*& zu, 
				       double*& lambda);
};

#endif
