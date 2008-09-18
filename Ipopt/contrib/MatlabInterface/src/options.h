#ifndef INCLUDE_OPTIONS
#define INCLUDE_OPTIONS

#include "mex.h"
#include "iterate.h"
#include "ipoptoptions.h"

// Class Options.
// -----------------------------------------------------------------
// This class processes the options input from MATLAB.
class Options {
public:

  // The constructor expects as input a point to a MATLAB array, in
  // particular a structure array with the appropriate fields. Note
  // that the Options object does not possess an independent copy of
  // some of the MATLAB data (such as the auxiliary data).
  Options (const Iterate& x, Ipopt::IpoptApplication& app, 
	   const mxArray* ptr);
  
  // The destructor.
  ~Options();

  // This function returns true if and only if both the lower and
  // upper bounds on the contraints have been specified.
  bool isConstrained() const { return cl && cu; };

  // Get the number of variables and the number of constraints.
  friend int numvars        (const Options& options) { return options.n; };
  friend int numconstraints (const Options& options) { return options.m; };

  // Access the lower and upper bounds on the variables and constraints. 
  const double* lowerbounds () const { return lb; };
  const double* upperbounds () const { return ub; };
  const double* constraintlb() const { return cl; };
  const double* constraintub() const { return cu; };

  // Access the auxiliary data.
  const mxArray* getAuxData() const { return auxdata; };

  // Access the IPOPT options object.
  const IpoptOptions ipoptOptions() const { return ipopt; };

protected:
  int            n;       // The number of optimization variables.
  int            m;       // The number of constraints.
  double*        lb;      // Lower bounds on the variables.
  double*        ub;      // Upper bounds on the variables.
  double*        cl;      // Lower bounds on constraints.
  double*        cu;      // Upper bounds on constraints.
  const mxArray* auxdata; // MATLAB array containing the auxiliary data.
  IpoptOptions   ipopt;   // The IPOPT options.

  // These are helper functions used by the class constructor.
  static double* loadLowerBounds      (int n, const mxArray* ptr, 
				       double neginfty);
  static double* loadUpperBounds      (int n, const mxArray* ptr, 
				       double posinfty);
  static int     loadConstraintBounds (const mxArray* ptr, double*& cl, 
				       double*& cu, double neginfty,
				       double posinfty);
};

#endif
