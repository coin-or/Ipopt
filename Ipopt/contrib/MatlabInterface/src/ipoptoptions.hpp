// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 15, 2008

#ifndef INCLUDE_IPOPTOPTIONS
#define INCLUDE_IPOPTOPTIONS

#include "mex.h"
#include "IpIpoptApplication.hpp"

// Class IpoptOptions.
// -----------------------------------------------------------------
// This class processes the IPOPT options as specified by a user in the 
// MATLAB environment. 
class IpoptOptions {
public:

  // The constructor accepts as input an IPOPT application object and
  // a MATLAB array. The latter input must be a structure array, with
  // field names corresponding to the names of options in IPOPT.
  IpoptOptions (Ipopt::IpoptApplication& app, const mxArray* ptr);

  // The destructor.
  ~IpoptOptions() { };

  // The first function returns true if and only if the user has
  // specified a quasi-Newton approximation to the Hessian instead of
  // the exact Hessian. The second function returns true if and only
  // if the user has activated the derivative checker. The third
  // function returns true if and only if a user-specified scaling of
  // the problem is activated. The fourth function returns the print
  // level for the IPOPT console. The remaining two functions return
  // the floating-point value for positive and negative infinity,
  // respectively.
  bool   useQuasiNewton () const;
  bool   useDerivChecker() const;
  bool   userScaling    () const;
  int    printLevel     () const;
  double getPosInfty    () const;
  double getNegInfty    () const;

protected:
  Ipopt::IpoptApplication& app;  // The IPOPT application object.

  // These three functions are used by the class constructor.
  void setOption        (const char* label, const mxArray* ptr);
  void setStringOption  (const char* label, const mxArray* ptr);
  void setIntegerOption (const char* label, const mxArray* ptr);
  void setNumberOption  (const char* label, const mxArray* ptr);
};

#endif
