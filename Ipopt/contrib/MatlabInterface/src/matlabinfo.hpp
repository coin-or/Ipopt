// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 25, 2008

#ifndef INCLUDE_MATLABINFO
#define INCLUDE_MATLABINFO

#include "mex.h"
#include "IpIpoptApplication.hpp"

using Ipopt::ApplicationReturnStatus;

// Class MatlabInfo.
// -----------------------------------------------------------------
// An object of this class stores all the information we will pass 
// back to MATLAB upon termination of IPOPT.
class MatlabInfo {
public:

  // Create a new info object and store the information in a MATLAB
  // array. The input pointer will point to the the newly created
  // MATLAB array. Since the user has an opportunity to modify the
  // MATLAB array pointed to by "ptr", we do not destroy the array
  // when the MatlabInfo object is destroyed. It is up to the user to
  // do that.
  explicit MatlabInfo (mxArray*& ptr);

  // The destructor.
  ~MatlabInfo() { };

  // Access and modify the exit status and solution statistics.
  ApplicationReturnStatus getExitStatus () const;
  void                    setExitStatus (ApplicationReturnStatus status);
  void                    setIterationCount (int iter);
  void                    setFuncEvals(int obj, int con, int grad, int jac, int hess);
  void                    setCpuTime (double cpu);

  // Access and modify the Lagrange multipliers.
  const double* getmultlb     () const;
  const double* getmultub     () const;
  const double* getmultconstr () const;
  void          setmultlb     (int n, const double* zl);
  void          setmultub     (int n, const double* zu);
  void          setmultconstr (int m, const double* lambda);

protected:
  mxArray* ptr;  // All the information is stored in a MATLAB array.
};

#endif
