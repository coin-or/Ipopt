#ifndef INCLUDE_MATLABINFO
#define INCLUDE_MATLABINFO

#include "mex.h"
#include "coin/IpIpoptApplication.hpp"

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

  // Access and modify the exit status.
  ApplicationReturnStatus getExitStatus () const;
  void                    setExitStatus (ApplicationReturnStatus status);

  // Access and modify the Lagrange multipliers.
  // TO DO.

protected:
  mxArray* ptr;  // All the information is stored in a MATLAB array.
};

#endif
