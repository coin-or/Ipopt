// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMuUpdate.hpp,v 1.1.1.1 2004/10/21 01:03:11 andreasw Exp $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpMuUpdate.hpp"

namespace Ipopt {

  DefineIpoptType(MuUpdate);

  MuUpdate::MuUpdate()
  {
    std::cout << "Bogus output" << std::endl;
  }

  void MuUpdate::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
     roptions->AddLowerBoundedNumberOption("kappa_epsilon", "???",
 					  0.0, true, 10.0);
     roptions->AddBoundedNumberOption("kappa_mu", "???",
 				     0.0, true, 1.0, true, 0.2);
     roptions->AddBoundedNumberOption("theta_mu", "???",
 				     1.0, true, 2.0, true, 1.5);
     roptions->AddBoundedNumberOption("tau_min", "???",
 				     0.0, true, 1.0, true, 0.99);
  }

} // namespace Ipopt
