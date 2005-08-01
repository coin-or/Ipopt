// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptApplication.hpp"
#include "MyNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

using namespace Ipopt;

int main(int argv, char* argc[])
{
  // Create an instance of your nlp...
  SmartPtr<TNLP> mynlp = new MyNLP();

  // Create an instance of the IpoptApplication
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  // To illustrate the use of ip_data and ip_cq, we will use the
  //  version of OptimzeTNLP that returns pointers to these objects
  // Usually, the standard IpoptApplication output is sufficient
  //  and you would only need to call OptimizeTNLP(mynlp);
  SmartPtr<IpoptData> ip_data = NULL;
  SmartPtr<IpoptCalculatedQuantities> ip_cq = NULL;
  ApplicationReturnStatus status = app->OptimizeTNLP(mynlp, ip_data, ip_cq);

  if (status == Solve_Succeeded) {
    // Retrieve some information from ip_data
    printf("\n\n*** The problem solved in %d iterations!\n", ip_data->iter_count());
    // Retrieve some information from ip_cq
    printf("\n\n*** The current value of the objective function is %g.\n", ip_cq->curr_f());
  }

  return (int) status;
}
