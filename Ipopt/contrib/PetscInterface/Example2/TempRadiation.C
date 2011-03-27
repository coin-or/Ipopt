// Copyright (C) 2010, 2011 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Johannes Huber, Andreas Waechter     IBM   2010-09-03
#include "mpi.h"
#include <iostream>
#include <stdexcept>

#include "IpIpoptApplication.hpp"
#include "IpPetscPDENLP.hpp"
#include "IpParTNLPAdapter.hpp"
#include "IpPetscPDETempRadiation.hpp"

using namespace libMesh;
using namespace Ipopt;

typedef libMesh::Number lm_Number;

int GetProcID()
{
  int proc_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  return proc_id;
}

// Begin the main program.
int main (int argc, char** argv)
{
  LibMeshInit init (argc, argv,MPI_COMM_WORLD);
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  PetscPDETempRadiation* pPDE = new PetscPDETempRadiation();
  pPDE->Init("Problem.dat");
  SmartPtr<ParTNLP> partnlp = new PetscPDENLP(*pPDE,*(app->Jnlst()));
  SmartPtr<NLP> nlp = new ParTNLPAdapter(partnlp, ConstPtr(app->Jnlst()));
  ApplicationReturnStatus status;

  Vec State,Control;
  pPDE->getControlVector(Control);
  VecSet(Control, 1.0);
  pPDE->getStateVector(State);
  VecSet(State, 0.5);

  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    delete pPDE;
    return (int) status;
  }
  // Ask Ipopt to solve the problem
  status = app->OptimizeNLP(nlp);
  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
  }
  else {
    printf("\n\n*** The problem FAILED!\n");
    delete pPDE;
    return (int)status;
  }
  delete pPDE;

  return 0;
}
