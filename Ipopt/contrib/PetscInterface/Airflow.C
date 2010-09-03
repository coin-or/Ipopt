// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors:  Johannes Huber, Andreas Waechter     IBM   2010-09-03
#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>

#include "libmesh.h"
#include "mesh.h"
#include "elem.h"
#include "boundary_info.h"

#include "IpIpoptApplication.hpp"
#include "IpLibMeshPDENLP.hpp"
#include "IpParTNLPAdapter.hpp"

#include "petsc.h"

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
  {
    LibMeshInit init (argc, argv,MPI_COMM_WORLD);

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    LibMeshPDEBase* pLibMeshPDE = new LibMeshPDEBase;
    std::ifstream f("ProblemGeometry.dat",std::ios::in);
    pLibMeshPDE->InitProblemData(f);
    pLibMeshPDE->reinit();
    SmartPtr<ParTNLP> partnlp = new LibMeshPDENLP(*pLibMeshPDE,*(app->Jnlst()));
    SmartPtr<NLP> nlp = new ParTNLPAdapter(partnlp, ConstPtr(app->Jnlst()));
    ApplicationReturnStatus status;
    { // Find initial point:
      // get bounds
      pLibMeshPDE->simulation_mode_=true;
      AutoPtr<  NumericVector< lm_Number > > state_l = pLibMeshPDE->getStateVector().clone();
      AutoPtr<  NumericVector< lm_Number > > state_u = pLibMeshPDE->getStateVector().clone();
      AutoPtr<  NumericVector< lm_Number > > control_l = pLibMeshPDE->getControlVector().clone();
      AutoPtr<  NumericVector< lm_Number > > control_u = pLibMeshPDE->getControlVector().clone();
      AutoPtr<  NumericVector< lm_Number > > aux_constr_l = pLibMeshPDE->getAuxConstrVector().clone();
      AutoPtr<  NumericVector< lm_Number > > aux_constr_u = pLibMeshPDE->getAuxConstrVector().clone();
      pLibMeshPDE->getControlVector() = 1.0;
      pLibMeshPDE->get_bounds(*state_l,*state_u,*control_l,*control_u,*aux_constr_l,*aux_constr_u);
      lm_Number MaxLowBd = aux_constr_l->max();
      // solve simulation (fixed controls)
      // Intialize the IpoptApplication and process the options
      status = app->Initialize();
      if (status != Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        return (int) status;
      }
      // Ask Ipopt to solve the problem
      status = app->OptimizeNLP(nlp);
      if (status == Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
      }
      else {
        printf("\n\n*** The problem FAILED!\n");
	return (int)status;
      }

      libMesh::NumericVector<lm_Number>* aux_constr;
      pLibMeshPDE->calcAux_constr(aux_constr);
      lm_Number MinAuxConstr = sqrt(aux_constr->min());
      //lm_Number fact = MaxLowBd/MinAuxConstr;
      lm_Number fact = 1.0;
      if(fact>0.9)
      { 
        pLibMeshPDE->getControlVector() *= fact;
        pLibMeshPDE->getStateVector() *= fact;
        // Now the inital point should be feasible
        #if 0
        // solve a second time to check 
        status = app->Initialize();
        if (status != Solve_Succeeded) {
          printf("\n\n*** Error during initialization!\n");
          return (int) status;
        }
        // Ask Ipopt to solve the problem
        status = app->OptimizeNLP(nlp);
        if (status == Solve_Succeeded) {
          printf("\n\n*** The problem solved!\n");
        }
        else {
          printf("\n\n*** The problem FAILED!\n");
        }
        #endif
      }

      // solve a second time to check 
      #if 0
      status = app->Initialize();
      if (status != Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        return (int) status;
      }
      // Ask Ipopt to solve the problem
      status = app->OptimizeNLP(nlp);
      if (status == Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
      }
      else {
        printf("\n\n*** The problem FAILED!\n");
	return (int)status;
      }
      #endif
    }
    // solve optimization problem
    pLibMeshPDE->simulation_mode_=false;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
      printf("\n\n*** Error during initialization!\n");
      return (int) status;
    }
    // Ask Ipopt to solve the problem
    status = app->OptimizeNLP(nlp);
    if (status == Solve_Succeeded) {
      printf("\n\n*** The problem solved!\n");
    }
    else {
      printf("\n\n*** The problem FAILED!\n");
    }
    delete pLibMeshPDE;
  }

  return 0;
}
