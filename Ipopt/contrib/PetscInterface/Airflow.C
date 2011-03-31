// Copyright (C) 2010, 2011 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Johannes Huber, Andreas Waechter     IBM   2010-09-03
#include "mpi.h"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <stdexcept>

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
    std::ifstream f("Problem.dat",std::ios::in);
    if (f.fail()) {
      std::cout << "Can't open file Problem.dat" << std::endl;
      exit(1);
    }
    try {
      pLibMeshPDE->InitProblemData(f);
    }
    catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      delete pLibMeshPDE;
      exit(1);
    }
    pLibMeshPDE->reinit();
    SmartPtr<ParTNLP> partnlp = new LibMeshPDENLP(*pLibMeshPDE,*(app->Jnlst()));
    SmartPtr<NLP> nlp = new ParTNLPAdapter(partnlp, ConstPtr(app->Jnlst()));
    ApplicationReturnStatus status;

    {
      // Find initial point:
      // get bounds
      pLibMeshPDE->simulation_mode_=true;
      double initContr = 20.;
      printf("Simulation with initial AC setting %e\n", initContr);
      pLibMeshPDE->getControlVector() = initContr;
      // solve simulation (fixed controls)
      // Intialize the IpoptApplication and process the options
      status = app->Initialize("simu.opt");
      if (status != Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        delete pLibMeshPDE;
        return (int) status;
      }
      // Ask Ipopt to solve the problem
      status = app->OptimizeNLP(nlp);
      if (status == Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
      }
      else {
        printf("\n\n*** The problem FAILED!\n");
        delete pLibMeshPDE;
        return (int)status;
      }

//#define SIMULATION_CONVERGENCE
#ifdef SIMULATION_CONVERGENCE
      std::cout << "************************************************" << std::endl;
      std::cout << "*                                              *" << std::endl;
      std::cout << "*               Refine and solve               *" << std::endl;
      std::cout << "*                                              *" << std::endl;
      std::cout << "************************************************" << std::endl;
      
      pLibMeshPDE->RefineMesh(1);
      AutoPtr<  NumericVector< lm_Number > > state1= pLibMeshPDE->getStateVector().clone();
      status = app->Initialize("simu.opt");
      if (status != Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        delete pLibMeshPDE;
        return (int) status;
      }
      // Ask Ipopt to solve the problem
      status = app->OptimizeNLP(nlp);
      if (status == Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
      }
      else {
        printf("\n\n*** The problem FAILED!\n");
        delete pLibMeshPDE;
        return (int)status;
      }
      
      std::cout << "L2-Diff: " << pLibMeshPDE->CalcL2Diff(state1.get()) << std::endl;
      exit(1);
#endif

// #define SCALE_INIT
#ifdef SCALE_INIT
      // scale the simulation solution to satisfy inequality constraints
      pLibMeshPDE->simulation_mode_=false;
      AutoPtr<  NumericVector< lm_Number > > state_l = pLibMeshPDE->getStateVector().clone();
      AutoPtr<  NumericVector< lm_Number > > state_u = pLibMeshPDE->getStateVector().clone();
      AutoPtr<  NumericVector< lm_Number > > control_l = pLibMeshPDE->getControlVector().clone();
      AutoPtr<  NumericVector< lm_Number > > control_u = pLibMeshPDE->getControlVector().clone();
      AutoPtr<  NumericVector< lm_Number > > aux_constr_l = pLibMeshPDE->getAuxConstrVector().clone();
      AutoPtr<  NumericVector< lm_Number > > aux_constr_u = pLibMeshPDE->getAuxConstrVector().clone();
      pLibMeshPDE->get_bounds(*state_l,*state_u,*control_l,*control_u,*aux_constr_l,*aux_constr_u);

      libMesh::NumericVector<lm_Number>* aux_constr;
      pLibMeshPDE->calcAux_constr(aux_constr);
      double LocMaxFact(0.0), CurFact, GlobMaxFact(0.0);
      for( int iAuxConstr=0; iAuxConstr<aux_constr->size(); iAuxConstr++) {
        if(pLibMeshPDE->IsLocalIneqIdx(iAuxConstr))
        {
          CurFact = aux_constr_l->el(iAuxConstr)/aux_constr->el(iAuxConstr);
          if(CurFact>LocMaxFact)
            LocMaxFact = CurFact;
        }
      }
      MPI_Allreduce(&LocMaxFact,&GlobMaxFact,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      double fact = sqrt(GlobMaxFact);
      fact = 2*fact;
      printf("Scaling simulation solution with fact = %e\n", fact);
      if (1) //fact>0.9)
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
#endif

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
    
    bool warmstart = false;
    for (int i=0; i<1; i++) {
      // solve optimization problem
      pLibMeshPDE->simulation_mode_=false;
      app = IpoptApplicationFactory();
      status = app->Initialize();
      if (status != Solve_Succeeded) {
	    printf("\n\n*** Error during initialization!\n");
	    return (int) status;
    }

    if (warmstart) {
      //app->Options()->SetStringValue("warm_start_init_point", "yes");
      app->Options()->SetStringValue("bound_mult_init_method", "mu-based");
      app->Options()->SetNumericValue("bound_push", 1e-6);
      app->Options()->SetNumericValue("bound_frac", 1e-6);
      app->Options()->SetNumericValue("mu_init", 1e-4);
    }
      // Ask Ipopt to solve the problem
      status = app->OptimizeNLP(nlp);
      if (status == Solve_Succeeded) {
	printf("\n\n*** The problem solved!\n");
      }
      else {
	printf("\n\n*** The problem FAILED!\n");
      }

      // Try one set of refinements
      //pLibMeshPDE->RefineMesh(i);
      warmstart = true;
    }
    delete pLibMeshPDE;
  }

  return 0;
}
