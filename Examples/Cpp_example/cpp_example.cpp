// Copyright (C) 2004, 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "IpIpoptApplication.hpp"
#include "MyNLP.hpp"

#include "IpDenseSymMatrix.hpp"

using namespace Ipopt;

int main(int argv, char* argc[])
{
  // Create an instance of your nlp...
  SmartPtr<TNLP> mynlp = new MyNLP();

  // Create an instance of the IpoptApplication
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  // Testing new stuff
  SmartPtr<const Journalist> jnlst = ConstPtr(app->Jnlst());

  Index dim = 3;
  SmartPtr<DenseSymMatrixSpace> dsm_space = new DenseSymMatrixSpace(dim);
  SmartPtr<DenseSymMatrix> dsm = dsm_space->MakeNewDenseSymMatrix();
  Number* dsm_values = dsm->Values();
  for (Index j=0; j<dim; j++) {
    for (Index i=j; i<dim; i++) {
      dsm_values[i+dim*j] = (double)(i+dim*j)+2;
    }
  }
  /*
  for (Index j=0; j<dim; j++) {
    for (Index i=j; i<dim; i++) {
      if (i==j) {
	dsm_values[i+dim*j] = i+1;
      }
      else {
	dsm_values[i+dim*j] = 0.;
      }
    }
  }
  */
  dsm->Print(jnlst, J_SUMMARY, J_MAIN, "dsm");

  SmartPtr<DenseGenMatrixSpace> dgm1_space = new DenseGenMatrixSpace(dim,dim);
  SmartPtr<DenseGenMatrix> dgm1 = dgm1_space->MakeNewDenseGenMatrix();
  SmartPtr<DenseVectorSpace> dv1_space = new DenseVectorSpace(dim);
  SmartPtr<DenseVector> dv1 = dv1_space->MakeNewDenseVector();
  dsm->ComputeEigenVectors(*dgm1, *dv1);
  dv1->Print(jnlst, J_SUMMARY, J_MAIN, "dv1");
  dgm1->Print(jnlst, J_SUMMARY, J_MAIN, "dgm1");

  dsm->ComputeCholeskyFactor(*dgm1);
  dgm1->Print(jnlst, J_SUMMARY, J_MAIN, "dgm1");

  SmartPtr<DenseSymMatrix> dsm2 = dsm_space->MakeNewDenseSymMatrix();
  dsm2->HighRankUpdate(false, 1., *dgm1, 0.);
  dsm2->Print(jnlst, J_SUMMARY, J_MAIN, "dsm2");
  

  Number* dv1_values = dv1->Values();
  for (Index i=0; i<dim; i++) {
    dv1_values[i] = 3+i;
  }
  dv1->Print(jnlst, J_SUMMARY, J_MAIN, "dv1 - rhs");
  dgm1->CholeskySolveVector(*dv1);
  dv1->Print(jnlst, J_SUMMARY, J_MAIN, "dv1 - sol");
  SmartPtr<DenseVector> dv2 = dv1_space->MakeNewDenseVector();
  dsm->MultVector(1., *dv1, 0., *dv2);
  dv2->Print(jnlst, J_SUMMARY, J_MAIN, "dv2 - rhs");

  ApplicationReturnStatus status = app->OptimizeTNLP(mynlp);

  if (status == Solve_Succeeded) {
    // Retrieve some statistics about the solve
    Index iter_count = app->Statistics()->IterationCount();
    printf("\n\n*** The problem solved in %d iterations!\n", iter_count);

    Number final_obj = app->Statistics()->FinalObjective();
    printf("\n\n*** The final value of the objective function is %e.\n", final_obj);
  }

  return (int) status;
}
