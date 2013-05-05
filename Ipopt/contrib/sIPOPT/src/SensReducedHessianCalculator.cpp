// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-08-01

#include "SensReducedHessianCalculator.hpp"
#include "IpDenseGenMatrix.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  ReducedHessianCalculator::ReducedHessianCalculator(SmartPtr<SchurData> hess_data,
						     SmartPtr<PCalculator> pcalc)
    :
    hess_data_(hess_data),
    pcalc_(pcalc)
  {
    DBG_START_METH("ReducedHessianCalculator::ReducedHessianCalculator", dbg_verbosity);
  }

  ReducedHessianCalculator::~ReducedHessianCalculator()
  {
    DBG_START_METH("ReducedHessianCalculator::~ReducedHessianCalculator", dbg_verbosity);
  }

  bool ReducedHessianCalculator::InitializeImpl(const OptionsList& options,
						const std::string& prefix)
  {
    DBG_START_METH("ReducedHessianCalculator::InitializeImpl", dbg_verbosity);
    options.GetBoolValue("rh_eigendecomp", compute_eigenvalues_, prefix);
    return true;
  }

  bool ReducedHessianCalculator::ComputeReducedHessian()
  {
    DBG_START_METH("ReducedHessianCalculator::ComputeReducedHessian", dbg_verbosity);

    Index dim_S = hess_data_->GetNRowsAdded();
    //SmartPtr<DenseGenMatrixSpace> S_space = new DenseGenMatrixSpace(dim_S, dim_S);
    //SmartPtr<DenseGenMatrix> S = new DenseGenMatrix(GetRawPtr(S_space));
    SmartPtr<Matrix> S;
    bool retval = pcalc_->GetSchurMatrix(GetRawPtr(hess_data_), S);

    SmartPtr<DenseSymMatrix> S_sym = dynamic_cast<DenseSymMatrix*>(GetRawPtr(S));
    if (!IsValid(S_sym)) {
      std::exception exc;
      throw (exc);
    }

    bool have_x_scaling, have_c_scaling, have_d_scaling;
    have_x_scaling = IpNLP().NLP_scaling()->have_x_scaling();
    have_c_scaling = IpNLP().NLP_scaling()->have_c_scaling();
    have_d_scaling = IpNLP().NLP_scaling()->have_d_scaling();

    if (have_x_scaling || have_c_scaling || have_d_scaling) {
      Jnlst().Printf(J_WARNING, J_MAIN,
		     "\n"
		     "-------------------------------------------------------------------------------\n"
		     "                              *** WARNING ***\n"
		     "-------------------------------------------------------------------------------\n"
		     "You are using the reduced hessian feature with scaling of\n");
      if (have_x_scaling) {
	Jnlst().Printf(J_WARNING, J_MAIN, "*** variables\n");
      }
      if (have_c_scaling) {
	Jnlst().Printf(J_WARNING, J_MAIN, "*** equality constraints\n");
      }
      if (have_d_scaling) {
	Jnlst().Printf(J_WARNING, J_MAIN, "*** inequality constraints\n");
      }
      Jnlst().Printf(J_WARNING, J_MAIN,
		     "enabled.\n"
		     "A correct unscaled solution of the reduced hessian cannot be guaranteed in this\n"
		     "case. Please consider rerunning with scaling turned off.\n"
		     "-------------------------------------------------------------------------------\n\n");

    }

    // Unscale by objective factor and multiply by (-1)
    Number obj_scal = IpNLP().NLP_scaling()->apply_obj_scaling(1.0);
    DBG_PRINT((dbg_verbosity, "Objective scaling = %f\n", obj_scal));
    Number* s_val = S_sym->Values();
    for (Index k=0; k<(S->NRows())*(S->NCols()); ++k) {
      s_val[k] *= -obj_scal;
    }

    S->Print(Jnlst(),J_INSUPPRESSIBLE,J_USER1,"RedHessian unscaled");

    if ( compute_eigenvalues_ ) {
      SmartPtr<DenseGenMatrixSpace> eigenvectorspace = new DenseGenMatrixSpace(dim_S, dim_S);
      SmartPtr<DenseGenMatrix> eigenvectors = new DenseGenMatrix(GetRawPtr(eigenvectorspace));
      SmartPtr<DenseVectorSpace> eigenvaluesspace = new DenseVectorSpace(dim_S);
      SmartPtr<DenseVector> eigenvalues = new DenseVector(GetRawPtr(eigenvaluesspace));

      eigenvectors->ComputeEigenVectors(*S_sym, *eigenvalues);
      eigenvalues->Print(Jnlst(),J_INSUPPRESSIBLE,J_USER1,"Eigenvalues of reduced hessian matrix");
    }

    return retval;
  }


}
