// Copyright 2009, 2010 Hans Pirnay
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Date   : 2009-08-01

#include "AsReducedHessianCalculator.hpp"
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
    return true;
  }

  bool ReducedHessianCalculator::ComputeReducedHessian()
  {
    DBG_START_METH("ReducedHessianCalculator::ComputeReducedHessian", dbg_verbosity);

    Index dim_S = hess_data_->GetNRowsAdded();
    SmartPtr<DenseGenMatrixSpace> S_space = new DenseGenMatrixSpace(dim_S, dim_S);
    SmartPtr<DenseGenMatrix> S = new DenseGenMatrix(GetRawPtr(S_space));
    bool retval = pcalc_->GetSchurMatrix(*hess_data_, *S);
    
    S->Print(Jnlst(),J_VECTOR,J_USER1,"RedHessian");

    // columns have to be unscaled. All elements are from x, so we only have to determine the x-scaling
    SmartPtr<IteratesVector> scaling_factors;
    scaling_factors = IpData().curr()->MakeNewIteratesVector();
    scaling_factors->Set(1.0);

    if (IpNLP().NLP_scaling()->have_x_scaling()) {
      SmartPtr<const Vector> x_scaling = IpNLP().NLP_scaling()->unapply_vector_scaling_x(scaling_factors->x());

      scaling_factors->Set(0.0);
      scaling_factors->Set_x(*x_scaling);

      scaling_factors->Print(Jnlst(),J_VECTOR,J_USER1,"scal_fac");
    }

    SmartPtr<DenseVectorSpace> red_hess_space = new DenseVectorSpace(hess_data_->GetNRowsAdded());
    SmartPtr<DenseVector> matrix_scaling = new DenseVector(GetRawPtr(red_hess_space));

    hess_data_->Multiply(*scaling_factors, *matrix_scaling);
    matrix_scaling->Scal(-1.0);

    S->ScaleColumns(*matrix_scaling);
    
    S->Print(Jnlst(),J_INSUPPRESSIBLE,J_USER1,"RedHessian");
    
    return retval;
  }


}
