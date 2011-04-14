// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-11-19

#include "SensDenseGenSchurDriver.hpp"
#include "SensIndexSchurData.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  DenseGenSchurDriver::DenseGenSchurDriver(SmartPtr<SensBacksolver> backsolver,
					   SmartPtr<PCalculator> pcalc,
					   SmartPtr<SchurData> data_B)
    :
    SchurDriver(pcalc,new IndexSchurData()),
    backsolver_(backsolver),
    S_(NULL)
  {
    DBG_START_METH("DenseGenSchurDriver::DenseGenSchurDriver", dbg_verbosity);
  }

  DenseGenSchurDriver::~DenseGenSchurDriver()
  {
    DBG_START_METH("DenseGenSchurDriver::~DenseGenSchurDriver", dbg_verbosity);
  }

  bool DenseGenSchurDriver::SchurBuild()
  {
    DBG_START_METH("DenseGenSchurDriver::SchurBuild", dbg_verbosity);
    /* This function is the very same as the one in DenseGenSchurDriver */
    bool retval = true;
    Index dim_S = 0;
    if (IsValid(data_B())) {
      dim_S = data_B()->GetNRowsAdded();
    }
    if (dim_S>0) {
      S_ = NULL;
      SmartPtr<DenseGenMatrixSpace> S_space = new DenseGenMatrixSpace(dim_S, dim_S);
      S_ = new DenseGenMatrix(GetRawPtr(S_space));
      SmartPtr<Matrix> S2 = GetRawPtr(S_);
      //retval = pcalc_nonconst()->GetSchurMatrix(GetRawPtr(data_B()), dynamic_cast<Matrix*>(GetRawPtr(S_)));
      retval = pcalc_nonconst()->GetSchurMatrix(data_B(), S2);
      S_->Print(Jnlst(),J_VECTOR,J_USER1,"S_");
    }
    return retval;
  }

  bool DenseGenSchurDriver::SchurFactorize()
  {
    DBG_START_METH("DenseGenSchurDriver::SchurFactorize", dbg_verbosity);
    /* This function is the very same as the one in DenseGenSchurDriver */
    bool retval;
    if (IsValid(S_)) {
      retval = S_->ComputeLUFactorInPlace();
      return retval;
    }
    return true;
  }

  /* The functions SchurSolve do IFT step, if S_==NULL, and DenseGenSchurDriver otherwise. */
  bool DenseGenSchurDriver::SchurSolve(SmartPtr<IteratesVector> lhs, // new left hand side will be stored here
				       SmartPtr<const IteratesVector> rhs, // rhs r_s
				       SmartPtr<Vector> delta_u,  // should be (u_p - u_0) WATCH OUT FOR THE SIGN! I like it this way, so that u_0+delta_u = u_p, but victor always used it the other way round, so be careful. At the end, delta_nu is saved in here.
				       SmartPtr<IteratesVector> sol) // the vector K^(-1)*r_s which usually should have been computed before.

  {
    DBG_START_METH("DenseGenSchurDriver::SchurSolve", dbg_verbosity);
    DBG_ASSERT(IsValid(S_));
    bool retval;

    // set up rhs of equation (3.48a)
    SmartPtr<Vector> delta_rhs = delta_u->MakeNew();
    data_B()->Multiply(*sol, *delta_rhs);
    delta_rhs->Print(Jnlst(),J_VECTOR,J_USER1,"delta_rhs");
    delta_rhs->Scal(-1.0);
    delta_rhs->Axpy(1.0, *delta_u);
    delta_rhs->Print(Jnlst(),J_VECTOR,J_USER1,"rhs 3.48a");

    // solve equation (3.48a) for delta_nu
    SmartPtr<DenseVector> delta_nu = dynamic_cast<DenseVector*>(GetRawPtr(delta_rhs))->MakeNewDenseVector();
    delta_nu->Copy(*delta_rhs);
    S_->LUSolveVector(*delta_nu); // why is LUSolveVector not bool??
    delta_nu->Print(Jnlst(),J_VECTOR,J_USER1,"delta_nu");

    // solve equation (3.48b) for lhs (=delta_s)
    SmartPtr<IteratesVector> new_rhs = lhs->MakeNewIteratesVector();
    data_A()->TransMultiply(*delta_nu, *new_rhs);
    new_rhs->Axpy(-1.0, *rhs);
    new_rhs->Scal(-1.0);
    new_rhs->Print(Jnlst(),J_VECTOR,J_USER1,"new_rhs");
    retval = backsolver_->Solve(lhs, ConstPtr(new_rhs));

    return retval;
  }
}
