// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-11-19

#include "AsIFTSchurDriver.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  IFTSchurDriver::IFTSchurDriver(SmartPtr<AsBacksolver> backsolver,
				 SmartPtr<SchurData> data_B)
    :
    SchurDriver(NULL,data_B),
    backsolver_(backsolver)
  {
    DBG_START_METH("IFTSchurDriver::IFTSchurDriver", dbg_verbosity);
  }

  IFTSchurDriver::~IFTSchurDriver()
  {
    DBG_START_METH("IFTSchurDriver::~IFTSchurDriver", dbg_verbosity);
  }

  bool IFTSchurDriver::SchurBuild()
  {
    DBG_START_METH("IFTSchurDriver::SchurBuild", dbg_verbosity);
    
    return true;
  }

  bool IFTSchurDriver::SchurFactorize()
  {
    DBG_START_METH("IFTSchurDriver::SchurFactorize", dbg_verbosity);

    return true;
  }

  bool IFTSchurDriver::SchurSolve(SmartPtr<IteratesVector> lhs, // new left hand side will be stored here 
				  SmartPtr<const IteratesVector> rhs, // rhs r_s
				  SmartPtr<IteratesVector> sol, // the vector K^(-1)*r_s which usually should have been computed before.
				  SmartPtr<Vector> delta_u)     // should be (u_p - u_0) WATCH OUT FOR THE SIGN! I like it this way, so that u_0+delta_u = u_p, but victor always used it the other way round, so be careful. At the end, delta_nu is saved in here.
  {
    DBG_START_METH("IFTSchurDriver::SchurSolve", dbg_verbosity);
    
    SmartPtr<IteratesVector> delta_rhs = rhs->MakeNewIteratesVector();
    data_B()->TransMultiply(*delta_u, *delta_rhs);
    delta_rhs->Scal(1.0);
    delta_rhs->Axpy(-1.0, *rhs); 
    delta_rhs->Print(Jnlst(),J_VECTOR,J_USER1,"delta_rhs");

    return backsolver_->Solve(lhs, ConstPtr(delta_rhs));
  }

  bool IFTSchurDriver::SchurSolve(SmartPtr<IteratesVector> lhs, // new left hand side will be stored here 
				  SmartPtr<const IteratesVector> rhs, // rhs r_s
				  SmartPtr<Vector> delta_u)     // should be (u_p - u_0) WATCH OUT FOR THE SIGN! I like it this way, so that u_0+delta_u = u_p, but victor always used it the other way round, so be careful.
  {
    DBG_START_METH("IFTSchurDriver::SchurSolve(3inputs)", dbg_verbosity);
    SmartPtr<IteratesVector> sol = rhs->MakeNewIteratesVector();
    /*    retval = backsolver_->Solve(sol, rhs);
    if (!retval) {
      return false;
      }*/
    return SchurSolve(lhs,rhs,sol,delta_u);
  }
}
