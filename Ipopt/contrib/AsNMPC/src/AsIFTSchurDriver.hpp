// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Date   : 2009-11-19

#ifndef __ASIFTSCHURDRIVER_HPP__
#define __ASIFTSCHURDRIVER_HPP__

#include "AsSchurDriver.hpp"
#include "AsAsBacksolver.hpp"

namespace Ipopt
{

  class IFTSchurDriver: public SchurDriver
  {

  public: 
    
    IFTSchurDriver(SmartPtr<AsBacksolver> backsolver,
		   SmartPtr<SchurData> data_B);

    virtual ~IFTSchurDriver();

    /** Creates the SchurMatrix from B and P */
    virtual bool SchurBuild();

    /** Calls the factorization routine for the SchurMatrix */
    virtual bool SchurFactorize();

    /** Performs a backsolve on S and K */
    virtual bool SchurSolve(SmartPtr<IteratesVector> lhs, 
			    SmartPtr<const IteratesVector> rhs,
			    SmartPtr<IteratesVector> sol,
			    SmartPtr<Vector> delta_u);

    /** Performs a backsolve on S and K */
    virtual bool SchurSolve(SmartPtr<IteratesVector> lhs, 
			    SmartPtr<const IteratesVector> rhs,
			    SmartPtr<Vector> delta_u);

  private:
    
    SmartPtr<AsBacksolver> backsolver_;

  };
}

#endif
