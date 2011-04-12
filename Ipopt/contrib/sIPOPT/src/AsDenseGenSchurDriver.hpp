// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-08

#ifndef __ASDENSEGENSCHURDRIVER_HPP__
#define __ASDENSEGENSCHURDRIVER_HPP__

#include "AsSchurDriver.hpp"
#include "AsAsBacksolver.hpp"
#include "IpDenseGenMatrix.hpp"

namespace Ipopt
{

  class DenseGenSchurDriver : public SchurDriver
  {
    /** This is the most basic of all possible implementations of the
     *  SchurDriver interface. It uses a simple backsolver as an interface
     *  to the KKT solver, a DenseGenMatrix as Schurmatrix, and LU factorization
     *  from LAPACK for the DenseGenMatrix (DGETRF) */

  public: 
    
    DenseGenSchurDriver(SmartPtr<AsBacksolver> backsolver,
			SmartPtr<PCalculator> pcalc,
			SmartPtr<SchurData> data_B);

    virtual ~DenseGenSchurDriver();

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

    SmartPtr<DenseGenMatrix> S_;
    
  };
}

#endif
