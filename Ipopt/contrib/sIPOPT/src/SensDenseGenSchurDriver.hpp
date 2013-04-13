// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-11-19

#ifndef __ASIFTSCHURDRIVER_HPP__
#define __ASIFTSCHURDRIVER_HPP__

#include "SensSchurDriver.hpp"
#include "SensBacksolver.hpp"
#include "IpDenseGenMatrix.hpp"

namespace Ipopt
{

  class DenseGenSchurDriver: public SchurDriver
  {

  public:

    DenseGenSchurDriver(SmartPtr<SensBacksolver> backsolver,
			SmartPtr<PCalculator> pcalc,
			SmartPtr<SchurData> data_B);

    virtual ~DenseGenSchurDriver();

    /** Creates the SchurMatrix from B and P */
    virtual bool SchurBuild();

    /** Calls the factorization routine for the SchurMatrix */
    virtual bool SchurFactorize();

    /** Performs a backsolve on S and : Solves the system
     *
     *  \f$\left[\begin{array}{c|c}
     *  K & E\\\hline
     *  E^T & 0
     *  \end{array}
     *  \right]
     *  \left[\begin{array}{c}x\\y\end{array}\right] =
     *  \left[\begin{array}{c}f\\g\end{array}\right]\f$
     *
     *  y will be stored in g at exit.
     *  Kf should hold
     *
     *  \f$K^{-1}f\f$
     *
     *  if it has been computed previously. If it is not available, just
     *  pass in Kf=NULL and it will be computed internally.
     */
    virtual bool SchurSolve(SmartPtr<IteratesVector> x,
			    SmartPtr<const IteratesVector> f,
			    SmartPtr<Vector> g,
			    SmartPtr<IteratesVector> Kf=NULL);

    /** DEPRECATED Performs a backsolve on S and K
	virtual bool SchurSolve(SmartPtr<IteratesVector> lhs,
	SmartPtr<const IteratesVector> rhs,
	SmartPtr<Vector> delta_u);
    */
  private:
    SmartPtr<SchurData> ift_data_;
    SmartPtr<SensBacksolver> backsolver_;
    SmartPtr<DenseGenMatrix> S_;

  };
}

#endif
