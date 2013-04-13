// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-06

#ifndef __ASSCHURDRIVER_HPP__
#define __ASSCHURDRIVER_HPP__

#include "SensSchurData.hpp"
#include "SensPCalculator.hpp"
#include "IpVector.hpp"
#include "IpIteratesVector.hpp"

namespace Ipopt
{

  class SchurDriver : public AlgorithmStrategyObject
  {
    /** This class is the interface for any class that deals with the Schur matrix
     *  from the point when it is constructed by the PCalculator to the solution
     *  against one vector. Specific implementations may also incorporate the
     *  treatment of adding rows/cols (like QPSchur).
     *
     *  The computations done by this class are
     *  1. Solve \f$S \Delta\nu = r_s\f$
     *  2. Solve \f$K\Delta s = ... - \Delta nu\f$ (really?)*/

  public:

    SchurDriver(SmartPtr<PCalculator> pcalc,
		SmartPtr<SchurData> data_B)
      :
      pcalc_(pcalc),
      data_B_(data_B)
    {
    }

    virtual ~SchurDriver()
    {
    }

    /** Overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix)
    {
      return true;
    }

    /** Const accessor methods to the SchurData for for the derived classes */
    virtual SmartPtr<const SchurData> data_A() const
    {
      return pcalc_->data_A();
    }

    virtual SmartPtr<const SchurData> data_B() const
    {
      return ConstPtr(data_B_);
    }

    virtual SmartPtr<SchurData> data_A_nonconst()
    {
      return pcalc_->data_A_nonconst();
    }

    virtual SmartPtr<SchurData> data_B_nonconst()
    {
      return data_B_;
    }

    virtual SmartPtr<const PCalculator> pcalc() const
    {
      return ConstPtr(pcalc_);
    }

    virtual SmartPtr<PCalculator> pcalc_nonconst()
    {
      return pcalc_;
    }

    /** Sets the Data for which this SchurMatrix will be built. */

    /** Creates the SchurMatrix from B and P */
    virtual bool SchurBuild() =0;

    /** Calls the factorization routine for the SchurMatrix */
    virtual bool SchurFactorize() =0;

    /** Performs a backsolve on S and K */
    virtual bool SchurSolve(SmartPtr<IteratesVector> lhs,
			    SmartPtr<const IteratesVector> rhs,
			    SmartPtr<Vector> delta_u,
			    SmartPtr<IteratesVector> sol=NULL)=0; // the vector K^(-1)*r_s which usually should have been computed before.


    /** Performs a backsolve on S and K; calls the latter with sol=K^(-1)*r_s=0
	virtual bool SchurSolve(SmartPtr<IteratesVector> lhs,
	SmartPtr<const IteratesVector> rhs,
	SmartPtr<Vector> delta_u) =0;
    */
  private:
    SchurDriver()
    {
    }

    SmartPtr<PCalculator> pcalc_;

    SmartPtr<SchurData> data_B_;
  };
}

#endif
