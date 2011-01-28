// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-08-31

#ifndef __IPINEXACTSEARCHDIRCALC_HPP__
#define __IPINEXACTSEARCHDIRCALC_HPP__

#include "IpSearchDirCalculator.hpp"
#include "IpInexactCq.hpp"
#include "IpInexactNormalStepCalc.hpp"
#include "IpInexactPDSolver.hpp"

namespace Ipopt
{
  /** Implementation of the search direction calculator that computes
   *  the search direction using iterative linear solvers.  Those
   *  steps do not necessarily satisfy the linearized KKT conditions
   *  with high accuracy.
   */
  class InexactSearchDirCalculator : public SearchDirectionCalculator
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    InexactSearchDirCalculator(SmartPtr<InexactNormalStepCalculator> normal_step_calculator,
                               SmartPtr<InexactPDSolver> inexact_pd_solver);

    /** Default destructor */
    virtual ~InexactSearchDirCalculator();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the search direction.  In this version,
     *  we compute a normal and a tangential component, which are
     *  stored in the InexactData object.  The overall step is still
     *  stored in the IpoptData object. */
    virtual bool ComputeSearchDirection();

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    InexactSearchDirCalculator();

    /** Copy Constructor */
    InexactSearchDirCalculator(const InexactSearchDirCalculator&);

    /** Overloaded Equals Operator */
    void operator=(const InexactSearchDirCalculator&);
    //@}

    /** Method to easily access Inexact data */
    InexactData& InexData()
    {
      InexactData& inexact_data =
        static_cast<InexactData&>(IpData().AdditionalData());
      DBG_ASSERT(dynamic_cast<InexactData*>(&IpData().AdditionalData()));
      return inexact_data;
    }

    /** Method to easily access Inexact calculated quantities */
    InexactCq& InexCq()
    {
      InexactCq& inexact_cq =
        static_cast<InexactCq&>(IpCq().AdditionalCq());
      DBG_ASSERT(dynamic_cast<InexactCq*>(&IpCq().AdditionalCq()));
      return inexact_cq;
    }

    /** @name Algorithmic options */
    //@{
    /** termination tolerance for local infeasibility */
    Number local_inf_Ac_tol_;
    //@}

    /** @name Strategy objects */
    //@{
    SmartPtr<InexactNormalStepCalculator> normal_step_calculator_;
    SmartPtr<InexactPDSolver> inexact_pd_solver_;
    //@}

    /** enumeration for decomposition options */
    enum DecompositionTypeEnum
    {
      ALWAYS=0,
      ADAPTIVE,
      SWITCH_ONCE
    };
    /** Type of decomposition */
    DecompositionTypeEnum decomposition_type_;
  };

} // namespace Ipopt

#endif
