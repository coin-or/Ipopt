// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPRESTOFILTERCONVCHECK_HPP__
#define __IPRESTOFILTERCONVCHECK_HPP__

#include "IpUtils.hpp"
#include "IpOptErrorConvCheck.hpp"
#include "IpFilterLineSearch.hpp"

namespace Ipopt
{

  /** Convergence check for the restoration phase as called by the
   *  filter.  This inherits from the OptimalityErrorConvergenceCheck
   *  so that the method for the regular optimality error convergence
   *  criterion can be checked as well.  In addition, this convergence
   *  check returns the CONVERGED message, if the current iteration is
   *  acceptable to the original filter.
   *
   *  Since this object needs to know about the original NLP, it also
   *  inherits from RestoProblemCoupler, so that the restoration phase
   *  object can call the SetObjs method to set the corresponding
   *  pointers before the Initilize for the restoration phase
   *  algorithm is called.
   */
  class RestoFilterConvergenceCheck :
        public OptimalityErrorConvergenceCheck
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    RestoFilterConvergenceCheck();

    /** Default destructor */
    virtual ~RestoFilterConvergenceCheck()
    {}
    //@}

    /** Set the object for the original filter line search. Here,
     *  filter_line_search must be the same strategy object to which
     *  the restoration phase object with this object is given.  This
     *  method must be called to finish the definition of the
     *  algorithm, before Initialize is called. */
    void SetOrigFilterLineSearch(const FilterLineSearch& orig_filter_line_search);

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** overloaded from ConvergenceCheck */
    virtual ConvergenceStatus CheckConvergence();

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    RestoFilterConvergenceCheck(const RestoFilterConvergenceCheck&);

    /** Overloaded Equals Operator */
    void operator=(const RestoFilterConvergenceCheck&);
    //@}

    /** @name Algorithmic parameters */
    //@{
    Number kappa_resto_;
    //@}

    /** Strategy object for the filter line search method for the
     *  original NLP */
    SmartPtr<const FilterLineSearch> orig_filter_line_search_;
  };

} // namespace Ipopt

#endif
