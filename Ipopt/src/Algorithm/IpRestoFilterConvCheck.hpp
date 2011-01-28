// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//
//           A Waechter: moved most code to IpRestoConvCheck.cpp 2008-06-24

#ifndef __IPRESTOFILTERCONVCHECK_HPP__
#define __IPRESTOFILTERCONVCHECK_HPP__

#include "IpRestoConvCheck.hpp"
#include "IpFilterLSAcceptor.hpp"

namespace Ipopt
{

  /** This is the implementation of the restoration convergence check
   *  is the original algorithm used the filter globalization
   *  mechanism.
   */
  class RestoFilterConvergenceCheck :
        public RestoConvergenceCheck
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    RestoFilterConvergenceCheck();

    /** Default destructor */
    virtual ~RestoFilterConvergenceCheck();
    //@}

    /** Set the object for the original filter line search. Here,
     *  orig_filter_ls_acceptor must be the same strategy object to
     *  which the restoration phase object with this object is given.
     *  This method must be called to finish the definition of the
     *  algorithm, before Initialize is called. */
    void SetOrigLSAcceptor(const BacktrackingLSAcceptor& orig_ls_acceptor);

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Methods used by IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}
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

    /** Method for checking progress with original filter
     *  globalization mechanism.  Overloaded from
     *  RestoConvergenceCheck. */
    virtual ConvergenceStatus
    TestOrigProgress(Number orig_trial_barr, Number orig_trial_theta);

    /** Strategy object for the filter line search method for the
     *  original NLP.  CAREFUL: We must not hold on to this object
     *  with a SmartPtr, because have otherwise circular references
     *  that prevent the destructor of the line search object to be
     *  called! */
    const FilterLSAcceptor* orig_filter_ls_acceptor_;
  };

} // namespace Ipopt

#endif
