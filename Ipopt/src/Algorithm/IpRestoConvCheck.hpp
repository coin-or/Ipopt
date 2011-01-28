// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//
//           was originally IpRestoFilterConvCheck.hpp (rev 781)
//             separated by A Waechter IBM  2008-06-24

#ifndef __IPRESTOCONVCHECK_HPP__
#define __IPRESTOCONVCHECK_HPP__

#include "IpOptErrorConvCheck.hpp"
#include "IpBacktrackingLSAcceptor.hpp"

namespace Ipopt
{

  /** Convergence check for the restoration phase.  This inherits from
   *  the OptimalityErrorConvergenceCheck so that the method for the
   *  regular optimality error convergence criterion can be checked as
   *  well.  In addition, this convergence check returns the CONVERGED
   *  message, if the current iteration is acceptable to the original
   *  globalization scheme.
   *
   */
  class RestoConvergenceCheck :
        public OptimalityErrorConvergenceCheck
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    RestoConvergenceCheck();

    /** Default destructor */
    virtual ~RestoConvergenceCheck();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** overloaded from ConvergenceCheck */
    virtual ConvergenceStatus CheckConvergence(bool call_intermediate_callback = true);

    /** Method for setting the LS acceptor from the main algorithm */
    virtual void SetOrigLSAcceptor(const BacktrackingLSAcceptor& orig_ls_acceptor) = 0;

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
    RestoConvergenceCheck(const RestoConvergenceCheck&);

    /** Overloaded Equals Operator */
    void operator=(const RestoConvergenceCheck&);
    //@}

    /** Method for checking progress with original globalization
     *  mechanism.  This needs to be overloaded */
    virtual ConvergenceStatus
    TestOrigProgress(Number orig_trial_barr, Number orig_trial_theta) = 0;

    /** @name Algorithmic parameters */
    //@{
    /** Fraction of required reduction in infeasibility before problem
     *  is considered to be solved. */
    Number kappa_resto_;
    /** Maximum number of iterations in restoration phase */
    Index maximum_iters_;
    /** Maximum number of succesive iterations in restoration phase */
    Index maximum_resto_iters_;
    /** Constraint violation tolerance for original algorithm */
    Number orig_constr_viol_tol_;
    //@}

    /** Flag indicating that this is the first call.  We don't want to
     *  leave the restoration phase without taking at least one step,
     *  so this flag is used to ensure this. */
    bool first_resto_iter_;

    /** Counter for successive iterations in restoration phase */
    Index successive_resto_iter_;
  };

} // namespace Ipopt

#endif
