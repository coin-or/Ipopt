// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPIPOPTALG_HPP__
#define __IPIPOPTALG_HPP__

#include "IpUtils.hpp"
#include "IpSmartPtr.hpp"
#include "IpAlgStrategy.hpp"
#include "IpPDSystemSolver.hpp"
#include "IpLineSearch.hpp"
#include "IpMuUpdate.hpp"
#include "IpConvCheck.hpp"
#include "IpOptionsList.hpp"
#include "IpIterateInitializer.hpp"
#include "IpIterationOutput.hpp"

namespace Ipopt
{

  /** forward declarations */

  /** The main ipopt algorithm class.
   *  Main Ipopt algorithm class, contains the main optimize method,
   *  handles the execution of the optimization.
   *  The constructor initializes the data structures through the nlp,
   *  and the Optimize method then assumes that everything is 
   *  initialized and ready to go.
   *  After an optimization is complete, the user can access the 
   *  solution through the passed in ip_data structure.
   *  Multiple calls to the Optimize method are allowed as long as the
   *  structure of the problem remains the same (i.e. starting point
   *  or nlp parameter changes only).
   */
  class IpoptAlgorithm : public AlgorithmStrategyObject
  {
  public:
    /**@name Enumerations */
    //@{
    /** enum for the return from the optimize algorithm
     *  (obviously we need to add more) */
    enum SolverReturn {
      SUCCESS,
      MAXITER_EXCEEDED,
      FAILED
    };
    //@}

    /**@name Constructors/Destructors */
    //@{
    /** Constructor. (The IpoptAlgorithm uses smart pointers for these
     *  passed-in pieces to make sure that a user of IpoptAlgoroithm
     *  cannot pass in an object created on the stack!)
     */
    IpoptAlgorithm(const SmartPtr<PDSystemSolver>& pd_solver,
                   const SmartPtr<LineSearch>& line_search,
                   const SmartPtr<MuUpdate>& mu_update,
                   const SmartPtr<ConvergenceCheck>& conv_check,
                   const SmartPtr<IterateInitializer>& iterate_initializer,
                   const SmartPtr<IterationOutput>& iter_output);

    /** Default destructor */
    virtual ~IpoptAlgorithm();
    //@}


    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Main solve method. */
    SolverReturn Optimize();

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
    IpoptAlgorithm();

    /** Copy Constructor */
    IpoptAlgorithm(const IpoptAlgorithm&);

    /** Overloaded Equals Operator */
    void operator=(const IpoptAlgorithm&);
    //@}

    /** @name Strategy objects */
    //@{
    SmartPtr<PDSystemSolver> pd_solver_;
    SmartPtr<LineSearch> line_search_;
    SmartPtr<MuUpdate> mu_update_;
    SmartPtr<ConvergenceCheck> conv_check_;
    SmartPtr<IterateInitializer> iterate_initializer_;
    SmartPtr<IterationOutput> iter_output_;
    //@}

    /** @name Main steps of the algorthim */
    //@{
    /** Method for updating the current Hessian.  This can either just
     *  evaluate the exact Hessian (based on the current iterate), or
     *  perform a quasi-Newton update.
     */
    void ActualizeHessian();

    /** Method to update the barrier parameter
     * ( this may later be made a strategy object
     *   and passed in ) */
    void UpdateBarrierParameter();

    /** Method to setup the call to the PDSystemSolver */
    void ComputeSearchDirection();

    /** Method computing the new iterate (usually vialine search).
     *  The acceptable point is the one in trial after return.
     */
    void ComputeAcceptableTrialPoint();

    /** Method for accepting the trial point as the new iteration,
     *  possibly after adjusting the variable bounds in the NLP. */
    void AcceptTrialPoint();

    /** Do all the output for one iteration */
    void OutputIteration();

    /** Sets up initial values for the iterates,
     * Corrects the initial values for x and s (force in bounds) 
     */
    void InitializeIterates();

    /** Print the problem size statistics */
    void PrintProblemStatistics();
    //@}

    /** @name internal flags */
    //@{
    /** Flag indicating if the statistic should not be printed */
    bool skip_print_problem_stats_;
    //@}

    /** @name auxilliary functions */
    //@{
    void calc_number_of_bounds(
      const Vector& x,
      const Vector& x_L,
      const Vector& x_U,
      const Matrix& Px_L,
      const Matrix& Px_U,
      Index& n_tot,
      Index& n_only_lower,
      Index& n_both,
      Index& n_only_upper);
    //@}
  };

} // namespace Ipopt

#endif
