// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IP_AUGSYSTEMSOLVER_HPP__
#define __IP_AUGSYSTEMSOLVER_HPP__

#include "IpUtils.hpp"
#include "IpSmartPtr.hpp"
#include "IpSymMatrix.hpp"
#include "IpSymLinearSolver.hpp"
#include "IpAlgStrategy.hpp"

namespace Ipopt
{

  /** Base class for Solver for the augmented system.  This is the
   *  base class for linear solvers that solve the augmented system,
   *  which is defined as
   *
   *  \f$\left[\begin{array}{cccc}
   *  W + D_x + \delta_xI & 0 & J_c^T & J_d^T\\
   *  0 & D_s + \delta_sI & 0 & -I \\
   *  J_c & 0 & D_c - \delta_cI & 0\\
   *  J_d & -I & 0 & D_d - \delta_dI
   *  \end{array}\right]
   *  \left(\begin{array}{c}sol_x\\sol_s\\sol_c\\sol_d\end{array}\right)=
   *  \left(\begin{array}{c}rhs_x\\rhs_s\\rhs_c\\rhs_d\end{array}\right)\f$
   *
   *  Since this system might be solved repeatedly for different right
   *  hand sides, it is desirable to step the factorization of a
   *  direct linear solver if possible.
   */
  class AugSystemSolver: public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default constructor. */
    AugSystemSolver()
    {}
    /** Default destructor */
    virtual ~AugSystemSolver()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** Set up the augmented system and solve it for a given right hand
     *  side.  If desired (i.e. if check_NegEVals is true), then the
     *  solution is only computed if the number of negative eigenvalues
     *  matches numberOfNegEVals.
     *
     *  The return value is the return value of the linear solver object.
     */
    virtual SymLinearSolver::ESolveStatus Solve(
      const SymMatrix* W,
      const Vector* D_x,
      double delta_x,
      const Vector* D_s,
      double delta_s,
      const Matrix* J_c,
      const Vector* D_c,
      double delta_c,
      const Matrix* J_d,
      const Vector* D_d,
      double delta_d,
      const Vector& rhs_x,
      const Vector& rhs_s,
      const Vector& rhs_c,
      const Vector& rhs_d,
      Vector& sol_x,
      Vector& sol_s,
      Vector& sol_c,
      Vector& sol_d,
      bool check_NegEVals,
      Index numberOfNegEVals) =0;

    /** Number of negative eigenvalues detected during last
     * solve.  Returns the number of negative eigenvalues of
     * the most recent factorized matrix.  This must not be called if
     * the linear solver does not compute this quantities (see
     * ProvidesInertia).
     */
    virtual Index NumberOfNegEVals() const =0;

    /** Query whether inertia is computed by linear solver.
     *  Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const =0;

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    AugSystemSolver(const AugSystemSolver&);

    /** Overloaded Equals Operator */
    void operator=(const AugSystemSolver&);
    //@}

  };

} // namespace Ipopt

#endif
