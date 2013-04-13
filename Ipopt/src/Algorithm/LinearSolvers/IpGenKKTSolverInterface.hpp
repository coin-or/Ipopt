// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter     IBM    2007-03-01

#ifndef __IPGENKKTSOLVERINTERFACE_HPP__
#define __IPGENKKTSOLVERINTERFACE_HPP__

#include "IpUtils.hpp"
#include "IpAlgStrategy.hpp"
#include "IpSymLinearSolver.hpp"

namespace Ipopt
{
  /** Base class for interfaces to symmetric indefinite linear solvers
   *  for generic matrices   */
  class GenKKTSolverInterface: public AlgorithmStrategyObject
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    GenKKTSolverInterface()
    {}

    virtual ~GenKKTSolverInterface()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Solve operation for multiple right hand sides.  The linear
     *  system is of the form
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
     * (see also AugSystemSolver).
     *
     *  The return code is SYMSOLV_SUCCESS if the factorization and
     *  solves were successful, SYMSOLV_SINGULAR if the linear system
     *  is singular, and SYMSOLV_WRONG_INERTIA if check_NegEVals is
     *  true and the number of negative eigenvalues in the matrix does
     *  not match numberOfNegEVals.  If SYMSOLV_CALL_AGAIN is
     *  returned, then the calling function will request the pointer
     *  for the array for storing a again (with GetValuesPtr), write
     *  the values of the nonzero elements into it, and call this
     *  MultiSolve method again with the same right-hand sides.  (This
     *  can be done, for example, if the linear solver realized it
     *  does not have sufficient memory and needs to redo the
     *  factorization; e.g., for MA27.)
     *
     *  The number of right-hand sides is given by nrhs, the values of
     *  the right-hand sides are given in rhs_vals (one full right-hand
     *  side stored immediately after the other), and solutions are
     *  to be returned in the same array.
     *
     *  check_NegEVals will not be chosen true, if ProvidesInertia()
     *  returns false.
     */
    virtual ESymSolverStatus MultiSolve(
      bool new_matrix /** If this flag is false, the same matrix as in the most recent call is given to the solver again */
      , Index n_x  /** Dimension of D_x */
      , Index n_c /** Dimension of D_s and D_c */
      , Index n_d /** Dimension of D_d */
      , SmartPtr<const SymMatrix> W /** Hessian of Lagrangian (as given by NLP) */
      , SmartPtr<const Matrix> Jac_c /** Jacobian of equality constraints (as given by NLP) */
      , SmartPtr<const Matrix> Jac_d /** Jacobian of inequality constraints (as given by NLP) */
      , const Number* D_x /** Array with the elements D_x (if NULL, assume all zero) */
      , const Number* D_s /** Array with the elements D_s (if NULL, assume all zero) */
      , const Number* D_c /** Array with the elements D_c (if NULL, assume all zero) */
      , const Number* D_d /** Array with the elements D_d (if NULL, assume all zero) */
      , Number delta_x /** \f$ \delta_x\f$ */
      , Number delta_s /** \f$ \delta_s\f$ */
      , Number delta_c /** \f$ \delta_c\f$ */
      , Number delta_d /** \f$ \delta_d\f$ */
      , Index n_rhs  /** Number of right hand sides */
      , Number* rhssol /** On input, this containts the right hand sides, and on successful termination of the solver, the solutions are expected in there on return. At the moment, the order is x,d,c,s, but this can be made flexible and chosen according to an option. */
      , bool check_NegEVals /** if true, we want to ensure that the inertia is correct */
      , Index numberOfNegEVals /** Required number of negative eigenvalues if check_NegEVals is true */
    )=0;

    /** Number of negative eigenvalues detected during last
     *  factorization.  Returns the number of negative eigenvalues of
     *  the most recent factorized matrix.  This must not be called if
     *  the linear solver does not compute this quantities (see
     *  ProvidesInertia).
     */
    virtual Index NumberOfNegEVals() const =0;
    //@}

    //* @name Options of Linear solver */
    //@{
    /** Request to increase quality of solution for next solve.  The
     *  calling class asks linear solver to increase quality of
     *  solution for the next solve (e.g. increase pivot tolerance).
     *  Returns false, if this is not possible (e.g. maximal pivot
     *  tolerance already used.)
     */
    virtual bool IncreaseQuality() =0;

    /** Query whether inertia is computed by linear solver.  Returns
     *  true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const =0;
    //@}
  };

} // namespace Ipopt

#endif
