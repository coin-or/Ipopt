// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPPDFULLSPACESOLVER_HPP__
#define __IPPDFULLSPACESOLVER_HPP__

#include "IpPDSystemSolver.hpp"
#include "IpAugSystemSolver.hpp"

namespace Ipopt
{

  /** This is the implemetation of the Primal-Dual System, using the
   *  full space approach with a direct linear solver.
   */
  class PDFullSpaceSolver: public PDSystemSolver
  {
  public:
    /** @name /Destructor */
    //@{
    /** Constructor that takes in the Augmented System solver that
     *  is to be used inside
     */
    PDFullSpaceSolver(AugSystemSolver& augSysSolver);

    /** Default destructor */
    virtual ~PDFullSpaceSolver();
    //@}

    /* overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    /** Solve the primal dual system, given one right hand side.
     */
    virtual void Solve(Number alpha,
                       Number beta,
                       const Vector& rhs_x,
                       const Vector& rhs_s,
                       const Vector& rhs_c,
                       const Vector& rhs_d,
                       const Vector& rhs_zL,
                       const Vector& rhs_zU,
                       const Vector& rhs_vL,
                       const Vector& rhs_vU,
                       Vector& res_x,
                       Vector& res_s,
                       Vector& res_c,
                       Vector& res_d,
                       Vector& res_zL,
                       Vector& res_zU,
                       Vector& res_vL,
                       Vector& res_vU,
                       bool allow_inexact=false);

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
    PDFullSpaceSolver();
    /** Overloaded Equals Operator */
    PDFullSpaceSolver& operator=(const PDFullSpaceSolver&);
    //@}

    /** Pointer to the Solver for the augmented system */
    SmartPtr<AugSystemSolver> augSysSolver_;

    /**@name Data about the correction made to the system */
    //@{
    /** A dummy cache to figure out if the deltas are still up to date*/
    CachedResults<void*> dummy_cache_;
    /** The current value for delta_x */
    Number delta_x_curr_;
    /** The current value for delta_s */
    Number delta_s_curr_;
    /** The current value for delta_c */
    Number delta_c_curr_;
    /** The current value for delta_d */
    Number delta_d_curr_;
    /** The last nonzero value for delta_x */
    Number delta_x_last_;
    /** The last nonzero value for delta_s */
    Number delta_s_last_;
    /** The last nonzero value for delta_c */
    Number delta_c_last_;
    /** The last nonzero value for delta_d */
    Number delta_d_last_;
    //@}

    /** @name Parameters */
    //@{
    /** Minimal number of iterative refinement performed per backsolve */
    Index num_min_iter_ref_;
    /** Maximal value for the regularization. */
    Number delta_regu_max_;
    //@}

    /** Internal function for a single backsolve (which will be used
     *  for iterative refinement on the outside).  This method returns
     *  false, if for some reason the linear system could not be
     *  solved (e.g. when the regularization parameter becomes too
     *  large.)
     */
    bool SolveOnce(const SymMatrix& W,
                   const Matrix& J_c,
                   const Matrix& J_d,
                   const Matrix& Px_L,
                   const Matrix& Px_U,
                   const Matrix& Pd_L,
                   const Matrix& Pd_U,
                   const Vector& z_L,
                   const Vector& z_U,
                   const Vector& v_L,
                   const Vector& v_U,
                   const Vector& slack_x_L,
                   const Vector& slack_x_U,
                   const Vector& slack_s_L,
                   const Vector& slack_s_U,
                   const Vector& sigma_x,
                   const Vector& sigma_s,
                   Number alpha,
                   Number beta,
                   const Vector& rhs_x,
                   const Vector& rhs_s,
                   const Vector& rhs_c,
                   const Vector& rhs_d,
                   const Vector& rhs_zL,
                   const Vector& rhs_zU,
                   const Vector& rhs_vL,
                   const Vector& rhs_vU,
                   Vector& res_x,
                   Vector& res_s,
                   Vector& res_c,
                   Vector& res_d,
                   Vector& res_zL,
                   Vector& res_zU,
                   Vector& res_vL,
                   Vector& res_vU);

    /** Internal function for computing the residual (resid) given the
     * right hand side (rhs) and the solution of the system (res).
     */
    void ComputeResiduals(const SymMatrix& W,
                          const Matrix& J_c,
                          const Matrix& J_d,
                          const Matrix& Px_L,
                          const Matrix& Px_U,
                          const Matrix& Pd_L,
                          const Matrix& Pd_U,
                          const Vector& z_L,
                          const Vector& z_U,
                          const Vector& v_L,
                          const Vector& v_U,
                          const Vector& slack_x_L,
                          const Vector& slack_x_U,
                          const Vector& slack_s_L,
                          const Vector& slack_s_U,
                          const Vector& sigma_x,
                          const Vector& sigma_s,
                          Number alpha,
                          Number beta,
                          const Vector& rhs_x,
                          const Vector& rhs_s,
                          const Vector& rhs_c,
                          const Vector& rhs_d,
                          const Vector& rhs_zL,
                          const Vector& rhs_zU,
                          const Vector& rhs_vL,
                          const Vector& rhs_vU,
                          const Vector& res_x,
                          const Vector& res_s,
                          const Vector& res_c,
                          const Vector& res_d,
                          const Vector& res_zL,
                          const Vector& res_zU,
                          const Vector& res_vL,
                          const Vector& res_vU,
                          Vector& resid_x,
                          Vector& resid_s,
                          Vector& resid_c,
                          Vector& resid_d,
                          Vector& resid_zL,
                          Vector& resid_zU,
                          Vector& resid_vL,
                          Vector& resid_vU);

    /** @name Auxilliary functions */
    //@{
    /** Compute
     * \f$ x = \alpha P S^{-1} z + \beta x \f$.
     */
    void AddPSinvZ(Number alpha, const Matrix& P,
                   const Vector& S, const Vector& Z,
                   Number beta, Vector& X);
    /** Compute \f$ x = S^{-1}(r + \alpha Z P^T d)\f$ */
    void SinvBlrmZPTdBr(Number alpha, const Vector& S,
                        const Vector& R, const Vector& Z,
                        const Matrix& P, const Vector& D, Vector& X);
    /** Compute \f$ y = \alpha* x + \beta * y \f$ */
    void AxpBy(Number alpha, const Vector& X, Number beta, Vector& Y);
    //@}
  };

} // namespace Ipopt

#endif
