// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-09

#ifndef __IPINEXACTPDSOLVER_HPP__
#define __IPINEXACTPDSOLVER_HPP__

#include "IpAlgStrategy.hpp"
#include "IpAugSystemSolver.hpp"
#include "IpPDPerturbationHandler.hpp"
#include "IpInexactCq.hpp"

namespace Ipopt
{

  /** This is the implemetation of the Primal-Dual System, allowing
   *  the usage of an inexact linear solver.  The step computed is
   *  usually for the tangential step.
   */
  class InexactPDSolver: public AlgorithmStrategyObject
  {
  public:
    /** @name /Destructor */
    //@{
    /** Constructor that takes in the Augmented System solver that
     *  is to be used inside
     */
    InexactPDSolver(AugSystemSolver& augSysSolver,
                    PDPerturbationHandler& perturbHandler);

    /** Default destructor */
    virtual ~InexactPDSolver();
    //@}

    /* overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    /** Solve the primal dual system, given one right hand side.
     */
    virtual bool Solve(const IteratesVector& rhs,
                       IteratesVector& sol);

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
    InexactPDSolver();
    /** Overloaded Equals Operator */
    InexactPDSolver& operator=(const InexactPDSolver&);
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

    /** @name Strategy objects to hold on to. */
    //@{
    /** Pointer to the Solver for the augmented system */
    SmartPtr<AugSystemSolver> augSysSolver_;
    /** Pointer to the Perturbation Handler. */
    SmartPtr<PDPerturbationHandler> perturbHandler_;
    //@}

    /** Internal function for computing the residual (resid) given the
     * right hand side (rhs) and the solution of the system (res).
     */
    void ComputeResiduals(const SymMatrix& W,
                          const Matrix& J_c,
                          const Matrix& J_d,
                          const Matrix& Pd_L,
                          const Matrix& Pd_U,
                          const Vector& v_L,
                          const Vector& v_U,
                          const Vector& slack_s_L,
                          const Vector& slack_s_U,
                          const Vector& sigma_s,
                          const IteratesVector& rhs,
                          const IteratesVector& res,
                          IteratesVector& resid);

    /** Method for checking if the Hessian matrix has to be modified.
     *  All required data is obtained from the Data objects, so those
     *  values have to be set before this is called. */
    bool HessianRequiresChange();

    /** @name Algorithmic options */
    //@{
    /** Psi factor in the tangential component condition */
    Number tcc_psi_;
    /** theta factor in the tangential component condition */
    Number tcc_theta_;
    /** mu exponent when multiplied to theta in the tangential
     *  component condition */
    Number tcc_theta_mu_exponent_;
    /** flag indicating if the Hessian for the (s,s) part should be
     *  modified with the slacks instead of the identity matrix */
    bool modify_hessian_with_slacks_;
    /** Threshold on line search evaluation count to trigger Hessia
     * modification */
    Index inexact_regularization_ls_count_trigger_;
    //@}

    /** flag indicating if we are dealing with the Pardiso solver
     *  (temporary) */
    bool is_pardiso_;

    Index last_info_ls_count_;
  };


} // namespace Ipopt

#endif
