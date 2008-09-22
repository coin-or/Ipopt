// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-19

#ifndef __IPITERATIVESOLVERTERMINATIONTESTER_HPP__
#define __IPITERATIVESOLVERTERMINATIONTESTER_HPP__

#include "IpAlgStrategy.hpp"
#include "IpInexactCq.hpp"

namespace Ipopt
{

  /** This object implements the termination tests for the iterative
   *  linear solver in the inexact version of Ipopt.
   */
  class IterativeSolverTerminationTester: public AlgorithmStrategyObject
  {
  public:
    /** Enum to report result of termination test */
    enum ETerminationTest {
      /** The current solution is not yet good enough */
      CONTINUE,
      /** Termination Test 1 is satisfied */
      TEST_1_SATISFIED,
      /** Termination Test 2 is satisfied */
      TEST_2_SATISFIED,
      /** Termination Test 3 is satisfied */
      TEST_3_SATISFIED,
      /** Hessian matrix should be modified */
      MODIFY_HESSIAN
    };

    /** @name /Destructor */
    //@{
    /** Default constructor
     */
    IterativeSolverTerminationTester();

    /** Default destructor */
    virtual ~IterativeSolverTerminationTester();
    //@}

    /* overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

    /** Method for initializing for the next iterative solve.  This
     *  must be call before the test methods are called. */
    bool InitializeSolve();

    /** This method checks if the current soltion of the iterative
     *  linear solver is good enough (by returning the corresponding
     *  satisfied termination test), or if the Hessian should be
     *  modified.  The input is the dimension of the augmented system,
     *  the current solution vector of the augmented system, the
     *  current residual vector. */
    ETerminationTest TestTerminaion(Index ndim, const Number* sol,
                                    const Number* resid);

    /** This method can be called after the Solve is over and we can
     *  delete anything that has been allocated to free memory. */
    void Clear();

    /** An easy way to get the journalist if accessed from the outside */
    const Journalist& GetJnlst() const
    {
      return Jnlst();
    }

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Overloaded Equals Operator */
    IterativeSolverTerminationTester& operator=(const IterativeSolverTerminationTester&);
    //@}

    /** Method for copying a long augmented system array into Vectors
     *  in Ipopt notation */
    void GetVectors(Index ndim, const Number* array,
                    SmartPtr<const Vector>& comp_x,
                    SmartPtr<const Vector>& comp_s,
                    SmartPtr<const Vector>& comp_c,
                    SmartPtr<const Vector>& comp_d);

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
    /** Psi factor in the tangential component condition */
    Number tcc_psi_;
    /** theta factor in the tangential component condition */
    Number tcc_theta_;
    /** zeta factor in the tangential component condition */
    Number tcc_zeta_;
    /** kappa_1 factor in termination test 1 */
    Number tt_kappa1_;
    /** kappa_2 factor in termination test 2 */
    Number tt_kappa2_;
    /** kappa_3 factor in termination test 3 */
    Number tt_kappa3_;
    /** eps_2 constant in termination test 2 */
    Number tt_eps2_;
    /** eps_3 constant in termination test 3 */
    Number tt_eps3_;
    /** rho constant from penalty parameter update.  This is called
     *  \tau_{\pi} in MIPS paper */
    Number rho_;
    //@}

    /** @name Quantities that are identical for all tests and can be
     *  precomputed */
    //@{
    SmartPtr<const Vector> curr_Av_c_;
    SmartPtr<const Vector> curr_Av_d_;
    Number c_norm_;
    SmartPtr<const Vector> curr_c_plus_Av_c_;
    SmartPtr<const Vector> curr_c_plus_Av_d_;
    Number c_plus_Av_norm_;
    Number v_norm_scaled_;
    SmartPtr<const Vector> curr_grad_barrier_obj_x_;
    SmartPtr<const Vector> curr_grad_barrier_obj_s_;
    SmartPtr<const Matrix> curr_jac_c_;
    SmartPtr<const Matrix> curr_jac_d_;
    SmartPtr<const Vector> curr_scaling_slacks_;
    Number curr_Av_norm_;
    Number curr_tt1_norm_;
    Number curr_tt2_norm_;
    SmartPtr<const Vector> curr_Wv_x_;
    SmartPtr<const Vector> curr_Wv_s_;
    //@}

    /** @name Quantities from previous iteration required in the
    tests */
    //@{
    SmartPtr<const Vector> last_grad_barrier_obj_x_;
    SmartPtr<const Vector> last_grad_barrier_obj_s_;
    SmartPtr<const Matrix> last_jac_c_;
    SmartPtr<const Matrix> last_jac_d_;
    SmartPtr<const Vector> last_scaling_slacks_;
    Number last_Av_norm_;
    Number last_tt1_norm_;
    //@}
  };

} // namespace Ipopt

#endif
