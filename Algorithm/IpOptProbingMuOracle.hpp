// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2004-11-12

#ifndef __IPOPTPROBINGMUORACLE_HPP__
#define __IPOPTPROBINGMUORACLE_HPP__

#include "IpMuOracle.hpp"
#include "IpPDSystemSolver.hpp"

namespace Ipopt
{

  /** Implementation of the probing strategy for computing the
   *  barrier parameter.
   */
  class OptProbingMuOracle : public MuOracle
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    OptProbingMuOracle(const SmartPtr<PDSystemSolver>& pd_solver);
    /** Default destructor */
    virtual ~OptProbingMuOracle();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the value of the barrier parameter that
     *  could be used in the current iteration (using the LOQO formula).
     */
    virtual Number CalculateMu();

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
    OptProbingMuOracle();
    /** Copy Constructor */
    OptProbingMuOracle(const OptProbingMuOracle&);

    /** Overloaded Equals Operator */
    void operator=(const OptProbingMuOracle&);
    //@}

    /** Pointer to the object that should be used to solve the
     *  primal-dual system.
     */
    SmartPtr<PDSystemSolver> pd_solver_;

    /** Auxilliary function for computing the average complementarity
     *  at a point, given step sizes and step
     */
    // ToDo Replace pointers by references
    Number CalculateQualityFunction(Number sigma,
                                    const Vector& step_aff_x_L,
                                    const Vector& step_aff_x_U,
                                    const Vector& step_aff_s_L,
                                    const Vector& step_aff_s_U,
				    const Vector& step_aff_y_c,
				    const Vector& step_aff_y_d,
                                    const Vector& step_aff_z_L,
                                    const Vector& step_aff_z_U,
                                    const Vector& step_aff_v_L,
                                    const Vector& step_aff_v_U,
                                    const Vector& step_cen_x_L,
                                    const Vector& step_cen_x_U,
                                    const Vector& step_cen_s_L,
                                    const Vector& step_cen_s_U,
				    const Vector& step_cen_y_c,
				    const Vector& step_cen_y_d,
                                    const Vector& step_cen_z_L,
                                    const Vector& step_cen_z_U,
                                    const Vector& step_cen_v_L,
                                    const Vector& step_cen_v_U,
				    SmartPtr<const Vector> jac_cT_times_step_aff_y_c,
				    SmartPtr<const Vector> jac_dT_times_step_aff_y_d,
				    SmartPtr<const Vector> jac_cT_times_step_cen_y_c,
				    SmartPtr<const Vector> jac_dT_times_step_cen_y_d);

    /** Auxilliary function performing the golden bisection */
    Number PerformGoldenBisection(Number sigma_up,
                                  Number sigma_lo,
                                  Number tol,
                                  const Vector& step_aff_x_L,
                                  const Vector& step_aff_x_U,
                                  const Vector& step_aff_s_L,
                                  const Vector& step_aff_s_U,
				  const Vector& step_aff_y_c,
				  const Vector& step_aff_y_d,
                                  const Vector& step_aff_z_L,
                                  const Vector& step_aff_z_U,
                                  const Vector& step_aff_v_L,
                                  const Vector& step_aff_v_U,
                                  const Vector& step_cen_x_L,
                                  const Vector& step_cen_x_U,
                                  const Vector& step_cen_s_L,
                                  const Vector& step_cen_s_U,
				  const Vector& step_cen_y_c,
				  const Vector& step_cen_y_d,
                                  const Vector& step_cen_z_L,
                                  const Vector& step_cen_z_U,
                                  const Vector& step_cen_v_L,
                                  const Vector& step_cen_v_U,
				  SmartPtr<const Vector> jac_cT_times_step_aff_y_c,
				  SmartPtr<const Vector> jac_dT_times_step_aff_y_d,
				  SmartPtr<const Vector> jac_cT_times_step_cen_y_c,
				  SmartPtr<const Vector> jac_dT_times_step_cen_y_d);

    /** @name Algorithmic parameters */
    //@{
    /** Upper bound on centering parameter sigma */
    Number sigma_max_;
    /** Norm to be used for the quality function. 1: 1-norm, 2: 2-norm */
    Index quality_function_norm_;
    /** Flag indicating whether the components of the quality function
     *  should be normalized. */
    bool quality_function_normalized_;
    /** Flag indicating how centrality should be involved in the
     *  quality function */
    Index quality_function_centrality_;
    /** Flag indicating what term is to be used for the dual
     *  infeasibility in the quality function. */
    Index quality_function_dual_inf_;
    /** Relative tolerance for golden bi-section algorithm. */
    Number bisection_tol_;
    /** Maximal number of bi-section steps in the golden bisection
     *  search for sigma. */
    Index max_bisection_steps_;
    /** Flag indicating whether the dual step size is to be used for
     *  the equality constraint multipliers. Note, this must be the
     *  same as in any line search option. */
    bool dual_alpha_for_y_;
    //@}

    /** @name scaling values */
    //@{
    Number dual_inf_scal_;
    Number primal_inf_scal_;
    Number compl_inf_scal_;
    //@}
  };

} // namespace Ipopt

#endif
