// Copyright (C) 2004,2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2004-11-12

#ifndef __IPQUALITYFUNCTIONMUORACLE_HPP__
#define __IPQUALITYFUNCTIONMUORACLE_HPP__

#include "IpMuOracle.hpp"
#include "IpPDSystemSolver.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

namespace Ipopt
{

  /** Implementation of the probing strategy for computing the
   *  barrier parameter.
   */
  class QualityFunctionMuOracle : public MuOracle
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    QualityFunctionMuOracle(const SmartPtr<PDSystemSolver>& pd_solver);
    /** Default destructor */
    virtual ~QualityFunctionMuOracle();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method for computing the value of the barrier parameter that
     *  could be used in the current iteration (using the LOQO formula).
     */
    virtual Number CalculateMu();

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

    /** @name Public enums.  Some of those are also used for the
     *  quality function */
    //@{
    /** enum for norm type */
    enum NormEnum
    {
      NM_NORM_1=0,
      NM_NORM_2_SQUARED,
      NM_NORM_MAX,
      NM_NORM_2
    };
    /** enum for centrality type */
    enum CentralityEnum
    {
      CEN_NONE=0,
      CEN_LOG,
      CEN_RECIPROCAL,
      CEN_CUBED_RECIPROCAL
    };
    /** enum for the quality function balancing term type */
    enum BalancingTermEnum
    {
      BT_NONE=0,
      BT_CUBIC
    };
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
    QualityFunctionMuOracle();
    /** Copy Constructor */
    QualityFunctionMuOracle(const QualityFunctionMuOracle&);

    /** Overloaded Equals Operator */
    void operator=(const QualityFunctionMuOracle&);
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
                                    const Vector& step_cen_v_U);

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
                                  const Vector& step_cen_v_U);

    /** Auxilliary function performing the golden bisection in the
     *  logarithmic scale */
    /* This doesn't seem to work well, so I took it out for now (AW)
    Number PerformGoldenBisectionLog(Number sigma_up,
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
                                     const Vector& step_cen_v_U);
    */

    /** @name Algorithmic parameters */
    //@{
    /** Upper bound on centering parameter sigma */
    Number sigma_max_;
    /** Norm to be used for the quality function. */
    NormEnum quality_function_norm_;
    /** Flag indicating how centrality should be involved in the
     *  quality function */
    CentralityEnum quality_function_centrality_;
    /** Flag indicating whether we use a balancing term in the quality
     *  function.
     */
    BalancingTermEnum quality_function_balancing_term_;
    /** Relative tolerance for golden bi-section algorithm. */
    Number bisection_tol_;
    /** Maximal number of bi-section steps in the golden bisection
     *  search for sigma. */
    Index max_bisection_steps_;
    //@}

    /** @name Temporary work space vectors.  We use those to avoid
     *  repeated reallocation in CalculateQualityFunction. */
    //@{
    SmartPtr<Vector> tmp_step_x_L_;
    SmartPtr<Vector> tmp_step_x_U_;
    SmartPtr<Vector> tmp_step_s_L_;
    SmartPtr<Vector> tmp_step_s_U_;
    SmartPtr<Vector> tmp_step_z_L_;
    SmartPtr<Vector> tmp_step_z_U_;
    SmartPtr<Vector> tmp_step_v_L_;
    SmartPtr<Vector> tmp_step_v_U_;

    SmartPtr<Vector> tmp_slack_x_L_;
    SmartPtr<Vector> tmp_slack_x_U_;
    SmartPtr<Vector> tmp_slack_s_L_;
    SmartPtr<Vector> tmp_slack_s_U_;
    SmartPtr<Vector> tmp_z_L_;
    SmartPtr<Vector> tmp_z_U_;
    SmartPtr<Vector> tmp_v_L_;
    SmartPtr<Vector> tmp_v_U_;
    //@}
  };

} // namespace Ipopt

#endif
