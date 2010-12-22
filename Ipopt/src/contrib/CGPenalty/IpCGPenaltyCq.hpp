// Copyright (C) 2007, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter           IBM    2007-06-04
//               derived from IpIpoptCalculatedQuantities.hpp

#ifndef __IPCGPENALTYCQ_HPP__
#define __IPCGPENALTYCQ_HPP__

#include "IpIpoptCalculatedQuantities.hpp"
#include "IpCGPenaltyData.hpp"

namespace Ipopt
{

  /** Class for all Chen-Goldfarb penalty method specific calculated
   *  quantities.
   */
  class CGPenaltyCq : public IpoptAdditionalCq
  {
  public:

    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    CGPenaltyCq(IpoptNLP* ip_nlp,
                IpoptData* ip_data,
                IpoptCalculatedQuantities* ip_cg);

    /** Default destructor */
    virtual ~CGPenaltyCq();
    //@}

    /** This method must be called to initialize the global
     *  algorithmic parameters.  The parameters are taken from the
     *  OptionsList object. */
    bool Initialize(const Journalist& jnlst,
                    const OptionsList& options,
                    const std::string& prefix);

    /**@name Methods for the Chen-Goldfarb line search */
    //@{

    /** Compute ||delta_c, delta_d||_infty */
    Number curr_jac_cd_norm(Index nrm_type);
    /** Compute gradient scaling based y->Amax */
    Number curr_scaled_y_Amax();
    /** Compute the 2-norm of y plus delta y */
    Number curr_added_y_nrm2();

    /** Method for the penalty function at current point */
    Number curr_penalty_function();
    /** Method for the penalty function at trial point */
    Number trial_penalty_function();
    /** Method for the directional derivative of the penalty function
     *  at current point with current step in delta */
    Number curr_direct_deriv_penalty_function();
    /** Method for the directional derivative of the penalty function
     *  at current point with current "fast" step in delta_cgpen */
    Number curr_fast_direct_deriv_penalty_function();

    /** Quality of d^T Aug(H) d */
    Number dT_times_barH_times_d();


    /** Method for the current value for the perturbation factor for
     *  the Chen-Goldfarb method.  The factor is computed as 2-norm of
     *  the constraints devided by the current penbalty parameter */
    Number curr_cg_pert_fact();

    /** Method for choose line search penalty parameter */
    Number compute_curr_cg_penalty(const Number );

    /** Method for choose penalty parameters for scaling the KKT system  */
    Number compute_curr_cg_penalty_scale();
    //@}

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(const SmartPtr<RegisteredOptions>& roptions);
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
    CGPenaltyCq();

    /** Copy Constructor */
    CGPenaltyCq(const CGPenaltyCq&);

    /** Overloaded Equals Operator */
    void operator=(const CGPenaltyCq&);
    //@}

    /** @name Pointers for easy access to data and NLP information. To
     *  avoid circular references of Smart Pointers, we use a regular
     *  pointer here. */
    //@{
    IpoptNLP* ip_nlp_;
    IpoptData* ip_data_;
    IpoptCalculatedQuantities* ip_cq_;
    //@}

    /** Method to easily access CGPenalty data */
    CGPenaltyData& CGPenData()
    {
      CGPenaltyData& cg_pen_data =
        static_cast<CGPenaltyData&>(ip_data_->AdditionalData());
      DBG_ASSERT(dynamic_cast<CGPenaltyData*>(&ip_data_->AdditionalData()));
      return cg_pen_data;
    }

    /**@name Caches for the Chen-Goldfarb line search */
    //@{
    CachedResults<Number> curr_fast_direct_deriv_penalty_function_cache_;
    CachedResults<Number> curr_jac_cd_norm_cache_;
    CachedResults<Number> curr_scaled_y_Amax_cache_;
    CachedResults<Number> curr_added_y_nrm2_cache_;
    /** Cache for the penalty function at current point */
    CachedResults<Number> curr_penalty_function_cache_;
    /** Cache for the penalty function at trial point */
    CachedResults<Number> trial_penalty_function_cache_;
    /** Cache for the directional derivative of the penalty function
     *  at current point with step in delta */
    CachedResults<Number> curr_direct_deriv_penalty_function_cache_;
    /** Cache for Chen-Goldfarb perturbation factor. */
    CachedResults<Number> curr_cg_pert_fact_cache_;
    //@}

    /** Parameters for penalty method */
    //@{
    Number reference_infeasibility_;
    //@}

    /** flag indicating if Initialize method has been called (for
     *  debugging) */
    bool initialize_called_;
  };

} // namespace Ipopt

#endif
