// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter           IBM    2008-08-31
//               derived from IpIpoptCalculatedQuantities.hpp

#ifndef __IPINEXACTCQ_HPP__
#define __IPINEXACTCQ_HPP__

#include "IpIpoptCalculatedQuantities.hpp"
#include "IpInexactData.hpp"

namespace Ipopt
{

  /** Class for all Chen-Goldfarb penalty method specific calculated
   *  quantities.
   */
  class InexactCq : public IpoptAdditionalCq
  {
  public:

    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    InexactCq(IpoptNLP* ip_nlp,
              IpoptData* ip_data,
              IpoptCalculatedQuantities* ip_cq);

    /** Default destructor */
    virtual ~InexactCq();
    //@}

    /** This method must be called to initialize the global
     *  algorithmic parameters.  The parameters are taken from the
     *  OptionsList object. */
    bool Initialize(const Journalist& jnlst,
                    const OptionsList& options,
                    const std::string& prefix);

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(const SmartPtr<RegisteredOptions>& roptions);
    //@}

    /** Gradient of infeasibility w.r.t. x.  Jacobian of equality
     *  constraints transpose times the equality constraints plus
     *  Jacobian of the inequality constraints transpose times the
     *  inequality constraints (including slacks). */
    SmartPtr<const Vector> curr_jac_cdT_times_curr_cdminuss();

    /** Vector of all inequality slacks for doing the slack-based scaling */
    SmartPtr<const Vector> curr_scaling_slacks();

    /** Vector with the slack-scaled d minus s inequalities */
    SmartPtr<const Vector> curr_slack_scaled_d_minus_s();

    /** Scaled norm of Ac */
    Number curr_scaled_Ac_norm();

    /** Scaled, squared norm of A */
    Number curr_scaled_A_norm2();

    /** Compute the 2-norm of a slack-scaled vector with x and s
     *  component */
    Number slack_scaled_norm(const Vector& x, const Vector &s);

    /** Compute x component of the W*vec product for the current
     *  Hessian and a vector */
    SmartPtr<const Vector> curr_W_times_vec_x(const Vector& vec_x);

    /** Compute s component of the W*vec product for the current
     *  Hessian and a vector */
    SmartPtr<const Vector> curr_W_times_vec_s(const Vector& vec_s);

    /** Compute x component of the W*u product for the current values.
     *  u here is the tangential step. */
    SmartPtr<const Vector> curr_Wu_x();

    /** Compute s component of the W*u product for the current values.
     *  u here is the tangential step. */
    SmartPtr<const Vector> curr_Wu_s();

    /** Compute the u^T*W*u product for the current values.  u here is the
    tangential step. */
    Number curr_uWu();

    /** Compute the c-component of the product of the current
     *  constraint Jacobian with the current normal step */
    SmartPtr<const Vector> curr_jac_times_normal_c();

    /** Compute the d-component of the product of the current
     *  constraint Jacobian with the current normal step */
    SmartPtr<const Vector> curr_jac_times_normal_d();

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
    InexactCq();

    /** Copy Constructor */
    InexactCq(const InexactCq&);

    /** Overloaded Equals Operator */
    void operator=(const InexactCq&);
    //@}

    /** @name Pointers for easy access to data and NLP information. To
     *  avoid circular references of Smart Pointers, we use a regular
     *  pointer here. */
    //@{
    IpoptNLP* ip_nlp_;
    IpoptData* ip_data_;
    IpoptCalculatedQuantities* ip_cq_;
    //@}

    /** Method to easily access Inexact data */
    InexactData& InexData()
    {
      InexactData& inexact_data =
        static_cast<InexactData&>(ip_data_->AdditionalData());
      DBG_ASSERT(dynamic_cast<InexactData*>(&ip_data_->AdditionalData()));
      return inexact_data;
    }

    /** @name Caches */
    //@{
    CachedResults<SmartPtr<const Vector> > curr_jac_cdT_times_curr_cdminuss_cache_;
    CachedResults<SmartPtr<const Vector> > curr_scaling_slacks_cache_;
    CachedResults<SmartPtr<const Vector> > curr_slack_scaled_d_minus_s_cache_;
    CachedResults<Number> curr_scaled_Ac_norm_cache_;
    CachedResults<Number> slack_scaled_norm_cache_;
    CachedResults<SmartPtr<const Vector> > curr_W_times_vec_x_cache_;
    CachedResults<SmartPtr<const Vector> > curr_W_times_vec_s_cache_;
    CachedResults<SmartPtr<const Vector> > curr_Wu_x_cache_;
    CachedResults<SmartPtr<const Vector> > curr_Wu_s_cache_;
    CachedResults<Number> curr_uWu_cache_;
    CachedResults<SmartPtr<const Vector> > curr_jac_times_normal_c_cache_;
    CachedResults<SmartPtr<const Vector> > curr_jac_times_normal_d_cache_;
    //@}

    /** Upper bound on slack-based scaling factors */
    Number slack_scale_max_;
  };

} // namespace Ipopt

#endif
