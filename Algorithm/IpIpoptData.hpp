// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPIPOPTDATA_HPP__
#define __IPIPOPTDATA_HPP__

#include "IpUtils.hpp"
#include "IpSmartPtr.hpp"
#include "IpVector.hpp"
#include "IpSymMatrix.hpp"
#include "IpReferenced.hpp"
#include "IpOptionsList.hpp"

namespace Ipopt
{
  /* Forward declaration */
  class IpoptNLP;

  /** Class to organize all the data required by the algorithm.
   *  Internally, once this Data object has been initialized, all
   *  internal curr_ vectors must always be set (so that prototyes are
   *  available).  The current values can only be set from the trial
   *  values.  The trial values can be set by copying from a vector or
   *  by adding some fraction of a step to the current values.  This
   *  object also stores steps, which allows to easily communicate the
   *  step from the step computation object to the line search object.
   */
  class IpoptData : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    IpoptData();

    /** Default destructor */
    ~IpoptData();
    //@}

    /** Initialize Data Structures */
    bool InitializeDataStructures(IpoptNLP& ip_nlp,
                                  bool want_x,
                                  bool want_y_c,
                                  bool want_y_d,
                                  bool want_z_L,
                                  bool want_z_U,
                                  bool want_v_L,
                                  bool want_v_U);

    /** This method must be called to initialize the global
     *  algorithmic parameters.  The parameters are taken from the
     *  OptionsList object. */
    bool Initialize(const Journalist& jnlst,
                    const OptionsList& options,
                    const std::string& prefix);

    /** @name Get Methods for Iterates */
    //@{
    /** Main iteration variables
     * (current iteration) */
    SmartPtr<const Vector> curr_x()
    {
      DBG_ASSERT(IsValid(curr_x_));
      return curr_x_;
    }

    /** Main iteration variables
     *  (trial calculations) */
    SmartPtr<const Vector> trial_x()
    {
      DBG_ASSERT(IsValid(trial_x_));
      return trial_x_;
    }

    /** Slacks from inequality constraints
     *  (current iteration) */
    SmartPtr<const Vector> curr_s()
    {
      DBG_ASSERT(IsValid(curr_s_));
      return curr_s_;
    }

    /** Slacks from inequality constraints
     *  (trial calculations) */
    SmartPtr<const Vector> trial_s()
    {
      DBG_ASSERT(IsValid(trial_s_));
      return trial_s_;
    }

    /** Multipliers for equality constraints
     *  (current iteration) */
    SmartPtr<const Vector> curr_y_c()
    {
      DBG_ASSERT(IsValid(curr_y_c_));
      return curr_y_c_;
    }

    /** Multipliers for equality constraints
     *  (trial iteration) */
    SmartPtr<const Vector> trial_y_c()
    {
      DBG_ASSERT(IsValid(trial_y_c_));
      return trial_y_c_;
    }

    /** Multipliers for the inequality constraints
     *  - Note, inequalities made equality 
     *  by introduction of slacks)
     *  (current iteration) */
    SmartPtr<const Vector> curr_y_d()
    {
      DBG_ASSERT(IsValid(curr_y_d_));
      return curr_y_d_;
    }

    /** Multipliers for the inequality constraints
     *  - Note, inequalities made equality 
     *  by introduction of slacks)
     *  (trial calculations) */
    SmartPtr<const Vector> trial_y_d()
    {
      DBG_ASSERT(IsValid(trial_y_d_));
      return trial_y_d_;
    }

    /** Multipliers for the lower bound on x
     *  (current iteration) */
    SmartPtr<const Vector> curr_z_L()
    {
      DBG_ASSERT(IsValid(curr_z_L_));
      return curr_z_L_;
    }

    /** Multipliers for the lower bound on x
     *  (trial calculations) */
    SmartPtr<const Vector> trial_z_L()
    {
      DBG_ASSERT(IsValid(trial_z_L_));
      return trial_z_L_;
    }

    /** Multipliers for the upper bound on x
     *  (current iteration) */
    SmartPtr<const Vector> curr_z_U()
    {
      DBG_ASSERT(IsValid(curr_z_U_));
      return curr_z_U_;
    }

    /** Multipliers for the upper bound on x
     * (trial calculations) */
    SmartPtr<const Vector> trial_z_U()
    {
      DBG_ASSERT(IsValid(trial_z_L_));
      return trial_z_U_;
    }

    /** Multipliers for the lower bound on s
     *  (current iteration) */
    SmartPtr<const Vector> curr_v_L()
    {
      DBG_ASSERT(IsValid(curr_v_L_));
      return curr_v_L_;
    }

    /** Multipliers for the lower bound on s
     *  (trial calculations) */
    SmartPtr<const Vector> trial_v_L()
    {
      DBG_ASSERT(IsValid(trial_v_L_));
      return trial_v_L_;
    }

    /** Multipliers for the upper bound on s
     *  (current iteration) */
    SmartPtr<const Vector> curr_v_U()
    {
      DBG_ASSERT(IsValid(curr_v_U_));
      return curr_v_U_;
    }

    /** Multipliers for the upper bound on s
     * (trial calculations) */
    SmartPtr<const Vector> trial_v_U()
    {
      DBG_ASSERT(IsValid(trial_v_U_));
      return trial_v_U_;
    }

    /** Hessian or Hessian approximation (do not hold on to it, it might be changed) */
    SmartPtr<const SymMatrix> W()
    {
      DBG_ASSERT(IsValid(W_));
      return W_;
    }

    /** Set Hessian approximation */
    void Set_W(SmartPtr<const SymMatrix> W)
    {
      W_ = W;
      ;
    }

    /** @name ("Main") Search direction.  They are used to store the
     *  search directions computed from solving the primal-dual
     *  system, and can be used in the line search.  They are
     *  overwritten in every iteration, so do not hold on to the
     *  pointers (make copies instead)*/
    //@{
    SmartPtr<const Vector> delta_x()
    {
      return ConstPtr(delta_x_);
    }

    SmartPtr<const Vector> delta_s()
    {
      return ConstPtr(delta_s_);
    }

    SmartPtr<const Vector> delta_y_c()
    {
      return ConstPtr(delta_y_c_);
    }

    SmartPtr<const Vector> delta_y_d()
    {
      return ConstPtr(delta_y_d_);
    }

    SmartPtr<const Vector> delta_z_L()
    {
      return ConstPtr(delta_z_L_);
    }

    SmartPtr<const Vector> delta_z_U()
    {
      return ConstPtr(delta_z_U_);
    }

    SmartPtr<const Vector> delta_v_L()
    {
      return ConstPtr(delta_v_L_);
    }

    SmartPtr<const Vector> delta_v_U()
    {
      return ConstPtr(delta_v_U_);
    }

    SmartPtr<Vector> NonConst_delta_x()
    {
      return delta_x_;
    }

    SmartPtr<Vector> NonConst_delta_s()
    {
      return delta_s_;
    }

    SmartPtr<Vector> NonConst_delta_y_c()
    {
      return delta_y_c_;
    }

    SmartPtr<Vector> NonConst_delta_y_d()
    {
      return delta_y_d_;
    }

    SmartPtr<Vector> NonConst_delta_z_L()
    {
      return delta_z_L_;
    }

    SmartPtr<Vector> NonConst_delta_z_U()
    {
      return delta_z_U_;
    }

    SmartPtr<Vector> NonConst_delta_v_L()
    {
      return delta_v_L_;
    }

    SmartPtr<Vector> NonConst_delta_v_U()
    {
      return delta_v_U_;
    }

    void SetFromPtr_delta_x(SmartPtr<Vector>& delta_x)
    {
      delta_x_ = delta_x;
    }

    void SetFromPtr_delta_s(SmartPtr<Vector>& delta_s)
    {
      delta_s_ = delta_s;
    }

    void SetFromPtr_delta_y_c(SmartPtr<Vector>& delta_y_c)
    {
      delta_y_c_ = delta_y_c;
    }

    void SetFromPtr_delta_y_d(SmartPtr<Vector>& delta_y_d)
    {
      delta_y_d_ = delta_y_d;
    }

    void SetFromPtr_delta_z_L(SmartPtr<Vector>& delta_z_L)
    {
      delta_z_L_ = delta_z_L;
    }

    void SetFromPtr_delta_z_U(SmartPtr<Vector>& delta_z_U)
    {
      delta_z_U_ = delta_z_U;
    }

    void SetFromPtr_delta_v_L(SmartPtr<Vector>& delta_v_L)
    {
      delta_v_L_ = delta_v_L;
    }

    void SetFromPtr_delta_v_U(SmartPtr<Vector>& delta_v_U)
    {
      delta_v_U_ = delta_v_U;
    }

    /** ToDo This is currently a hack and should be handled differently */
    bool HaveDeltas() const
    {
      return have_deltas_;
    }

    void SetHaveDeltas(bool have_deltas)
    {
      have_deltas_ = have_deltas;
    }
    //@}

    /** @name Public Methods for updating iterates */
    //@{
    void SetTrialPrimalVariables(const Vector& x, const Vector& s);
    void SetTrialXVariables(const Vector& x);
    void SetTrialSVariables(const Vector& s);
    void SetTrialEqMultipliers(const Vector& y_c, const Vector& y_d);
    void SetTrialBoundMultipliers(const Vector& z_L, const Vector& z_U,
                                  const Vector& v_L, const Vector& v_U);
    /** Set the values of the primal trial variables (x and s) from
     *  provided Step with step length alpha.
     */
    void SetTrialPrimalVariablesFromStep(Number alpha,
                                         const Vector& delta_x,
                                         const Vector& delta_s);
    /** Set the values of the trial values for the equality constraint
     *  multipliers (y_c and y_d) from provided step with step length
     *  alpha.
     */
    void SetTrialEqMultipilersFromStep(Number alpha,
                                       const Vector& delta_y_c,
                                       const Vector& delta_y_d);
    /** Set the value of the trial values for the bound multipliers
     *  (z_L, z_U, v_L, v_U) from provided step with step length
     *  alpha.
     */
    void SetTrialBoundMutlipliersFromStep(Number alpha,
                                          const Vector& delta_z_L,
                                          const Vector& delta_z_U,
                                          const Vector& delta_v_L,
                                          const Vector& delta_v_U);

    /** Copy the trial values to the current values */
    void CopyTrialToCurrent();

    /** Set the current iterate values from the
     *  trial values. */
    void AcceptTrialPoint();

    /** Set the pointers for the primal trial variables.  Make sure
    that you are not modifying the incoming vectors afterwards! */
    void SetTrialPrimalVariablesFromPtr(const SmartPtr<const Vector>& xptr,
                                        const SmartPtr<const Vector>& sptr);
    /** Set the pointers for the trial bound multipliers.  Make sure
    that you are not modifying the incoming vectors afterwards! */
    void SetTrialBoundMultipliersFromPtr(const SmartPtr<const Vector>& z_Lptr,
                                         const SmartPtr<const Vector>& z_Uptr,
                                         const SmartPtr<const Vector>& v_Lptr,
                                         const SmartPtr<const Vector>& v_Uptr);
    //@}

    /** @name General algorithmic data */
    //@{
    Index iter_count() const
    {
      return iter_count_;
    }
    void Set_iter_count(Index iter_count)
    {
      iter_count_ = iter_count;
    }

    Number curr_mu() const
    {
      DBG_ASSERT(mu_initialized_);
      return curr_mu_;
    }
    void Set_mu(Number mu)
    {
      curr_mu_ = mu;
      mu_initialized_ = true;
    }
    bool MuInitialized() const
    {
      return mu_initialized_;
    }

    Number curr_tau() const
    {
      DBG_ASSERT(tau_initialized_);
      return curr_tau_;
    }
    void Set_tau(Number tau)
    {
      curr_tau_ = tau;
      tau_initialized_ = true;
    }
    bool TauInitialized() const
    {
      return tau_initialized_;
    }

    void SetFreeMuMode(bool free_mu_mode)
    {
      free_mu_mode_ = free_mu_mode;
    }
    bool FreeMuMode() const
    {
      return free_mu_mode_;
    }
    //@}

    /**@name Algorithm Parameters (these will later be
     *  moved to a specific IpoptParameters class, once
     *  I get time to write one) */
    //@{
    Number epsilon_tol() const
    {
      DBG_ASSERT(initialize_called_);
      return epsilon_tol_;
    }
    //@}

    /** @name Information gathered for iteration output */
    //@{
    Number info_regu_x() const
    {
      return info_regu_x_;
    }
    void Set_info_regu_x(Number regu_x)
    {
      info_regu_x_ = regu_x;
    }
    Number info_alpha_primal() const
    {
      return info_alpha_primal_;
    }
    void Set_info_alpha_primal(Number alpha_primal)
    {
      info_alpha_primal_ = alpha_primal;
    }
    char info_alpha_primal_char() const
    {
      return info_alpha_primal_char_;
    }
    void Set_info_alpha_primal_char(char info_alpha_primal_char)
    {
      info_alpha_primal_char_ = info_alpha_primal_char;
    }
    Number info_alpha_dual() const
    {
      return info_alpha_dual_;
    }
    void Set_info_alpha_dual(Number alpha_dual)
    {
      info_alpha_dual_ = alpha_dual;
    }
    Index info_ls_count() const
    {
      return info_ls_count_;
    }
    void Set_info_ls_count(Index ls_count)
    {
      info_ls_count_ = ls_count;
    }
    bool info_skip_output() const
    {
      return info_skip_output_;
    }
    void Append_info_string(const std::string& add_str)
    {
      info_string_ += add_str;
    }
    const std::string& info_string() const
    {
      return info_string_;
    }
    /** Set this to true, if the next time when output is written, the
     *  summary line should not be printed. */
    void Set_info_skip_output(bool info_skip_output)
    {
      info_skip_output_ = info_skip_output;
    }

    /** Reset all info fields */
    void ResetInfo()
    {
      info_regu_x_ = 0;
      info_alpha_primal_ = 0;
      info_alpha_dual_ = 0.;
      info_alpha_primal_char_ = ' ';
      info_ls_count_ = -1;
      info_skip_output_ = false;
      info_string_.clear();
    }
    //@}

  private:
    /** @name Iterates */
    //@{
    /** Main iteration variables
     * (current iteration) */
    SmartPtr<const Vector> curr_x_;

    /** Main iteration variables
     *  (trial calculations) */
    SmartPtr<const Vector> trial_x_;

    /** Slacks from inequality constraints
     *  (current iteration) */
    SmartPtr<const Vector> curr_s_;

    /** Slacks from inequality constraints
     *  (trial calculations) */
    SmartPtr<const Vector> trial_s_;

    /** Multipliers for equality constraints
     *  (current iteration) */
    SmartPtr<const Vector> curr_y_c_;

    /** Multipliers for equality constraints
     *  (trial iteration) */
    SmartPtr<const Vector> trial_y_c_;

    /** Multipliers for the inequality constraints
     *  - Note, inequalities made equality 
     *  by introduction of slacks)
     *  (current iteration) */
    SmartPtr<const Vector> curr_y_d_;

    /** Multipliers for the inequality constraints
     *  - Note, inequalities made equality 
     *  by introduction of slacks)
     *  (trial calculations) */
    SmartPtr<const Vector> trial_y_d_;

    /** Multipliers for the lower bound on x
     *  (current iteration) */
    SmartPtr<const Vector> curr_z_L_;

    /** Multipliers for the lower bound on x
     *  (trial calculations) */
    SmartPtr<const Vector> trial_z_L_;

    /** Multipliers for the upper bound on x
     *  (current iteration) */
    SmartPtr<const Vector> curr_z_U_;

    /** Multipliers for the upper bound on x
     * (trial calculations) */
    SmartPtr<const Vector> trial_z_U_;

    /** Multipliers for the lower bound on s
     *  (current iteration) */
    SmartPtr<const Vector> curr_v_L_;

    /** Multipliers for the lower bound on s
     *  (trial calculations) */
    SmartPtr<const Vector> trial_v_L_;

    /** Multipliers for the upper bound on s
     *  (current iteration) */
    SmartPtr<const Vector> curr_v_U_;

    /** Multipliers for the upper bound on s
     * (trial calculations) */
    SmartPtr<const Vector> trial_v_U_;

    /** Hessian (approximation) - might be changed elsewhere! */
    SmartPtr<const SymMatrix> W_;

    /** @name Delta Variables (steps) */
    //@{
    SmartPtr<Vector> delta_x_;
    SmartPtr<Vector> delta_s_;
    SmartPtr<Vector> delta_y_c_;
    SmartPtr<Vector> delta_y_d_;
    SmartPtr<Vector> delta_z_L_;
    SmartPtr<Vector> delta_z_U_;
    SmartPtr<Vector> delta_v_L_;
    SmartPtr<Vector> delta_v_U_;
    //@}

    /** iteration count */
    Index iter_count_;

    /** current barrier parameter */
    Number curr_mu_;
    bool mu_initialized_;

    /** current fraction to the boundary parameter */
    Number curr_tau_;
    bool tau_initialized_;

    /** flag indicating if Initialize method has been called (for
     *  debugging) */
    bool initialize_called_;

    /** flag for debugging whether we have already curr_ values
     *  available (from which new Vectors can be generated */
    bool have_prototypes_;

    /** The following flag is set to true, if some other part of the
     *  algorithm (like the optprobing heuristic) has already computed
     *  the search direction.  This flag is reset when the
     *  AcceptTrialPoint method is called. */
    bool have_deltas_;

    /** @name Global algorithm parameters.  Those are options that can
     *  be modified by the user and appear at different places in the
     *  algorithm.  They are set using an OptionsList object in the
     *  Initialize method.  */
    //@{
    /** Overall convergence tolerance */
    Number epsilon_tol_;
    //@}

    /** @name Gathered information for iteration output */
    //@{
    /** Size of regularization for the Hessian */
    Number info_regu_x_;
    /** Primal step size */
    Number info_alpha_primal_;
    /** Info character for primal step size */
    char info_alpha_primal_char_;
    /** Dual step size */
    Number info_alpha_dual_;
    /** Number of backtracking trial steps */
    Index info_ls_count_;
    /** true, if next summary output line should not be printed (eg
     *  after restoration phase. */
    bool info_skip_output_;
    /** any string of characters for the end of the output line */
    std::string info_string_;
    /** flag indicating whether the algorithm is in the free mu mode */
    bool free_mu_mode_;
    //@}

    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    IpoptData(const IpoptData&);

    /** Overloaded Equals Operator */
    void operator=(const IpoptData&);
    //@}

  };

} // namespace Ipopt

#endif
