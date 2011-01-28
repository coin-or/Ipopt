// Copyright (C) 2008, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter           IBM     2008-08-31

#ifndef __IPINEXACTDATA_HPP__
#define __IPINEXACTDATA_HPP__

#include "IpIpoptData.hpp"

namespace Ipopt
{

  /** Class to organize all the additional data required by the
   *  Chen-Goldfarb penalty function algorithm. */
  class InexactData : public IpoptAdditionalData
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    InexactData();

    /** Default destructor */
    ~InexactData();
    //@}

    /** @name Methods overloaded from IpoptAdditionalData */
    //@{
    /** This method must be called to initialize the global
     *  algorithmic parameters.  The parameters are taken from the
     *  OptionsList object. */
    bool Initialize(const Journalist& jnlst,
                    const OptionsList& options,
                    const std::string& prefix);

    /** Initialize Data Structures at the beginning. */
    bool InitializeDataStructures();

    /** Do whatever is necessary to accept a trial point as current
     *  iterate.  This is also used to finish an iteration, i.e., to
     *  release memory, and to reset any flags for a new iteration. */
    void AcceptTrialPoint();
    //@}

    /** @name Normal step set and accessor methods */
    //@{
    void set_normal_x(SmartPtr<Vector>& normal_x)
    {
      normal_x_ = ConstPtr(normal_x);
      normal_x = NULL;
    }
    void set_normal_s(SmartPtr<Vector>& normal_s)
    {
      normal_s_ = ConstPtr(normal_s);
      normal_s = NULL;
    }
    SmartPtr<const Vector> normal_x()
    {
      return normal_x_;
    }
    SmartPtr<const Vector> normal_s()
    {
      return normal_s_;
    }
    //@}

    /** @name Tangential step set and accessor methods */
    //@{
    void set_tangential_x(SmartPtr<const Vector>& tangential_x)
    {
      tangential_x_ = tangential_x;
      tangential_x = NULL;
    }
    void set_tangential_s(SmartPtr<const Vector>& tangential_s)
    {
      tangential_s_ = tangential_s;
      tangential_s = NULL;
    }
    SmartPtr<const Vector> tangential_x()
    {
      return tangential_x_;
    }
    SmartPtr<const Vector> tangential_s()
    {
      return tangential_s_;
    }
    //@}

    /** @name Flag indicating if most recent step has been fully
     *  accepted.  This is used to determine if the trust region
     *  radius should be increased. */
    //@{
    void set_full_step_accepted(bool full_step_accepted)
    {
      full_step_accepted_ = full_step_accepted;
    }
    bool full_step_accepted()
    {
      return full_step_accepted_;
    }
    //@}

    /** @name Current value of penalty parameter */
    //@{
    void set_curr_nu(Number nu)
    {
      curr_nu_ = nu;
    }
    Number curr_nu()
    {
      return curr_nu_;
    }
    //@}

    /** @name Current normal step computation flag */
    //@{
    void set_compute_normal(bool compute_normal)
    {
      compute_normal_ = compute_normal;
    }
    bool compute_normal()
    {
      return compute_normal_;
    }
    //@}

    /** @name Next iteration normal step computation flag */
    //@{
    void set_next_compute_normal(bool next_compute_normal)
    {
      next_compute_normal_ = next_compute_normal;
    }
    bool next_compute_normal()
    {
      return next_compute_normal_;
    }
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
    /** Copy Constructor */
    InexactData(const InexactData&);

    /** Overloaded Equals Operator */
    void operator=(const InexactData&);
    //@}

    /** @name Normal step */
    //@{
    SmartPtr<const Vector> normal_x_;
    SmartPtr<const Vector> normal_s_;
    //@}

    /** @name Tangential step */
    //@{
    SmartPtr<const Vector> tangential_x_;
    SmartPtr<const Vector> tangential_s_;
    //@}

    /** Flag indicating if most recent step has been fully accepted */
    bool full_step_accepted_;

    /** current value of penalty parameter */
    Number curr_nu_;

    /** current normal step computation flag */
    bool compute_normal_;

    /** next iteration normal step computation flag */
    bool next_compute_normal_;
  };

} // namespace Ipopt

#endif
