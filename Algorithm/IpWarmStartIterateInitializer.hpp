// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter              IBM    2004-09-24

#ifndef __IPWARMSTARTITERATEINITIALIZER_HPP__
#define __IPWARMSTARTITERATEINITIALIZER_HPP__

#include "IpUtils.hpp"
#include "IpIterateInitializer.hpp"
#include "IpEqMultCalculator.hpp"
#include "IpIpoptType.hpp"

namespace Ipopt
{

  DeclareIpoptType(WarmStartIterateInitializer);

  /** Class implementing an initialization procedure for warm starts.
   */
  class WarmStartIterateInitializer: public IterateInitializer
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor. */
    WarmStartIterateInitializer();

    /** Default destructor */
    virtual ~WarmStartIterateInitializer()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Compute the initial iterates and set the into the curr field
     *  of the ip_data object. */
    virtual bool SetInitialIterates();

    /** Methods used by IpoptType */
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
    /** Copy Constructor */
    WarmStartIterateInitializer(const WarmStartIterateInitializer&);

    /** Overloaded Equals Operator */
    void operator=(const WarmStartIterateInitializer&);
    //@}

    /**@name Algorithmic Parameters */
    //@{
    /** Parameters for bumping x0 in warm start mode */
    Number warm_start_bound_push_;
    /** Parameters for bumping x0 in warm start mode */
    Number warm_start_bound_frac_;
    /** Parameters for bumping initial bound multipliers */
    Number warm_start_mult_bound_push_;
    /** Maximal size of entries in bound and equality constraint
     *  multipliers in magnitute.  If chosen less of equal to zero, no
     *  upper limit is imposed.  Otherwise, the entries exceeding the
     *  given limit are set to the value closest to the limit. */
    Number warm_start_mult_init_max_;
    /** Target values for the barrier parameter in warm start option.
     */
    Number warm_start_target_mu_;
    //@}

    /** @name Auxilliary functions */
    //@{
    void push_variables(Number bound_bush,
                        Number bound_frac,
                        std::string name,
                        const Vector& orig_x,
                        SmartPtr<const Vector>& new_x,
                        const Vector& x_L,
                        const Vector& x_U,
                        const Matrix& Px_L,
                        const Matrix& Px_U);

    void process_target_mu(Number factor,
                           const Vector& curr_vars,
                           const Vector& curr_slacks,
                           const Vector& curr_mults,
                           const Matrix& P,
                           SmartPtr<const Vector>& ret_vars,
                           SmartPtr<const Vector>& ret_mults);

    void adapt_to_target_mu(Vector& new_s,
                            Vector& new_z,
                            Number target_mu);
    //@}
  };

} // namespace Ipopt

#endif
