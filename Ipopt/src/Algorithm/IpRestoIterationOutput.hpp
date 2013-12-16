// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter       IBM    2004-09-27

#ifndef __IPRESTOITERATIONOUTPUT_HPP__
#define __IPRESTOITERATIONOUTPUT_HPP__

#include "IpIterationOutput.hpp"
#include "IpOrigIterationOutput.hpp"

namespace Ipopt
{

  /** Class for the iteration summary output for the restoration
   *  phase.  This prints information for the ORIGINAL NLP (and
   *  possibly for the restoration phase NLP.
   */
  class RestoIterationOutput: public IterationOutput
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  If resto_orig_iteration_output is not NULL, the
     *  output will be done twice per iteration, first for the
     *  restoration phase problem, and secondly using the functions
     *  for the original NLP. */
    RestoIterationOutput(const SmartPtr<OrigIterationOutput>& resto_orig_iteration_output);

    /** Default destructor */
    virtual ~RestoIterationOutput();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method to do all the summary output per iteration.  This
     *  include the one-line summary output as well as writing the
     *  details about the iterates if desired */
    virtual void WriteOutput();

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    RestoIterationOutput();

    /** Copy Constructor */
    RestoIterationOutput(const RestoIterationOutput&);

    /** Overloaded Equals Operator */
    void operator=(const RestoIterationOutput&);
    //@}

    /** Pointer to output strategy object during regular iterations. */
    SmartPtr<OrigIterationOutput> resto_orig_iteration_output_;

    /** Flag indicating weather info string should be printed at end
     *  of iteration summary line. */
    bool print_info_string_;

    /** Option indication what should be printed in inf_pr column */
    InfPrOutput inf_pr_output_;

    /** Option indicating at which iteration frequency the summary line should be printed */
    int print_frequency_iter_;
    /** Option indicating at which time frequency the summary line should be printed */
    Number print_frequency_time_;
  };

} // namespace Ipopt

#endif
