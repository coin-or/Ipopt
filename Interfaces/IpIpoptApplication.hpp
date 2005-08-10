// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPIPOPTAPPLICATION_HPP__
#define __IPIPOPTAPPLICATION_HPP__

#include "IpJournalist.hpp"
#include "IpTNLP.hpp"
#include "IpRegOptions.hpp"
#include "IpOptionsList.hpp"

namespace Ipopt
{
  /** Return codes for the Optimize call for an application */
  enum ApplicationReturnStatus
  {
    Solve_Succeeded,
    Maximum_Iterations_Exceeded,
    Solved_To_Best_Possible_Precision,
    Solved_To_Acceptable_Level,
    Infeasible_Problem_Detected,
    Restoration_Failed,
    Not_Enough_Degrees_Of_Freedom,
    Unrecoverable_Exception,
    NonIpopt_Exception_Thrown,
    Insufficient_Memory,
    Internal_Error
  };

  DeclareOptionsRegistrar(IpoptApplication);

  /* forward declarations */
  class IpoptData;
  class IpoptCalculatedQuantities;
  class NLP;

  /** This is the main application class for making calls
  *     to Ipopt. */
  class IpoptApplication : public ReferencedObject
  {
  public:
    IpoptApplication(bool read_params_dat = true, bool create_console_out = true);

    virtual ~IpoptApplication();

    /**@name Solve methods */
    //@{
    /** Solve a problem that inherits from TNLP */
    ApplicationReturnStatus OptimizeTNLP(const SmartPtr<TNLP>& nlp);

    /** Solve a problem that inherits from TNLP
     *   use this method when you want access to ip_data and ip_cq after the solve
     */
    ApplicationReturnStatus OptimizeTNLP(const SmartPtr<TNLP>& nlp,
                                         SmartPtr<IpoptData>& ip_data,
                                         SmartPtr<IpoptCalculatedQuantities>& ip_cq);

    /** Solve a problem that inherits from NLP */
    ApplicationReturnStatus OptimizeNLP(const SmartPtr<NLP>& nlp);

    /** Solve a problem that inherits from NLP
     *   use this method when you want access to ip_data and ip_cq after the solve
     */
    ApplicationReturnStatus OptimizeNLP(const SmartPtr<NLP>& nlp,
                                        SmartPtr<IpoptData>& ip_data,
                                        SmartPtr<IpoptCalculatedQuantities>& ip_cq);
    //@}

    /**@name Accessor methods */
    //@{
    /** Get the Journalist for printing output */
    SmartPtr<Journalist> Jnlst()
    {
      return jnlst_;
    }

    /** Get the options list for setting options */
    SmartPtr<OptionsList> Options()
    {
      return options_;
    }
    //@}

    /** @name Methods for IpoptTypeInfo */
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
    // IpoptApplication();

    /** Copy Constructor */
    IpoptApplication(const IpoptApplication&);

    /** Overloaded Equals Operator */
    void operator=(const IpoptApplication&);
    //@}

    /**@name Variables that customize the application behavior */
    //@{
    /** Decide whether or not the PARAMS.DAT file should be read */
    bool read_params_dat_;
    //@}

    /** Journalist for reporting output */
    SmartPtr<Journalist> jnlst_;

    /** OptionsList used for the application */
    SmartPtr<OptionsList> options_;
  };

} // namespace Ipopt

#endif
