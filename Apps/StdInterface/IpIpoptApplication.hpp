// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPIPOPTAPPLICATION_HPP__
#define __IPIPOPTAPPLICATION_HPP__

#include "IpReferenced.hpp"
#include "IpSmartPtr.hpp"
#include "IpJournalist.hpp"
#include "IpTNLP.hpp"
#include "IpNLP.hpp"

namespace Ipopt
{

  /* forward declarations */
  class IpoptData;
  class IpoptCalculatedQuantities;

  /** This is the main application class for making calls
  *     to Ipopt. */
  class IpoptApplication : public ReferencedObject
  {
  public:
    IpoptApplication();

    virtual ~IpoptApplication();

    /**@name Solve methods */
    //@{
    /** Solve a problem that inherits from TNLP */
    ApplicationReturnStatus OptimizeTNLP(const SmartPtr<TNLP>& nlp);
  
    /** Solve a problem that inherits from TNLP
     *   use this method when you want access to ip_data nad ip_cq after the solve
     */
    ApplicationReturnStatus OptimizeTNLP(const SmartPtr<TNLP>& nlp, 
					 SmartPtr<IpoptData>& ip_data, 
					 SmartPtr<IpoptCalculatedQuantities>& ip_cq);

    /** Solve a problem that inherits from NLP */
    ApplicationReturnStatus OptimizeNLP(const SmartPtr<NLP>& nlp);

    /** Solve a problem that inherits from NLP
     *   use this method when you want access to ip_data nad ip_cq after the solve
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

    /**@name Methods to set the options for the application */
    //@{
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
    /** Name of the file for detailed output */
    /** Decide whether or not the PARAMS.DAT file should be read */
    bool read_params_dat_;
    /** Decide whether or not to report the solution to the journalist */
    bool report_solution_;
    /** Decide whether or not to force reporting of the solution to the console */
    bool force_report_solution_to_console_;
    /** Decide whether or not to report statistics */
    bool report_statistics_;
    //@}

    /** Journalist for reporting output */
    SmartPtr<Journalist> jnlst_;

    /** OptionsList used for the application */
    SmartPtr<OptionsList> options_;
  };
} // namespace Ipopt

#endif
