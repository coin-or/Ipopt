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
#include "IpTNLP.hpp"
#include "IpNLP.hpp"

namespace Ipopt
{
   enum ApplicationReturnStatus
   {
      Solve_Succeeded,
      Maximum_Iterations_Exceeded,
      Solve_Failed,
      NonIpopt_Exception_Thrown,
      Internal_Error
   };

   /** This is the main application class for making calls
   *     to Ipopt. */
   class IpoptApplication : public ReferencedObject
   {
   public:
      IpoptApplication();

      virtual ~IpoptApplication();


      ApplicationReturnStatus OptimizeTNLP(const SmartPtr<TNLP>& nlp, 
            const SmartPtr<OptionsList> additional_options = NULL);

      ApplicationReturnStatus OptimizeNLP(const SmartPtr<NLP>& nlp,
            const SmartPtr<OptionsList> additional_options = NULL);

      /**@name Methods to set the options for the application */
      //@{
      //@}


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
    std::string output_file_;
    /** Print level of the file for detailed output */
    EJournalLevel output_file_print_level_;
    /** Decide whether or not the PARAMS.DAT file should be read */
    bool read_params_dat_;
    /** Decide whether or not to report the algorithm return status to the journalist */
    bool report_solve_status_;
    /** Decide whether or not to report the solution to the journalist */
    bool report_solution_;
    /** Decide whether or not to force reporting of the solution to the console */
    bool force_report_solution_to_console_;
    /** Decide whether or not to report statistics */
    bool report_statistics_;
   //@}
   };
} // namespace Ipopt

#endif
