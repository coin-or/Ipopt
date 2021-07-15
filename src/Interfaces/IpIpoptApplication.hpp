// Copyright (C) 2004, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPIPOPTAPPLICATION_HPP__
#define __IPIPOPTAPPLICATION_HPP__

#include <iostream>

#include "IpJournalist.hpp"
#include "IpTNLP.hpp"
#include "IpNLP.hpp"
#include "IpReturnCodes.hpp"

namespace Ipopt
{
DECLARE_STD_EXCEPTION(IPOPT_APPLICATION_ERROR);

/* forward declarations */
class IpoptAlgorithm;
class IpoptNLP;
class IpoptData;
class IpoptCalculatedQuantities;
class AlgorithmBuilder;
class RegisteredOptions;
class OptionsList;
class SolveStatistics;

/** This is the main application class for making calls to Ipopt. */
class IPOPTLIB_EXPORT IpoptApplication: public ReferencedObject
{
public:
   IpoptApplication(
      bool create_console_out = true,
      bool create_empty = false
   );

   /** Another constructor that assumes that the code in the
    *  (default) constructor has already been executed
    */
   IpoptApplication(
      SmartPtr<RegisteredOptions> reg_options,
      SmartPtr<OptionsList>       options,
      SmartPtr<Journalist>        jnlst
   );

   virtual ~IpoptApplication();

   /** Method for creating a new IpoptApplication that uses the same
    *  journalist and registered options, and a copy of the options list.
    */
   virtual SmartPtr<IpoptApplication> clone();

   /** Initialization method.
    *
    *  This method reads options from the
    *  input stream and initializes the journalists.
    *
    *  @return Solve_Succeeded or something else if there was a
    *  problem in the initialization (such as an invalid option).
    *
    *  You should call one of the initialization methods at some
    *  point before the first optimize call.
    *  Set @par allow_clobber to true if you want to allow
    *  overwriting options that are set by the input stream.
    */
   virtual ApplicationReturnStatus Initialize(
      std::istream& is,
      bool          allow_clobber = false
   );

   /** Initialization method.
    *
    *  This method reads options from the
    *  params file and initializes the journalists.

    *  @return Solve_Succeeded or something else if there was a
    *  problem in the initialization (such as an invalid option).

    *  You should call one of the initialization methods at some
    *  point before the first optimize call.
    *
    *  @note You can skip the processing of a params file by
    *  setting params_file to "".
    *
    *  Set @par allow_clobber to true if you want to allow
    *  overwriting options that are set by the params file.
    */
   virtual ApplicationReturnStatus Initialize(
      std::string params_file,
      bool        allow_clobber = false
   );

   /** Initialization method.
    *
    *  This method reads options from the
    *  params file and initializes the journalists.

    *  @return Solve_Succeeded or something else if there was a
    *  problem in the initialization (such as an invalid option).
    *
    *  You should call one of the initialization methods at some
    *  point before the first optimize call.
    *
    *  @note You can skip the processing of a params file by
    *  setting params_file to "".
    *
    *  Set @par allow_clobber to true if you want to allow
    *  overwriting options that are set by the params file.
    */
   virtual ApplicationReturnStatus Initialize(
      const char* params_file,
      bool        allow_clobber = false
   )
   {
      return Initialize(std::string(params_file), allow_clobber);
   }

   /** Initialize method.
    *
    *  This method reads the options file specified
    *  by the option_file_name option and initializes the journalists.

    *  @return Solve_Succeeded or something else if there was a
    *  problem in the initialization (such as an invalid option).

    *  You should call this method at some point before the first optimize
    *  call.

    *  Set @par allow_clobber to true if you want to allow
    *  overwriting options that are set by the options file.
    */
   virtual ApplicationReturnStatus Initialize(
      bool allow_clobber = false
   );

   /**@name Solve methods */
   ///@{
   /** Solve a problem that inherits from TNLP */
   virtual ApplicationReturnStatus OptimizeTNLP(
      const SmartPtr<TNLP>& tnlp
   );

   /** Solve a problem that inherits from NLP */
   virtual ApplicationReturnStatus OptimizeNLP(
      const SmartPtr<NLP>& nlp
   );

   /** Solve a problem that inherits from NLP */
   virtual ApplicationReturnStatus OptimizeNLP(
      const SmartPtr<NLP>&        nlp,
      SmartPtr<AlgorithmBuilder>& alg_builder
   );

   /** Solve a problem (that inherits from TNLP) for a repeated time.
    *
    *  The OptimizeTNLP method must have been called before.  The
    *  TNLP must be the same object. The IpoptAlgorithm object from the
    *  previous solve will be reused.
    */
   virtual ApplicationReturnStatus ReOptimizeTNLP(
      const SmartPtr<TNLP>& tnlp
   );

   /** Solve a problem (that inherits from NLP) for a repeated time.
    *
    *  The OptimizeNLP method must have been called before.  The
    *  NLP must be the same object. The IpoptAlgorithm object from the
    *  previous solve will be reused.
    */
   virtual ApplicationReturnStatus ReOptimizeNLP(
      const SmartPtr<NLP>& nlp
   );
   ///@}

   /** Method for opening an output file with given print_level.
    *
    *  @return false if there was a problem
    */
   virtual bool OpenOutputFile(
      std::string  file_name,
      EJournalLevel print_level
   );

   /**@name Accessor methods */
   ///@{
   /** Get the Journalist for printing output */
   virtual SmartPtr<Journalist> Jnlst()
   {
      return jnlst_;
   }

   /** Get a pointer to RegisteredOptions object to add new options */
   virtual SmartPtr<RegisteredOptions> RegOptions()
   {
      return reg_options_;
   }

   /** Get the options list for setting options */
   virtual SmartPtr<OptionsList> Options()
   {
      return options_;
   }

   /** Get the options list for setting options (const version) */
   virtual SmartPtr<const OptionsList> Options() const
   {
      return ConstPtr(options_);
   }

   /** Get the object with the statistics about the most recent
    *  optimization run.
    *
    *  @note Statistics are not available if optimization terminated
    *  with a serious problem, that is, an ApplicationReturnStatus of
    *  Not_Enough_Degrees_Of_Freedom or lower.
    */
   virtual SmartPtr<SolveStatistics> Statistics();

   /** Get the IpoptNLP Object */
   virtual SmartPtr<IpoptNLP> IpoptNLPObject();

   /** Get the IpoptData Object */
   SmartPtr<IpoptData> IpoptDataObject();

   /** Get the IpoptCQ Object */
   virtual SmartPtr<IpoptCalculatedQuantities> IpoptCQObject();

   /** Get the Algorithm Object */
   SmartPtr<IpoptAlgorithm> AlgorithmObject();
   ///@}

   /** Method for printing Ipopt copyright message now instead of
    *  just before the optimization.
    *
    *  If you want to have the copy right message printed earlier
    *  than by default, call this method at the convenient time.
    */
   void PrintCopyrightMessage();

   /** Method to set whether non-ipopt non-bad_alloc non-overflow_error exceptions
    * are rethrown by Ipopt.
    *
    * By default, non-Ipopt and non-bad_alloc and non-overflow_error exceptions are
    * caught by Ipopts initialization and optimization methods
    * and the status NonIpopt_Exception_Thrown is returned.
    * This function allows to enable rethrowing of such exceptions.
    *
    * @return Returns whether non-ipopt exceptions were rethrown before.
    */
   bool RethrowNonIpoptException(
      bool dorethrow
   )
   {
      bool oldval = rethrow_nonipoptexception_;
      rethrow_nonipoptexception_ = dorethrow;
      return oldval;
   }

   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );

   /** Method to register all Ipopt options. */
   static void
   RegisterAllIpoptOptions(
      const SmartPtr<RegisteredOptions>& roptions
   );

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
   ///@{
   /** Copy Constructor */
   IpoptApplication(
      const IpoptApplication&
   );

   /** Default Assignment Operator */
   void operator=(
      const IpoptApplication&
   );
   ///@}

   /** Method for the actual optimize call of the Ipopt algorithm.
    *
    *  This is used both for Optimize and ReOptimize
    */
   ApplicationReturnStatus call_optimize();

   /**@name Variables that customize the application behavior */
   ///@{
   /** Decide whether or not the ipopt.opt file should be read */
   bool read_params_dat_;

   /** Decide whether non-ipopt non-bad_alloc non-overflow_error exceptions should be rethrown */
   bool rethrow_nonipoptexception_;
   ///@}

   /** Journalist for reporting output */
   SmartPtr<Journalist> jnlst_;

   /** RegisteredOptions */
   SmartPtr<RegisteredOptions> reg_options_;

   /** OptionsList used for the application */
   SmartPtr<OptionsList> options_;

   /** Object for storing statistics about the most recent
    *  optimization run.
    */
   SmartPtr<SolveStatistics> statistics_;

   /** Object with the algorithm skeleton.
    */
   SmartPtr<IpoptAlgorithm> alg_;

   /** IpoptNLP Object for the NLP.
    *
    *  We keep this around for a ReOptimize warm start.
    */
   SmartPtr<IpoptNLP> ip_nlp_;

   /** IpoptData Object for the NLP.
    *
    *  We keep this around for a ReOptimize warm start.
    */
   SmartPtr<IpoptData> ip_data_;

   /** IpoptCalculatedQuantities Object for the NLP.
    *
    *  We keep this around for a ReOptimize warm start.
    */
   SmartPtr<IpoptCalculatedQuantities> ip_cq_;

   /** Pointer to the TNLPAdapter used to convert the TNLP to an NLP.
    *
    *  We keep this around for the ReOptimizerTNLP call.
    */
   SmartPtr<NLP> nlp_adapter_;

   /** @name Algorithmic parameters */
   ///@{
   /** Flag indicating if we are to use the inexact linear solver option */
   bool inexact_algorithm_;

   /** Flag indicating if all bounds should be replaced by inequality
    *  constraints.
    *
    *  This is necessary for the inexact algorithm.
    */
   bool replace_bounds_;
   ///@}
};

} // namespace Ipopt

extern "C" IPOPTLIB_EXPORT class Ipopt::IpoptApplication* IPOPT_CALLCONV IpoptApplicationFactory();

#endif
