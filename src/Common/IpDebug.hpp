// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPDEBUG_HPP__
#define __IPDEBUG_HPP__

#include "IpoptConfig.h"
#include "IpTypes.hpp"

#ifndef IPOPT_CHECKLEVEL
#define IPOPT_CHECKLEVEL 0
#endif

#if IPOPT_CHECKLEVEL > 0
# ifdef NDEBUG
#  undef NDEBUG
# endif
# include <cassert>
# define DBG_ASSERT(test) assert(test)
# define DBG_ASSERT_EXCEPTION(__condition, __except_type, __msg) \
   ASSERT_EXCEPTION( (__condition), __except_type, __msg);
# define DBG_DO(__cmd) __cmd
#else
# define DBG_ASSERT(test)
# define DBG_ASSERT_EXCEPTION(__condition, __except_type, __msg)
# define DBG_DO(__cmd)
#endif

#ifndef IPOPT_VERBOSITY
#define IPOPT_VERBOSITY 0
#endif

#if IPOPT_VERBOSITY < 1
# define DBG_START_FUN(__func_name, __verbose_level)
# define DBG_START_METH(__func_name, __verbose_level)
# define DBG_PRINT(__printf_args)
# define DBG_PRINT_VECTOR(__verbose_level, __vec_name, __vec)
# define DBG_PRINT_MATRIX(__verbose_level, __mat_name, __mat)
# define DBG_EXEC(__verbosity, __cmd)
# define DBG_VERBOSITY() 0
#else
#include <string>

namespace Ipopt
{
// forward definition
class Journalist;

/** Class that lives throughout the execution of a method or
 *  function for which debug output is to be generated.
 *
 *  The output
 *  is sent to the unique debug journalist that is set with
 *  SetJournalist at the beginning of program execution.
 */
class IPOPTLIB_EXPORT DebugJournalistWrapper
{
public:
   /** @name Constructors/Destructors. */
   ///@{
   DebugJournalistWrapper(
      const std::string& func_name,
      Index              verbose_level
   );

   DebugJournalistWrapper(
      const std::string& func_name,
      Index              verbose_level,
      const void* const  method_owner
   );
   ~DebugJournalistWrapper();
   ///@}

   /** @name accessor methods */
   ///@{
   Index Verbosity()
   {
      return verbose_level_;
   }
   const Journalist* Jnlst()
   {
      return jrnl_;
   }
   Index IndentationLevel()
   {
      return indentation_level_;
   }
   ///@}

   /** Printing */
#ifdef __GNUC__
   __attribute__((format(printf, 3, 4)))
#endif
   void DebugPrintf(
      Index       verbosity,
      const char* pformat,
      ...
   );

private:
   friend class IpoptApplication;
   /* Method for initialization of the static GLOBAL journalist,
    * through with all debug printout is to be written.
    *
    * This needs to be set before any debug printout can be done.
    * It is expected that this is only called by the IpoptApplication constructor.
    */
   static void SetJournalist(
      Journalist* jrnl
   );

   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    *
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called.
    */
   ///@{
   /** default constructor */
   DebugJournalistWrapper();

   /** copy contructor */
   DebugJournalistWrapper(
      const DebugJournalistWrapper&
   );

   /** Default Assignment Operator */
   DebugJournalistWrapper& operator=(
      const DebugJournalistWrapper&
   );
   ///@}

   static Index indentation_level_;
   std::string func_name_;
   Index verbose_level_;
   const void* method_owner_;

   static Journalist* jrnl_;
};
}

# define DBG_START_FUN(__func_name, __verbose_level) \
  DebugJournalistWrapper dbg_jrnl((__func_name), (__verbose_level)); \

# define DBG_START_METH(__func_name, __verbose_level) \
  DebugJournalistWrapper dbg_jrnl((__func_name), (__verbose_level), this);

# define DBG_PRINT(__args) \
  dbg_jrnl.DebugPrintf __args;

# define DBG_EXEC(__verbose_level, __cmd) \
  if (dbg_jrnl.Verbosity() >= (__verbose_level)) { \
    (__cmd); \
  }

# define DBG_VERBOSITY() \
  dbg_jrnl.Verbosity()

#endif

#endif
