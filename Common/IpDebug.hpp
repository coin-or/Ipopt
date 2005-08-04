// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPDEBUG_HPP__
#define __IPDEBUG_HPP__

#include "config.h"
#include "IpTypes.hpp"
#ifdef OLD_C_HEADERS
# include <assert.h>
#else
# include <cassert>
#endif

#ifdef IP_DEBUG
# define DBG_ASSERT(test) assert(test)
# define DBG_ASSERT_EXCEPTION(__condition, __except_type, __msg) \
   ASSERT_EXCEPTION( (__condition), __except_type, __msg);
# define DBG_SET_VERBOSITY(__level) \
   static const Index dbg_verbosity = __level;
#else
# define DBG_ASSERT(test)
# define DBG_ASSERT_EXCEPTION(__condition, __except_type, __msg)
# define DBG_SET_VERBOSITY(__level)
#endif

#ifndef IP_DEBUG
# define DBG_START_FUN(__func_name, __verbose_level)
# define DBG_START_METH(__func_name, __verbose_level)
# define DBG_PRINT(__printf_args)
# define DBG_PRINT_VECTOR(__verbose_level, __vec_name, __vec)
# define DBG_PRINT_MATRIX(__verbose_level, __mat_name, __mat)
# define DBG_EXEC(__verbosity, __cmd)
# define DBG_VERBOSITY() 0
#else
#include <string>
//#include "IpVector.hpp"

namespace Ipopt
{

  // forward declarations
  class Journalist;
  class Vector;
  class Matrix;

  class DebugJournalistWrapper
  {
  public:
    DebugJournalistWrapper(std::string func_name, Index verbose_level);
    DebugJournalistWrapper(std::string func_name, Index verbose_level,
                           const void* const method_owner);
    ~DebugJournalistWrapper();

    Index Verbosity();

    void DebugPrintf(Index verbosity, const char* pformat, ...);

    void DebugPrintVector(const std::string& vec_name, const Vector& vec);

    void DebugPrintMatrix(const std::string& mat_name, const Matrix& mat);

    static void SetJournalist(Journalist* jrnl);

  private:

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

# define DBG_PRINT_VECTOR(__verbose_level, __vec_name, __vec) \
   if (dbg_jrnl.Verbosity() >= (__verbose_level)) { \
      dbg_jrnl.DebugPrintVector(__vec_name, __vec); \
   }

# define DBG_PRINT_MATRIX(__verbose_level, __mat_name, __mat) \
   if (dbg_jrnl.Verbosity() >= (__verbose_level)) { \
      dbg_jrnl.DebugPrintMatrix(__mat_name, __mat); \
   }

# define DBG_EXEC(__verbose_level, __cmd) \
  if (dbg_jrnl.Verbosity() >= (__verbose_level)) { \
    (__cmd); \
  }

# define DBG_VERBOSITY() \
  dbg_jrnl.Verbosity()

#endif


#endif
