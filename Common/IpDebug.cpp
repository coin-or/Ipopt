// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpDebug.hpp"
#include "IpJournalist.hpp"

#ifdef IP_DEBUG

namespace Ipopt
{

  Index DebugJournalistWrapper::indentation_level_ = 0;
  Journalist* DebugJournalistWrapper::jrnl_ = NULL;

  DebugJournalistWrapper::DebugJournalistWrapper(char* func_name,
      Index verbose_level)
      :
      verbose_level_(verbose_level),
      method_owner_(NULL)
  {
    strcpy(func_name_, func_name);

    // AW:  I changed the following check so that the journalist can interit
    //      from ReferencedObject.
    // DBG_ASSERT(jrnl_);
    if (jrnl_==NULL) {
      verbose_level_ = 0;
      return;
    }
    DebugPrintf(1, "-> Calling to: %s\n", func_name_);
    if (verbose_level_>0) {
      indentation_level_++;
    }
  }

  DebugJournalistWrapper::DebugJournalistWrapper(
    char* func_name, Index verbose_level, const void* const method_owner)
      :
      verbose_level_(verbose_level),
      method_owner_(method_owner)
  {
    strcpy(func_name_, func_name);

    // AW:  I changed the following check so that the journalist can interit
    //      from ReferencedObject.
    // DBG_ASSERT(jrnl_);
    if (jrnl_==NULL) {
      verbose_level_ = 0;
      return;
    }
    DebugPrintf(1, "-> Calling to: %s in obj: 0x%x\n", func_name_,
                method_owner_);
    if (verbose_level_>0) {
      indentation_level_++;
    }
  }

  DebugJournalistWrapper::~DebugJournalistWrapper()
  {
    if (verbose_level_>0) {
      indentation_level_--;
    }
    if (jrnl_) {
      if (method_owner_ == NULL) {
        DebugPrintf(1, "<- Returning from : %s\n", func_name_);
      }
      else {
        DebugPrintf(1, "<- Returning from : %s in obj: 0x%x\n",
                    func_name_, method_owner_);
      }
    }
  }


  void DebugJournalistWrapper::SetJournalist(Journalist* jrnl)
  {

    if (jrnl == NULL) {
      jrnl_->Printf(J_ERROR, J_DBG,
                    "# Setting Journalist to NULL in DebugJournalistWrapper.\n");
      jrnl_ = NULL;
    }
    else if (!jrnl_) {
      jrnl_ = jrnl;
    }
  }


  Index DebugJournalistWrapper::Verbosity()
  {
    return verbose_level_;
  }

  void DebugJournalistWrapper::DebugPrintf(Index verbosity, char* pformat, ...)
  {

    if (Verbosity() >= verbosity) {
      va_list(ap);

      va_start(ap, pformat);

      DBG_ASSERT(jrnl_);
      jrnl_->PrintfIndented(J_ERROR, J_DBG, indentation_level_, "# ");
      jrnl_->VPrintf(J_ERROR, J_DBG, pformat, ap);

      va_end(ap);
    }
  }

  void DebugJournalistWrapper::DebugPrintVector(const std::string& vec_name, const Vector& vec)
  {
    jrnl_->PrintVector(J_ERROR, J_DBG, vec_name, vec, indentation_level_*2, "# ");
  }

  void DebugJournalistWrapper::DebugPrintMatrix(const std::string& mat_name, const Matrix& mat)
  {
    jrnl_->PrintMatrix(J_ERROR, J_DBG, mat_name, mat, indentation_level_*2, "# ");
  }



} // namespace Ipopt

#endif // #ifdef IP_DEBUG
