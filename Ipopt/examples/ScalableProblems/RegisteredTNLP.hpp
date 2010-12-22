// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2005-10-20

#ifndef __REGISTEREDTNLPS_HPP__
#define __REGISTEREDTNLPS_HPP__

#include "IpTNLP.hpp"
#include <map>

using namespace Ipopt;

/** Class implemented the NLP discretization of
 *
 */
class RegisteredTNLP : public TNLP
{
public:
  /** Initialize internal parameters, where N is a parameter
   *  determining the problme size.  This returns false, if N has an
   *  invalid value. */
  virtual bool InitializeProblem(Index N) = 0;
};

class RegisteredTNLPs
{
public:
  RegisteredTNLPs(const SmartPtr<RegisteredTNLP>& tnlp, const std::string name)
  {
    RegisterTNLP(tnlp, name);
  }
  virtual ~RegisteredTNLPs()
  {}
  static SmartPtr<RegisteredTNLP> GetTNLP(const std::string name);
  static void PrintRegisteredProblems();
private:
  void RegisterTNLP(const SmartPtr<RegisteredTNLP>& tnlp,
                    const std::string name);
  SmartPtr<RegisteredTNLP> tnlp_;
};

#define REGISTER_TNLP(class_constructor, name) \
class RegisteredTNLP_Setup_ ## name : public RegisteredTNLPs \
{ \
public: \
  RegisteredTNLP_Setup_ ## name() \
    : \
    RegisteredTNLPs(new class_constructor, #name) \
  { } \
  RegisteredTNLP_Setup_ ## name* KeepCompilerFromRemovingThis(); \
}; \
 \
RegisteredTNLP_Setup_ ## name RegisteredTNLP_Setup_ ## name ## instance_; \
RegisteredTNLP_Setup_ ## name* \
RegisteredTNLP_Setup_ ## name::KeepCompilerFromRemovingThis() \
{ return &RegisteredTNLP_Setup_ ## name ## instance_; }


//static RegisteredTNLP_Setup_ ## name RegisteredTNLP_Setup_ ## name ## instance
#endif
