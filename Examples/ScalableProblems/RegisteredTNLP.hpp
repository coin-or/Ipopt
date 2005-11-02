// Copyright (C) 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
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

class RegisteredTNLPList
{
public:
  static void RegisterTNLP(const SmartPtr<RegisteredTNLP>& tnlp,
			   const std::string name);
  static SmartPtr<RegisteredTNLP> GetTNLP(const std::string name);
  static void PrintRegisteredProblems();
private:
};

#define REGISTER_TNLP(class_constructor, name) \
class RegisteredTNLP_Setup_ ## name \
{ \
public: \
  RegisteredTNLP_Setup_ ## name() \
  { \
    RegisteredTNLPList::RegisterTNLP(new class_constructor, #name); \
  } \
}; \
 \
static RegisteredTNLP_Setup_ ## name RegisteredTNLP_Setup_ ## name ## instance

#endif
