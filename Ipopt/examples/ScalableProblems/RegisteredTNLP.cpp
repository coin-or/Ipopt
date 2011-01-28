// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2005-10-20

#include "RegisteredTNLP.hpp"

#include <cstdio>

std::map<std::string, SmartPtr<RegisteredTNLP> >& RegisteredTNLPListMap()
{
  static std::map<std::string, SmartPtr<RegisteredTNLP> > tnlp_map_;
  return tnlp_map_;
}

void
RegisteredTNLPs::RegisterTNLP(const SmartPtr<RegisteredTNLP>& tnlp,
                              const std::string name)
{
  RegisteredTNLPListMap()[name] = GetRawPtr(tnlp);
}

SmartPtr<RegisteredTNLP>
RegisteredTNLPs::GetTNLP(const std::string name)
{
  SmartPtr<RegisteredTNLP> retval = NULL;
  std::map<std::string, SmartPtr<RegisteredTNLP> >::iterator it;
  it = RegisteredTNLPListMap().find(name);
  if (it != RegisteredTNLPListMap().end()) {
    retval = it->second;
  }
  return retval;
}

void
RegisteredTNLPs::PrintRegisteredProblems()
{
  for (std::map<std::string, SmartPtr<RegisteredTNLP> >::iterator it = RegisteredTNLPListMap().begin();
       it != RegisteredTNLPListMap().end(); it++) {
    printf("%s\n", it->first.c_str());
  }
}
