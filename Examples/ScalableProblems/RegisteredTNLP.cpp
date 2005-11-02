// Copyright (C) 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2005-10-20

#include "RegisteredTNLP.hpp"

static std::map<std::string, SmartPtr<RegisteredTNLP> > tnlp_map_;

void
RegisteredTNLPList::RegisterTNLP(const SmartPtr<RegisteredTNLP>& tnlp,
				 const std::string name)
{
  tnlp_map_[name] = GetRawPtr(tnlp);
}

SmartPtr<RegisteredTNLP>
RegisteredTNLPList::GetTNLP(const std::string name)
{
  SmartPtr<RegisteredTNLP> retval = NULL;
  std::map<std::string, SmartPtr<RegisteredTNLP> >::iterator it;
  it = tnlp_map_.find(name);
  if (it != tnlp_map_.end()) {
    retval = it->second;
  }
  return retval;
}

void
RegisteredTNLPList::PrintRegisteredProblems()
{
  for (std::map<std::string, SmartPtr<RegisteredTNLP> >::iterator it = tnlp_map_.begin();
       it != tnlp_map_.end(); it++) {
    printf("%s\n", it->first.c_str());
  }
}
