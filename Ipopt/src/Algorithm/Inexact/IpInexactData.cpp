// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-08-31

#include "IpInexactData.hpp"

namespace Ipopt
{
  InexactData::InexactData()
  {}

  InexactData::~InexactData()
  {}

  bool
  InexactData::Initialize(const Journalist& jnlst,
                          const OptionsList& options,
                          const std::string& prefix)
  {
    full_step_accepted_ = false;
    return true;
  }

  bool
  InexactData::InitializeDataStructures()
  {
    return true;
  }

  void
  InexactData::AcceptTrialPoint()
  {
    // delete data
    normal_x_ = NULL;
    normal_s_ = NULL;
  }

}
