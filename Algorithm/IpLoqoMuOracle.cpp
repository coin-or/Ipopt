// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpLoqoMuOracle.hpp"

#include <limits>

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  LoqoMuOracle::LoqoMuOracle()
      :
      MuOracle()
  {}

  LoqoMuOracle::~LoqoMuOracle()
  {}

  bool LoqoMuOracle::InitializeImpl(const OptionsList& options,
                                    const std::string& prefix)
  {
    // The following line is only here so that
    // IpoptCalculatedQuantities::CalculateSafeSlack and the first
    // output line have something to work with
    IpData().Set_mu(1.);

    return true;
  }

  Number LoqoMuOracle::CalculateMu()
  {
    DBG_START_METH("LoqoMuOracle::CalculateMu",
                   dbg_verbosity);

    Number avrg_compl = IpCq().curr_avrg_compl();
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  Average complemantarity is %lf\n", avrg_compl);

    Number xi = IpCq().curr_centrality_measure();
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  Xi (distance from uniformity) is %lf\n", xi);

    //Number factor = 1.-tau_min_;   //This is the original values
    Number factor = 0.05;   //This is the value I used otherwise
    Number sigma = 0.1*pow(Min(factor*(1.-xi)/xi,2),3.);

    Number mu = sigma*avrg_compl;
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  Barrier parameter proposed by LOQO rule is %lf\n", mu);

    // DELETEME
    char ssigma[40];
    sprintf(ssigma, " sigma=%8.2e", sigma);
    IpData().Append_info_string(ssigma);
    sprintf(ssigma, " xi=%8.2e ", IpCq().curr_centrality_measure());
    IpData().Append_info_string(ssigma);

    return mu;
  }

} // namespace Ipopt
