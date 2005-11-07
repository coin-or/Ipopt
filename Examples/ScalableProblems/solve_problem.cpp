// Copyright (C) 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2004-11-05

#include "IpIpoptApplication.hpp"
#include "RegisteredTNLP.hpp"

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

using namespace Ipopt;

// This could probably be done more elegant and automatically, but I
// can't get it to work right now.  For now, list explicitly the
// problems we want to include:
#include "LuksanVlcek1.hpp"
REGISTER_TNLP(LuksanVlcek1(0,0), LukVlE1);
REGISTER_TNLP(LuksanVlcek1(-1.,0.), LukVlI1);
#include "LuksanVlcek2.hpp"
REGISTER_TNLP(LuksanVlcek2(0,0), LukVlE2);
REGISTER_TNLP(LuksanVlcek2(-1.,0.), LukVlI2);
#include "LuksanVlcek3.hpp"
REGISTER_TNLP(LuksanVlcek3(0,0), LukVlE3);
REGISTER_TNLP(LuksanVlcek3(-1.,0.), LukVlI3);
#include "LuksanVlcek4.hpp"
REGISTER_TNLP(LuksanVlcek4(0,0), LukVlE4);
REGISTER_TNLP(LuksanVlcek4(-1.,0.), LukVlI4);
#include "LuksanVlcek5.hpp"
REGISTER_TNLP(LuksanVlcek5(0,0), LukVlE5);
REGISTER_TNLP(LuksanVlcek5(-1.,0.), LukVlI5);
#include "LuksanVlcek6.hpp"
REGISTER_TNLP(LuksanVlcek6(0,0), LukVlE6);
REGISTER_TNLP(LuksanVlcek6(-1.,0.), LukVlI6);
#include "LuksanVlcek7.hpp"
REGISTER_TNLP(LuksanVlcek7(0,0), LukVlE7);
REGISTER_TNLP(LuksanVlcek7(-1.,0.), LukVlI7);


#include "MittelmannBndryCntrlDiri.hpp"
REGISTER_TNLP(MittelmannBndryCntrlDiri1, MBndryCntrl1);
REGISTER_TNLP(MittelmannBndryCntrlDiri2, MBndryCntrl2);
REGISTER_TNLP(MittelmannBndryCntrlDiri3, MBndryCntrl3);
REGISTER_TNLP(MittelmannBndryCntrlDiri4, MBndryCntrl4);

#include "MittelmannBndryCntrlNeum.hpp"
REGISTER_TNLP(MittelmannBndryCntrlNeum1, MBndryCntrl5);
REGISTER_TNLP(MittelmannBndryCntrlNeum2, MBndryCntrl6);
REGISTER_TNLP(MittelmannBndryCntrlNeum3, MBndryCntrl7);
REGISTER_TNLP(MittelmannBndryCntrlNeum4, MBndryCntrl8);

#include "MittelmannDistCntrlDiri.hpp"
REGISTER_TNLP(MittelmannDistCntrlDiri1, MDistCntrl1);
REGISTER_TNLP(MittelmannDistCntrlDiri2, MDistCntrl2);
REGISTER_TNLP(MittelmannDistCntrlDiri3, MDistCntrl3);
REGISTER_TNLP(MittelmannDistCntrlDiri4, MDistCntrl3a);

#include "MittelmannDistCntrlNeumA.hpp"
REGISTER_TNLP(MittelmannDistCntrlNeumA1, MDistCntrl4);
REGISTER_TNLP(MittelmannDistCntrlNeumA2, MDistCntrl5);
REGISTER_TNLP(MittelmannDistCntrlNeumA3, MDistCntrl6a);

#include "MittelmannDistCntrlNeumB.hpp"
REGISTER_TNLP(MittelmannDistCntrlNeumB1, MDistCntrl4a);
REGISTER_TNLP(MittelmannDistCntrlNeumB2, MDistCntrl5a);
REGISTER_TNLP(MittelmannDistCntrlNeumB3, MDistCntrl6);

#include "MittelmannParaCntrl.hpp"
REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_1>, MPara5_1);
REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_2_1>, MPara5_2_1);
REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_2_2>, MPara5_2_2);
REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_2_3>, MPara5_2_3);
//REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_try>, MPara5_try);

static void print_problems()
{
  printf("\nList of all registered problems:\n\n");
  RegisteredTNLPs::PrintRegisteredProblems();
}

int main(int argv, char* argc[])
{
  if (argv==2 && !strcmp(argc[1],"list")) {
    print_problems();
    return 0;
  }

  if (argv!=3) {
    printf("Usage: %s ProblemName N\n", argc[0]);
    printf("          where N is a positive parameter determining problem size\n");
    printf("       %s list\n", argc[0]);
    printf("          to list all registered problems.\n");
    return -1;
  }

  // Create an instance of your nlp...
  SmartPtr<RegisteredTNLP> tnlp =
    RegisteredTNLPs::GetTNLP(argc[1]);
  if (!IsValid(tnlp)) {
    printf("Problem with name \"%s\" not known.\n", argc[1]);
    print_problems();
    return -2;
  }

  Index N = atoi(argc[2]);
  if (N <= 0) {
    printf("Given problem size is invalid.\n");
    return -3;
  }

  bool retval = tnlp->InitializeProblem(N);
  if (!retval) {
    printf("Cannot initialize problem.  Abort.\n");
    return -4;
  }

  // Create an instance of the IpoptApplication
  SmartPtr<IpoptApplication> app = new IpoptApplication();

  ApplicationReturnStatus status = app->OptimizeTNLP(GetRawPtr(tnlp));

  return (int) status;
}
