// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2004-11-05

#include "IpIpoptApplication.hpp"
#include "RegisteredTNLP.hpp"

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "configall_system.h"
#endif

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# else
#  error "don't have header file for stdlib"
# endif
#endif

#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif

#include <iostream>

//**********************************************************************
// Stuff for benchmarking
// #define TIME_LIMIT
#ifdef TIME_LIMIT

#include <unistd.h>
#include <pthread.h>

extern "C" void* killer_thread(void* arg)
{
  int runtime = *reinterpret_cast<int *>(arg);
  if (runtime <= 0) {
    printf("Invalid argument for run time (%d)\n", runtime);
    exit(-100);
  }
  printf("Limiting wall clock time to %d seconds.\n", runtime);
  sleep(runtime);
  printf("EXIT: Exceeding wall clock time limit of %d seconds.\n", runtime);
  exit(-999);
  return NULL;
}
#endif
//**********************************************************************


using namespace Ipopt;
using namespace std;

// This could probably be done more elegant and automatically, but I
// can't get it to work right now.  For now, list explicitly the
// problems we want to include:
#include "LuksanVlcek1.hpp"
REGISTER_TNLP(LuksanVlcek1(0,0), LukVlE1)
REGISTER_TNLP(LuksanVlcek1(-1.,0.), LukVlI1)
#include "LuksanVlcek2.hpp"
REGISTER_TNLP(LuksanVlcek2(0,0), LukVlE2)
REGISTER_TNLP(LuksanVlcek2(-1.,0.), LukVlI2)
#include "LuksanVlcek3.hpp"
REGISTER_TNLP(LuksanVlcek3(0,0), LukVlE3)
REGISTER_TNLP(LuksanVlcek3(-1.,0.), LukVlI3)
#include "LuksanVlcek4.hpp"
REGISTER_TNLP(LuksanVlcek4(0,0), LukVlE4)
REGISTER_TNLP(LuksanVlcek4(-1.,0.), LukVlI4)
#include "LuksanVlcek5.hpp"
REGISTER_TNLP(LuksanVlcek5(0,0), LukVlE5)
REGISTER_TNLP(LuksanVlcek5(-1.,0.), LukVlI5)
#include "LuksanVlcek6.hpp"
REGISTER_TNLP(LuksanVlcek6(0,0), LukVlE6)
REGISTER_TNLP(LuksanVlcek6(-1.,0.), LukVlI6)
#include "LuksanVlcek7.hpp"
REGISTER_TNLP(LuksanVlcek7(0,0), LukVlE7)
REGISTER_TNLP(LuksanVlcek7(-1.,0.), LukVlI7)


#include "MittelmannBndryCntrlDiri.hpp"
REGISTER_TNLP(MittelmannBndryCntrlDiri1, MBndryCntrl1)
REGISTER_TNLP(MittelmannBndryCntrlDiri2, MBndryCntrl2)
REGISTER_TNLP(MittelmannBndryCntrlDiri3, MBndryCntrl3)
REGISTER_TNLP(MittelmannBndryCntrlDiri4, MBndryCntrl4)

#include "MittelmannBndryCntrlDiri3D.hpp"
REGISTER_TNLP(MittelmannBndryCntrlDiri3D, MBndryCntrl_3D)

#include "MittelmannBndryCntrlDiri3D_27.hpp"
REGISTER_TNLP(MittelmannBndryCntrlDiri3D_27, MBndryCntrl_3D_27)
REGISTER_TNLP(MittelmannBndryCntrlDiri3D_27BT, MBndryCntrl_3D_27BT)

#include "MittelmannBndryCntrlDiri3Dsin.hpp"
REGISTER_TNLP(MittelmannBndryCntrlDiri3Dsin, MBndryCntrl_3Dsin)

#include "MittelmannBndryCntrlNeum.hpp"
REGISTER_TNLP(MittelmannBndryCntrlNeum1, MBndryCntrl5)
REGISTER_TNLP(MittelmannBndryCntrlNeum2, MBndryCntrl6)
REGISTER_TNLP(MittelmannBndryCntrlNeum3, MBndryCntrl7)
REGISTER_TNLP(MittelmannBndryCntrlNeum4, MBndryCntrl8)

#include "MittelmannDistCntrlDiri.hpp"
REGISTER_TNLP(MittelmannDistCntrlDiri1, MDistCntrl1)
REGISTER_TNLP(MittelmannDistCntrlDiri2, MDistCntrl2)
REGISTER_TNLP(MittelmannDistCntrlDiri3, MDistCntrl3)
REGISTER_TNLP(MittelmannDistCntrlDiri3a, MDistCntrl3a)

#include "MittelmannDistCntrlNeumA.hpp"
REGISTER_TNLP(MittelmannDistCntrlNeumA1, MDistCntrl4)
REGISTER_TNLP(MittelmannDistCntrlNeumA2, MDistCntrl5)
REGISTER_TNLP(MittelmannDistCntrlNeumA3, MDistCntrl6a)

#include "MittelmannDistCntrlNeumB.hpp"
REGISTER_TNLP(MittelmannDistCntrlNeumB1, MDistCntrl4a)
REGISTER_TNLP(MittelmannDistCntrlNeumB2, MDistCntrl5a)
REGISTER_TNLP(MittelmannDistCntrlNeumB3, MDistCntrl6)

#include "MittelmannParaCntrl.hpp"
REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_1>, MPara5_1)
REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_2_1>, MPara5_2_1)
REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_2_2>, MPara5_2_2)
REGISTER_TNLP(MittelmannParaCntrlBase<MittelmannParaCntrl5_2_3>, MPara5_2_3)

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

#ifdef TIME_LIMIT
  if (argv==4) {
    int runtime = atoi(argc[3]);
    pthread_t thread;
    pthread_create(&thread, NULL, killer_thread, &runtime);
  }
  else
#endif
    if (argv!=3 && argv!=1) {
      printf("Usage: %s (this will ask for problem name)\n", argc[0]);
      printf("       %s ProblemName N\n", argc[0]);
      printf("          where N is a positive parameter determining problem size\n");
      printf("       %s list\n", argc[0]);
      printf("          to list all registered problems.\n");
      return -1;
    }

  SmartPtr<RegisteredTNLP> tnlp;
  Index N;

  if (argv!=1) {
    // Create an instance of your nlp...
    tnlp = RegisteredTNLPs::GetTNLP(argc[1]);
    if (!IsValid(tnlp)) {
      printf("Problem with name \"%s\" not known.\n", argc[1]);
      print_problems();
      return -2;
    }

    N = atoi(argc[2]);
  }
  else {
    bool done = false;
    while (!done) {
      string inputword;
      cout << "Enter problem name (or \"list\" for all available names):\n";
      cin >> inputword;
      if (inputword=="list") {
        print_problems();
      }
      else {
        tnlp = RegisteredTNLPs::GetTNLP(inputword.c_str());
        if (!IsValid(tnlp)) {
          printf("Problem with name \"%s\" not known.\n", inputword.c_str());
        }
        else {
          done = true;
        }
      }
    }
    cout << "Enter problem size:\n";
    cin >> N;
  }

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
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }
  // Set option to use internal scaling
  // DOES NOT WORK FOR VLUKL* PROBLEMS:
  // app->Options()->SetStringValueIfUnset("nlp_scaling_method", "user-scaling");

  status = app->OptimizeTNLP(GetRawPtr(tnlp));

  return (int) status;
}
