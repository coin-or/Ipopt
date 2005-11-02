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

static void print_problems()
{
  printf("\nList of all registered problems:\n\n");
  RegisteredTNLPList::PrintRegisteredProblems();
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
    RegisteredTNLPList::GetTNLP(argc[1]);
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
