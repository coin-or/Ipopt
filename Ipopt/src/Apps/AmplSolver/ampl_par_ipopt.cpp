// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptApplication.hpp"
#include "IpParTNLPAdapter.hpp"
#include "AmplTNLP.hpp"
#include "AmplParTNLP.hpp"

#include "IpoptConfig.h"
#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif

// for printf
#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

// FIXME - proper header files
//extern "C" {
#define MPICH_SKIP_MPICXX
#include "mpi.h"
//}

int main(int argv, char**argc)
{
  using namespace Ipopt;

  int DebugWait = 0; // set this to 1 to attach debugger
  if (argv > 2) DebugWait = 1;
  while (DebugWait);

  MPI_Init(&argv, &argc);

  SmartPtr<IpoptApplication> app = new IpoptApplication();

  // Check if executable is run only to print out options documentation
  if (argv == 2) {
    bool print_options = false;
    bool print_latex_options = false;
    if (!strcmp(argc[1],"--print-options")) {
      print_options = true;
    }
    else if (!strcmp(argc[1],"--print-latex-options")) {
      print_options = true;
      print_latex_options = true;
    }
    if (print_options) {
      SmartPtr<OptionsList> options = app->Options();
      options->SetStringValue("print_options_documentation", "yes");
      if (print_latex_options) {
        options->SetStringValue("print_options_latex_mode", "yes");
      }
      app->Initialize("");
      return 0;
    }
  }

  // Call Initialize the first time to create a journalist, but ignore
  // any options file
  ApplicationReturnStatus retval;
  retval = app->Initialize("");
  if (retval != Solve_Succeeded) {
    printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
    exit(-100);
  }

  // Add the suffix handler for scaling
  SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);
  // Modified for warm-start from AMPL
  suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

  SmartPtr<AmplTNLP> ampl_tnlp = new AmplTNLP(ConstPtr(app->Jnlst()),
                                 app->Options(),
                                 argc, suffix_handler);

  // Call Initialize again to process output related options
  retval = app->Initialize();
  if (retval != Solve_Succeeded) {
    printf("ampl_ipopt.cpp: Error in second Initialize!!!!\n");
    exit(-101);
  }

  SmartPtr<ParTNLP> parampl = new AmplParTNLP(ampl_tnlp);
  SmartPtr<NLP> nlp = new ParTNLPAdapter(parampl, ConstPtr(app->Jnlst()));

  const int n_loops = 1; // make larger for profiling
  for (Index i=0; i<n_loops; i++) {
    bool status;
    status = app->OptimizeNLP(nlp);
  }

  // finalize_solution method in AmplTNLP writes the solution file

  // clean up memory
  app = NULL;
  nlp = NULL;
  ampl_tnlp = NULL;
  parampl = NULL;

  MPI_Finalize();

  return 0;
}


