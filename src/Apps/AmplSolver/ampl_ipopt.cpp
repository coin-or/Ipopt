// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "AmplTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "IpoptConfig.h"

#include <cstring>
#include <cstdio>

int main(
   int argc,
   char** args)
{
   using namespace Ipopt;

   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
   app->RethrowNonIpoptException(false);

   // Check if executable is run only to print out options documentation
   if( argc == 2 )
   {
      bool print_options = false;
      std::string print_options_mode("text");
      if( !strcmp(args[1], "--print-options=latex") )
      {
         print_options = true;
         print_options_mode = "latex";
      }
      else if( !strcmp(args[1], "--print-options=doxygen") )
      {
         print_options = true;
         print_options_mode = "doxygen";
      }
      else if( !strcmp(args[1], "--print-options") )
      {
         print_options = true;
      }
      else if( !strcmp(args[1], "--print-latex-options") )
      {
         fprintf(stderr, "ampl_ipopt.cpp: Options --print-latex-options has been replaced by --print-options=latex. Please adjust your call.\n");
         exit(-200);
      }
      if( print_options )
      {
         SmartPtr<OptionsList> options = app->Options();
         options->SetStringValue("print_options_documentation", "yes");
         options->SetStringValue("print_advanced_options", "yes");
         options->SetStringValue("print_options_mode", print_options_mode);
         app->Initialize("");
         return 0;
      }
   }

   // Call Initialize the first time to create a journalist, but ignore
   // any options file
   ApplicationReturnStatus retval;
   retval = app->Initialize("");
   if( retval != Solve_Succeeded )
   {
      printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
      exit(-100);
   }

   // Add the suffix handler for scaling
   SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
   suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source,
                                      AmplSuffixHandler::Number_Type);
   // Modified for warm-start from AMPL
   suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);

   SmartPtr<TNLP> ampl_tnlp = new AmplTNLP(ConstPtr(app->Jnlst()), app->RegOptions(), app->Options(), args, suffix_handler);

   // Call Initialize again to process output related options
   retval = app->Initialize();
   if( retval != Solve_Succeeded )
   {
      printf("ampl_ipopt.cpp: Error in second Initialize!!!!\n");
      exit(-101);
   }

   const int n_loops = 1; // make larger for profiling
   for( Index i = 0; i < n_loops; i++ )
   {
      app->OptimizeTNLP(ampl_tnlp);
   }

   // finalize_solution method in AmplTNLP writes the solution file

   return 0;
}

