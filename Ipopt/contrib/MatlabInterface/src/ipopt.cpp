// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 18, 2008

#include "mex.h"
#include "iterate.hpp"
#include "options.hpp"
#include "matlabexception.hpp"
#include "callbackfunctions.hpp"
#include "matlabinfo.hpp"
#include "matlabjournal.hpp"
#include "matlabprogram.hpp"
#include "IpRegOptions.hpp"
#include "IpJournalist.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"

using Ipopt::IsValid;
using Ipopt::RegisteredOption;
using Ipopt::EJournalLevel;
using Ipopt::Journal;
using Ipopt::MatlabJournal;
using Ipopt::IpoptApplication;
using Ipopt::SmartPtr;
using Ipopt::TNLP;
using Ipopt::ApplicationReturnStatus;
using Ipopt::SolveStatistics;

extern void _main();

// Function definitions.
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[]) 
  try {

    // Check to see if we have the correct number of input and output
    // arguments.
    if (nrhs != 3)
      throw MatlabException("Incorrect number of input arguments");
    if (nlhs != 2)
      throw MatlabException("Incorrect number of output arguments");

    // Get the first input which specifies the initial iterate.
    Iterate x0(mxDuplicateArray(prhs[0]));

    // Get the second input which specifies the callback functions.
    CallbackFunctions funcs(prhs[1]);

    // Create a new IPOPT application object and process the options.
    IpoptApplication app(false);
    Options          options(x0,app,prhs[2]);

    // The first output argument is the value of the optimization
    // variables obtained at the solution.
    plhs[0] = mxDuplicateArray(x0);
    Iterate x(plhs[0]);

    // The second output argument stores other information, such as
    // the exit status, the value of the Lagrange multipliers upon
    // termination, the final state of the auxiliary data, and so on.
    MatlabInfo info(plhs[1]);

    // Check to see whether the user provided a callback function for
    // computing the Hessian. This is not needed in the special case
    // when a quasi-Newton approximation to the Hessian is being used.
    if (!options.ipoptOptions().useQuasiNewton() && 
	!funcs.hessianFuncIsAvailable())
      throw MatlabException("You must supply a callback function for \
computing the Hessian unless you decide to use a quasi-Newton \
approximation to the Hessian");

    // If the user tried to use her own scaling, report an error.
    if (options.ipoptOptions().userScaling())
      throw MatlabException("The user-defined scaling option does not \
work in the MATLAB interface for IPOPT");

    // If the user supplied initial values for the Lagrange
    // multipliers, activate the "warm start" option in IPOPT.
    if (options.multlb() && options.multub() &&
	(numconstraints(options)==0 || options.multconstr()) )
      app.Options()->SetStringValue("warm_start_init_point","yes");

    // Set up the IPOPT console.
    EJournalLevel printLevel = (EJournalLevel) 
      options.ipoptOptions().printLevel();
    if(printLevel > 0) { //prevents IPOPT display if we don't want it
        SmartPtr<Journal> console = new MatlabJournal(printLevel);
        app.Jnlst()->AddJournal(console);
    }

    // Intialize the IpoptApplication object and process the options.
    ApplicationReturnStatus exitstatus;
    exitstatus = app.Initialize();
    if (exitstatus != Ipopt::Solve_Succeeded)
      throw MatlabException("IPOPT solver initialization failed");

    // Create a new instance of the constrained, nonlinear program.
    MatlabProgram* matlabProgram 
      = new MatlabProgram(x0,funcs,options,x,info);
    SmartPtr<TNLP> program = matlabProgram;

    // Ask Ipopt to solve the problem.
    exitstatus = app.OptimizeTNLP(program);
    info.setExitStatus(exitstatus);

    // Collect statistics about Ipopt run
    if (IsValid(app.Statistics())) {
      SmartPtr<SolveStatistics> stats = app.Statistics();
      info.setIterationCount(stats->IterationCount());
      //Get Function Calls
      int obj, con, grad, jac, hess;
      stats->NumberOfEvaluations(obj,con,grad,jac,hess);
      info.setFuncEvals(obj, con, grad, jac, hess);      
      //CPU Time
      info.setCpuTime(stats->TotalCpuTime());
    }

    // Free the dynamically allocated memory.
    mxDestroyArray(x0);

  } catch (std::exception& error) {
    mexErrMsgTxt(error.what());
  }
