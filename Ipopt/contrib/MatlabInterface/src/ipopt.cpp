// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabexception.h"
#include "array.h"
#include "matlabscalar.h"
#include "matlabstring.h"
#include "matlabmatrix.h"
#include "arrayofmatrices.h"
#include "matlabjournal.h"
#include "matlabprogram.h"
#include "matlaboption.h"
#include "mex.h"
#include "matrix.h"
#include "ipopt/IpRegOptions.hpp"
#include "ipopt/IpJournalist.hpp"
#include "ipopt/IpIpoptApplication.hpp"

using Ipopt::IsValid;
using Ipopt::RegisteredOption;
using Ipopt::EJournalLevel;
using Ipopt::Journal;
using Ipopt::MatlabJournal;
using Ipopt::IpoptApplication;
using Ipopt::SmartPtr;
using Ipopt::TNLP;
using Ipopt::ApplicationReturnStatus;

extern void _main();

// Constants.
// -----------------------------------------------------------------
const int           minNumInputArgs   = 10;
const EJournalLevel defaultPrintLevel = Ipopt::J_ITERSUMMARY;

// Function declarations.
// -----------------------------------------------------------------
void convertMatlabConstraints (Array<double>& constraints, double infty);

// Function definitions.
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[]) 
  try {

    // Check to see if we have the correct number of input and output
    // arguments.
    if (nrhs < minNumInputArgs)
      throw MatlabException("Incorrect number of input arguments");

    // Get the starting point for the variables. This is specified in
    // the first input argument. The variables must be either a single
    // matrix or a cell array of matrices.
    int k = 0;  // The index of the current input argument.
    ArrayOfMatrices x0(prhs[k++]);

    // Create the output, which stores the solution obtained from
    // running IPOPT. There should be as many output arguments as cell
    // entries in X.
    if (nlhs != x0.length())
      throw MatlabException("Incorrect number of output arguments");
    ArrayOfMatrices x(plhs,x0);

    // Load the lower and upper bounds on the variables as
    // ArrayOfMatrices objects. They should have the same structure as
    // the ArrayOfMatrices object "x".
    ArrayOfMatrices lb(prhs[k++]);
    ArrayOfMatrices ub(prhs[k++]);

    // Check to make sure the bounds make sense.
    if (lb != x || ub != x)
      throw MatlabException("Input arguments LB and UB must have the same \
structure as X");

    // Load the lower and upper bounds on the constraints. Each of
    // these is a Matlab array. They should have the same length.
    Matrix constraintlb(prhs[k++]);
    Matrix constraintub(prhs[k++]);
    if (constraintlb.length() != constraintub.length())
      throw MatlabException("Input arguments CONSTRAINTLB and CONSTRAINTUB \
should have the same number of elements");

    // Get the Matlab callback functions.
    MatlabString objFunc(prhs[k++]);
    MatlabString gradFunc(prhs[k++]);
    MatlabString constraintFunc(prhs[k++]);
    MatlabString jacobianFunc(prhs[k++]);
    MatlabString hessianFunc(prhs[k++]);

    // Get the auxiliary data.
    const mxArray* auxData;
    const mxArray* ptr = prhs[k++];
    if (nrhs > 10) {
      if (mxIsEmpty(ptr))
	auxData = 0;
      else
	auxData = ptr;
    }
    else
      auxData = 0;

    // Get the iterative callback function.
    MatlabString* iterFunc;
    ptr = prhs[k++];
    if (nrhs > 11) {
      if (mxIsEmpty(ptr))
	iterFunc = new MatlabString("");
      else
	iterFunc = new MatlabString(ptr);
    }
    else
      iterFunc = new MatlabString("");

    // Create a new instance of IpoptApplication.
    SmartPtr<Journal> console = new MatlabJournal(defaultPrintLevel);
    IpoptApplication  app(false);
    app.Jnlst()->AddJournal(console);
    if (!iterFunc->isempty())
      app.Options()->SetIntegerValue("print_level",Ipopt::J_NONE);
    
    // If an upper/lower bound is set to infinity, then it means that
    // the variable is not upper/lower bounded.
    double lower_infty;
    double upper_infty;
    app.Options()->GetNumericValue("nlp_lower_bound_inf",lower_infty,"");
    app.Options()->GetNumericValue("nlp_upper_bound_inf",upper_infty,"");

    for (int i = 0; i < lb.length(); i++)
      convertMatlabConstraints(*lb[i],lower_infty);
    for (int i = 0; i < ub.length(); i++)
      convertMatlabConstraints(*ub[i],upper_infty);
    convertMatlabConstraints(constraintlb,lower_infty);
    convertMatlabConstraints(constraintub,upper_infty);

    // Process the remaining input arguments, which set options for
    // the IPOPT algorithm.
    bool useQuasiNewton = false;
    while (k < nrhs) {

      // Get the option label from the Matlab input argument.
      MatlabString optionLabel(prhs[k++]);

      if (k < nrhs) {

	// Get the option value from the Matlab input argument.
	MatlabOption optionValue(prhs[k++]);

	// Next, check to make sure we have a valid option.
	SmartPtr<const RegisteredOption> option
	  = app.RegOptions()->GetOption(optionLabel);
	if (!IsValid(option)) {
	  delete iterFunc;
	  throw MatlabException("Nonexistent IPOPT option");
	}

	if (optionValue.isString()) {

	  // Set the string option.
	  if (!app.Options()->SetStringValue(optionLabel,optionValue)) {
	    delete iterFunc;
	    throw MatlabException("Invalid IPOPT option value");
	  }
	  if ((strcmp(optionLabel,"hessian_approximation") == 0) && 
	      (strcmp(optionValue,"limited-memory")) == 0)
	    useQuasiNewton = true;
	} else {

	  // Set the numeric option.
	  if (option->Type() == Ipopt::OT_Integer) {
	    if (!app.Options()->SetIntegerValue(optionLabel,optionValue)) {
	      delete iterFunc;
	      throw MatlabException("Invalid IPOPT option value");
	    }
	  }
	  else if (!app.Options()->SetNumericValue(optionLabel,optionValue)) {
	    delete iterFunc;
	    throw MatlabException("Invalid IPOPT option value");
	  }
	}
      }
    }

    // Intialize the IpoptApplication object and process the options.
    app.Initialize();

    // Create a new instance of the constrained, nonlinear program.
    SmartPtr<TNLP> program 
      = new MatlabProgram(x0,lb,ub,constraintlb,constraintub,objFunc,
			  gradFunc,constraintFunc,jacobianFunc,hessianFunc,
			  *iterFunc,auxData,x,useQuasiNewton);

    // Ask Ipopt to solve the problem.
    ApplicationReturnStatus exitstatus = app.OptimizeTNLP(program);

    // Get rid of the dynamically allocated memory.
    delete iterFunc;

    // Throw an exception if the solver terminates before finding a
    // local solution.
    switch (exitstatus) {
    case Ipopt::Solve_Succeeded:
    case Ipopt::Solved_To_Acceptable_Level:
    case Ipopt::User_Requested_Stop:
    case Ipopt::Maximum_Iterations_Exceeded:
      break;
    default:
      throw MatlabException("IPOPT solver terminated unexpectedly");
    }
  } catch (std::exception& error) {
    mexErrMsgTxt(error.what());
  }
  
void convertMatlabConstraints (Array<double>& constraints, double infty) {
  int n = constraints.length();
  for (int i = 0; i < n; i++)
    if (mxIsInf(constraints[i]))
      constraints[i] = infty;
}

