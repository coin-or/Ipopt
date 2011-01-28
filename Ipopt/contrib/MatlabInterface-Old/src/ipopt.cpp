// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include "matlabexception.h"
#include "array.h"
#include "matlabscalar.h"
#include "matlabstring.h"
#include "matlabfunctionhandle.h"
#include "matlabmatrix.h"
#include "arrayofmatrices.h"
#include "matlabjournal.h"
#include "matlabprogram.h"
#include "matlaboption.h"
#include "multipliers.h"
#include "mex.h"
#include "matrix.h"
#include "IpRegOptions.hpp"
#include "IpJournalist.hpp"
#include "IpIpoptApplication.hpp"

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

    // Check to see if we have the correct number of input arguments.
    if (nrhs < minNumInputArgs)
      throw MatlabException("Incorrect number of input arguments");

    // Get the starting point for the variables. This is specified in
    // the first input argument. The variables must be either a single
    // matrix or a cell array of matrices.
    int k = 0;  // The index of the current input argument.
    int l = 0;  // The index of the current output argument.
    ArrayOfMatrices x0(prhs[k++]);

    // Check to see if we have the correct number of output arguments.
    if ((nlhs < x0.length() + 1) || nlhs > x0.length() + 3)
      throw MatlabException("Incorrect number of output arguments");

    // Create the first output, which will store the termination
    // status of the optimization algorithm.
    MatlabScalar status(plhs[l++],0);

    // Create the second set of outputs, which will store the solution
    // obtained from running IPOPT. There should be at least as many
    // output arguments as cell entries in X.
    ArrayOfMatrices x(&plhs[l],x0);
    l += x0.length();

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

    // If requested, create the output that will store the values of
    // the multipliers obtained when IPOPT converges to a stationary
    // point. This output is a MATLAB structure with three fields, two
    // for the final values of the upper and lower bound multipliers,
    // and one for the constraint multipliers.
    Multipliers* multipliers = 0;
    if (nlhs > l)
      multipliers = new Multipliers(plhs[l++],x.numelems(),
				    constraintlb.length());

    // If requested, create the output that will store the number if
    // iterations needed to attain convergence to a stationary point.
    MatlabScalar* numiter = 0;
    if (nlhs > l)
      numiter = new MatlabScalar(plhs[l++],0);

    // Get the Matlab callback functions.
    MatlabFunctionHandle objFunc(prhs[k++]);
    MatlabFunctionHandle gradFunc(prhs[k++]);
    MatlabFunctionHandle constraintFunc(prhs[k++]);
    MatlabFunctionHandle jacobianFunc(prhs[k++]);
    MatlabFunctionHandle hessianFunc(prhs[k++]);

    // Get the auxiliary data.
    const mxArray* auxData = 0;
    const mxArray* ptr     = prhs[k++];
    if (nrhs > 10) 
      if (!mxIsEmpty(ptr))
	auxData = ptr;

    // Get the iterative callback function.
    MatlabFunctionHandle* iterFunc = new MatlabFunctionHandle();
    ptr = prhs[k++];
    if (nrhs > 11)
      if (!mxIsEmpty(ptr)) {
	delete iterFunc;
	iterFunc = new MatlabFunctionHandle(ptr);
      }

    // Create a new instance of IpoptApplication.
    EJournalLevel printLevel = defaultPrintLevel;
    if (*iterFunc)
      printLevel = Ipopt::J_NONE;
    SmartPtr<Journal> console = new MatlabJournal(printLevel);
    IpoptApplication  app(false);
    app.Jnlst()->AddJournal(console);

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

    // Get the initial Lagrange multipliers, if provided.
    Multipliers* initialMultipliers = 0;
    ptr = prhs[k++];
    app.Options()->SetStringValue("warm_start_init_point","no");
    if (nrhs > 12) 
      if (!mxIsEmpty(ptr)) {
	initialMultipliers = new Multipliers(ptr);

	// Notify the IPOPT algorithm that we will provide our own
	// values for the initial Lagrange multipliers.
	app.Options()->SetStringValue("warm_start_init_point","yes");
      }


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

	  // Special treatment is needed for the limited-memory
	  // quasi-Newton approximation to the Hessian.
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

	    // Special treatment is needed if the print level is set
	    // by the user.
	    if (strcmp(optionLabel,"print_level") == 0) {
	      int printLevel = optionValue;
	      console->SetAllPrintLevels((EJournalLevel) printLevel);
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
    ApplicationReturnStatus exitstatus = app.Initialize();
    if (exitstatus != Ipopt::Solve_Succeeded) {
      throw MatlabException("IPOPT solver initialization failed");
    }

    // Create a new instance of the constrained, nonlinear program.
    MatlabProgram* matlabprogram
      = new MatlabProgram(x0,lb,ub,constraintlb,constraintub,objFunc,
			  gradFunc,constraintFunc,jacobianFunc,hessianFunc,
			  *iterFunc,auxData,x,useQuasiNewton,
			  initialMultipliers,multipliers);
    SmartPtr<TNLP> program = matlabprogram;

    // Ask Ipopt to solve the problem.
    exitstatus = app.OptimizeTNLP(program);
    status     = (double) exitstatus;

    // Return the number of IPOPT iterations, if requested.
    if (numiter)
      *numiter = matlabprogram->getnumiterations();

    // Get rid of the dynamically allocated memory.
    if (multipliers)        delete multipliers;
    if (numiter)            delete numiter;
    if (iterFunc)           delete iterFunc;
    if (initialMultipliers) delete initialMultipliers;

  } catch (std::exception& error) {
    mexErrMsgTxt(error.what());
  }
  
void convertMatlabConstraints (Array<double>& constraints, double infty) {
  int n = constraints.length();
  for (int i = 0; i < n; i++)
    if (mxIsInf(constraints[i]))
      constraints[i] = infty;
}

