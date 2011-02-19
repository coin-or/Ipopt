/*
 * Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * File:   ipoptr.cpp
 * Author: Jelmer Ypma
 * Date:   18 April 2010
 *
 * This file defines the main function IpoptRSolve that provides
 * an interface to Ipopt from R.
 * The function converts and R object containing objective function,
 * constraints, etc. into an IpoptApplication, solves the problem,
 * and returns the result.
 *
 * Financial support of the UK Economic and Social Research Council 
 * through a grant (RES-589-28-0001) to the ESRC Centre for Microdata 
 * Methods and Practice (CeMMAP) is gratefully acknowledged.
 *
 * 30/01/2011: added IpoptRJournal to correctly direct output to R terminal.
 */

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "IpoptRNLP.hpp"
#include "IpoptRJournal.hpp"

#include <R.h>
#include <Rdefines.h>
// Rdefines.h is somewhat more higher level then Rinternal.h, and is preferred if the code might be shared with S at any stage.

#include <string>


//
// Extracts element with name 'str' from R object 'list'
// and returns that element.
//
SEXP
getListElement(SEXP list, std::string str)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < length(list); i++)
	if(str.compare(CHAR(STRING_ELT(names, i))) == 0) {
		elmt = VECTOR_ELT(list, i);
		break;
	}
	return elmt;
}


//
// Set options specified in R object opts in IpoptApplcation app.
// If we set the option to approximate the Hessian, then IpoptRNLP
// needs to know that, so it can return false when calling eval_h,
// instead of trying to find an R function that evaluates the Hessian.
//
// SEXP opts is an R list containing three sub-lists, with names
// integer, string and numeric. These sub-lists contain the actual
// options and there values that were specified by the user.
// Passing the options in this way makes it easier to call different
// SetValue function in IpoptApplication of the different types.
//
void setApplicationOptions( Ipopt::SmartPtr<Ipopt::IpoptApplication> app, 
                            Ipopt::SmartPtr<IpoptRNLP> nlp,
							SEXP opts ) {
	
	// extract the sub-lists with options of the different types into separate lists
	SEXP opts_integer = getListElement(opts, "integer");
	SEXP opts_numeric = getListElement(opts, "numeric");
	SEXP opts_string = getListElement(opts, "string");
	
	// loop over the integer options and set them
	SEXP opts_integer_names;
	opts_integer_names = getAttrib(opts_integer, R_NamesSymbol);
	for (int list_cnt=0;list_cnt<length( opts_integer );list_cnt++) {
		
		SEXP opt_value;
		PROTECT(opt_value = AS_INTEGER(VECTOR_ELT(opts_integer, list_cnt)));
		
		app->Options()->SetIntegerValue(CHAR(STRING_ELT(opts_integer_names, list_cnt)), INTEGER(opt_value)[0]);
		UNPROTECT(1);	
	}
	
	// loop over the numeric options and set them
	SEXP opts_numeric_names;
	opts_numeric_names = getAttrib(opts_numeric, R_NamesSymbol);
	for (int list_cnt=0;list_cnt<length( opts_numeric );list_cnt++) {
		
		SEXP opt_value;
		PROTECT(opt_value = VECTOR_ELT(opts_numeric, list_cnt));
		
		app->Options()->SetNumericValue(CHAR(STRING_ELT(opts_numeric_names, list_cnt)), REAL(opt_value)[0]);
		UNPROTECT(1);	
	}
	
	// loop over the string options and set them
	SEXP opts_string_names;
	opts_string_names = getAttrib(opts_string, R_NamesSymbol);
	for (int list_cnt=0;list_cnt<length( opts_string );list_cnt++) {
		
		// opt_value will contain the first (should be the only one) element of the list
		SEXP opt_value;
		PROTECT(opt_value = STRING_ELT(VECTOR_ELT(opts_string, list_cnt),0));
		
		app->Options()->SetStringValue(CHAR(STRING_ELT(opts_string_names, list_cnt)), CHAR(opt_value));
		
		// change the setting to approximate the hessian in nlp if this is part of the options
		if ( std::string( CHAR(STRING_ELT(opts_string_names, list_cnt)) ) == "hessian_approximation" &&
		     std::string( CHAR(opt_value) ) == "limited-memory" ) 
		{
			nlp->set_hessian_approximation( true );	
		}
  		     
		UNPROTECT(1);	
	}

}


// we want this function to be available in R, so we put extern around it.
extern "C" {

SEXP IpoptRSolve( SEXP args )
{
	// Create an instance of your nlp...
	//  (use a SmartPtr, not raw)
  	Ipopt::SmartPtr<IpoptRNLP> ipoptr_nlp;
  	Ipopt::SmartPtr<Ipopt::TNLP> ipoptr_tnlp;
  
  	ipoptr_nlp = new IpoptRNLP();
  	ipoptr_tnlp = GetRawPtr(ipoptr_nlp);
  
	// Set initial values, bounds, functions, etc.
	// These are taking from args, which is a list containing the needed elements.
	// Checking whether all elements are there and have correct values is done in R.
	ipoptr_nlp->set_R_environment( getListElement(args, "environment") );
	ipoptr_nlp->set_R_init_values( getListElement(args, "x0" ) );
	ipoptr_nlp->set_R_lower_bounds( getListElement(args, "lower_bounds" ) );
	ipoptr_nlp->set_R_upper_bounds( getListElement(args, "upper_bounds" ) );
	ipoptr_nlp->set_R_eval_f( getListElement(args, "eval_f") );
	ipoptr_nlp->set_R_eval_grad_f( getListElement(args, "eval_grad_f") );
  
	ipoptr_nlp->set_R_eval_g( getListElement(args, "eval_g") );
	ipoptr_nlp->set_R_eval_jac_g( getListElement(args, "eval_jac_g") );
	ipoptr_nlp->set_R_eval_jac_g_structure( getListElement(args, "eval_jac_g_structure") );
	ipoptr_nlp->set_R_constraint_lower_bounds( getListElement(args, "constraint_lower_bounds") );
	ipoptr_nlp->set_R_constraint_upper_bounds( getListElement(args, "constraint_upper_bounds") );
  
	ipoptr_nlp->set_R_eval_h( getListElement(args, "eval_h") );
	ipoptr_nlp->set_R_eval_h_structure( getListElement(args, "eval_h_structure") );
  
	// Create an instance of the IpoptApplication
	Ipopt::SmartPtr<Ipopt::IpoptApplication> app = new Ipopt::IpoptApplication();

	// Set options that were passed from R
	setApplicationOptions( app, ipoptr_nlp, getListElement(args, "options") );
  
    // Set up the IPOPT console.
    //
    // Get print_level from options
    Ipopt::Index print_level;
    app->Options()->GetIntegerValue( "print_level", print_level, "" );
    
    // Set print_level to 0 for default console (to avoid double output under Linux)
    app->Options()->SetIntegerValue( "print_level", 0 );

    // Add new journal with user-supplied print_level to print output to R console
    Ipopt::SmartPtr<Ipopt::Journal> console = new IpoptRJournal( static_cast<Ipopt::EJournalLevel>( print_level ) );
    app->Jnlst()->AddJournal(console);
    
	// Initialize the IpoptApplication and process the options
	Ipopt::ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Ipopt::Solve_Succeeded) {
		printf("\n\n*** Error during initialization!\n");
    
		SEXP answer;
		PROTECT(answer = allocVector(INTSXP,1));
		INTEGER(answer)[0] = (int) status;
		UNPROTECT(1);
		return(answer);
	}

	// Solve the optimization problem
	status = app->OptimizeTNLP(ipoptr_tnlp);

	// Retrieve results that were saved in IpoptRNLP::finalize_solution()
	SEXP R_result_list;
    PROTECT(R_result_list = ipoptr_nlp->get_R_result_list());
	
	// convert number of iterations to an R object and add to the result_list
	// memory in R_result_list has already been allocated in IpoptRNLP::finalize_solution()
	SEXP R_num_iterations;
	PROTECT(R_num_iterations = allocVector(INTSXP,1));
	INTEGER(R_num_iterations)[0] = (int) app->Statistics()->IterationCount();
	
	SET_VECTOR_ELT(R_result_list, 2, R_num_iterations);
	UNPROTECT(2);
	
	return(R_result_list);
}

} // extern "C"

