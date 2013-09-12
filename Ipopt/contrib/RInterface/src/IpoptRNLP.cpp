/*
 * Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * file:   IpoptRNLP.cpp
 * author: Jelmer Ypma
 * date:   18 April 2010
 *
 * This file defines a C++ class that derives from Ipopt::TNLP. The class
 * takes care of interaction between Ipopt and user-defined functions in R.
 *
 * Financial support of the UK Economic and Social Research Council 
 * through a grant (RES-589-28-0001) to the ESRC Centre for Microdata 
 * Methods and Practice (CeMMAP) is gratefully acknowledged. 
 *
 * Changelog:
 *   09/03/2012: added outputs in finalize_solution; z_L, z_U, constraints, lambda (thanks to Michael Schedl)
 */

#include "IpoptRNLP.hpp"

/* Constructor. */
IpoptRNLP::IpoptRNLP() : 
    d_hessian_approximation( false ),
    d_num_protected_members( 0 )
{}

IpoptRNLP::~IpoptRNLP()
{
    // UNPROTECT all SEXP members that we PROTECT
    UNPROTECT( d_num_protected_members );
}

// 
// Functions to load R Objects into IpoptRProblem
//
void IpoptRNLP::set_R_environment( SEXP env ) 
{
	PROTECT(R_environment = env);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_eval_f( SEXP f ) 
{
	PROTECT(R_eval_f = f);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_eval_grad_f( SEXP f ) 
{
	PROTECT(R_eval_grad_f = f);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_init_values( SEXP x0 )
{
	PROTECT(R_init_values = x0);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_lower_bounds( SEXP lb ) 
{
	PROTECT(R_lower_bounds = lb);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_upper_bounds( SEXP ub ) 
{
	PROTECT(R_upper_bounds = ub);
    d_num_protected_members++;
}
 
void IpoptRNLP::set_R_eval_g( SEXP g )
{
	PROTECT(R_eval_g = g);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_eval_jac_g( SEXP g )
{
	PROTECT(R_eval_jac_g = g);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_eval_jac_g_structure( SEXP s )
{
	PROTECT(R_eval_jac_g_structure = s);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_constraint_lower_bounds( SEXP lb )
{
	PROTECT(R_constraint_lower_bounds = lb);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_constraint_upper_bounds( SEXP ub )
{
	PROTECT(R_constraint_upper_bounds = ub);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_eval_h( SEXP h )
{
	PROTECT(R_eval_h = h);
    d_num_protected_members++;
}

void IpoptRNLP::set_R_eval_h_structure( SEXP s )
{
	PROTECT(R_eval_h_structure = s);
    d_num_protected_members++;
}

void IpoptRNLP::set_hessian_approximation( bool b )
{
    d_hessian_approximation = b;
}

SEXP IpoptRNLP::get_R_result_list() {
	return R_result_list;
}

bool IpoptRNLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                         Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	// Check for user interruption from R
	R_CheckUserInterrupt();
	
	// number of control variables
	n = length( R_init_values );

	// number of constraints
	m = length( R_constraint_lower_bounds );
	
	// Loop over the elements in R_eval_jac_g_structure and count the number of non-zero indices
	// in the Jacobian. As far as I know unlist() does not exist in C, so we cannot call that directly.
	nnz_jac_g = 0;
	for (int list_cnt=0;list_cnt<length( R_eval_jac_g_structure );list_cnt++) {
		
		SEXP R_list_element;
		PROTECT(R_list_element = AS_INTEGER(VECTOR_ELT(R_eval_jac_g_structure, list_cnt)));

		nnz_jac_g += length( R_list_element );
		UNPROTECT(1);	
	}
	
	// Loop over the elements in R_eval_h_structure and count the number of non-zero indices
	// in the hessian of the lagrangian (combined hessian of the objective and hessian of the constraints).
	
	nnz_h_lag = 0;
	for (int list_cnt=0;list_cnt<length( R_eval_h_structure );list_cnt++) {
		
		SEXP R_list_element;
		PROTECT(R_list_element = AS_INTEGER(VECTOR_ELT(R_eval_h_structure, list_cnt)));
		
		nnz_h_lag+=length( R_list_element );
		UNPROTECT(1);	
	}	
	
	// We use the standard Fortran Ipopt::Index style for row/col entries,
	// This is the same as R, start counting indices in the structure matrices at 1
	index_style = FORTRAN_STYLE;

	return true;
}

bool IpoptRNLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                            Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
	// Check that the number of controls, n, and the number of constraints, m
    // are of the same length as the R variables that were passed.
	assert(n == length( R_init_values ));
	assert(n == length( R_lower_bounds ));
	assert(n == length( R_upper_bounds ));
	assert(m == length( R_constraint_lower_bounds ));
	assert(m == length( R_constraint_upper_bounds ));
	
	// Check for user interruption from R
	R_CheckUserInterrupt();
	
    // set the upper and lower bounds of the control
	for (Ipopt::Index i=0;i<n;i++) {
		x_l[i] = REAL(R_lower_bounds)[i];				// lower bound
		x_u[i] = REAL(R_upper_bounds)[i];				// upper bound
	}


	// set the upper and lower bounds of the inequality constraints
	for (Ipopt::Index i=0;i<m;i++) {
		g_l[i] = REAL(R_constraint_lower_bounds)[i];	// lower bound
		g_u[i] = REAL(R_constraint_upper_bounds)[i];	// upper bound
	}
	
	return true;
}

bool IpoptRNLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                               bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                               Ipopt::Index m, bool init_lambda,
                               Ipopt::Number* lambda)
{
	// We have starting values for the control, x, only.
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	// Check for user interruption from R
	R_CheckUserInterrupt();
	
	// set initial values of the controls
	for (Ipopt::Index i=0;i<n;i++) {
		x[i] = REAL(R_init_values)[i];
	}
	
	return true;
}

bool IpoptRNLP::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value)
{
	// Calculate and return the value of the objective function
  
	// Check for user interruption from R
	R_CheckUserInterrupt();
  
	SEXP rargs,Rcall,result;
  
	// Allocate memory for a vector of reals.
	// This vector will contain the elements of x,
	// x is the argument to the R function R_eval_f
	PROTECT(rargs = allocVector(REALSXP,n));
	for (Ipopt::Index i=0;i<n;i++) {
		REAL(rargs)[i] = x[i];
	}
  
	// evaluate R function R_eval_f with the control x as an argument
	PROTECT(Rcall = lang2(R_eval_f,rargs));
	PROTECT(result = eval(Rcall,R_environment));
  
	// recode the return value from SEXP to Number
	obj_value = REAL(result)[0];
 
	UNPROTECT(3);
  
	return true;
}

bool IpoptRNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f)
{
	// Calculate and return the gradient of the objective function grad_{x} f(x)

	// if we have two controls, x1 and x2:
	// grad_f[0] = grad_{x1} f(x)
	// grad_f[1] = grad_{x2} f(x)
	
	// Check for user interruption from R
	R_CheckUserInterrupt();
	
	SEXP rargs,Rcall,result;
  
	// allocate memory for a vector of reals
	// this vector will contain the elements of x
	// x is the argument to the R function R_eval_grad_f
	PROTECT(rargs = allocVector(REALSXP,n));
	for (Ipopt::Index i=0;i<n;i++) {
		REAL(rargs)[i] = x[i];
	}
  
	// evaluate R function R_eval_grad_f with the control x as an argument
	PROTECT(Rcall = lang2(R_eval_grad_f,rargs));
	PROTECT(result = eval(Rcall,R_environment));
  
  	// recode the return values from SEXP to Numbers
	for (Ipopt::Index i=0;i<n;i++) {
		grad_f[i] = REAL(result)[i];
	}
 
	UNPROTECT(3);
  
	return true;
}

bool IpoptRNLP::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
	// Calculate and return the value of the constraints: g(x)
  
	// Check for user interruption from R
	R_CheckUserInterrupt();
  
	SEXP rargs,Rcall,result;
  
	// Allocate memory for a vector of reals
	// this vector will contain the elements of x
	// x is the argument to the R function R_eval_g
	PROTECT(rargs = allocVector(REALSXP,n));
	for (Ipopt::Index i=0;i<n;i++) {
		REAL(rargs)[i] = x[i];
	}
  
	PROTECT(Rcall = lang2(R_eval_g,rargs));
	PROTECT(result = eval(Rcall,R_environment));
  
	for (Ipopt::Index i=0;i<m;i++) {
		g[i] = REAL(result)[i];
	}
 
	UNPROTECT(3);
  
	return true;
}

bool IpoptRNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                       Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
                       Ipopt::Number* values)
{
	// These use Fortran indexing style and start counting at 1

	// Check for user interruption from R
	R_CheckUserInterrupt();
	
	if (values == NULL) {
		// return the structure of the jacobian of the constraints

		// element at 1,1: grad_{x1} g_{1}(x)
		//iRow[0] = 1;
		//jCol[0] = 1;

		// element at 1,2: grad_{x2} g_{1}(x)
		//iRow[1] = 1;
		//jCol[1] = 2;
		
		Ipopt::Index total_cnt = 0;
		for (int list_cnt=0;list_cnt<length( R_eval_jac_g_structure );list_cnt++) {
			
			SEXP R_list_element;
			PROTECT(R_list_element = AS_INTEGER(VECTOR_ELT(R_eval_jac_g_structure, list_cnt)));
			for (int vector_cnt=0;vector_cnt< length(R_list_element); vector_cnt++) {
				iRow[ total_cnt ] = list_cnt+1;		// we have to add 1 to turn it into Fortran styl indexing
				jCol[ total_cnt ] = INTEGER(R_list_element)[vector_cnt];
				total_cnt++;
			}
			UNPROTECT(1);	
		}
	
	}
	else {
		// return the values of the jacobian of the constraints
		
		SEXP rargs,Rcall,result;
	  
		// allocate memory for a vector of reals
		// this vector will contain the elements of x
		// x is the argument to the R function R_eval_g_jac
		PROTECT(rargs = allocVector(REALSXP,n));
		for (Ipopt::Index i=0;i<n;i++) {
			REAL(rargs)[i] = x[i];
		}
	  
		PROTECT(Rcall = lang2(R_eval_jac_g,rargs));
		PROTECT(result = eval(Rcall,R_environment));
	  
		for (Ipopt::Index i=0;i<nele_jac;i++) {
			values[i] = REAL(result)[i];
		}
	 
		UNPROTECT(3);
		
	}

	return true;
}

bool IpoptRNLP::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                   Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                   bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                   Ipopt::Index* jCol, Ipopt::Number* values)
{

	// Check for user interruption from R
	R_CheckUserInterrupt();

	if ( d_hessian_approximation ) {
		return false;
	}
	else {
	
		if (values == NULL) {
			// return the structure. This is a symmetric matrix, fill the lower left
			// triangle only.
			// Note: off-diagonal elements are zero for this problem
			// element at 1,1: grad^2_{x1,x1} L(x,lambda)
			// iRow[0] = 1;
			// jCol[0] = 1;

			// element at 2,2: grad^2_{x2,x2} L(x,lambda)
			// iRow[1] = 2;
			// jCol[1] = 2;

			Ipopt::Index total_cnt = 0;
			for (int list_cnt=0;list_cnt<length( R_eval_h_structure );list_cnt++) {
				
				SEXP R_list_element;
				PROTECT(R_list_element = AS_INTEGER(VECTOR_ELT(R_eval_h_structure, list_cnt)));
				for (int vector_cnt=0;vector_cnt< length(R_list_element); vector_cnt++) {
					iRow[ total_cnt ] = list_cnt+1;		// we have to add 1 to turn it into Fortran styl indexing
					jCol[ total_cnt ] = INTEGER(R_list_element)[vector_cnt];
					total_cnt++;
				}
				UNPROTECT(1);	
			}	
			
		}
		else {
			// return the values

			// element at 1,1: grad^2_{x1,x1} L(x,lambda)
			// values[0] = -2.0 * lambda[0];

			// element at 2,2: grad^2_{x2,x2} L(x,lambda)
			// values[1] = -2.0 * obj_factor;

			SEXP rargs_x;
            PROTECT(rargs_x = allocVector(REALSXP,n));
			for (Ipopt::Index i=0;i<n;i++) {
				REAL(rargs_x)[i] = x[i];
			}
			
            SEXP rargs_obj_factor;
			PROTECT(rargs_obj_factor = allocVector(REALSXP,1));
			REAL(rargs_obj_factor)[0] = obj_factor;

			SEXP rargs_lambda;
            PROTECT(rargs_lambda = allocVector(REALSXP,m));
			for (Ipopt::Index i=0;i<m;i++) {
				REAL(rargs_lambda)[i] = lambda[i];
			}
			
			SEXP Rcall, result;
			PROTECT(Rcall = lang4(R_eval_h, rargs_x, rargs_obj_factor, rargs_lambda));
			PROTECT(result = eval(Rcall, R_environment));
			
			for (Ipopt::Index i=0;i<nele_hess;i++) {
				values[i] = REAL(result)[i];
			}
		 
			UNPROTECT(5);
		}

		return true;
	}
}

void IpoptRNLP::finalize_solution(Ipopt::SolverReturn status,
                              Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
                              Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                              Ipopt::Number obj_value,
			                  const Ipopt::IpoptData* ip_data,
			                  Ipopt::IpoptCalculatedQuantities* ip_cq)
{
	// Here we convert the results from c++ to an SEXP list with elements
    // 0. 	status;      integer with convergence status
	// 1.   message;     string with convergence status
	// 2.   iterations;  number of iterations
	// 3.   objective;   final value of the objective function
    // 4.   solution;    final values for the control variables
	// 5.   z_L;         final values for the lower bound multipliers
	// 6.   z_U;         final values for the upper bound multipliers
	// 7.   constraints; final values for the constraints
	// 8.   lambda;      final values for the Lagrange mulipliers
	int num_return_elements = 9;
  
    // R_result_list is a member object, which has been protected in the constructor
    // and will be unprotected in the destructor.
	PROTECT(R_result_list = allocVector(VECSXP, num_return_elements));
	d_num_protected_members++;
	
	// attach names to the return list
	SEXP names;
	PROTECT(names = allocVector(STRSXP, num_return_elements));
	
	SET_STRING_ELT(names, 0, mkChar("status"));
	SET_STRING_ELT(names, 1, mkChar("message"));
	SET_STRING_ELT(names, 2, mkChar("iterations"));
	SET_STRING_ELT(names, 3, mkChar("objective"));
	SET_STRING_ELT(names, 4, mkChar("solution"));
	SET_STRING_ELT(names, 5, mkChar("z_L"));
	SET_STRING_ELT(names, 6, mkChar("z_U"));
	SET_STRING_ELT(names, 7, mkChar("constraints"));
	SET_STRING_ELT(names, 8, mkChar("lambda"));
	setAttrib(R_result_list, R_NamesSymbol, names);
	
	// convert status to an R object
	SEXP R_status;
	PROTECT(R_status = allocVector(INTSXP,1));
	INTEGER(R_status)[0] = (int) status;
	
	
	SEXP R_status_message;
	PROTECT(R_status_message = allocVector(STRSXP, 1));
	switch ( status )
    {
        case Ipopt::SUCCESS: 
            SET_STRING_ELT(R_status_message, 0, mkChar("SUCCESS: Algorithm terminated successfully at a locally optimal point, satisfying the convergence tolerances (can be specified by options)."));
                break;
        case Ipopt::MAXITER_EXCEEDED: 
            SET_STRING_ELT(R_status_message, 0, mkChar("MAXITER_EXCEEDED: Maximum number of iterations exceeded (can be specified by an option)."));
                break;
        case Ipopt::STOP_AT_TINY_STEP:  
            SET_STRING_ELT(R_status_message, 0, mkChar("STOP_AT_TINY_STEP: Algorithm proceeds with very little progress."));
                break;
        case Ipopt::STOP_AT_ACCEPTABLE_POINT:  
            SET_STRING_ELT(R_status_message, 0, mkChar("STOP_AT_ACCEPTABLE_POINT: Algorithm stopped at a point that was converged, not to ``desired'' tolerances, but to ``acceptable'' tolerances (see the acceptable-... options)."));
                break;
        case Ipopt::LOCAL_INFEASIBILITY:  
            SET_STRING_ELT(R_status_message, 0, mkChar("LOCAL_INFEASIBILITY: Algorithm converged to a point of local infeasibility. Problem may be infeasible."));
                break;
        case Ipopt::USER_REQUESTED_STOP:  
            SET_STRING_ELT(R_status_message, 0, mkChar("USER_REQUESTED_STOP: The user call-back function intermediate_callback (see Section 3.3.4) returned false, i.e., the user code requested a premature termination of the optimization."));
                break;
        case Ipopt::DIVERGING_ITERATES:  
            SET_STRING_ELT(R_status_message, 0, mkChar("DIVERGING_ITERATES: It seems that the iterates diverge."));
                break;
        case Ipopt::RESTORATION_FAILURE:  
            SET_STRING_ELT(R_status_message, 0, mkChar("RESTORATION_FAILURE: Restoration phase failed, algorithm doesn't know how to proceed."));
                break;
        case Ipopt::ERROR_IN_STEP_COMPUTATION:  
            SET_STRING_ELT(R_status_message, 0, mkChar("ERROR_IN_STEP_COMPUTATION: An unrecoverable error occurred while IPOPT tried to compute the search direction."));
                break;
        case Ipopt::INVALID_NUMBER_DETECTED:  
            SET_STRING_ELT(R_status_message, 0, mkChar("INVALID_NUMBER_DETECTED: Algorithm received an invalid number (such as NaN or Inf) from the NLP; see also option check_derivatives_for_naninf."));
                break;
        case Ipopt::INTERNAL_ERROR:  
            SET_STRING_ELT(R_status_message, 0, mkChar("INTERNAL_ERROR: An unknown internal error occurred. Please contact the IPOPT authors through the mailing list."));
                break;
		default:
			SET_STRING_ELT(R_status_message, 0, mkChar("Return status not recognized."));
	
	}

	
	// !!! we add number of iterations in the main program
	
	// convert value of objective function to an R object
	SEXP R_objective;
	PROTECT(R_objective = allocVector(REALSXP,1));
	REAL(R_objective)[0] = obj_value;
	
	// convert the value of the controls to an R object
	SEXP R_solution;
	PROTECT(R_solution = allocVector(REALSXP,n));
	for (Ipopt::Index i=0;i<n;i++) {
		REAL(R_solution)[i] = x[i];
	}
    
    SEXP R_z_L;
    PROTECT(R_z_L = allocVector(REALSXP,n));
    for (Ipopt::Index i=0;i<n;i++) {
        REAL(R_z_L)[i] = z_L[i];
    }
    
    SEXP R_z_U;
    PROTECT(R_z_U = allocVector(REALSXP,n));
    for (Ipopt::Index i=0;i<n;i++) {
        REAL(R_z_U)[i] = z_U[i];
    }
    
	SEXP R_constraints;
	PROTECT(R_constraints = allocVector(REALSXP,m));
	for (Ipopt::Index i=0;i<m;i++) {
		REAL(R_constraints)[i] = g[i];
	}
    
    SEXP R_lambda;
    PROTECT(R_lambda = allocVector(REALSXP,m));
    for (Ipopt::Index i=0;i<m;i++) {
        REAL(R_lambda)[i] = lambda[i];
    }

	// add elements to the list
	SET_VECTOR_ELT(R_result_list, 0, R_status);
	SET_VECTOR_ELT(R_result_list, 1, R_status_message);
	SET_VECTOR_ELT(R_result_list, 3, R_objective);
	SET_VECTOR_ELT(R_result_list, 4, R_solution);
    SET_VECTOR_ELT(R_result_list, 5, R_z_L);
	SET_VECTOR_ELT(R_result_list, 6, R_z_U);
	SET_VECTOR_ELT(R_result_list, 7, R_constraints);
	SET_VECTOR_ELT(R_result_list, 8, R_lambda);
	
	UNPROTECT(num_return_elements);
}
