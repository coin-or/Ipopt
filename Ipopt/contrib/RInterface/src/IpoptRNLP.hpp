/*
 * Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * file:   IpoptRNLP.hpp
 * author: Jelmer Ypma
 * date:   18 April 2010
 *
 * This file defines a C++ class that derives from Ipopt::TNLP. The class
 * takes care of interaction between Ipopt and user-defined functions in R.
 *
 * Financial support of the UK Economic and Social Research Council 
 * through a grant (RES-589-28-0001) to the ESRC Centre for Microdata 
 * Methods and Practice (CeMMAP) is gratefully acknowledged.
 */

#ifndef __IpoptRNLP_HPP__
#define __IpoptRNLP_HPP__

#include "IpTNLP.hpp"               // ISA TNLP

#include <assert.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
// Rdefines.h is somewhat more higher level then Rinternal.h, and is preferred if the code might be shared with S at any stage.
// Utils.h defines void R_CheckUserInterrupt(void); to allow user interuption from R


class IpoptRNLP : public Ipopt::TNLP
{
	SEXP R_environment;	            // this is the environment that the function gets called in. This environment can be used
                                    // to pass common data to the R functions.
	SEXP R_eval_f;			        // objective function
	SEXP R_eval_grad_f;		        // gradient of objective function
	
	SEXP R_init_values;		        // vector with initial values, we get the Ipopt::Number of controls from the length of this vector
	
	SEXP R_lower_bounds;	        // lower bounds of the control x
	SEXP R_upper_bounds;	        // upper bounds of the control x
	
	SEXP R_eval_g;			        // function to evaluate constraints
	SEXP R_eval_jac_g;			    // function to evaluate jacobian of constraints
	SEXP R_eval_jac_g_structure;	// list with non-zero elements in the Jacobian, this defines the sparse structure
	
	SEXP R_constraint_lower_bounds;	// lower bounds of the contraint function g()
	SEXP R_constraint_upper_bounds;	// upper bounds of the contraint function g()
	
	SEXP R_eval_h;                  // function to evaluate Hessian
	SEXP R_eval_h_structure;        // list with non-zero elements of the Hessian, this defines the sparse structure
	
	SEXP R_result_list;				// structure that will contain the return values
	
	bool d_hessian_approximation;   // should we approximate the Hessian? default: false
	
    int d_num_protected_members;    // counter of the number of PROTECT calls of the SEXPs above
public:
  /** default constructor */
  IpoptRNLP();

  /** default destructor */
  virtual ~IpoptRNLP();

  void set_R_environment( SEXP env );
  
  
  void set_R_eval_f( SEXP f );
  void set_R_eval_grad_f( SEXP f );
  
  void set_R_init_values( SEXP x0 );
  void set_R_lower_bounds( SEXP lb );
  void set_R_upper_bounds( SEXP ub );
  
  void set_R_eval_g( SEXP g );
  void set_R_eval_jac_g( SEXP g );
  void set_R_eval_jac_g_structure( SEXP s );
  
  void set_R_constraint_lower_bounds( SEXP lb );
  void set_R_constraint_upper_bounds( SEXP ub );
  
  void set_R_eval_h( SEXP h );
  void set_R_eval_h_structure( SEXP s );
  
  void set_hessian_approximation( bool b );
  
  SEXP get_R_result_list();
  
  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                            Ipopt::Index& nnz_h_lag, IndexStyleEnum& Index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                               Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                  bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                  Ipopt::Index m, bool init_lambda,
                                  Ipopt::Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                          Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
                          Ipopt::Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                      bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                      Ipopt::Index* jCol, Ipopt::Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(Ipopt::SolverReturn status,
                                 Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
                                 Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                                 Ipopt::Number obj_value,
				                 const Ipopt::IpoptData* ip_data,
				                 Ipopt::IpoptCalculatedQuantities* ip_cq);
  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  IpoptRNLP();
  IpoptRNLP(const IpoptRNLP&);
  IpoptRNLP& operator=(const IpoptRNLP&);
  //@}
};


#endif
