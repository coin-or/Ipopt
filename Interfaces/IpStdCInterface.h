/*************************************************************************
   Copyright (C) 2004, International Business Machines and others.
   All Rights Reserved.
   This code is published under the Common Public License.

   $Id$

   Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02
 *************************************************************************/

#ifndef __IPSTDCINTERFACE_H__
#define __IPSTDCINTERFACE_H__

#ifdef __cplusplus
extern "C" {
#endif

  /** This has to be same as in IpTypes.hpp */
  typedef int Index;
  /** This has to be same as in IpTypes.hpp */
  typedef int Int;
  /** This has to be same as in IpTypes.hpp */
  typedef double Number;
  typedef int Bool;
  typedef void * UserDataPtr;
  typedef void * OptionsPtr;

  /** Type defining the callback function for evaluating the value of
   *  the objective function. */
  typedef Bool (*Eval_F_CB)(Index n, Number* x, Bool new_x,
			    Number* obj_value, UserDataPtr user_data);

  /** Type defining the callback function for evaluating the gradient of
   *  the objective function. */
  typedef Bool (*Eval_Grad_F_CB)(Index n, Number* x, Bool new_x,
				 Number* grad_f, UserDataPtr user_data);

  /** Type defining the callback function for evaluating the value of
   *  the constraint functions. */
  typedef Bool (*Eval_G_CB)(Index n, Number* x, Bool new_x,
			    Index m, Number* g, UserDataPtr user_data);

  /** Type defining the callback function for evaluating the Jacobian of
   *  the constrant functions */
  typedef Bool (*Eval_Jac_G_CB)(Index n, Number *x, Bool new_x, Index nele_jac,
				Index *iRow, Index *jCol, Number *values,
				UserDataPtr user_data);

  /** Type defining the callback function for evaluating the Hessian of
   *  the Lagrangian function */
  typedef Bool (*Eval_H_CB)(Index n, Number *x, Bool new_x, Number obj_factor,
			    Index m, Number *lambda, Bool new_lambda, 
			    Index nele_hess, Index *iRow, Index *jCol,
			    Number *values, UserDataPtr user_data);

  /** Function for creating a new Options Object */
  OptionsPtr Ipopt_NewOptions(void);

  /** Function for deleting an Options Object */
  void Ipopt_DeleteOptions(OptionsPtr options);

  /** Function for adding a string option */
  Bool Ipopt_AddOption(OptionsPtr options, char* keyword, char* val);

  /** Function for adding a Number option */
  Bool Ipopt_AddNumOption(OptionsPtr options, char* keyword, Number val);

  /** Function for adding an Int option */
  Bool Ipopt_AddIntOption(OptionsPtr options, char* keyword, Int val);

  /** Main IPOPT function.  This is the main entrance point for the
   *  standard C interface.  Returns error code (0 = success)
   */
  Int IpoptSolve(Index n             /** Number of optimization variables */
               , Number* x           /** Input:  Starting point
                                         Output: Optimal solution */
	       , Number* x_L         /** Lower bounds on variables */
	       , Number* x_U         /** Upper bounds on variables */
	       , Index m             /** Number of constraints */
	       , Number* g           /** Values of constraint at final point
                                         (output only - ignored if set to NULL) */
	       , Number* g_L         /** Lower bounds on constraints */
               , Number* g_U         /** Upper bounds on constraints */
	       , Index nele_jac      /** Number of non-zero elements in
					 constraint Jacobian */
	       , Index nele_hess     /** Number of non-zero elements in
					 Hessian of Lagrangian */
               , Number* obj_val     /** Final value of objective function
                                         (output only - ignored if set to NULL) */
	       , Number* mult_g      /** Final multipliers for constraints
                                         (output only - ignored if set to NULL) */
	       , Number* mult_x_L    /** Final multipliers for lower variable bounds
                                         (output only - ignored if set to NULL) */
	       , Number* mult_x_U    /** Final multipliers for upper variable bounds
                                         (output only - ignored if set to NULL) */
	       , Eval_F_CB eval_f    /** Callback function for evaluating
					 objective function */
	       , Eval_G_CB eval_g    /** Callback function for evaluating
					 constraint functions */
	       , Eval_Grad_F_CB eval_grad_f
                                     /** Callback function for evaluating
				         gradient of objective function */
               , Eval_Jac_G_CB eval_jac_g
                                     /** Callback function for evaluating
				         Jacobian of constraint functions */
               , Eval_H_CB eval_h    /** Callback function for evaluating
                                         Hessian of Lagrangian function */
               , OptionsPtr options  /** Options for algorithm (NULL: default) */
               , UserDataPtr user_data /** Pointer to user data */
	         );

#ifdef __cplusplus
} /* extern "C" { */
#endif

#endif
