/* Copyright (C) 2004, 2010 International Business Machines and others.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02
 */

#ifndef __IPSTDCINTERFACE_H__
#define __IPSTDCINTERFACE_H__

#include <stdbool.h>

#include "IpoptConfig.h"
#include "IpTypes.h"
#include "IpReturnCodes.h"

#ifndef IPOPT_EXPORT
/** @deprecated Use IPOPT_CALLCONV instead. */
#define IPOPT_EXPORT(type) type IPOPT_CALLCONV
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/** Type for all number.
 * @deprecated Use ipnumber instead.
 */
IPOPT_DEPRECATED
typedef ipnumber Number;

/** Type for all indices.
 * @deprecated Use ipindex instead.
 */
IPOPT_DEPRECATED
typedef int Index;

/** Type for all integers.
 * @deprecated Use int instead.
 */
IPOPT_DEPRECATED
typedef int Int;

/** Structure collecting all information about the problem definition and solve statistics etc. */
struct IpoptProblemInfo;

/** Pointer to an Ipopt Problem. */
typedef struct IpoptProblemInfo* IpoptProblem;

/** define a boolean type for C
 * @deprecated Use bool instead.
 */
typedef bool Bool;
#ifndef TRUE
/* @deprecated Use true instead. */
# define TRUE (1)
#endif
/* @deprecated Use false instead. */
#ifndef FALSE
# define FALSE (0)
#endif

/** A pointer for anything that is to be passed between the called and individual callback function */
typedef void* UserDataPtr;

/** Type defining the callback function for evaluating the value of the objective function.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also Ipopt::TNLP::eval_f.
 */
typedef bool (*Eval_F_CB)(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   obj_value,
   UserDataPtr user_data
);

/** Type defining the callback function for evaluating the gradient of the objective function.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also Ipopt::TNLP::eval_grad_f.
 */
typedef bool (*Eval_Grad_F_CB)(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber*   grad_f,
   UserDataPtr user_data
);

/** Type defining the callback function for evaluating the value of the constraint functions.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also Ipopt::TNLP::eval_g.
 */
typedef bool (*Eval_G_CB)(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipindex     m,
   ipnumber*   g,
   UserDataPtr user_data
);

/** Type defining the callback function for evaluating the Jacobian of the constrant functions.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also Ipopt::TNLP::eval_jac_g.
 */
typedef bool (*Eval_Jac_G_CB)(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipindex     m,
   ipindex     nele_jac,
   ipindex*    iRow,
   ipindex*    jCol,
   ipnumber*   values,
   UserDataPtr user_data
);

/** Type defining the callback function for evaluating the Hessian of the Lagrangian function.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also Ipopt::TNLP::eval_h.
 */
typedef bool (*Eval_H_CB)(
   ipindex     n,
   ipnumber*   x,
   bool        new_x,
   ipnumber    obj_factor,
   ipindex     m,
   ipnumber*   lambda,
   bool        new_lambda,
   ipindex     nele_hess,
   ipindex*    iRow,
   ipindex*    jCol,
   ipnumber*   values,
   UserDataPtr user_data
);

/** Type defining the callback function for giving intermediate execution control to the user.
 *
 *  If set, it is called once per iteration, providing the user
 *  with some information on the state of the optimization.
 *  This can be used to print some user-defined output.
 *  It also gives the user a way to terminate the optimization prematurely.
 *  If this method returns false, Ipopt will terminate the optimization.
 *
 *  See also Ipopt::TNLP::intermediate_callback.
 */
typedef bool (*Intermediate_CB)(
   ipindex     alg_mod,   /**< algorithm mode: 0 is regular, 1 is restoration */
   ipindex     iter_count,/**< iteration number */
   ipnumber    obj_value, /**< objective function value */
   ipnumber    inf_pr,    /**< primal infeasibility */
   ipnumber    inf_du,    /**< dual infeasibility */
   ipnumber    mu,        /**< barrier parameter */
   ipnumber    d_norm,    /**< infinity-norm of primal step */
   ipnumber    regularization_size,  /**< size of regularization of Hessian of Lagrangian */
   ipnumber    alpha_du,  /**< step length for dual variables */
   ipnumber    alpha_pr,  /**< step length for primal variables */
   ipindex     ls_trials, /**< number of backtracking line search steps */
   UserDataPtr user_data  /**< user data */
);

/** Function for creating a new Ipopt Problem object.
 *
 *  This function returns an object that can be passed to the IpoptSolve call.
 *  It contains the basic definition of the optimization problem, such
 *  as number of variables and constraints, bounds on variables and
 *  constraints, information about the derivatives, and the callback
 *  function for the computation of the optimization problem
 *  functions and derivatives.  During this call, the options file
 *  PARAMS.DAT is read as well.
 *
 *  If NULL is returned, there was a problem with one of the inputs
 *  or reading the options file.
 *
 *  See also Ipopt::TNLP::get_nlp_info and Ipopt::TNLP::get_bounds_info.
 */
IPOPTLIB_EXPORT IpoptProblem IPOPT_CALLCONV CreateIpoptProblem(
   ipindex       n,           /**< Number of optimization variables */
   ipnumber*     x_L,         /**< Lower bounds on variables
                               *
                               * This array of size n is copied internally, so that the
                               * caller can change the incoming data after
                               * return without that IpoptProblem is modified.
                               * Any value less or equal than the number specified by option
                               * 'nlp_lower_bound_inf' is interpreted to be minus infinity.
                               */
   ipnumber*      x_U,         /**< Upper bounds on variables
                               *
                               * This array of size n is copied internally,
                               * so that the caller can change the incoming data after
                               * return without that IpoptProblem is modified.
                               * Any value greater or equal than the number specified by option
                               * 'nlp_upper_bound_inf' is interpreted to be plus infinity.
                               */
   ipindex        m,           /**< Number of constraints */
   ipnumber*      g_L,         /**< Lower bounds on constraints
                               *
                               * This array of size m is copied internally,
                               * so that the caller can change the incoming data after
                               * return without that IpoptProblem is modified.
                               * Any value less or equal than the number specified by option
                               * 'nlp_lower_bound_inf' is interpreted to be minus infinity.
                               */
   ipnumber*      g_U,         /**< Upper bounds on constraints
                               *
                               * This array of size m is copied internally,
                               * so that the caller can change the incoming data after
                               * return without that IpoptProblem is modified.
                               * Any value greater or equal than the number specified by option
                               * 'nlp_upper_bound_inf' is interpreted to be plus infinity.
                               */
   ipindex        nele_jac,    /**< Number of non-zero elements in constraint Jacobian */
   ipindex        nele_hess,   /**< Number of non-zero elements in Hessian of Lagrangian */
   ipindex        index_style, /**< Indexing style for iRow & jCol, 0 for C style, 1 for Fortran style */
   Eval_F_CB      eval_f,      /**< Callback function for evaluating objective function */
   Eval_G_CB      eval_g,      /**< Callback function for evaluating constraint functions */
   Eval_Grad_F_CB eval_grad_f, /**< Callback function for evaluating gradient of objective function */
   Eval_Jac_G_CB  eval_jac_g,  /**< Callback function for evaluating Jacobian of constraint functions */
   Eval_H_CB      eval_h       /**< Callback function for evaluating Hessian of Lagrangian function */
);

/** Method for freeing a previously created IpoptProblem.
 *
 * After freeing an IpoptProblem, it cannot be used anymore.
 */
IPOPTLIB_EXPORT void IPOPT_CALLCONV FreeIpoptProblem(
   IpoptProblem ipopt_problem
);

/** Function for adding a string option.
 *
 * @return false, if the option could not be set (e.g., if keyword is unknown)
 */
IPOPTLIB_EXPORT bool IPOPT_CALLCONV AddIpoptStrOption(
   IpoptProblem ipopt_problem,
   char*        keyword,
   char*        val
);

/** Function for adding a Number option.
 *
 * @return false, if the option could not be set (e.g., if keyword is unknown)
 */
IPOPTLIB_EXPORT bool IPOPT_CALLCONV AddIpoptNumOption(
   IpoptProblem ipopt_problem,
   char*        keyword,
   ipnumber     val
);

/** Function for adding an Integer option.
 *
 * @return false, if the option  could not be set (e.g., if keyword is unknown)
 @*/
IPOPTLIB_EXPORT bool IPOPT_CALLCONV AddIpoptIntOption(
   IpoptProblem ipopt_problem,
   char*        keyword,
   ipindex      val
);

/** Function for opening an output file for a given name with given printlevel.
 *
 * @return false, if there was a problem opening the file.
 */
IPOPTLIB_EXPORT bool IPOPT_CALLCONV OpenIpoptOutputFile(
   IpoptProblem ipopt_problem,
   char*        file_name,
   int          print_level
);

/** Optional function for setting scaling parameter for the NLP.
 *
 *  This corresponds to the TNLP::get_scaling_parameters method.
 *  If the pointers x_scaling or g_scaling are NULL, then no scaling
 *  for x resp. g is done.
 */
IPOPTLIB_EXPORT bool IPOPT_CALLCONV SetIpoptProblemScaling(
   IpoptProblem ipopt_problem,
   ipnumber     obj_scaling,
   ipnumber*    x_scaling,
   ipnumber*    g_scaling
);

/** Setting a callback function for the "intermediate callback" method in the TNLP.
 *
 *  This gives control back to the user once
 *  per iteration.  If set, it provides the user with some
 *  information on the state of the optimization.  This can be used
 *  to print some user-defined output.  It also gives the user a way
 *  to terminate the optimization prematurely.  If the callback
 *  method returns false, Ipopt will terminate the optimization.
 *  Calling this set method to set the CB pointer to NULL disables
 *  the intermediate callback functionality. */
IPOPTLIB_EXPORT bool IPOPT_CALLCONV SetIntermediateCallback(
   IpoptProblem    ipopt_problem,
   Intermediate_CB intermediate_cb
);

/** Function calling the Ipopt optimization algorithm for a problem
 * previously defined with CreateIpoptProblem.
 *
 * @return outcome of the optimization procedure (e.g., success, failure etc).
 */
IPOPTLIB_EXPORT enum ApplicationReturnStatus IPOPT_CALLCONV IpoptSolve(
   IpoptProblem ipopt_problem, /**< Problem that is to be optimized.
                                *
                                * Ipopt will use the options previously specified with
                                * AddIpoptOption (etc) for this problem.
                                */
   ipnumber*    x,             /**< Input: Starting point; Output: Optimal solution */
   ipnumber*    g,             /**< Values of constraint at final point (output only; ignored if set to NULL) */
   ipnumber*    obj_val,       /**< Final value of objective function (output only; ignored if set to NULL) */
   ipnumber*    mult_g,        /**< Input: Initial values for the constraint multipliers (only if warm start option is chosen);
                                *  Output: Final multipliers for constraints (ignored if set to NULL)
                                */
   ipnumber*    mult_x_L,      /**< Input: Initial values for the multipliers for lower variable bounds (only if warm start option is chosen);
                                *  Output: Final multipliers for lower variable bounds (ignored if set to NULL)
                                */
   ipnumber*    mult_x_U,      /**< Input: Initial values for the multipliers for upper variable bounds (only if warm start option is chosen);
                                *  Output: Final multipliers for upper variable bounds (ignored if set to NULL)
                                */
   UserDataPtr  user_data      /**< Pointer to user data.
                                *
                                * This will be passed unmodified to the callback functions.
                                */
);

/** Get primal and dual variable values of the current iterate.
 *
 * This method can be used to get the values of the current iterate during the intermediate callback set by SetIntermediateCallback().
 * The method expects the number of variables (dimension of x), number of constraints (dimension of g(x)),
 * and allocated arrays of appropriate lengths as input.
 *
 * The method translates the x(), c(), d(), y_c(), y_d(), z_L(), and z_U() vectors from Ipopt::IpoptData::curr()
 * of the internal NLP representation into the form used by the TNLP.
 * For the correspondence between scaled and unscaled solutions, see the detailed description of Ipopt::OrigIpoptNLP.
 * If %Ipopt is in restoration mode, it maps the current iterate of restoration %NLP (see Ipopt::RestoIpoptNLP) back to the original TNLP.
 *
 * If there are fixed variables and fixed_variable_treatment=make_parameter, then requesting z_L and z_U can trigger a reevaluation of
 * the Gradient of the objective function and the Jacobian of the constraint functions.
 *
 * @param ipopt_problem (in) Problem that is currently optimized.
 * @param n       (in)  the number of variables \f$x\f$ in the problem; can be arbitrary if skipping x, z_L, and z_U
 * @param scaled  (in)  whether to retrieve scaled or unscaled iterate
 * @param x       (out) buffer to store value of primal variables \f$x\f$, must have length at least n; pass NULL to skip retrieving x
 * @param z_L     (out) buffer to store the lower bound multipliers \f$z_L\f$, must have length at least n; pass NULL to skip retrieving z_L and Z_U
 * @param z_U     (out) buffer to store the upper bound multipliers \f$z_U\f$, must have length at least n; pass NULL to skip retrieving z_L and Z_U
 * @param m       (in)  the number of constraints \f$g(x)\f$; can be arbitrary if skipping g and lambda
 * @param g       (out) buffer to store the constraint values \f$g(x)\f$, must have length at least m; pass NULL to skip retrieving g
 * @param lambda  (out) buffer to store the constraint multipliers \f$\lambda\f$, must have length at least m; pass NULL to skip retrieving lambda
 *
 * @return Whether Ipopt has successfully filled the given arrays
 * @since 3.14.0
 */
IPOPTLIB_EXPORT bool IPOPT_CALLCONV GetIpoptCurrentIterate(
   IpoptProblem    ipopt_problem,
   bool            scaled,
   ipindex         n,
   ipnumber*       x,
   ipnumber*       z_L,
   ipnumber*       z_U,
   ipindex         m,
   ipnumber*       g,
   ipnumber*       lambda
);

/** Get primal and dual infeasibility of the current iterate.
 *
 * This method can be used to get the violations of constraints and optimality conditions
 * at the current iterate during the intermediate callback set by SetIntermediateCallback().
 * The method expects the number of variables (dimension of x), number of constraints (dimension of g(x)),
 * and allocated arrays of appropriate lengths as input.
 *
 * The method makes the vectors behind (unscaled_)curr_orig_bounds_violation(), (unscaled_)curr_nlp_constraint_violation(), (unscaled_)curr_dual_infeasibility(),
 * (unscaled_)curr_complementarity() from Ipopt::IpoptCalculatedQuantities of the internal NLP representation available into the form used by the TNLP.
 * If %Ipopt is in restoration mode, it maps the current iterate of restoration %NLP (see Ipopt::RestoIpoptNLP) back to the original TNLP.
 *
 * @note If in restoration phase, then requesting grad_lag_x can trigger a call to Eval_F_CB.
 *
 * @note Ipopt by default relaxes variable bounds (option bound_relax_factor > 0.0).
 *   x_L_violation and x_U_violation report the violation of a solution w.r.t. the original unrelaxed bounds.
 *   However, compl_x_L and compl_x_U use the relaxed variable bounds to calculate the complementarity.
 *
 * @param ipopt_problem (in) Problem that is currently optimized.
 * @param scaled     (in)  whether to retrieve scaled or unscaled violations
 * @param n          (in)  the number of variables \f$x\f$ in the problem; can be arbitrary if skipping compl_x_L, compl_x_U, and grad_lag_x
 * @param x_L_violation (out) buffer to store violation of original lower bounds on variables (max(orig_x_L-x,0)), must have length at least n; pass NULL to skip retrieving orig_x_L
 * @param x_U_violation (out) buffer to store violation of original upper bounds on variables (max(x-orig_x_U,0)), must have length at least n; pass NULL to skip retrieving orig_x_U
 * @param compl_x_L  (out) buffer to store violation of complementarity for lower bounds on variables (\f$(x-x_L)z_L\f$), must have length at least n; pass NULL to skip retrieving compl_x_L
 * @param compl_x_U  (out) buffer to store violation of complementarity for upper bounds on variables (\f$(x_U-x)z_U\f$), must have length at least n; pass NULL to skip retrieving compl_x_U
 * @param grad_lag_x (out) buffer to store gradient of Lagrangian w.r.t. variables \f$x\f$, must have length at least n; pass NULL to skip retrieving grad_lag_x
 * @param m          (in)  the number of constraints \f$g(x)\f$; can be arbitrary if skipping lambda
 * @param nlp_constraint_violation (out) buffer to store violation of constraints \f$max(g_l-g(x),g(x)-g_u,0)\f$, must have length at least m; pass NULL to skip retrieving constraint_violation
 * @param compl_g    (out) buffer to store violation of complementarity of constraint (\f$(g(x)-g_l)*\lambda^+ + (g_l-g(x))*\lambda^-\f$, where \f$\lambda^+=max(0,\lambda)\f$ and \f$\lambda^-=max(0,-\lambda)\f$ (componentwise)), must have length at least m; pass NULL to skip retrieving compl_g
 *
 * @return Whether Ipopt has successfully filled the given arrays
 * @since 3.14.0
 */
IPOPTLIB_EXPORT bool IPOPT_CALLCONV GetIpoptCurrentViolations(
   IpoptProblem  ipopt_problem,
   bool          scaled,
   ipindex       n,
   ipnumber*     x_L_violation,
   ipnumber*     x_U_violation,
   ipnumber*     compl_x_L,
   ipnumber*     compl_x_U,
   ipnumber*     grad_lag_x,
   ipindex       m,
   ipnumber*     nlp_constraint_violation,
   ipnumber*     compl_g
);

#ifdef __cplusplus
} /* extern "C" { */
#endif

#endif
