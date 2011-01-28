// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-04-14

#ifndef __IPPARTNLP_HPP__
#define __IPPARTNLP_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
#include "IpException.hpp"
#include "IpAlgTypes.hpp"
#include "IpReturnCodes.hpp"

#include <string>
#include <map>

namespace Ipopt
{
  // forward declarations
  class IpoptData;
  class IpoptCalculatedQuantities;
  class IteratesVector;

  /** Base class for all parallel NLP's that use standard triplet
   *  matrix format.  The class ParTNLPAdapter then converts this
   *  interface to an NLP class that can be used directly by ipopt.
   *
   *  This interface presents the problem form:
   *  
   *     min f(x)
   *
   *     s.t. gL <= g(x) <= gU
   *
   *          xL <=  x   <= xU
   *
   *  In order to specify an equality constraint, set gL_i = gU_i =
   *  rhs.  The value that indicates "infinity" for the bounds
   *  (i.e. the variable or constraint has no lower bound (-infinity)
   *  or upper bound (+infinity)) is set through the option
   *  nlp_lower_bound_inf and nlp_upper_bound_inf.  To indicate that a
   *  variable has no upper or lower bound, set the bound to
   *  -ipopt_inf or +ipopt_inf respectively.
   *
   *  In the MPI-based parallel version, the methods of this class are
   *  executed on each processor.  Each method receives the total
   *  number of processors (num_proc) and its particular processor ID
   *  (0 <= proc_id < num_proc).
   *
   *  The constraints are partitioned across the processors, and the
   *  partition is specified as a collection of disjoint num_proc
   *  sub-intervals of [0..m-1] (or [1..m], depending on index_style,
   *  see get_nlp_info) where m is the total number of constraints.
   *  The user determines the partition of constraints and returns it
   *  to Ipopt. Each processor sets m_first to the first constraint
   *  index it is responsible for, and m_last to the last
   *  index. Subsequently this partition is fixed throughout the
   *  optimization. For the constraints with indices in [m_first
   *  .. m_last] on a processor, the user returns their bounds,
   *  computes their values and Jacobians, and only the Hessian
   *  entries for these constraints.
   *
   *  The Hessian of the Lagrangian also includes entries from the
   *  objective function.  The user can choose to partition these
   *  entries across all processors, so that the sum of the Hessian
   *  matrices computed on all processors is the overall Hessian of
   *  the Lagrangian.
   *
   *  The user also provides a partition of the primal x variables as
   *  sub-intervals [n_first .. n_last] of [0 .. n-1] (or [1..n],
   *  depending on index_style).  When asked for certain vector
   *  components (bounds or objective gradient), each processor is
   *  responsible for only setting its part.  The computation of the
   *  objective function is also partitioned; the objective function
   *  value is computed as the sum of the eval_f value returned on all
   *  processors.  However, all evaluation methods receive a vector x
   *  with ALL components.
   *
   */
  class ParTNLP : public ReferencedObject
  {
  public:
    /** Type of the constraints*/
    enum LinearityType
    {
      LINEAR/** Constraint/Variable is linear.*/,
      NON_LINEAR/**Constraint/Varaible is non-linear.*/
    };

    /**@name Constructors/Destructors */
    //@{
    ParTNLP()
    {}

    /** Default destructor */
    virtual ~ParTNLP()
    {}
    //@}

    DECLARE_STD_EXCEPTION(INVALID_PARTNLP);

    /**@name methods to gather information about the NLP */
    //@{
    /** Basic NLP information. Overload this method to return the
     *  total number of variables (n) and constraints (m), the
     *  sub-interval of variable indices this processor is responsible
     *  for (all variables with indices in [n_first .. n_last]), and the
     *  sub-interval of constraint indices this processor is responsible
     *  for ([m_first .. m_last]).  The numbering starts according to the
     *  choice of index_style.
     *
     *  The number of nonzeros for the constraint Jacobian (nnz_jac_g)
     *  and Hessian of the Lagrangian (nnz_h_lag) is only for those
     *  constraints that this processor is responsible for.
     *
     *  The index_style parameter lets you specify C or Fortran style
     *  indexing for the sparse matrix iRow and jCol parameters, as
     *  well as the constraint and variable indices in the
     *  partitioning.  C_STYLE is 0-based (counting starts at 0), and
     *  FORTRAN_STYLE is 1-based (counting starts at 1).
     */
    enum IndexStyleEnum { C_STYLE=0, FORTRAN_STYLE=1 };
    virtual bool get_nlp_info(Index num_proc, Index proc_id,
                              Index& n, Index& n_first, Index& n_last,
                              Index& m, Index& m_first, Index& m_last,
                              Index& nnz_jac_g_part, Index& nnz_h_lag_part,
                              IndexStyleEnum& index_style)=0;

    typedef std::map<std::string, std::vector<std::string> > StringMetaDataMapType;
    typedef std::map<std::string, std::vector<Index> > IntegerMetaDataMapType;
    typedef std::map<std::string, std::vector<Number> > NumericMetaDataMapType;

    /** overload this method to return any meta data for
     *  the variables and the constraints */
    virtual bool get_var_con_metadata(Index num_proc, Index proc_id,
                                      Index n, Index n_first, Index n_last,
                                      StringMetaDataMapType& var_string_md,
                                      IntegerMetaDataMapType& var_integer_md,
                                      NumericMetaDataMapType& var_numeric_md,
                                      Index m, Index m_first, Index m_last,
                                      StringMetaDataMapType& con_string_md,
                                      IntegerMetaDataMapType& con_integer_md,
                                      NumericMetaDataMapType& con_numeric_md)

    {
      return false;
    }

    /** Return bounds information.  Overload this method to return the
     *  information about the bound on the variables and
     *  constraints. The value that indicates that a bound does not
     *  exist is specified in the parameters nlp_lower_bound_inf and
     *  nlp_upper_bound_inf.  By default, nlp_lower_bound_inf is -1e19
     *  and nlp_upper_bound_inf is 1e19.
     *
     *  Each processor receives only the section of the bound array
     *  for the variables and constraints that it is responsible for,
     *  e.g., x_l and x_u have the size 'n_last-n_first+1'. */
    virtual bool get_bounds_info(Index num_proc, Index proc_id,
                                 Index n, Index n_first, Index n_last,
                                 Number* x_l_part, Number* x_u_part,
                                 Index m, Index m_first, Index m_last,
                                 Number* g_l_part, Number* g_u_part)=0;

    /** overload this method to return scaling parameters. This is
     *  only called if the options are set to retrieve user scaling.
     *  There, use_x_scaling (or use_g_scaling) should get set to true
     *  only if the variables (or constraints) are to be scaled.  This
     *  method should return true only if the scaling parameters could
     *  be provided.
     *
     *  Each processor receives only the section of the bound array
     *  for the variables and constraints that it is responsible for.
     *  The objective scaling factor is taken from the prociess with
     *  proc_id 0.  use_x_scaling and use_g_scaling must be set to the
     *  same on each process.
     */
    virtual bool get_scaling_parameters(Index num_proc, Index proc_id,
                                        Number& obj_scaling,
                                        bool& use_x_scaling,
                                        Index n, Index n_first, Index n_last,
                                        Number* x_scaling_part,
                                        bool& use_g_scaling,
                                        Index m, Index m_first, Index m_last,
                                        Number* g_scaling_part)
    {
      return false;
    }

    /** Compute the starting point.  Overload this method to return
     *  the starting point. The bool variables indicate whether the
     *  algorithm wants you to initialize x, z_L/z_u, and lambda,
     *  respectively.  If, for some reason, the algorithm wants you to
     *  initialize these and you cannot, return false, which will
     *  cause Ipopt to stop.  You will have to run Ipopt with
     *  different options then.
     *
     *  The x, z_L and z_U variables are partitioned according to
     *  n_first and n_last, and the lambda variables according to
     *  m_first and m_last.
     */
    virtual bool
    get_starting_point(Index num_proc, Index proc_id,
                       Index n, Index n_first, Index n_last,
                       bool init_x, Number* x_part,
                       bool init_z, Number* z_L_part, Number* z_U_part,
                       Index m, Index m_first, Index m_last,
                       bool init_lambda, Number* lambda_part)=0;

    /** Compute the objective function. Overload this method to return
     *  a "component" of the value of the objective function; the
     *  actual objective function value is the sum of the values
     *  returned in obj_value across all processors
     */
    virtual bool eval_f(Index num_proc, Index proc_id,
                        Index n, Index n_first, Index n_last,
                        const Number* x, bool new_x, Number& obj_value)=0;

    /** Compute the gradient of the objective function.  Overload this
     *  method to return the vector of the gradient of the objective
     *  w.r.t. the x variables with indices in [n_first .. n_last];
     *  grad_f[0] (if C_STYLE indexing is used, grad_f[1] otherwise)
     *  contains the partial derivative of f with respect to
     *  x_{n_first}.
     */
    virtual bool eval_grad_f(Index num_proc, Index proc_id,
                             Index n,  Index n_first, Index n_last,
                             const Number* x, bool new_x,
                             Number* grad_f_part)=0;

    /** Compute constraint values.  Overload this method to return the
     *  vector of constraint values with indices in [m_first
     *  .. m_last]; g[0] (if C_STYLE indexing is used, g[1] otherwise)
     *  contains the value of the constraint g_{m_first} evaluated at
     *  x = (x_0, .. x_{n-1}).
     */
    virtual bool eval_g(Index num_proc, Index proc_id,
                        Index n, const Number* x, bool new_x,
                        Index m, Index m_first, Index m_last,
                        Number* g_part)=0;

    /** Compute structure and values of the constraint Jacobian.
     *  Overload this method to return the jacobian of the constraints
     *  with indices in [m_first .. m_last].  The vectors iRow and
     *  jCol only need to be set once. The first call is used to set
     *  the structure only (iRow and jCol will be non-NULL, and values
     *  will be NULL) For subsequent calls, iRow and jCol will be
     *  NULL.  iRow must only contain values in [0..(m_last-m_first)],
     *  the sub-interval of constraints the processor proc_num is
     *  responsible for, while jCol will contain values from 0 to n-1
     *  (assuming C_STYLE indexing).  This means, the row counting
     *  starts at 0 (or 1) for the constraints for each individual
     *  processor.
     */
    virtual bool eval_jac_g(Index num_proc, Index proc_id,
                            Index n, const Number* x, bool new_x,
                            Index m, Index m_first, Index m_last,
                            Index nele_jac_part, Index* iRow_part,
                            Index *jCol_part, Number* values_part)=0;

    /** Compute structure and values of the Hessian of the Lagrangian
     *  function.  Overload this method to return the part of the
     *  Hessian of the Lagrangian the processor proc_id is responsible
     *  for.  The vectors iRow and jCol only need to be set once
     *  (during the first call). The first call is used to set the
     *  structure only (iRow and jCol will be non-NULL, and values
     *  will be NULL) For subsequent calls, iRow and jCol will be
     *  NULL. This matrix is symmetric - specify the lower diagonal
     *  (or upper diagonal) part only.  A default implementation is
     *  provided, in case the user wants to use quasi-Newton
     *  approximations to estimate the second derivatives and doesn't
     *  not need to implement this method.  iRow and jCol will contain
     *  values from 0 to n-1 (assuming C_STYLE indexing).
     *
     *  The vector lambda contains the multipliers for ALL constraints
     *  (not only those corresponding to the constraints that this
     *  processor is responsible for).  The user can partition the
     *  computation in any way, as long as the sum of all matrices
     *  from all processors give the overal Hessian matrix.  Note that
     *  non-zero positions do not have to be disjunct among
     *  processors.
     */
    virtual bool eval_h(Index num_proc, Index proc_id,
                        Index n, Index n_first, Index n_last,
                        const Number* x, bool new_x, Number obj_factor,
                        Index m, Index m_first, Index m_last,
                        const Number* lambda,
                        bool new_lambda, Index nele_hess_part,
                        Index* iRow_part, Index* jCol_part,
                        Number* values_part)
    {
      return false;
    }
    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the
     *  ParTNLP can store/write the solution.  All vectors (x, g, etc)
     *  are provided as full arrays, not just with portitions
     *  according to the partitioning. */
    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq)=0;

    /** Intermediate Callback method for the user.  Providing dummy
     *  default implementation.  For details see IntermediateCallBack
     *  in IpNLP.hpp. */
    virtual bool intermediate_callback(AlgorithmMode mode,
                                       Index iter, Number obj_value,
                                       Number inf_pr, Number inf_du,
                                       Number mu, Number d_norm,
                                       Number regularization_size,
                                       Number alpha_du, Number alpha_pr,
                                       Index ls_trials,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq)
    {
      return true;
    }
    //@}

    /** @name Methods for quasi-Newton approximation.  If the second
     *  derivatives are approximated by Ipopt, it is better to do this
     *  only in the space of nonlinear variables.  The following
     *  methods are call by Ipopt if the quasi-Newton approximation is
     *  selected.  If -1 is returned as number of nonlinear variables,
     *  Ipopt assumes that all variables are nonlinear.  Otherwise, it
     *  calls get_list_of_nonlinear_variables with an array into which
     *  the indices of the nonlinear variables should be written - the
     *  array has the lengths num_nonlin_vars, which is identical with
     *  the return value of get_number_of_nonlinear_variables().  It
     *  is assumed that the indices are counted starting with 1 in the
     *  FORTRAN_STYLE, and 0 for the C_STYLE. */
    //@{
    virtual Index get_number_of_nonlinear_variables()
    {
      return -1;
    }

    virtual bool get_list_of_nonlinear_variables(Index num_nonlin_vars,
        Index* pos_nonlin_vars)
    {
      return false;
    }
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    //ParTNLP();

    /** Copy Constructor */
    ParTNLP(const ParTNLP&);

    /** Overloaded Equals Operator */
    void operator=(const ParTNLP&);
    //@}
  };

} // namespace Ipopt

#endif
