// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Johannes Huber, Andreas Waechter     IBM    2010-09-03

#ifndef __IPLIBMESHPDENLP_HPP__
#define __IPLIBMESHPDENLP_HPP__

#include "IpParTNLP.hpp"
#include "IpLibMeshPDE.hpp"

namespace Ipopt
{
  // forward declarations
  class Journalist;

  /** Instantiation of a ParTNLP to formulate PDE constrained
      optimization problems where the PDE is provided by libMesh.*/
  class LibMeshPDENLP : public ParTNLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    LibMeshPDENLP(LibMeshPDEBase& libmeshPDE,
                  Journalist& jnlst);

    /** Default destructor */
    virtual ~LibMeshPDENLP()
    {}
    //@}

    /**@name methods to gather information about the NLP. Overloaded
       from ParTNLP. */
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
    virtual bool get_nlp_info(Index num_proc, Index proc_id,
                              Index& n, Index& n_first, Index& n_last,
                              Index& m, Index& m_first, Index& m_last,
                              Index& nnz_jac_g_part, Index& nnz_h_lag_part,
                              IndexStyleEnum& index_style);

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
                                 Number* g_l_part, Number* g_u_part);

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
                       bool init_lambda, Number* lambda_part);

    /** Compute the objective function. Overload this method to return
     *  a "component" of the value of the objective function; the
     *  actual objective function value is the sum of the values
     *  returned in obj_value across all processors
     */
    virtual bool eval_f(Index num_proc, Index proc_id,
                        Index n, Index n_first, Index n_last,
                        const Number* x, bool new_x, Number& obj_value);

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
                             Number* grad_f_part);

    /** Compute constraint values.  Overload this method to return the
     *  vector of constraint values with indices in [m_first
     *  .. m_last]; g[0] (if C_STYLE indexing is used, g[1] otherwise)
     *  contains the value of the constraint g_{m_first} evaluated at
     *  x = (x_0, .. x_{n-1}).
     */
    virtual bool eval_g(Index num_proc, Index proc_id,
                        Index n, const Number* x, bool new_x,
                        Index m, Index m_first, Index m_last,
                        Number* g_part);

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
                            Index *jCol_part, Number* values_part);

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
                        Number* values_part);
    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the
     *  LibMeshPDENLP can store/write the solution.  All vectors (x, g, etc)
     *  are provided as full arrays, not just with portitions
     *  according to the partitioning. */
    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq);

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
    LibMeshPDENLP();

    /** Copy Constructor */
    LibMeshPDENLP(const LibMeshPDENLP&);

    /** Overloaded Equals Operator */
    void operator=(const LibMeshPDENLP&);
    //@}

    /** Object containing all information regarding the PDE */
    LibMeshPDEBase* libmeshPDE_;
    int control_first_;
    int control_last_;
    int state_first_;
    int state_last_;
    int n_first_;
    int n_last_;
    int m_first_;
    int m_last_;
    int pde_first_;
    int pde_last_;
    int aux_first_;
    int aux_last_;
    SmartPtr<Journalist> jnlst_;

    void optim_var_libMesh2local(libMesh::NumericVector<libMesh::Number>& state, libMesh::NumericVector<libMesh::Number>& control, Number* plocalDest);

    void optim_var_global2libMesh(const Number* pglobal, libMesh::NumericVector<libMesh::Number>& state, libMesh::NumericVector<libMesh::Number>& control);
    void update_x(const Number* x);
    int GlobStateControlIdx2OptVarIdx(int i_var,int num_of_procs,
                                      const int* proc_state_first_arr, const int* proc_state_last_arr,
                                      const int* proc_var_first, const int* proc_loc_offset, int& prev_proc);


  };

} // namespace Ipopt

#endif
