// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPTNLP_HPP__
#define __IPTNLP_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
#include "IpException.hpp"

namespace Ipopt
{

  /** Base class for all NLP's that use standard triplet matrix form
   *  and dense vectors.  This is the standard base class for all
   *  NLP's that use the standard triplet matrix form (as for Harwell
   *  routines) and dense vectors. The class TNLPAdapter then converts
   *  this interface to an interface that can be used directly by
   *  ipopt.
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
   *  -ipopt_inf or +ipopt_inf respectively
   */
  class TNLP : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    TNLP()
    {}
    ;

    /** Default destructor */
    virtual ~TNLP()
    {}
    ;
    //@}

    DECLARE_STD_EXCEPTION(INVALID_TNLP);

    /**@name methods to gather information about the NLP */
    //@{
    /** overload this method to return the number of variables
     *  and constraints, and the number of non-zeros in the jacobian and
     *  the hessian. */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag)=0;

    /** overload this method to return the information about the bound
     *  on the variables and constraints. The value that indicates
     *  that a bound does not exist is specified in the parameters
     *  nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
     *  nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is
     *  1e19. (see TNLPAdapter) */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u)=0;

    /** overload this method to return the starting point. The bools
     *  init_x and init_lambda are both inputs and outputs. As inputs,
     *  they indicate whether or not the algorithm wants you to
     *  initialize x and lambda respectively. If, for some reason, the
     *  algorithm wants you to initialize these and you cannot, set
     *  the respective bool to false.
     */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda)=0;

    /** overload this method to return the value of the objective function */
    virtual bool eval_f(Index n, const Number* x, bool new_x,
                        Number& obj_value)=0;

    /** overload this method to return the vector of the gradient of
     *  the objective w.r.t. x */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
                             Number* grad_f)=0;

    /** overload this method to return the vector of constraint values */
    virtual bool eval_g(Index n, const Number* x, bool new_x,
                        Index m, Number* g)=0;
    /** overload this method to return the jacobian of the
     *  constraints. The vectors iRow and jCol only need to be set
     *  once. The first call is used to set the structure only (iRow
     *  and jCol will be non-NULL, and values will be NULL) For
     *  subsequent calls, iRow and jCol will be NULL. */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow,
                            Index *jCol, Number* values)=0;

    /** overload this method to return the hessian of the
     *  lagrangian. The vectors iRow and jCol only need to be set once
     *  (during the first call). The first call is used to set the
     *  structure only (iRow and jCol will be non-NULL, and values
     *  will be NULL) For subsequent calls, iRow and jCol will be
     *  NULL. This matrix is symmetric - specify the lower diagonal
     *  only */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess,
                        Index* iRow, Index* jCol, Number* values)=0;
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
    //TNLP();

    /** Copy Constructor */
    TNLP(const TNLP&);

    /** Overloaded Equals Operator */
    void operator=(const TNLP&);
    //@}

    // CARL: I think we don't need this
    //Number nlp_lower_bound_inf_;
    //Number nlp_upper_bound_inf_;

  };

} // namespace Ipopt

#endif
