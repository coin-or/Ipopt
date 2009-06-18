// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-06-17

#ifndef __IPPARTNLPWRAPPER_HPP__
#define __IPPARTNLPWRAPPER_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
#include "IpException.hpp"
#include "IpAlgTypes.hpp"
#include "IpReturnCodes.hpp"
#include "IpTNLP.hpp"
#include "IpParTNLP.hpp"

namespace Ipopt
{
  // forward declarations
  class IpoptData;
  class IpoptCalculatedQuantities;
  class IteratesVector;

  /** wrapper which allows creation of an ParTNLP object from a NLP object*/
  class ParTNLPWrapper : ParTNLP
  {
  public:

    /**@name Constructors/Destructors */
    //@{
    ParTNLPWrapper(SmartPtr<TNLP> tnlpobj)
    { tnlpobj_ = tnlpobj; }

    /** Default destructor */
    ~ParTNLPWrapper()
    {}
    //@}

    DECLARE_STD_EXCEPTION(INVALID_PARParTNLP);

    /**@name methods to gather information about the NLP */
    virtual bool get_nlp_info(Index num_proc, Index proc_id,
			      Index& n, Index& n_first, Index& n_last,
			      Index& m, Index& m_first, Index& m_last,
			      Index& nnz_jac_g_part, Index& nnz_h_lag_part,
			      TNLP::IndexStyleEnum& index_style);

    virtual bool get_bounds_info(Index num_proc, Index proc_id,
				 Index n, Index n_first, Index n_last,
				 Number* x_l_part, Number* x_u_part,
				 Index m, Index m_first, Index m_last,
                                 Number* g_l_part, Number* g_u_part);

    virtual bool
    get_starting_point(Index num_proc, Index proc_id,
		       Index n, Index n_first, Index n_last,
		       bool init_x, Number* x_part,
		       bool init_z, Number* z_L_part, Number* z_U_part,
		       Index m, Index m_first, Index m_last,
		       bool init_lambda, Number* lambda_part);

    virtual bool eval_f(Index num_proc, Index proc_id,
			Index n, Index n_first, Index n_last,
			const Number* x, bool new_x, Number& obj_value);

    virtual bool eval_grad_f(Index num_proc, Index proc_id,
			     Index n,  Index n_first, Index n_last,
			     const Number* x, bool new_x,
			     Number* grad_f_part);

    virtual bool eval_g(Index num_proc, Index proc_id,
			Index n, const Number* x, bool new_x,
                        Index m, Index m_first, Index m_last,
			Number* g_part);

    virtual bool eval_jac_g(Index num_proc, Index proc_id,
			    Index n, const Number* x, bool new_x,
                            Index m, Index m_first, Index m_last,
			    Index nele_jac_part, Index* iRow_part,
			    Index *jCol_part, Number* values_part);


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
    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq);

    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */

    /** Copy Constructor */
    ParTNLPWrapper(const ParTNLPWrapper&);

    /** Overloaded Equals Operator */
    void operator=(const ParTNLPWrapper&);
    //@}

    SmartPtr<TNLP> tnlpobj_;
  };

} // namespace Ipopt

#endif
