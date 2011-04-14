// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-11

#ifndef __SENSAMPLTNLP_HPP__
#define __SENSAMPLTNLP_HPP__

#include "AmplTNLP.hpp"

namespace Ipopt 
{

  DECLARE_STD_EXCEPTION(SUFFIX_EMPTY);

  class SensAmplTNLP : public AmplTNLP
  {
    /** This class is the Sens-wrapper for the ampltnlp, adapts 
     *  the get bounds function and some others to our needs. */

  public:
    /** constructor */
    SensAmplTNLP(const SmartPtr<const Journalist>& jnlst,
		 const SmartPtr<OptionsList> options,
		 char**& argv,
		 SmartPtr<AmplSuffixHandler> suffix_handler= NULL,
		 bool allow_discrete = false,
		 SmartPtr<AmplOptionsList> ampl_options_list= NULL ,
		 const char* ampl_option_string = NULL ,
		 const char* ampl_invokation_string  = NULL,
		 const char* ampl_banner_string  = NULL,
		 std::string* nl_file_content  = NULL);

    virtual ~SensAmplTNLP();

    /** returns bounds of the nlp. Overloaded from AmplTNLP */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);

    void set_sens_solution(Index idx, SmartPtr<IteratesVector> sens_sol);

    virtual void finalize_metadata(Index n,
				   const StringMetaDataMapType& var_string_md,
				   const IntegerMetaDataMapType& var_integer_md,
				   const NumericMetaDataMapType& var_numeric_md,
				   Index m,
				   const StringMetaDataMapType& con_string_md,
				   const IntegerMetaDataMapType& con_integer_md,
				   const NumericMetaDataMapType& con_numeric_md);

    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq);    

    const Index* get_index_suffix(const char* suffix_name);

    std::vector<Index> get_index_suffix_vec(const char* suffix_name);

    const Number* get_number_suffix(const char* suffix_name);

    std::vector<Number> get_number_suffix_vec(const char* suffix_name);

    const Index* get_index_suffix_constr(const char* suffix_name);

    std::vector<Index> get_index_suffix_constr_vec(const char* suffix_name);

    const Number* get_number_suffix_constr(const char* suffix_name);

    virtual bool get_var_con_metadata(Index n,
				      StringMetaDataMapType& var_string_md,
				      IntegerMetaDataMapType& var_integer_md,
				      NumericMetaDataMapType& var_numeric_md,
				      Index m,
				      StringMetaDataMapType& con_string_md,
				      IntegerMetaDataMapType& con_integer_md,
				      NumericMetaDataMapType& con_numeric_md);

  private:
    
    /** local copy of current lower and upper bounds - needed for parameter change */
    //    Number* x_L;
    //Number* x_U;

    SmartPtr<const Journalist> jnlst_;
    SmartPtr<OptionsList> options_;

    bool have_parameters_;

    Index* parameter_flags_;
    Number* parameter_values_;

    /** important Options */
    Index n_sens_steps_;
    bool run_sens_;
    bool compute_red_hessian_;

    std::vector< SmartPtr<IteratesVector> > sens_sol_;
    
  };

}

#endif
