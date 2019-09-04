// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-11


#include "SensAmplTNLP.hpp"
#include "SensUtils.hpp"
#include "IpDenseVector.hpp"
#include "IpIteratesVector.hpp"
#include "IpBlas.hpp"
#include "IpIpoptData.hpp"

/* AMPL includes */
#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  SensAmplTNLP::SensAmplTNLP(const SmartPtr<const Journalist>& jnlst,
			     const SmartPtr<OptionsList> options,
			     char**& argv,
			     SmartPtr<AmplSuffixHandler> suffix_handler /* = NULL */,
			     bool allow_discrete /* = false */,
			     SmartPtr<AmplOptionsList> ampl_options_list /* = NULL */,
			     const char* ampl_option_string /* = NULL */,
			     const char* ampl_invokation_string /* = NULL */,
			     const char* ampl_banner_string /* = NULL */,
			     std::string* nl_file_content /* = NULL */)
    :
    AmplTNLP(jnlst, // gotta call constructor of base class properly
	     options,
	     argv,
	     suffix_handler /* = NULL */,
	     allow_discrete /* = false */,
	     ampl_options_list /* = NULL */,
	     ampl_option_string /* = NULL */,
	     ampl_invokation_string /* = NULL */,
	     ampl_banner_string /* = NULL */,
	     nl_file_content /* = NULL */),
    jnlst_(jnlst),
    options_(options),
    have_parameters_(false),
    parameter_flags_(NULL),
    parameter_values_(NULL)
  {
    DBG_START_METH("SensAmplTNLP::SensAmplTNLP", dbg_verbosity);

    SmartPtr<AmplSuffixHandler> suff_handler = get_suffix_handler();
    ASL_pfgh* asl = AmplSolverObject();

    const Index* parameter_flags = suff_handler->GetIntegerSuffixValues("parameter", AmplSuffixHandler::Variable_Source);

    if (parameter_flags_!=NULL) {
      have_parameters_ = true;
      for (Index i=0; i<n_var; ++i) {
	parameter_flags_[i] = parameter_flags[i];
      }
      parameter_values_ = new Number[n_var];
      const Number* nominal_values =
	suff_handler->GetNumberSuffixValues("nominal_value", AmplSuffixHandler::Variable_Source);
      if (nominal_values==NULL) {
	for (Index i=0; i<n_var; ++i) {
	  parameter_values_[i] = 0;
	}
      }
      else {
	for (Index i=0; i<n_var; ++i) {
	  parameter_values_[i] = nominal_values[i];
	}
      }
    }
    std::string prefix = "";
    options->GetIntegerValue("n_sens_steps",n_sens_steps_,prefix);
    sens_sol_.resize(n_sens_steps_, NULL);
    if ( n_sens_steps_==0 ) {
      options->SetStringValue("run_sens","no");
      run_sens_ = false;
    }

    options->GetBoolValue("run_sens", run_sens_, "");
    options->GetBoolValue("compute_red_hessian", compute_red_hessian_, "");
  }

  SensAmplTNLP::~SensAmplTNLP()
  {
    DBG_START_METH("SensAmplTNLP::~SensAmplTNLP", dbg_verbosity);

    delete[] parameter_values_;
  }


  bool SensAmplTNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
				     Index m, Number* g_l, Number* g_u)
  {
    DBG_START_METH("SensAmplTNLP::get_bounds_info", dbg_verbosity);

    ASL_pfgh* asl = AmplSolverObject();

    DBG_ASSERT(n == n_var);
    DBG_ASSERT(m == n_con);

    for (Index i=0; i<n; i++) {
      x_l[i] = LUv[2*i];
      x_u[i] = LUv[2*i+1];
    }

    for (Index i=0; i<m; i++) {
      g_l[i] = LUrhs[2*i];
      g_u[i] = LUrhs[2*i+1];
    }

    /*  Sensitivity: Fix parameters by bounds */
    if (have_parameters_) {
      // parameters are set in the model. Set Bounds to current parameters
      for (Index i=0; i<n; i++) {
	if (parameter_flags_[i]) {
	  x_l[i] = parameter_values_[i];
	  x_u[i] = parameter_values_[i];
	}
      }
    }
    return true;
  }

  void SensAmplTNLP::set_sens_solution(Index idx, SmartPtr<IteratesVector> sens_sol)
  {
    DBG_START_METH("SensAmplTNLP::set_sens_solution", dbg_verbosity);

    DBG_PRINT((dbg_verbosity, "n_sens_steps=%d\n", sens_sol_.size()));
    DBG_ASSERT(idx>0);
    DBG_ASSERT(idx<=(Index)sens_sol_.size());

    sens_sol_[idx-1] = sens_sol;
  }

  void SensAmplTNLP::finalize_metadata(Index n,
				       const StringMetaDataMapType& var_string_md,
				       const IntegerMetaDataMapType& var_integer_md,
				       const NumericMetaDataMapType& var_numeric_md,
				       Index m,
				       const StringMetaDataMapType& con_string_md,
				       const IntegerMetaDataMapType& con_integer_md,
				       const NumericMetaDataMapType& con_numeric_md)
  {
    DBG_START_METH("SensAmplTNLP::finalize_metadata", dbg_verbosity);
    ASL_pfgh* asl = AmplSolverObject();

    if (run_sens_) {
      for (Index step=1; step<=n_sens_steps_; ++step) {
	std::string sol_state_id = "sens_sol_state_";
	append_Index(sol_state_id, step);
	NumericMetaDataMapType::const_iterator num_it;
	num_it = var_numeric_md.find(sol_state_id);
	if (num_it!=var_numeric_md.end()) {
	  suf_rput(sol_state_id.c_str(), ASL_Sufkind_var, const_cast<Number*>(&num_it->second[0]));
	}
	std::string sol_state_z_L_id = sol_state_id + "_z_L";
	num_it = var_numeric_md.find(sol_state_z_L_id);
	if (num_it!=var_numeric_md.end()) {
	  suf_rput(sol_state_z_L_id.c_str(), ASL_Sufkind_var, const_cast<Number*>(&num_it->second[0]));
	}
	std::string sol_state_z_U_id = sol_state_id + "_z_U";
	num_it = var_numeric_md.find(sol_state_z_U_id);
	if (num_it!=var_numeric_md.end()) {
	  suf_rput(sol_state_z_U_id.c_str(), ASL_Sufkind_var, const_cast<Number*>(&num_it->second[0]));
	}
	num_it = con_numeric_md.find(sol_state_id);
	if (num_it!=con_numeric_md.end()) {
	  suf_rput(sol_state_id.c_str(), ASL_Sufkind_con, const_cast<Number*>(&num_it->second[0]));
	}
      }
    }
  }


  void SensAmplTNLP::finalize_solution(SolverReturn status,
				       Index n, const Number* x, const Number* z_L, const Number* z_U,
				       Index m, const Number* g, const Number* lambda,
				       Number obj_value,
				       const IpoptData* ip_data,
				       IpoptCalculatedQuantities* ip_cq)
  {
    DBG_START_METH("SensAmplTNLP::finalize_solution", dbg_verbosity);


    AmplTNLP::finalize_solution(status,
				n, x, z_L,  z_U,
				m, g, lambda,
				obj_value,
				ip_data,
				ip_cq);
  }

  bool SensAmplTNLP::get_var_con_metadata(Index n,
					  StringMetaDataMapType& var_string_md,
					  IntegerMetaDataMapType& var_integer_md,
					  NumericMetaDataMapType& var_numeric_md,
					  Index m,
					  StringMetaDataMapType& con_string_md,
					  IntegerMetaDataMapType& con_integer_md,
					  NumericMetaDataMapType& con_numeric_md)
  {
    DBG_START_METH("SensAmplTNLP::get_var_con_metadata", dbg_verbosity);

    try {
      if (run_sens_) {
	// Get Sens Suffixes
	std::string sens_state = "sens_state_";
	std::vector<Index> state;
	for (Index i=1; i<=n_sens_steps_; ++i) {
	  append_Index(sens_state,i);
	  state = get_index_suffix_vec(sens_state.c_str());
	  set_integer_metadata_for_var(sens_state, state);
	  sens_state = "sens_state_";
	}
	std::string sens_state_value = "sens_state_value_";
	std::vector<Number> state_val;
	for (Index i=1; i<=n_sens_steps_; ++i) {
	  append_Index(sens_state_value,i);
	  state_val = get_number_suffix_vec(sens_state_value.c_str());
	  set_numeric_metadata_for_var(sens_state_value, state_val);
	  sens_state_value = "sens_state_value_";
	}
	std::string init_constr = "sens_init_constr";
	if (n_sens_steps_ >0) {
	  std::vector<Index> init_idx = get_index_suffix_constr_vec(init_constr.c_str());
	  set_integer_metadata_for_con(init_constr,init_idx);
	}
      }
    }
    catch ( SUFFIX_EMPTY& exc ) {
      //exc.ReportException(*jnlst_);
      // const std::string exc_mess = exc.Message();
      const std::string exc_mess = exc.Message();
      jnlst_->Printf(J_WARNING, J_INITIALIZATION,
		     "    WARNING: Will not run sIPOPT "
		     "because of incorrect AMPL suffix!\n"
		     "      Message: %s\n\n", exc_mess.c_str() );
      options_->SetStringValue("sens_internal_abort", "yes");
      bool ignore_suffix_error;
      options_->GetBoolValue("ignore_suffix_error", ignore_suffix_error, "");
      if ( !ignore_suffix_error ) {
	THROW_EXCEPTION(SUFFIX_EMPTY,
			"Encountered Suffix Error");
      }
    }

    try {
      if (compute_red_hessian_) {
	std::string red_hess_str = "red_hessian";
	std::vector<Index> red_hess_idx = get_index_suffix_vec(red_hess_str.c_str());
	set_integer_metadata_for_var(red_hess_str.c_str(), red_hess_idx);
      }
    }
    catch ( SUFFIX_EMPTY& exc ) {
      const std::string exc_mess = exc.Message();
      jnlst_->Printf(J_WARNING, J_INITIALIZATION,
		     "    WARNING: Will not run reduced hessian computation "
		     "because of incorrect AMPL suffix!\n"
		     "      Message: %s\n\n", exc_mess.c_str() );
      options_->SetStringValue("sens_internal_abort", "yes");
      bool ignore_suffix_error;
      options_->GetBoolValue("ignore_suffix_error", ignore_suffix_error, "");
      if ( !ignore_suffix_error ) {
	THROW_EXCEPTION(SUFFIX_EMPTY,
			"Encountered Suffix Error");
      }
    }

    bool retval = AmplTNLP::get_var_con_metadata(n,
						 var_string_md,
						 var_integer_md,
						 var_numeric_md,
						 m,
						 con_string_md,
						 con_integer_md,
						 con_numeric_md);
    return retval;
  }


  const Index* SensAmplTNLP::get_index_suffix(const char* suffix_name)
  {
    DBG_START_METH("SensAmplTNLP::get_index_suffix", dbg_verbosity);

    SmartPtr<AmplSuffixHandler> suffix_handler = get_suffix_handler();

    const Index* index_suffix=
      suffix_handler->GetIntegerSuffixValues(suffix_name, AmplSuffixHandler::Variable_Source);

    return index_suffix;
  }

  std::vector<Index> SensAmplTNLP::get_index_suffix_vec(const char* suffix_name)
  {
    DBG_START_METH("SensAmplTNLP::get_index_suffix_vec", dbg_verbosity);

    ASL_pfgh* asl = AmplSolverObject();
    SmartPtr<AmplSuffixHandler> suffix_handler = get_suffix_handler();
    DBG_ASSERT(IsValid(suffix_handler));

    std::vector<Index> index_suffix=
      suffix_handler->GetIntegerSuffixValues(n_var, suffix_name, AmplSuffixHandler::Variable_Source);
    if ( index_suffix.size()==0 ) {
      index_suffix.resize(n_var, 0);
    }
    return index_suffix;
  }

  const Number* SensAmplTNLP::get_number_suffix(const char* suffix_name)
  {
    DBG_START_METH("SensAmplTNLP::get_number_suffix", dbg_verbosity);

    SmartPtr<AmplSuffixHandler> suffix_handler = get_suffix_handler();

    const Number* number_suffix=
      suffix_handler->GetNumberSuffixValues(suffix_name, AmplSuffixHandler::Variable_Source);

    if (number_suffix==NULL) {
      // suffix invalid
      std::string except = suffix_name;
      except.append(" is empty");
      THROW_EXCEPTION(SUFFIX_EMPTY, except);
    }

    return number_suffix;
  }

  std::vector<Number> SensAmplTNLP::get_number_suffix_vec(const char* suffix_name)
  {
    DBG_START_METH("SensAmplTNLP::get_number_suffix_vec", dbg_verbosity);
    ASL_pfgh* asl = AmplSolverObject();
    SmartPtr<AmplSuffixHandler> suffix_handler = get_suffix_handler();

    std::vector<Number> number_suffix=
      suffix_handler->GetNumberSuffixValues(n_var, suffix_name, AmplSuffixHandler::Variable_Source);

    if (number_suffix.empty()) {
      // suffix invalid
      std::string except = suffix_name;
      except.append(" is empty");
      THROW_EXCEPTION(SUFFIX_EMPTY, except);
    }

    return number_suffix;
  }

  const Index* SensAmplTNLP::get_index_suffix_constr(const char* suffix_name)
  {
    DBG_START_METH("SensAmplTNLP::get_index_suffix_constr", dbg_verbosity);

    SmartPtr<AmplSuffixHandler> suffix_handler = get_suffix_handler();

    const Index* index_suffix=
      suffix_handler->GetIntegerSuffixValues(suffix_name, AmplSuffixHandler::Constraint_Source);

    if (index_suffix==NULL) {
      // suffix invalid
      std::string except = suffix_name;
      except.append(" is empty");
      THROW_EXCEPTION(SUFFIX_EMPTY, except);
    }

    return index_suffix;
  }

  std::vector<Index> SensAmplTNLP::get_index_suffix_constr_vec(const char* suffix_name)
  {
    DBG_START_METH("SensAmplTNLP::get_index_suffix_constr_vec", dbg_verbosity);
    ASL_pfgh* asl = AmplSolverObject();
    SmartPtr<AmplSuffixHandler> suffix_handler = get_suffix_handler();

    std::vector<Index> index_suffix=
      suffix_handler->GetIntegerSuffixValues(n_con, suffix_name, AmplSuffixHandler::Constraint_Source);

    if (index_suffix.empty()) {
      // suffix invalid
      std::string except = suffix_name;
      except.append(" is empty");
      THROW_EXCEPTION(SUFFIX_EMPTY, except);
    }

    return index_suffix;
  }

  const Number* SensAmplTNLP::get_number_suffix_constr(const char* suffix_name)
  {
    DBG_START_METH("SensAmplTNLP::get_number_suffix_constr", dbg_verbosity);

    SmartPtr<AmplSuffixHandler> suffix_handler = get_suffix_handler();

    const Number* number_suffix=
      suffix_handler->GetNumberSuffixValues(suffix_name, AmplSuffixHandler::Constraint_Source);

    if (number_suffix==NULL) {
      // suffix invalid
      std::string except = suffix_name;
      except.append(" is empty");
      THROW_EXCEPTION(SUFFIX_EMPTY, except);
    }

    return number_suffix;
  }


}
