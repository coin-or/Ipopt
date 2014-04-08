// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"

#include "AmplTNLP.hpp"
#include "IpDenseVector.hpp"
#include "IpGenTMatrix.hpp"
#include "IpSymTMatrix.hpp"
#include "IpBlas.hpp"

#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif

/* AMPL includes */
#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  AmplTNLP::AmplTNLP(const SmartPtr<const Journalist>& jnlst,
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
      TNLP(),
      jnlst_(jnlst),
      asl_(NULL),
      obj_sign_(1),
      nz_h_full_(-1),
      x_sol_(NULL),
      z_L_sol_(NULL),
      z_U_sol_(NULL),
      g_sol_(NULL),
      lambda_sol_(NULL),
      obj_sol_(0.0),
      objval_called_with_current_x_(false),
      conval_called_with_current_x_(false),
      hesset_called_(false),
      set_active_objective_called_(false),
      Oinfo_ptr_(NULL),
      suffix_handler_(suffix_handler)
  {
    DBG_START_METH("AmplTNLP::AmplTNLP",
                   dbg_verbosity);

    // The ASL include files #define certain
    // variables that they expect you to work with.
    // These variables then appear as though they are
    // global variables when, in fact, they are not
    // Most of them are data members of an asl object

    // Create the ASL structure
    ASL_pfgh* asl = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh);
    DBG_ASSERT(asl);
    asl_ = asl; // keep the pointer for ourselves to use later...

    // First assume that we don't want to halt on error (default)
    fint* fint_nerror = new fint;
    *fint_nerror = 0;
    nerror_ = (void*) fint_nerror;

    // Read the options and stub
    char* stub = get_options(options, ampl_options_list,
                             ampl_option_string, ampl_invokation_string,
                             ampl_banner_string, argv);
    FILE*nl = NULL;
    if (nl_file_content) {
      nl = jac0dim(const_cast<char*>(nl_file_content->c_str()),
                   -(ftnlen)nl_file_content->length());
    }
    else {
      if (!stub) {
        jnlst_->Printf(J_ERROR, J_MAIN, "No .nl file given!\n");
        THROW_EXCEPTION(INVALID_TNLP, "No .nl file given!\n");
      }
      nl = jac0dim(stub, (fint)strlen(stub));
      DBG_ASSERT(nl);
    }
    jnlst_->Printf(J_SUMMARY, J_MAIN, "\n");

    // check the problem statistics (see Table 1 in AMPL doc)
    DBG_ASSERT(n_var > 0); // need some continuous variables
    if (!allow_discrete && (nbv>0 || niv>0 || nlvbi>0 || nlvci>0 || nlvoi>0) ) {
      jnlst_->Printf(J_WARNING, J_MAIN, "==> Warning: Treating %d binary and %d integer variables as continous.\n\n", nbv, niv+nlvbi+nlvci+nlvoi);
      allow_discrete = true;
    }
    allow_discrete = true;
    ASSERT_EXCEPTION(allow_discrete || (nbv == 0 && niv == 0 && nlvbi == 0 && nlvci == 0 && nlvoi == 0),
                     IpoptException,
                     "Discrete variables not allowed when the allow_discrete flag is false, "
                     "Either remove the integer variables, or change the flag in the constructor of AmplTNLP"
                    );

    if (n_cc != 0) {
      jnlst_->Printf(J_ERROR, J_MAIN, "\n\n***** Currently, the AMPL interface to Ipopt does not support\n      complementarity constraints.  Aborting. *****\n\n");
      ASSERT_EXCEPTION(n_cc == 0, IpoptException,
		       "Currently, the AMPL interface to Ipopt does not support complementarity constraints.");
    }

    DBG_ASSERT(nlo == 0 || nlo == 1); // Can handle nonlinear obj.
    DBG_ASSERT(nwv == 0); // Don't know what "linear arc" variables are
    DBG_ASSERT(nlnc == 0); // Don't know what "nonlinear network"constraints are
    DBG_ASSERT(lnc == 0); // Don't know what "linear network" constraints are

    // Set options in the asl structure
    want_xpi0 = 1 | 2;  // allocate initial values for primal and dual if available
    obj_no = 0;
    DBG_ASSERT((want_xpi0 & 1) == 1 && (want_xpi0 & 2) == 2);

    // allocate space for initial values
    X0 = new real[n_var];
    havex0 = new char[n_var];
    pi0 = new real[n_con];
    havepi0 = new char[n_con];

    // prepare for suffixes
    if (IsValid(suffix_handler)) {
      suffix_handler->PrepareAmplForSuffixes(asl_);
    }

    // read the rest of the nl file
    int retcode = pfgh_read(nl, ASL_return_read_err | ASL_findgroups);

    switch (retcode) {
    case ASL_readerr_none : {}
      break;
    case ASL_readerr_nofile : {
        jnlst_->Printf(J_ERROR, J_MAIN, "Cannot open .nl file\n");
        THROW_EXCEPTION(INVALID_TNLP, "Cannot open .nl file");
      }
      break;
    case ASL_readerr_nonlin : {
        DBG_ASSERT(false); // this better not be an error!
        jnlst_->Printf(J_ERROR, J_MAIN, "model involves nonlinearities (ed0read)\n");
        THROW_EXCEPTION(INVALID_TNLP, "model involves nonlinearities (ed0read)");
      }
      break;
    case  ASL_readerr_argerr : {
        jnlst_->Printf(J_ERROR, J_MAIN, "user-defined function with bad args\n");
        THROW_EXCEPTION(INVALID_TNLP, "user-defined function with bad args");
      }
      break;
    case ASL_readerr_unavail : {
        jnlst_->Printf(J_ERROR, J_MAIN, "user-defined function not available\n");
        THROW_EXCEPTION(INVALID_TNLP, "user-defined function not available");
      }
      break;
    case ASL_readerr_corrupt : {
        jnlst_->Printf(J_ERROR, J_MAIN, "corrupt .nl file\n");
        THROW_EXCEPTION(INVALID_TNLP, "corrupt .nl file");
      }
      break;
    case ASL_readerr_bug : {
        jnlst_->Printf(J_ERROR, J_MAIN, "bug in .nl reader\n");
        THROW_EXCEPTION(INVALID_TNLP, "bug in .nl reader");
      }
      break;
    case ASL_readerr_CLP : {
        jnlst_->Printf(J_ERROR, J_MAIN, "Ampl model contains a constraint without \"=\", \">=\", or \"<=\".\n");
        THROW_EXCEPTION(INVALID_TNLP, "Ampl model contains a constraint without \"=\", \">=\", or \"<=\".");
      }
      break;
    default: {
        jnlst_->Printf(J_ERROR, J_MAIN, "Unknown error in stub file read. retcode = %d\n", retcode);
        THROW_EXCEPTION(INVALID_TNLP, "Unknown error in stub file read");
      }
      break;
    }
  }

  void AmplTNLP::set_active_objective(Index in_obj_no)
  {
    if (hesset_called_) {
      jnlst_->Printf(J_ERROR, J_MAIN, "Internal error: AmplTNLP::set_active_objective called after AmplTNLP::call_hesset.\n");
      THROW_EXCEPTION(INVALID_TNLP, "Internal error: AmplTNLP::set_active_objective called after AmplTNLP::call_hesset.");
    }
    ASL_pfgh* asl = asl_;
    obj_no = in_obj_no;
    set_active_objective_called_ = true;
  }

  void AmplTNLP::call_hesset()
  {
    if (hesset_called_) {
      jnlst_->Printf(J_ERROR, J_MAIN, "Internal error: AmplTNLP::call_hesset is called twice.\n");
      THROW_EXCEPTION(INVALID_TNLP, "Internal error: AmplTNLP::call_hesset is called twice.");
    }

    ASL_pfgh* asl = asl_;

    if (n_obj == 0) {
      hesset(1,0,0,0,nlc);
    }
    else {
      if (n_obj>1 && !set_active_objective_called_) {
        jnlst_->Printf(J_ERROR, J_MAIN,
                       "There is more than one objective function in the AMPL model, but AmplTNLP::set_active_objective has not been called.\n");
        THROW_EXCEPTION(INVALID_TNLP, "There is more than one objective function in the AMPL model, but AmplTNLP::set_active_objective has not been called");
      }
      // see "changes" in solvers directory of ampl code...
      hesset(1,obj_no,1,0,nlc);
    }

    obj_sign_ = 1; // minimization
    if (objtype[obj_no] != 0) {
      obj_sign_ = -1;
    }

    // find the nonzero structure for the hessian parameters to
    // sphsetup:
    int coeff_obj = 1;
    int mult_supplied = 1; // multipliers will be supplied
    int uptri = 1; // only need the upper triangular part
    nz_h_full_ = sphsetup(-1, coeff_obj, mult_supplied, uptri);

    hesset_called_ = true;
  }

  AmplTNLP::~AmplTNLP()
  {
    ASL_pfgh* asl = asl_;

    if (asl) {
      if (X0) {
        delete [] X0;
        X0 = NULL;
      }
      if (havex0) {
        delete [] havex0;
        havex0 = NULL;
      }
      if (pi0) {
        delete [] pi0;
        pi0 = NULL;
      }
      if (havepi0) {
        delete [] havepi0;
        havepi0 = NULL;
      }
      ASL* asl_to_free = (ASL*)asl_;
      ASL_free(&asl_to_free);
      asl_ = NULL;
    }

    delete [] x_sol_;
    x_sol_ = NULL;
    delete [] z_L_sol_;
    z_L_sol_ = NULL;
    delete [] z_U_sol_;
    z_U_sol_ = NULL;
    delete [] g_sol_;
    g_sol_ = NULL;
    delete [] lambda_sol_;
    lambda_sol_ = NULL;

    if (Oinfo_ptr_) {
      Option_Info* Oinfo = (Option_Info*) Oinfo_ptr_;
      delete [] Oinfo->sname;
      delete [] Oinfo->bsname;
      delete [] Oinfo->opname;
      delete Oinfo;
    }

    delete (fint*) nerror_;
  }

  bool AmplTNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, IndexStyleEnum& index_style)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);

    if (!hesset_called_) {
      call_hesset();
    }

    n = n_var; // # of variables (variable types have been asserted in the constructor
    m = n_con; // # of constraints
    nnz_jac_g = nzc; // # of non-zeros in the jacobian
    nnz_h_lag = nz_h_full_; // # of non-zeros in the hessian

    index_style = TNLP::FORTRAN_STYLE;

    return true;
  }

  bool AmplTNLP::get_var_con_metadata(Index n,
                                      StringMetaDataMapType& var_string_md,
                                      IntegerMetaDataMapType& var_integer_md,
                                      NumericMetaDataMapType& var_numeric_md,
                                      Index m,
                                      StringMetaDataMapType& con_string_md,
                                      IntegerMetaDataMapType& con_integer_md,
                                      NumericMetaDataMapType& con_numeric_md)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);

    // pick up the variable and constraints names if available
    Index rlen = maxrownamelen;
    Index clen = maxcolnamelen;

    if (rlen > 0) {
      std::vector<std::string> var_names(n);
      for (Index i=0; i<n; i++) {
        var_names[i] = var_name(i);
      }
      var_string_md_["idx_names"] = var_names;
    }

    if (clen > 0) {
      std::vector<std::string> con_names(m);
      for (Index i=0; i<m; i++) {
        con_names[i] = con_name(i);
      }
      con_string_md_["idx_names"] = con_names;
    }

    if (var_string_md_.size() > 0 || var_integer_md_.size() > 0 || var_numeric_md_.size() > 0
        || con_string_md_.size() > 0 || con_integer_md_.size() > 0 || con_numeric_md_.size() > 0) {
      var_string_md = var_string_md_;
      var_integer_md = var_integer_md_;
      var_numeric_md = var_numeric_md_;
      con_string_md = con_string_md_;
      con_integer_md = con_integer_md_;
      con_numeric_md = con_numeric_md_;
      return true;
    }

    return false;
  }

  bool AmplTNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);

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

    return true;
  }


  bool AmplTNLP::get_constraints_linearity(Index n,
      LinearityType* const_types)
  {
    ASL_pfgh* asl = AmplSolverObject();
    //check that n is good
    DBG_ASSERT(n == n_con);
    // check that there are no network constraints
    DBG_ASSERT(nlnc == 0 && lnc == 0);
    //the first nlc constraints are non linear the rest is linear
    for (Index i=0; i<nlc; i++) {
      const_types[i]=NON_LINEAR;
    }
    // the rest is linear
    for (Index i=nlc; i<n_con; i++)
      const_types[i]=LINEAR;
    return true;
  }

  bool AmplTNLP::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    DBG_ASSERT(n == n_var);
    DBG_ASSERT(m == n_con);

    if (init_x)
      for (Index i=0; i<n; i++)
        x[i] = havex0[i] ? X0[i] : Max(LUv[2*i], Min(LUv[2*i+1], 0.0));

    if (init_z) {
      // Modified for warm-start from AMPL
      DBG_ASSERT(IsValid(suffix_handler_));
      const double* zL_init = suffix_handler_->GetNumberSuffixValues("ipopt_zL_in", AmplSuffixHandler::Variable_Source);
      const double* zU_init = suffix_handler_->GetNumberSuffixValues("ipopt_zU_in", AmplSuffixHandler::Variable_Source);
      for (Index i=0; i<n; i++) {
        z_L[i] = zL_init != NULL ?  obj_sign_ * zL_init[i] : 1.0;
        z_U[i] = zU_init != NULL ? -obj_sign_ * zU_init[i] : 1.0;
      }
    }

    if (init_lambda)
      for (Index i=0; i<m; i++)
        lambda[i] = havepi0[i] ? -obj_sign_ * pi0[i] : 0.0;

    return true;
  }

  bool AmplTNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
  {
    DBG_START_METH("AmplTNLP::eval_f",
                   dbg_verbosity);
    if (!apply_new_x(new_x, n, x)) {
      return false;
    }

    return internal_objval(x, obj_value);
  }

  bool AmplTNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
  {
    DBG_START_METH("AmplTNLP::eval_grad_f",
                   dbg_verbosity);
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);

    if (!apply_new_x(new_x, n, x)) {
      return false;
    }

    if (n_obj==0) {
      for (Index i=0; i<n; i++) {
        grad_f[i] = 0.;
      }
    }
    else {
      objgrd(obj_no, const_cast<Number*>(x), grad_f, (fint*)nerror_);
      if (!nerror_ok(nerror_)) {
        return false;
      }

      if (obj_sign_==-1) {
        for (Index i=0; i<n; i++) {
          grad_f[i] *= -1.;
        }
      }
    }
    return true;
  }

  bool AmplTNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
  {
    DBG_START_METH("AmplTNLP::eval_g", dbg_verbosity);

    DBG_DO(ASL_pfgh* asl = asl_);
    DBG_ASSERT(n == n_var);
    DBG_ASSERT(m == n_con);

    if (!apply_new_x(new_x, n, x)) {
      return false;
    }

    return internal_conval(x, m, g);
  }

  bool AmplTNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow,
                            Index *jCol, Number* values)
  {
    DBG_START_METH("AmplTNLP::eval_jac_g",
                   dbg_verbosity);
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    DBG_ASSERT(n == n_var);
    DBG_ASSERT(m == n_con);

    if (iRow && jCol && !values) {
      // setup the structure
      Index current_nz = 0;
      for (Index i=0; i<n_con; i++) {
        for (cgrad* cg=Cgrad[i]; cg; cg = cg->next) {
          iRow[cg->goff] = i + 1;
          jCol[cg->goff] = cg->varno + 1;
          //    iRow[current_nz] = i + 1;
          //    jCol[current_nz] = cg->varno+1;
          current_nz++;
        }
      }
      DBG_ASSERT(current_nz == nele_jac);
      return true;
    }
    else if (!iRow && !jCol && values) {
      if (!apply_new_x(new_x, n, x)) {
        return false;
      }

      jacval(const_cast<Number*>(x), values, (fint*)nerror_);
      if (nerror_ok(nerror_)) {
        return true;
      }
    }
    else {
      DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
    }

    return false;
  }

  bool AmplTNLP::eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values)
  {
    DBG_START_METH("AmplTNLP::eval_h",
                   dbg_verbosity);
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    DBG_ASSERT(n == n_var);
    DBG_ASSERT(m == n_con);

    if (iRow && jCol && !values) {
      // setup the structure
      int k=0;
      for (int i=0; i<n; i++) {
        for (int j=sputinfo->hcolstarts[i]; j<sputinfo->hcolstarts[i+1]; j++) {
          iRow[k] = i + 1;
          jCol[k] = sputinfo->hrownos[j]+1;
          k++;
        }
      }
      DBG_ASSERT(k==nele_hess);
      return true;
    }
    else if (!iRow && !jCol && values) {
      if (!apply_new_x(new_x, n, x)) {
        return false;
      }
      if (!objval_called_with_current_x_) {
        Number dummy;
        internal_objval(x, dummy);
        internal_conval(x, m);
      }
      if (!conval_called_with_current_x_) {
        internal_conval(x, m);
      }

      real* OW = new real[Max(1,n_obj)];
      if (n_obj>0) {
        for (Index i=0; i<n_obj; i++) {
          OW[i] = 0.;
        }
        OW[obj_no] = obj_sign_*obj_factor;
      }
      sphes(values, -1, OW, const_cast<Number*>(lambda));
      delete [] OW;
      return true;
    }
    else {
      DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
    }

    return false;
  }

  void AmplTNLP::finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq)
  {
    ASL_pfgh* asl = asl_;

    if (!x_sol_) {
      x_sol_ = new Number[n];
    }
    if (!z_L_sol_) {
      z_L_sol_ = new Number[n];
    }
    if (!z_U_sol_) {
      z_U_sol_ = new Number[n];
    }
    if (!g_sol_) {
      g_sol_ = new Number[m];
    }
    if (!lambda_sol_) {
      lambda_sol_ = new Number[m];
    }

    IpBlasDcopy(n, x, 1, x_sol_, 1);
    IpBlasDcopy(m, g, 1, g_sol_, 1);
    if( obj_sign_ == -1.0 ) // maximization
    {
      for( int i = 0; i < n; ++i )
        z_L_sol_[i] = -z_L[i];
      IpBlasDcopy(n, z_U, 1, z_U_sol_, 1);
      IpBlasDcopy(m, lambda, 1, lambda_sol_, 1);
    }
    else
    {
      IpBlasDcopy(n, z_L, 1, z_L_sol_, 1);
      for( int i = 0; i < n; ++i )
        z_U_sol_[i] = -z_U[i];
      for( int i = 0; i < m; ++i )
        lambda_sol_[i] = -lambda[i];
    }
    obj_sol_ = obj_value;

    std::string message;
    if (status == SUCCESS) {
      message = "Optimal Solution Found";
      solve_result_num = 0;
    }
    else if (status == MAXITER_EXCEEDED) {
      message = "Maximum Number of Iterations Exceeded.";
      solve_result_num = 400;
    }
    else if (status == CPUTIME_EXCEEDED) {
      message = "Maximum CPU Time Exceeded.";
      solve_result_num = 401;
    }
    else if (status == STOP_AT_TINY_STEP) {
      message = "Search Direction becomes Too Small.";
      solve_result_num = 500;
    }
    else if (status == STOP_AT_ACCEPTABLE_POINT) {
      message = "Solved To Acceptable Level.";
      solve_result_num = 1;
    }
    else if (status == FEASIBLE_POINT_FOUND) {
      message = "Found feasible point for square problem.";
      solve_result_num = 2;
    }
    else if (status == LOCAL_INFEASIBILITY) {
      message = "Converged to a locally infeasible point. Problem may be infeasible.";
      solve_result_num = 200;
    }
    else if (status == RESTORATION_FAILURE) {
      message = "Restoration Phase Failed.";
      solve_result_num = 501;
    }
    else if (status == DIVERGING_ITERATES) {
      message = "Iterates diverging; problem might be unbounded.";
      solve_result_num = 300;
    }
    else {
      message = "Unknown Error";
      solve_result_num = 502;
    }

    if (IsValid(suffix_handler_)) {
      // Modified for warm-start from AMPL. Assign Bound Multipliers as Suffixes
      suf_rput("ipopt_zL_out", ASL_Sufkind_var,  z_L_sol_);
      suf_rput("ipopt_zU_out", ASL_Sufkind_var,  z_U_sol_);
    }

    // Write the .sol file
    message = " \nIpopt " IPOPT_VERSION ": " + message;
    write_solution_file(message.c_str());
  }

  bool AmplTNLP::internal_objval(const Number* x, Number& obj_val)
  {
    DBG_START_METH("AmplTNLP::internal_objval",
                   dbg_verbosity);
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    objval_called_with_current_x_ = false; // in case the call below fails

    if (n_obj==0) {
      obj_val = 0;
      objval_called_with_current_x_ = true;
      return true;
    }
    else {
      Number retval = objval(obj_no, const_cast<Number*>(x), (fint*)nerror_);
      if (nerror_ok(nerror_)) {
        obj_val = obj_sign_*retval;
        objval_called_with_current_x_ = true;
        return true;
      }
    }

    return false;
  }

  bool AmplTNLP::internal_conval(const Number* x, Index m, Number* g)
  {
    DBG_START_METH("AmplTNLP::internal_conval",
                   dbg_verbosity);
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    DBG_ASSERT(m == n_con);
    conval_called_with_current_x_ = false; // in case the call below fails

    bool allocated = false;
    if (!g) {
      g = new double[m];
      allocated = true;
    }
    conval(const_cast<Number*>(x), g, (fint*)nerror_);

    if (allocated) {
      delete [] g;
      g = NULL;
    }

    if (nerror_ok(nerror_)) {
      conval_called_with_current_x_ = true;
      return true;
    }
    return false;
  }


  bool AmplTNLP::apply_new_x(bool new_x, Index n, const Number* x)
  {
    DBG_START_METH("AmplTNLP::apply_new_x",
                   dbg_verbosity);

    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);

    if (new_x) {
      if (!hesset_called_) {
        call_hesset();
      }

      DBG_PRINT((1, "Set new x.\n"));
      // update the flags so these methods are called
      // before evaluating the hessian
      conval_called_with_current_x_ = false;
      objval_called_with_current_x_ = false;

      // tell ampl that we have a new x
      xknowne(const_cast<Number*>(x), (fint*)nerror_);
      return nerror_ok(nerror_);
    }

    return true;
  }

  void AmplTNLP::write_solution_file(const std::string& message) const
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl);
    DBG_ASSERT(x_sol_ && lambda_sol_);

    // We need to copy the message into a non-const char array to make
    // it work with the AMPL C function.
    char* cmessage = new char[message.length()+1];
    strcpy(cmessage, message.c_str());

    write_sol(cmessage, x_sol_, lambda_sol_, (Option_Info*)Oinfo_ptr_);

    delete [] cmessage;
  }

  void AmplTNLP::get_discrete_info(Index& nlvb_,
                                   Index& nlvbi_,
                                   Index& nlvc_,
                                   Index& nlvci_,
                                   Index& nlvo_,
                                   Index& nlvoi_,
                                   Index& nbv_,
                                   Index& niv_) const
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl);

    nlvb_ = nlvb;
    nlvbi_ = nlvbi;
    nlvc_ = nlvc;
    nlvci_ = nlvci;
    nlvo_ = nlvo;
    nlvoi_ = nlvoi;
    nbv_ = nbv;
    niv_ = niv;
  }

  bool AmplTNLP::get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling)
  {
    DBG_ASSERT(IsValid(suffix_handler_));
    const double* obj = suffix_handler_->GetNumberSuffixValues("scaling_factor", AmplSuffixHandler::Objective_Source);
    obj_scaling = (obj) ? obj[0] : 1.0;

    const double* x = suffix_handler_->GetNumberSuffixValues("scaling_factor", AmplSuffixHandler::Variable_Source);
    if (x) {
      use_x_scaling = true;
      for (int i=0; i < n; i++) {
        if (x[i] > 0.0) {
          x_scaling[i] = x[i];
        }
        else {
          x_scaling[i] = 1.0;
        }
      }
    }
    else {
      use_x_scaling = false;
    }

    const double* g = suffix_handler_->GetNumberSuffixValues("scaling_factor", AmplSuffixHandler::Constraint_Source);
    if (g) {
      use_g_scaling = true;
      for (int i=0; i < m; i++) {
        if (g[i] > 0) {
          g_scaling[i] = g[i];
        }
        else {
          g_scaling[i] = 1.0;
        }
      }
    }
    else {
      use_g_scaling = false;
    }

    return true;
  }

  Index AmplTNLP::get_number_of_nonlinear_variables()
  {
    ASL_pfgh* asl = asl_;
    return Max(nlvo,nlvc);
  }

  bool AmplTNLP::get_list_of_nonlinear_variables(Index num_nonlin_vars,
      Index* pos_nonlin_vars)
  {
    DBG_DO(ASL_pfgh* asl = asl_;)
    DBG_ASSERT(num_nonlin_vars == Max(nlvo,nlvc));

    // The first variables are the nonlinear ones (using Fortran
    // numbering)
    for (Index i=0; i<num_nonlin_vars; i++) {
      pos_nonlin_vars[i] = i+1;
    }
    return true;
  }

  extern "C"
  {
    static char* get_num_opt(Option_Info *oi, keyword *kw, char *value) {
      AmplOptionsList::PrivatInfo*
      pinfo = (AmplOptionsList::PrivatInfo*) kw->info;

      real real_val;
      kw->info = &real_val;
      char* retval = D_val(oi, kw, value);
      kw->info = (void*) pinfo;

      if (!pinfo->Options()->SetNumericValue(pinfo->IpoptName().c_str(), real_val)) {
        pinfo->Jnlst()->Printf(J_ERROR, J_MAIN,
                               "\nInvalid value \"%s\" for option %s.\n", value, kw->name);
        THROW_EXCEPTION(OPTION_INVALID, "Invalid numeric option");
      }

      return retval;
    }

    static char* get_int_opt(Option_Info *oi, keyword *kw, char *value) {
      AmplOptionsList::PrivatInfo*
      pinfo = (AmplOptionsList::PrivatInfo*) kw->info;

      int int_val;
      kw->info = &int_val;
      char* retval = I_val(oi, kw, value);
      kw->info = (void*) pinfo;

      if (!pinfo->Options()->SetIntegerValue(pinfo->IpoptName().c_str(), int_val)) {
        pinfo->Jnlst()->Printf(J_ERROR, J_MAIN,
                               "\nInvalid value \"%s\" for option %s.\n", value, kw->name);
        THROW_EXCEPTION(OPTION_INVALID, "Invalid integer option");
      }

      return retval;
    }

    static char* get_str_opt(Option_Info *oi, keyword *kw, char *value) {
      AmplOptionsList::PrivatInfo*
      pinfo = (AmplOptionsList::PrivatInfo*) kw->info;

      char* str_val;
      kw->info = &str_val;
      char* retval = C_val(oi, kw, value);
      kw->info = (void*) pinfo;

      if (!pinfo->Options()->SetStringValue(pinfo->IpoptName().c_str(), str_val)) {
        pinfo->Jnlst()->Printf(J_ERROR, J_MAIN,
                               "\nInvalid value \"%s\" for option %s.\n", value, kw->name);
        THROW_EXCEPTION(OPTION_INVALID, "Invalid string option");
      }

      return retval;
    }

    static char* get_haltonerror_opt(Option_Info *oi, keyword *kw, char *value) {
      AmplOptionsList::PrivatInfo*
      pinfo = (AmplOptionsList::PrivatInfo*) kw->info;

      char* str_val;
      kw->info = &str_val;
      char* retval = C_val(oi, kw, value);
      kw->info = (void*) pinfo;

      fint** nerror = (fint**) pinfo->NError();

      if (strcmp(str_val, "yes")==0) {
        delete *nerror;
        *nerror = NULL;
      }
      else if (strcmp(str_val, "no")==0) {
        delete *nerror;
        *nerror = new fint;
        **nerror = 0;
      }
      else {
        pinfo->Jnlst()->Printf(J_ERROR, J_MAIN,
                               "\nInvalid value \"%s\" for option %s.\n", value, kw->name);
        THROW_EXCEPTION(OPTION_INVALID, "Invalid option");
      }

      return retval;
    }
  }


  AmplOptionsList::AmplOption::AmplOption(const std::string ipopt_option_name,
                                          AmplOptionType type,
                                          const std::string description)
      :
      ipopt_option_name_(ipopt_option_name),
      type_(type)
  {
    description_ = new char[description.size()+1];
    strcpy(description_, description.c_str());
  }

  AmplOptionsList::~AmplOptionsList()
  {
    if (keywds_) {
      DBG_ASSERT(nkeywds_>0);
      keyword* keywords = (keyword*) keywds_;
      for (Index i=0; i<nkeywds_; i++) {
        PrivatInfo* pinfo = (PrivatInfo*) keywords[i].info;
        delete pinfo;
        delete [] keywords[i].name;
      }
      delete [] keywords;
    }
  }

  void* AmplOptionsList::Keywords(const SmartPtr<OptionsList>& options,
                                  SmartPtr<const Journalist> jnlst,
                                  void** nerror)
  {
    if (keywds_) {
      DBG_ASSERT(nkeywds_>0);
      keyword* keywords = (keyword*) keywds_;
      for (Index i=0; i<nkeywds_; i++) {
        PrivatInfo* pinfo = (PrivatInfo*) keywords[i].info;
        delete pinfo;
        delete [] keywords[i].name;
      }
      delete [] keywords;
      nkeywds_ = 0;
    }

    /* PrintLatex(jnlst); */

    Index n_options = NumberOfAmplOptions();
    keyword* keywords = new keyword[n_options];

    Index ioption = 0;
    for (std::map<std::string, SmartPtr<const AmplOption> >::iterator
         iter = ampl_options_map_.begin();
         iter != ampl_options_map_.end(); iter++) {
      keywords[ioption].name = new char[iter->first.size()+1];
      strcpy(keywords[ioption].name, iter->first.c_str());
      keywords[ioption].desc = iter->second->Description();
      switch (iter->second->Type()) {
      case String_Option: {
          PrivatInfo* pinfo = new PrivatInfo(iter->second->IpoptOptionName(), options, jnlst);
          keywords[ioption].info = (void*) pinfo;
          keywords[ioption].kf = get_str_opt;
        }
        break;
      case Number_Option: {
          PrivatInfo* pinfo = new PrivatInfo(iter->second->IpoptOptionName(), options, jnlst);
          keywords[ioption].info = (void*) pinfo;
          keywords[ioption].kf = get_num_opt;
        }
        break;
      case Integer_Option: {
          PrivatInfo* pinfo = new PrivatInfo(iter->second->IpoptOptionName(), options, jnlst);
          keywords[ioption].info = (void*) pinfo;
          keywords[ioption].kf = get_int_opt;
        }
        break;
      case WS_Option:
        keywords[ioption].info = NULL;
        keywords[ioption].kf = WS_val;
        break;
      case HaltOnError_Option:
        PrivatInfo* pinfo = new PrivatInfo(iter->second->IpoptOptionName(), options, jnlst, nerror);
        keywords[ioption].info = (void*) pinfo;
        keywords[ioption].kf = get_haltonerror_opt;
        break;
      }
      ioption++;
    }

    DBG_ASSERT(ioption==n_options);
    nkeywds_ = n_options;
    keywds_ = (void*) keywords;
    return keywds_;
  }

  void AmplOptionsList::MakeValidLatexString(std::string source, std::string& dest) const
  {
    std::string::iterator c;
    for (c=source.begin(); c!=source.end(); c++) {
      if (*c == '_') {
        dest.append("\\_");
      }
      else if (*c == '^') {
        dest.append("\\^");
      }
      else if (*c == '>' || *c == '<') {
        ;
      }
      else if (*c == '\n') {
        dest += ", ";
      }
      else {
        dest += *c;
      }
    }
  }

  void AmplOptionsList::PrintLatex(SmartPtr<const Journalist> jnlst)
  {
     jnlst->Printf(J_SUMMARY, J_DOCUMENTATION, "\\begin{description}\n");
     for( std::map<std::string, SmartPtr<const AmplOption> >::iterator iter = ampl_options_map_.begin();
          iter != ampl_options_map_.end(); ++iter)
     {
        std::string amplname;
        std::string ipoptname;
        std::string descr;

        MakeValidLatexString(iter->first.c_str(), amplname);
        MakeValidLatexString(iter->second->IpoptOptionName(), ipoptname);
        MakeValidLatexString(iter->second->Description(), descr);

        if( ipoptname.length() > 0 )
        {
           jnlst->Printf(J_SUMMARY, J_DOCUMENTATION, "\\item[{\\htmlref{%s}{opt:%s}}]",
              amplname.c_str(), iter->second->IpoptOptionName().c_str());
           if( amplname != ipoptname )
              jnlst->Printf(J_SUMMARY, J_DOCUMENTATION, " (Ipopt name: \\htmlref{%s}{opt:%s})",
                 ipoptname.c_str(), iter->second->IpoptOptionName().c_str());
        }
        else
        {
           jnlst->Printf(J_SUMMARY, J_DOCUMENTATION, "\\item[%s]", amplname.c_str());
        }
        jnlst->Printf(J_SUMMARY, J_DOCUMENTATION, " %s\n", descr.c_str());
     }
     jnlst->Printf(J_SUMMARY, J_DOCUMENTATION, "\\end{description}\n");
  }

  char*
  AmplTNLP::get_options(const SmartPtr<OptionsList>& options,
                        SmartPtr<AmplOptionsList>& ampl_options_list,
                        const char* ampl_option_string,
                        const char* ampl_invokation_string,
                        const char* ampl_banner_string, char**& argv)
  {
    ASL_pfgh* asl = asl_;

    if (!IsValid(ampl_options_list)) {
      ampl_options_list = new AmplOptionsList();
    }

    // Output
    ampl_options_list->AddAmplOption("print_level",
                                     "print_level",
                                     AmplOptionsList::Integer_Option,
                                     "Verbosity level");
    ampl_options_list->AddAmplOption("outlev",
                                     "print_level",
                                     AmplOptionsList::Integer_Option,
                                     "Verbosity level");
    ampl_options_list->AddAmplOption("print_user_options",
                                     "print_user_options",
                                     AmplOptionsList::String_Option,
                                     "Toggle printing of user options");
    ampl_options_list->AddAmplOption("print_options_documentation",
                                     "print_options_documentation",
                                     AmplOptionsList::String_Option,
                                     "Print all available options (for ipopt.opt)");
    ampl_options_list->AddAmplOption("output_file",
                                     "output_file",
                                     AmplOptionsList::String_Option,
                                     "File name of an output file (leave unset for no file output)");
    ampl_options_list->AddAmplOption("file_print_level",
                                     "file_print_level",
                                     AmplOptionsList::Integer_Option,
                                     "Verbosity level for output file");
    ampl_options_list->AddAmplOption("option_file_name",
                                     "option_file_name",
                                     AmplOptionsList::String_Option,
                                     "File name of options file (default: ipopt.opt)");

    // Termination
    ampl_options_list->AddAmplOption("tol",
                                     "tol",
                                     AmplOptionsList::Number_Option,
                                     "Desired convergence tolerance (relative)");
    ampl_options_list->AddAmplOption("max_iter",
                                     "max_iter",
                                     AmplOptionsList::Integer_Option,
                                     "Maximum number of iterations");
    ampl_options_list->AddAmplOption("maxit",
                                     "max_iter",
                                     AmplOptionsList::Integer_Option,
                                     "Maximum number of iterations");
    ampl_options_list->AddAmplOption("max_cpu_time",
                                     "max_cpu_time",
                                     AmplOptionsList::Number_Option,
                                     "CPU time limit");
    ampl_options_list->AddAmplOption("compl_inf_tol",
                                     "compl_inf_tol",
                                     AmplOptionsList::Number_Option,
                                     "Acceptance threshold for the complementarity conditions");
    ampl_options_list->AddAmplOption("dual_inf_tol",
                                     "dual_inf_tol",
                                     AmplOptionsList::Number_Option,
                                     "Desired threshold for the dual infeasibility");
    ampl_options_list->AddAmplOption("constr_viol_tol",
                                     "constr_viol_tol",
                                     AmplOptionsList::Number_Option,
                                     "Desired threshold for the constraint violation");
    ampl_options_list->AddAmplOption("acceptable_tol",
                                     "acceptable_tol",
                                     AmplOptionsList::Number_Option,
                                     "Acceptable convergence tolerance (relative)");
    ampl_options_list->AddAmplOption("acceptable_compl_inf_tol",
                                     "acceptable_compl_inf_tol",
                                     AmplOptionsList::Number_Option,
                                     "Acceptance threshold for the complementarity conditions");
    ampl_options_list->AddAmplOption("acceptable_dual_inf_tol",
                                     "acceptable_dual_inf_tol",
                                     AmplOptionsList::Number_Option,
                                     "Acceptance threshold for the dual infeasibility");
    ampl_options_list->AddAmplOption("acceptable_constr_viol_tol",
                                     "acceptable_constr_viol_tol",
                                     AmplOptionsList::Number_Option,
                                     "Acceptance threshold for the constraint violation");

    ampl_options_list->AddAmplOption("diverging_iterates_tol",
                                     "diverging_iterates_tol",
                                     AmplOptionsList::Number_Option,
                                     "Threshold for maximal value of primal iterates");

    // NLP scaling
    ampl_options_list->AddAmplOption("obj_scaling_factor",
                                     "obj_scaling_factor",
                                     AmplOptionsList::Number_Option,
                                     "Scaling factor for the objective function");
    ampl_options_list->AddAmplOption("nlp_scaling_method",
                                     "nlp_scaling_method",
                                     AmplOptionsList::String_Option,
                                     "Select the technique used for scaling the NLP");
    ampl_options_list->AddAmplOption("nlp_scaling_max_gradient",
                                     "nlp_scaling_max_gradient",
                                     AmplOptionsList::Number_Option,
                                     "Maximum gradient after scaling");

    // NLP corrections
    ampl_options_list->AddAmplOption("bound_relax_factor",
                                     "bound_relax_factor",
                                     AmplOptionsList::Number_Option,
                                     "Factor for initial relaxation of the bounds");
    ampl_options_list->AddAmplOption("honor_original_bounds",
                                     "honor_original_bounds",
                                     AmplOptionsList::String_Option,
                                     "If no, solution might slightly violate bounds");

    // Barrier parameter
    ampl_options_list->AddAmplOption("mu_strategy",
                                     "mu_strategy",
                                     AmplOptionsList::String_Option,
                                     "Update strategy for barrier parameter");
    ampl_options_list->AddAmplOption("mu_oracle",
                                     "mu_oracle",
                                     AmplOptionsList::String_Option,
                                     "Oracle for a new barrier parameter in the adaptive strategy");
    // Barrier parameter
    ampl_options_list->AddAmplOption("mu_max",
                                     "mu_max",
                                     AmplOptionsList::Number_Option,
                                     "Maximal value for barrier parameter for adaptive strategy");
    ampl_options_list->AddAmplOption("mu_init",
                                     "mu_init",
                                     AmplOptionsList::Number_Option,
                                     "Initial value for the barrier parameter");

    // Initialization
    ampl_options_list->AddAmplOption("bound_frac",
                                     "bound_frac",
                                     AmplOptionsList::Number_Option,
                                     "Desired minimal relative distance of initial point to bound");
    ampl_options_list->AddAmplOption("bound_push",
                                     "bound_push",
                                     AmplOptionsList::Number_Option,
                                     "Desired minimal absolute distance of initial point to bound");
    ampl_options_list->AddAmplOption("slack_bound_frac",
                                     "slack_bound_frac",
                                     AmplOptionsList::Number_Option,
                                     "Desired minimal relative distance of initial slack to bound");
    ampl_options_list->AddAmplOption("slack_bound_push",
                                     "slack_bound_push",
                                     AmplOptionsList::Number_Option,
                                     "Desired minimal absolute distance of initial slack to bound");
    ampl_options_list->AddAmplOption("bound_mult_init_val",
                                     "bound_mult_init_val",
                                     AmplOptionsList::Number_Option,
                                     "Initial value for the bound multipliers");
    ampl_options_list->AddAmplOption("constr_mult_init_max",
                                     "constr_mult_init_max",
                                     AmplOptionsList::Number_Option,
                                     "Maximal allowed least-square guess of constraint multipliers");

    // Multiplier updates
    ampl_options_list->AddAmplOption("alpha_for_y",
                                     "alpha_for_y",
                                     AmplOptionsList::String_Option,
                                     "Step size for constraint multipliers");

    // Line search
    ampl_options_list->AddAmplOption("max_soc",
                                     "max_soc",
                                     AmplOptionsList::Integer_Option,
                                     "Maximal number of second order correction trial steps");
    ampl_options_list->AddAmplOption("watchdog_shortened_iter_trigger",
                                     "watchdog_shortened_iter_trigger",
                                     AmplOptionsList::Integer_Option,
                                     "Trigger counter for watchdog procedure");

    // Restoration phase
    ampl_options_list->AddAmplOption("expect_infeasible_problem",
                                     "expect_infeasible_problem",
                                     AmplOptionsList::String_Option,
                                     "Enable heuristics to quickly detect an infeasible problem");
    ampl_options_list->AddAmplOption("required_infeasibility_reduction",
                                     "required_infeasibility_reduction",
                                     AmplOptionsList::Number_Option,
                                     "Required infeasibility reduction in restoration phase");

    // Added for Warm-Start
    ampl_options_list->AddAmplOption("warm_start_init_point",
                                     "warm_start_init_point",
                                     AmplOptionsList::String_Option,
                                     "Enables to specify bound multiplier values");
    ampl_options_list->AddAmplOption("warm_start_bound_push",
                                     "warm_start_bound_push",
                                     AmplOptionsList::Number_Option,
                                     "Enables to specify how much should variables should be pushed inside the feasible region");
    ampl_options_list->AddAmplOption("warm_start_mult_bound_push",
                                     "warm_start_mult_bound_push",
                                     AmplOptionsList::Number_Option,
                                     "Enables to specify how much should bound multipliers should be pushed inside the feasible region");

    // Quasi-Newton
    ampl_options_list->AddAmplOption("hessian_approximation",
                                     "hessian_approximation",
                                     AmplOptionsList::String_Option,
                                     "Can enable Quasi-Newton approximation of hessian");
    // Linear solver
    ampl_options_list->AddAmplOption("linear_solver",
                                     "linear_solver",
                                     AmplOptionsList::String_Option,
                                     "Linear solver to be used for step calculation");
    ampl_options_list->AddAmplOption("linear_system_scaling",
                                     "linear_system_scaling",
                                     AmplOptionsList::String_Option,
                                     "Method for scaling the linear systems");
    ampl_options_list->AddAmplOption("linear_scaling_on_demand",
                                     "linear_scaling_on_demand",
                                     AmplOptionsList::String_Option,
                                     "Enables heuristic for scaling only when seems required");
    ampl_options_list->AddAmplOption("max_refinement_steps",
                                     "max_refinement_steps",
                                     AmplOptionsList::Integer_Option,
                                     "Maximal number of iterative refinement steps per linear system solve");
    ampl_options_list->AddAmplOption("min_refinement_steps",
                                     "min_refinement_steps",
                                     AmplOptionsList::Integer_Option,
                                     "Minimum number of iterative refinement steps per linear system solve");


    // Quasi-Newton
    ampl_options_list->AddAmplOption("hessian_approximation",
                                     "hessian_approximation",
                                     AmplOptionsList::String_Option,
                                     "Can enable Quasi-Newton approximation of hessian");

    // Special linear solver options
    ampl_options_list->AddAmplOption("ma27_pivtol",
                                     "ma27_pivtol",
                                     AmplOptionsList::Number_Option,
                                     "Pivot tolerance for the linear solver MA27");
    ampl_options_list->AddAmplOption("ma27_pivtolmax",
                                     "ma27_pivtolmax",
                                     AmplOptionsList::Number_Option,
                                     "Maximal pivot tolerance for the linear solver MA27");

    ampl_options_list->AddAmplOption("ma57_pivtol",
                                     "ma57_pivtol",
                                     AmplOptionsList::Number_Option,
                                     "Pivot tolerance for the linear solver MA57");
    ampl_options_list->AddAmplOption("ma57_pivtolmax",
                                     "ma57_pivtolmax",
                                     AmplOptionsList::Number_Option,
                                     "Maximal pivot tolerance for the linear solver MA57");
    ampl_options_list->AddAmplOption("ma57_pivot_order",
                                     "ma57_pivot_order",
                                     AmplOptionsList::Integer_Option,
                                     "Controls pivot order in MA57");

    ampl_options_list->AddAmplOption("pardiso_matching_strategy",
                                     "pardiso_matching_strategy",
                                     AmplOptionsList::String_Option,
                                     "Matching strategy for linear solver Pardiso");
    //ampl_options_list->AddAmplOption("pardiso_out_of_core_power",
    //                                 "pardiso_out_of_core_power",
    //                                 AmplOptionsList::Integer_Option,
    //                                 "Enables out-of-core version of linear solver Pardiso");

#ifdef HAVE_WSMP

    ampl_options_list->AddAmplOption("wsmp_num_threads",
                                     "wsmp_num_threads",
                                     AmplOptionsList::Integer_Option,
                                     "Number of threads to be used in WSMP");
    ampl_options_list->AddAmplOption("wsmp_pivtol",
                                     "wsmp_pivtol",
                                     AmplOptionsList::Number_Option,
                                     "Pivot tolerance for the linear solver WSMP");
    ampl_options_list->AddAmplOption("wsmp_pivtolmax",
                                     "wsmp_pivtolmax",
                                     AmplOptionsList::Number_Option,
                                     "Maximum pivot tolerance for the linear solver WSMP");
    ampl_options_list->AddAmplOption("wsmp_scaling",
                                     "wsmp_scaling",
                                     AmplOptionsList::Integer_Option,
                                     "Determines how the matrix is scaled by WSMP");
#endif

    // AMPL's wantsol option
    ampl_options_list->AddAmplOption("wantsol", "",
                                     AmplOptionsList::WS_Option,
                                     WS_desc_ASL+5);

    // special AMPL option to exit when there is in error in the
    // function evaluation
    ampl_options_list->AddAmplOption("halt_on_ampl_error", "",
                                     AmplOptionsList::HaltOnError_Option,
                                     "Exit with message on evaluation error");

    int n_options = ampl_options_list->NumberOfAmplOptions();

    keyword* keywds =
      (keyword*) ampl_options_list->Keywords(options, jnlst_,
                                             (void**)&nerror_);

    static const char sname_default[] = "ipopt";
    static const char bsname_default[] = "Ipopt " IPOPT_VERSION;
    static const char opname_default[] = "ipopt_options";
    const char* sname;
    const char* bsname;
    const char* opname;
    if (ampl_option_string) {
      opname = ampl_option_string;
    }
    else {
      opname = opname_default;
    }
    if (ampl_invokation_string) {
      sname = ampl_invokation_string;
    }
    else {
      sname = sname_default;
    }
    if (ampl_banner_string) {
      bsname = ampl_banner_string;
    }
    else {
      bsname = bsname_default;
    }

    DBG_ASSERT(!Oinfo_ptr_);
    Option_Info* Oinfo = new Option_Info;
    Oinfo->sname = new char[strlen(sname)+1];
    strcpy(Oinfo->sname, sname);
    Oinfo->bsname = new char[strlen(bsname)+1];
    strcpy(Oinfo->bsname, bsname);
    Oinfo->opname = new char[strlen(opname)+1];
    strcpy(Oinfo->opname, opname);
    Oinfo->keywds = keywds;
    Oinfo->n_keywds = n_options;
    // Set the default for the remaining entries
    Oinfo->flags = 0;
    Oinfo->version = NULL;
    Oinfo->usage = NULL;
    Oinfo->kwf = NULL;
    Oinfo->feq = NULL;
    Oinfo->options = NULL;
    Oinfo->n_options = 0;
    Oinfo->driver_date = 0;
    Oinfo->wantsol = 0;
    Oinfo->nS = 0;
    Oinfo->S = NULL;
    Oinfo->uinfo = NULL;
    Oinfo->asl = NULL;
    Oinfo->eqsign = NULL;
    Oinfo->n_badopts = 0;
    Oinfo->option_echo = 0;
    Oinfo->nnl = 0;

    Oinfo_ptr_ = Oinfo;

    char* stub = getstops(argv, Oinfo);

    return stub;
  }

  bool AmplTNLP::nerror_ok(void* nerror)
  {
    DBG_START_METH("AmplTNLP::nerror_ok",
                   dbg_verbosity);

    if (nerror == NULL || *((fint*)nerror) == 0) {
      return true;
    }
    jnlst_->Printf(J_ERROR, J_MAIN, "Error in an AMPL evaluation. Run with \"halt_on_ampl_error yes\" to see details.\n");
    DBG_PRINT((1, "nerror = %d\n", *((fint*)nerror)));
    return false;
  }

  AmplSuffixHandler::AmplSuffixHandler()
      :
      asl_(NULL),
      suftab_ (NULL)
  {}

  AmplSuffixHandler::~AmplSuffixHandler()
  {
    if (suftab_) {
      Index n = (Index)suffix_ids_.size();
      for (Index i=0; i<n; i++) {
        delete [] suftab_[i].name;
        suftab_[i].name = NULL;
      }
    }
    delete [] suftab_;
    suftab_ = NULL;
  }

  void AmplSuffixHandler::PrepareAmplForSuffixes(ASL_pfgh* asl)
  {
    DBG_ASSERT(asl);
    asl_ = asl;

    Index n = (Index)suffix_ids_.size();
    suftab_ = new SufDecl[n];
    for (Index i=0; i<n; i++) {
      Index id_len = (Index)strlen(suffix_ids_[i].c_str());
      suftab_[i].name = new char[id_len + 1];
      strcpy(suftab_[i].name, suffix_ids_[i].c_str());

      suftab_[i].table = 0;

      if (suffix_sources_[i] == Variable_Source) {
        suftab_[i].kind = ASL_Sufkind_var;
      }
      else if (suffix_sources_[i]  == Constraint_Source) {
        suftab_[i].kind = ASL_Sufkind_con;
      }
      else if (suffix_sources_[i] == Objective_Source) {
        suftab_[i].kind = ASL_Sufkind_obj;
      }
      else if (suffix_sources_[i] == Problem_Source) {
        suftab_[i].kind = ASL_Sufkind_prob;
      }
      else {
        DBG_ASSERT(false && "Unknown suffix source in PrepareAmplForSuffixes");
      }

      if (suffix_types_[i] == Number_Type) {
        suftab_[i].kind = suftab_[i].kind | ASL_Sufkind_real;
      }

      suftab_[i].nextra = 0;
    }

    suf_declare(suftab_, n);
  }

  const Index*
  AmplSuffixHandler::GetIntegerSuffixValues(std::string suffix_string,
      Suffix_Source source) const
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl);

    int kind;
    if (source == Variable_Source) {
      kind = ASL_Sufkind_var;
    }
    else if (source == Constraint_Source) {
      kind = ASL_Sufkind_con;
    }
    else if (source == Objective_Source) {
      kind = ASL_Sufkind_obj;
    }
    else if (source == Problem_Source) {
      kind = ASL_Sufkind_prob;
    }
    else {
      kind = 0;
      DBG_ASSERT(false && "Unknown suffix source in GetIntegerSuffixValues");
    }
    SufDesc* dp = suf_get(suffix_string.c_str(), kind);
    return dp->u.i;
  }

  std::vector<Index>
  AmplSuffixHandler::GetIntegerSuffixValues(Index n, std::string suffix_string,
      Suffix_Source source) const
  {
    std::vector<Index> ret;
    const Index* ptr = GetIntegerSuffixValues(suffix_string, source);
    if (ptr) {
      ret.reserve(n);
      for (Index i=0; i<n; i++) {
        ret.push_back(ptr[i]);
      }
    }
    return ret;
  }

  const Number*
  AmplSuffixHandler::GetNumberSuffixValues(std::string suffix_string,
      Suffix_Source source) const
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl);

    int kind;
    if (source == Variable_Source) {
      kind = ASL_Sufkind_var;
    }
    else if (source == Constraint_Source) {
      kind = ASL_Sufkind_con;
    }
    else if (source == Objective_Source) {
      kind = ASL_Sufkind_obj;
    }
    else if (source == Problem_Source) {
      kind = ASL_Sufkind_prob;
    }
    else {
      kind = 0;
      DBG_ASSERT(false && "Unknown suffix source in GetNumberSuffixValues");
    }
    SufDesc* dp = suf_get(suffix_string.c_str(), kind);
    return dp->u.r;
  }

  std::vector<Number>
  AmplSuffixHandler::GetNumberSuffixValues(Index n, std::string suffix_string,
      Suffix_Source source) const
  {
    std::vector<Number> ret;
    const Number* ptr = GetNumberSuffixValues(suffix_string, source);
    if (ptr) {
      ret.reserve(n);
      for (Index i=0; i<n; i++) {
        ret.push_back(ptr[i]);
      }
    }
    return ret;
  }

} // namespace Ipopt



