// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

// TODO:
// - clean up the boolean for initialization
// - pass in a tag so we can know when x has changed
// - look closer at the interface - it should pass in non-zeros as well as m in the jacobian stuff (maybe hessian stuff)

/* AMPL includes */
extern "C"
{
#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"
}

#include "AmplTNLP.hpp"
#include "IpDenseVector.hpp"
#include "IpGenTMatrix.hpp"
#include "IpSymTMatrix.hpp"
#include "IpBlas.hpp"

namespace Ipopt
{

  AmplTNLP::AmplTNLP(const SmartPtr<Journalist>& jnlst, char**& argv,
                     SmartPtr<AmplSuffixHandler> suffix_handler /* = NULL */,
                     bool allow_discrete /* = false */)
      :
      TNLP(),
      jnlst_(ConstPtr(jnlst)),
      asl_(NULL),
      obj_sign_(1),
      nz_h_full_(-1),
      non_const_x_(NULL),
      x_sol_(NULL),
      z_L_sol_(NULL),
      z_U_sol_(NULL),
      g_sol_(NULL),
      obj_sol_(0.0),
      lambda_sol_(NULL),
      objval_called_with_current_x_(false),
      conval_called_with_current_x_(false)
  {

    // The ASL include files #define certain
    // variables that they expect you to work with.
    // These variables then appear as though they are
    // global variables when, in fact, they are not
    // Most of them are data members of an asl object

    // Create the ASL structure
    ASL_pfgh* asl = (ASL_pfgh*)ASL_alloc(ASL_read_pfgh);
    DBG_ASSERT(asl);
    asl_ = asl; // keep the pointer for ourselves to use later...

    // Read the options and stub
    // ToDo: Figure out the options stuff
    char* stub = getstub(&argv, NULL); // need to deal with options here
    if (!stub) {
      printf("No .nl file given!\n");
      exit(1);
    }

    jnlst_->Printf(J_SUMMARY, J_MAIN, "Ampl Model: %s\n", stub);

    // Parse the first part of the nl file
    //char* stub = argv[0];
    FILE* nl = jac0dim(stub, (fint)strlen(stub));
    DBG_ASSERT(nl);
    // check the problem statistics (see Table 1 in AMPL doc)
    DBG_ASSERT(n_var > 0); // need some continuous variables
    //    DBG_ASSERT(nbv == 0); // Cannot handle binary variables
    //    DBG_ASSERT(niv == 0); // Cannot handle integer variables
    ASSERT_EXCEPTION(allow_discrete || (nbv == 0 && niv == 0),
                     IpoptException,
                     "Discrete variables not allowed when the allow_discrete flag is false, "
                     "Either remove the integer variables, or change the flag in the constructor of AmplTNLP"
                    );
    // n_con can be >= 0
    DBG_ASSERT(n_obj == 1); // Currently can handle only 1 objective
    DBG_ASSERT(nlo == 0 || nlo == 1); // Can handle nonlinear obj.
    DBG_ASSERT(nwv == 0); // Don't know what "linear arc" variables are
    DBG_ASSERT(nlnc == 0); // Don't know what "nonlinear network"constraints are
    DBG_ASSERT(lnc == 0); // Don't know what "linear network" constraints are

    // Set options in the asl structure
    want_xpi0 = 1 | 2;  // allocate initial values for primal and dual if available
    DBG_ASSERT((want_xpi0 & 1) == 1 && (want_xpi0 & 2) == 2);
    obj_no = 0; // always want to work with the first (and only?) objective

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

    // see "changes" in solvers directory of ampl code...
    hesset(1,0,1,0,nlc);

    switch (retcode) {
      case ASL_readerr_none : {}
      break;
      case ASL_readerr_nofile : {
        printf("Cannot open .nl file\n");
        exit(1);
      }
      break;
      case ASL_readerr_nonlin : {
        DBG_ASSERT(false); // this better not be an error!
        printf("model involves nonlinearities (ed0read)\n");
        exit(1);
      }
      break;
      case  ASL_readerr_argerr : {
        printf("user-defined function with bad args\n");
        exit(1);
      }
      break;
      case ASL_readerr_unavail : {
        printf("user-defined function not available\n");
        exit(1);
      }
      break;
      case ASL_readerr_corrupt : {
        printf("corrupt .nl file\n");
        exit(1);
      }
      break;
      case ASL_readerr_bug : {
        printf("bug in .nl reader\n");
        exit(1);
      }
      break;
      default: {
        printf("Unknown error in stub file read\n");
        exit(1);
      }
      break;
    }

    obj_sign_ = 1; // minimization
    if (objtype[obj_no] != 0) {
      obj_sign_ = -1;
    }

    // find the nonzero structure for the hessian
    // parameters to sphsetup:
    int coeff_obj = 1; // coefficient of the objective fn ???
    int mult_supplied = 1; // multipliers will be supplied
    int uptri = 1; // only need the upper triangular part
    nz_h_full_ = sphsetup(-1, coeff_obj, mult_supplied, uptri);
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
      ASL_free((ASL**)&asl_);
      asl_ = NULL;
    }

    delete [] non_const_x_;
    non_const_x_ = NULL;
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
  }

  bool AmplTNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);

    n = n_var; // # of variables (variable types have been asserted in the constructor
    m = n_con; // # of constraints
    nnz_jac_g = nzc; // # of non-zeros in the jacobian
    nnz_h_lag = nz_h_full_; // # of non-zeros in the hessian

    return true;
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

  bool AmplTNLP::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    DBG_ASSERT(n == n_var);
    DBG_ASSERT(m == n_con);

    if (init_x) {
      for (Index i=0; i<n; i++) {
        if (havex0[i]) {
          x[i] = X0[i];
        }
        else {
          x[i] = 0.0;
        }
      }
    }

    if (init_z) {
      for (Index i=0; i<n; i++) {
        z_L[i] = z_U[i] = 1.0;
      }
    }

    if (init_lambda) {
      for (Index i=0; i<m; i++) {
        if (havepi0[i]) {
          lambda[i] = pi0[i];
        }
        else {
          lambda[i] = 0.0;
        }
      }
    }

    return true;
  }

  bool AmplTNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
  {
    apply_new_x(new_x, n, x);

    return internal_objval(obj_value);
  }

  bool AmplTNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);

    apply_new_x(new_x, n, x);

    fint nerror = 0;
    objgrd(0, non_const_x_, grad_f, &nerror);
    if (obj_sign_==-1) {
      for (Index i=0; i<n; i++) {
        grad_f[i] *= -1.;
      }
    }
    if (nerror == 0) {
      return true;
    }
    return false;
  }

  bool AmplTNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    DBG_ASSERT(n == n_var);
    DBG_ASSERT(m == n_con);

    apply_new_x(new_x, n, x);

    fint nerror = 0;
    return internal_conval(m, g);
  }

  bool AmplTNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow,
                            Index *jCol, Number* values)
  {
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
          //				iRow[current_nz] = i + 1;
          //				jCol[current_nz] = cg->varno+1;
          current_nz++;
        }
      }
      DBG_ASSERT(current_nz == nele_jac);
      return true;
    }
    else if (!iRow && !jCol && values) {
      apply_new_x(new_x, n, x);

      fint nerror = 0;
      jacval(non_const_x_, values, &nerror);

      if (nerror == 0) {
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
    else if (!iRow & !jCol && values) {
      apply_new_x(new_x, n, x);
      if (!objval_called_with_current_x_) {
        Number dummy;
        internal_objval(dummy);
      }
      if (!conval_called_with_current_x_) {
        internal_conval(m);
      }
      // copy lambda to non_const_lambda - note, we do not store a copy like
      // we do with x since lambda is only used here and not in other calls
      Number* non_const_lambda = new Number[m];
      for (Index i=0; i<m; i++) {
        non_const_lambda[i] = lambda[i];
      }

      fint nerror = 0;
      real ow=obj_sign_*obj_factor;
      sphes(values, -1, &ow, non_const_lambda);

      delete [] non_const_lambda;
      return true;
    }
    else {
      DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
    }

    return false;
  }

  void AmplTNLP::finalize_solution(ApplicationReturnStatus status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value)
  {
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
    IpBlasDcopy(n, z_L, 1, z_L_sol_, 1);
    IpBlasDcopy(n, z_U, 1, z_U_sol_, 1);
    IpBlasDcopy(m, g, 1, g_sol_, 1);
    IpBlasDcopy(m, lambda, 1, lambda_sol_, 1);
    obj_sol_ = obj_value;

    std::string message;
    if (status == Solve_Succeeded) {
      message = "Optimal Solution Found";
    }
    else if (status == Maximum_Iterations_Exceeded) {
      message = "Maximum Iterations Exceeded";
    }
    else if (status == Solve_Failed) {
      message = "Solve Failed";
    }
    else if (status == NonIpopt_Exception_Thrown) {
      message = "Non-Ipopt Exception Thrown";
    }
    else if (status == Internal_Error) {
      message = "Ipopt Internal Error";
    }
    else {
      message = "Unkown Error";
    }

    // Write the .sol file
    message = " \nEXIT: " + message;
    write_solution_file(message.c_str());
  }

  bool AmplTNLP::internal_objval(Number& obj_val)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    objval_called_with_current_x_ = false; // in case the call below fails

    fint nerror = 0;
    Number retval = objval(0, non_const_x_, &nerror);
    if (nerror == 0) {
      obj_val = obj_sign_*retval;
      objval_called_with_current_x_ = true;
      return true;
    }

    //DBG_ASSERT(false && "Error evaluating AMPL objective.\n");
    return false;
  }

  bool AmplTNLP::internal_conval(Index m, Number* g)
  {
    DBG_START_METH("AmplTNLP::internal_conval", 0);
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);
    DBG_ASSERT(m == n_con);
    conval_called_with_current_x_ = false; // in case the call below fails

    bool allocated = false;
    if (!g) {
      g = new double[m];
      allocated = true;
    }

    fint nerror = 0;
    conval(non_const_x_, g, &nerror);

    if (allocated) {
      delete [] g;
      g = NULL;
    }

    if (nerror == 0) {
      conval_called_with_current_x_ = true;
      return true;
    }
    return false;
  }


  void AmplTNLP::apply_new_x(bool new_x, Index n, const Number* x)
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl_);

    if (new_x) {
      // update the flags so these methods are called
      // before evaluating the hessian
      conval_called_with_current_x_ = false;
      objval_called_with_current_x_ = false;

      //copy the data to the non_const_x_
      if (!non_const_x_) {
        non_const_x_ = new Number[n];
      }

      for (Index i=0; i<n; i++) {
        non_const_x_[i] = x[i];
      }

      // tell ampl that we have a new x
      xknown(non_const_x_);
    }
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

    write_sol(cmessage, x_sol_, lambda_sol_, NULL);

    delete [] cmessage;
  }

  void AmplTNLP::get_discrete_info(Index& n_binaries, Index& n_integers) const
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl);

    n_binaries = nbv;
    n_integers = niv;
  }

  AmplSuffixHandler::AmplSuffixHandler()
      :
      asl_(NULL),
      suftab_ (NULL)
  {}

  AmplSuffixHandler::~AmplSuffixHandler()
  {
    if (suftab_) {
      Index n = suffix_ids_.size();
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

    Index n = suffix_ids_.size();
    suftab_ = new SufDecl[n];
    for (Index i=0; i<n; i++) {
      Index id_len = strlen(suffix_ids_[i].c_str());
      suftab_[i].name = new char[id_len + 1];
      strcpy(suftab_[i].name, suffix_ids_[i].c_str());

      suftab_[i].table = 0;

      if (suffix_sources_[i] == Variable_Source) {
        suftab_[i].kind = ASL_Sufkind_var;
      }
      else if (suffix_sources_[i]  == Constraint_Source) {
        suftab_[i].kind = ASL_Sufkind_con;
      }
      else {
        DBG_ASSERT(false && "Unknown suffix source in PrepareAmplForSuffixes");
      }

      suftab_[i].nextra = 0;
    }

    suf_declare(suftab_, n);
  }

  const Index* AmplSuffixHandler::GetIntegerSuffixValues(std::string suffix_string, Suffix_Source source) const
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl);

    int kind = (source == Variable_Source) ? ASL_Sufkind_var : ASL_Sufkind_con;
    SufDesc* dp = suf_get(suffix_string.c_str(), kind);
    return dp->u.i;
  }

  const Number* AmplSuffixHandler::GetNumberSuffixValues(std::string suffix_string, Suffix_Source source) const
  {
    ASL_pfgh* asl = asl_;
    DBG_ASSERT(asl);

    int kind = (source == Variable_Source) ? ASL_Sufkind_var : ASL_Sufkind_con;
    SufDesc* dp = suf_get(suffix_string.c_str(), kind);
    return dp->u.r;
  }


} // namespace Ipopt



