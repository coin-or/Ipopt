// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-11


#include "AsNMPCApplication.hpp"
#include "AsSchurBuilder.hpp"
#include "AsNmpcUtils.hpp"
#include "AsAsNMPCRegOp.hpp"

// Ipopt includes
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif


  NmpcApplication::NmpcApplication(SmartPtr<Journalist> jnlst,
				   SmartPtr<OptionsList> options,
				   SmartPtr<RegisteredOptions> reg_options)
    :
    jnlst_(jnlst),
    options_(options),
    reg_options_(reg_options)
  {
    DBG_START_METH("NmpcApplication::NmpcApplication", dbg_verbosity);
    
    // Initialize Journalist
    DBG_DO(SmartPtr<Journal> sens_jrnl = jnlst_->AddFileJournal("Sensitivity","sensdebug.out",J_ITERSUMMARY);
	   sens_jrnl->SetPrintLevel(J_USER1,J_ALL));

  }

  NmpcApplication::~NmpcApplication()
  {
    DBG_START_METH("NmpcApplication::~NmpcApplication", dbg_verbosity);
  }

  void NmpcApplication::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    // Options for NMPC Parameter Sensitivity
    roptions->SetRegisteringCategory("AsNmpc");
    roptions->AddLowerBoundedIntegerOption(
					   "n_sens_steps", "Number of steps computed by sIPOPT",
					   0, 0,
					   "");
    roptions->AddStringOption2(
			       "sens_boundcheck",
			       "Activate boundcheck and re-solve for sIPOPT",
			       "no",
			       "no", "don't check bounds and do another SchurSolve",
			       "yes", "check bounds and resolve Schur decomposition",
			       "If this option is activated, the algorithm will check the iterate after an initilal Schursolve and will resolve the decomposition if any bounds are not satisfied");
    roptions->AddStringOption2(
			       "nmpc_calc_style",
			       "Determines which implementation of SchurData and PCalculator will be used",
			       "index",
			       "std", "use StdPCalculator and StdSchurData",
			       "index", "use IndexPCalculator and IndexSchurData",
			       "Determines which implementation of SchurData and PCalculator will be used. Only the index style can do boundchecking.");
    roptions->AddLowerBoundedNumberOption(
					  "sens_bound_eps",
					  "Bound accuracy within which a bound still is considered to be valid",
					  0, true, 1e-3,
					  "The schur complement solution cannot make sure that variables stay inside bounds."
					  "I cannot use the primal-frac-to-the-bound step because I don't know if the initial iterate is feasible."
					  "To make things easier for me I have decided to make bounds not so strict.");
    roptions->AddStringOption2(
			       "compute_red_hessian",
			       "Determines if reduced hessian should be computed",
			       "no",
			       "yes", "compute reduced hessian",
			       "no", "don't compute reduced hessian",
			       "");
    roptions->AddStringOption4(
			       "select_step",
			       "Choose by which formula the step is computed",
			       "iftsensitivity",
			       "advanced","use advanced step based on KKT",
			       "sensitivity","use sensitivity step",
			       "ift","use one-parametric step with multiplier correction",
			       "iftsensitivity","use one-parametric step without multiplier correction"
			       "see original paper for more information on each step computation");
    // This option must be in IpInterfacesRegOp.cpp
    roptions->AddStringOption2(
			       "run_sens",
			       "Determines if sIPOPT alg runs",
			       "no",
			       "yes", "run sIPOPT",
			       "no", "don't run sIPOPT",
			       "");
    roptions->AddStringOption2(
			       "nmpc_internal_abort",
			       "Internal option - if set (internally), nmpc algorithm is not conducted",
			       "no",
			       "yes", "abort sIPOPT",
			       "no", "run sIPOPT",
			       "");
    roptions->AddStringOption2(
			       "redhess_internal_abort",
			       "Internal option - if set (internally), reduced hessian computation is not conducted",
			       "no",
			       "yes", "abort redhess computation",
			       "no", "run redhess computation",
			       "");
    roptions->AddStringOption2(
			       "ignore_suffix_error",
			       "If set, IPOPT runs even if there are errors in the suffixes",
			       "no",
			       "yes", "don't abort on suffix error",
			       "no", "abort on suffix error",
			       "");
    roptions->AddLowerBoundedNumberOption(
					  "sens_max_pdpert",
					  "Maximum perturbation of primal dual system, for that the sIPOPT algorithm will not abort",
					  0.0, true, 1e-3, 
					  "For certain problems, IPOPT uses inertia correction of the primal dual matrix to achieve"
					  "better convergence properties. This inertia correction changes the matrix and renders it"
					  "useless for the use with sIPOPT. This option sets an upper bound, which the inertia correction"
					  "may have. If any of the inertia correction values is above this bound, the sIPOPT algorithm"
					  "is aborted.");
    roptions->AddStringOption2(
			       "rh_eigendecomp",
			       "If yes, the eigenvalue decomposition of the reduced hessian matrix is computed",
			       "no",
			       "yes", "compute eigenvalue decomposition of reduced hessian",
			       "no", "don't compute eigenvalue decomposition of reduced hessian",
			       "The eigenvalue decomposition of the reduced hessian has different meanings depending on the specific problem. For parameter estimation problems, the eigenvalues are linked to the confidence interval of the parameters. See for example Victor Zavala's Phd thesis, chapter 4 for details.");
  }

  NmpControllerExitStatus NmpcApplication::Run()
  {
    DBG_START_METH("NmpcApplication::Run", dbg_verbosity);

    NmpControllerExitStatus retval = SOLVE_SUCCESS;

    bool nmpc_internal_abort, redhess_internal_abort;
    Options()->GetBoolValue("nmpc_internal_abort", nmpc_internal_abort, "");
    Options()->GetBoolValue("redhess_internal_abort", redhess_internal_abort, "");

    // Check for perturbation of primal dual system
    Number max_pdpert;
    Options()->GetNumericValue("sens_max_pdpert", max_pdpert, "");
    Number pdpert_x, pdpert_s, pdpert_c, pdpert_d;
    ip_data_->getPDPert(pdpert_x, pdpert_s, pdpert_c, pdpert_d);
    if (Max(pdpert_x, pdpert_s, pdpert_c, pdpert_d)>max_pdpert) {
      jnlst_->Printf(J_WARNING, J_MAIN, "\n\t--------------= Warning =--------------\nInertia correction of primal dual system is too large for meaningful AsNMPC results.\n"
		    "\t... aborting computation.\n"
		     "Set option sens_max_pdpert to a higher value (current: %f) to run AsNMPC algorithm anyway\n", max_pdpert);
      nmpc_internal_abort = true;
      redhess_internal_abort = true;
    }


    if (compute_red_hessian_ && !redhess_internal_abort) {
      SmartPtr<SchurBuilder> schur_builder = new SchurBuilder();
      const std::string prefix = ""; // I should be getting this somewhere else...
      SmartPtr<ReducedHessianCalculator> red_hess_calc = schur_builder->BuildRedHessCalc(*jnlst_,
											 *options_,
											 prefix,
											 *ip_nlp_,
											 *ip_data_,
											 *ip_cq_,
											 *pd_solver_);
      
      red_hess_calc->ComputeReducedHessian();
    }

    if (run_sens_ && n_sens_steps_>0 && !nmpc_internal_abort) {
      SmartPtr<SchurBuilder> schur_builder = new SchurBuilder();
      const std::string prefix = ""; // I should be getting this somewhere else...
      SmartPtr<AsNmpController> controller = schur_builder->BuildNmpc(*jnlst_,
								      *options_,
								      prefix,
								      *ip_nlp_,
								      *ip_data_,
								      *ip_cq_,
								      *pd_solver_);

      retval = controller->Run();
    }
    else if (run_sens_) {
      if (n_sens_steps_<=0) {
	jnlst_->Printf(J_WARNING, J_MAIN, "\n"
		       "The run_sens option was set to true, but the specified\n"
		       "number of advanced steps was set to zero.\n"
		       "Computation is aborted.\n\n");
      }
    }

    SolverReturn status = SUCCESS;

    if (IsValid(ip_data_->curr()) && IsValid(ip_data_->curr()->x())) {
      SmartPtr<const Vector> c;
      SmartPtr<const Vector> d;
      SmartPtr<const Vector> zL;
      SmartPtr<const Vector> zU;
      SmartPtr<const Vector> yc;
      SmartPtr<const Vector> yd;
      Number obj = 0.;

      switch (status) {
      case SUCCESS:
	c = ip_cq_->curr_c();
	d = ip_cq_->curr_d();
	obj = ip_cq_->curr_f();
	zL = ip_data_->curr()->z_L();
	zU = ip_data_->curr()->z_U();
	yc = ip_data_->curr()->y_c();
	yd = ip_data_->curr()->y_d();
      case MAXITER_EXCEEDED:
      case STOP_AT_TINY_STEP:
      case STOP_AT_ACCEPTABLE_POINT:
      case LOCAL_INFEASIBILITY:
      case USER_REQUESTED_STOP:
      case FEASIBLE_POINT_FOUND:
      case DIVERGING_ITERATES:
      case RESTORATION_FAILURE:
      case ERROR_IN_STEP_COMPUTATION:
	c = ip_cq_->curr_c();
	d = ip_cq_->curr_d();
	obj = ip_cq_->curr_f();
	zL = ip_data_->curr()->z_L();
	zU = ip_data_->curr()->z_U();
	yc = ip_data_->curr()->y_c();
	yd = ip_data_->curr()->y_d();
	break;
      default: {
	SmartPtr<Vector> tmp = ip_data_->curr()->y_c()->MakeNew();
	tmp->Set(0.);
	c = ConstPtr(tmp);
	yc = ConstPtr(tmp);
	tmp = ip_data_->curr()->y_d()->MakeNew();
	tmp->Set(0.);
	d = ConstPtr(tmp);
	yd = ConstPtr(tmp);
	tmp = ip_data_->curr()->z_L()->MakeNew();
	tmp->Set(0.);
	zL = ConstPtr(tmp);
	tmp = ip_data_->curr()->z_U()->MakeNew();
	tmp->Set(0.);
	zU = ConstPtr(tmp);
      }
      }

      if (redhess_internal_abort) {
	jnlst_->Printf(J_WARNING, J_MAIN, "\nReduced hessian was not computed "
		       "because an error occured.\n"
		       "See exception message above for details.\n\n");
      }
      if (nmpc_internal_abort) {
	jnlst_->Printf(J_WARNING, J_MAIN, "\nNMPC controller was not called "
		       "because an error occured.\n"
		       "See exception message above for details.\n\n");
      }

      ip_nlp_->FinalizeSolution(status,
				*ip_data_->curr()->x(),
				*zL, *zU, *c, *d, *yc, *yd,
				obj, GetRawPtr(ip_data_), GetRawPtr(ip_cq_));
    }
    return retval;
  }

  void NmpcApplication::Initialize() 
  {
    DBG_START_METH("NmpcApplication::Initialize", dbg_verbosity);

    const std::string prefix = ""; // I should be getting this somewhere else...

    Options()->GetIntegerValue("n_sens_steps",n_sens_steps_, prefix.c_str());
    Options()->GetBoolValue("run_sens", run_sens_, prefix.c_str());
    Options()->GetBoolValue("compute_red_hessian", compute_red_hessian_, prefix.c_str());

    // make sure run_sens and skip_finalize_solution_call are consistent
    if (run_sens_ || compute_red_hessian_) {
      Options()->SetStringValue("skip_finalize_solution_call", "yes");
    }
    else {
      Options()->SetStringValue("skip_finalize_solution_call", "no");
    }

  }

  void NmpcApplication::SetIpoptAlgorithmObjects(SmartPtr<IpoptApplication> app_ipopt,
						 ApplicationReturnStatus ipopt_retval)
  {
    DBG_START_METH("NmpcApplication::SetIpoptAlgorithmObjects", dbg_verbosity);

    // get optionsList and Journalist
    options_ = app_ipopt->Options();
    jnlst_ = app_ipopt->Jnlst();

    // Check whether Ipopt solved to optimality - if not, end computation.
    if ( ipopt_retval != Solve_Succeeded ) {
      jnlst_->Printf(J_ERROR, J_MAIN, "ASNMPC: Aborting AsNMPC computation, because IPOPT did not succeed\n\n");
      options_->SetStringValue("nmpc_internal_abort", "yes");
      options_->SetStringValue("redhess_internal_abort", "yes");
    }

    // get pointers from IpoptApplication assessor methods
    SmartPtr<IpoptAlgorithm> alg = app_ipopt->AlgorithmObject();

    SmartPtr<PDSearchDirCalculator> pd_search;
    pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));

    // get PD_Solver
    pd_solver_ = pd_search->PDSolver();
  

    // get data 
    ip_data_ = app_ipopt->IpoptDataObject();

    // get calulated quantities
    ip_cq_ = app_ipopt->IpoptCQObject();

    // get NLP
    ip_nlp_ = app_ipopt->IpoptNLPObject();

    options_->GetIntegerValue("n_sens_steps",n_sens_steps_,"");

    // This checking should be rewritten
    /*    if (false && run_sens_) {
    // check suffixes
    std::string state;
    std::string state_value;
    const Index* index;
    const Number* number;
    Index n_nmpc_indices, n_this_nmpc_indices;
    // collect information from suffixes
    state = "sens_state_1";
    //index = ampl_tnlp_->get_index_suffix(state.c_str());
    if (index==NULL) {
    THROW_EXCEPTION(NMPC_SUFFIX_ERROR, "Suffix sens_state_1 is not set");
    }
    n_nmpc_indices = AsIndexSum(ip_data_->curr()->x()->Dim(), index, 1);
    for (Index i=1; i<=n_sens_steps_; ++i) {
    state = "sens_state_";
    state_value = "sens_state_value_";
    append_Index(state, i);
    append_Index(state_value, i);
    //index = ampl_tnlp_->get_index_suffix(state.c_str());
    if (index==NULL) {
    std::string msg = "Suffix " + state + " is not set";
    THROW_EXCEPTION(NMPC_SUFFIX_ERROR, msg);
    }
    n_this_nmpc_indices = AsIndexSum(ip_data_->curr()->x()->Dim(), index, 1);
    if (n_this_nmpc_indices!=n_nmpc_indices) {
    std::string msg = "Suffix" + state + "does not have the correct number of flags";
    THROW_EXCEPTION(NMPC_SUFFIX_ERROR, msg);
    }
    //number = ampl_tnlp_->get_number_suffix(state_value.c_str());
    if (number==NULL) {
    std::string msg = "Suffix " + state_value + " is not set";
    THROW_EXCEPTION(NMPC_SUFFIX_ERROR, msg);
    }
    }
    } */
  }

}
