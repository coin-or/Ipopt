// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-16


#include "SensStdStepCalc.hpp"
#include "IpDenseVector.hpp"
#include "IpIteratesVector.hpp"
#include "IpBlas.hpp"
#include "SensIndexSchurData.hpp"


namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 1;
#endif

  StdStepCalculator::StdStepCalculator(SmartPtr<SchurData> ift_data,
				       SmartPtr<SensBacksolver> backsolver)
    :
    ift_data_(ift_data),
    backsolver_(backsolver),
    bound_eps_(1e-3),
    kkt_residuals_(true)
  {
    DBG_START_METH("StdStepCalculator::StdStepCalculator", dbg_verbosity);
  }

  StdStepCalculator::~StdStepCalculator()
  {
    DBG_START_METH("StdStepCalculator::~StdStepCalculator", dbg_verbosity);
  }

  bool StdStepCalculator::InitializeImpl(const OptionsList& options,
					 const std::string& prefix)
  {
    options.GetNumericValue("sens_bound_eps", bound_eps_, prefix);
    options.GetBoolValue("sens_kkt_residuals", kkt_residuals_, prefix);
    SensitivityStepCalculator::InitializeImpl(options,
					      prefix);
    return true;
  }


  bool StdStepCalculator::Step(DenseVector& delta_u,
			       IteratesVector& sol)
  {
    DBG_START_METH("StdStepCalculator::Step", dbg_verbosity);

    bool retval;

    SmartPtr<IteratesVector> delta_u_long = IpData().trial()->MakeNewIteratesVector();
    ift_data_->TransMultiply(delta_u, *delta_u_long);

    SmartPtr<IteratesVector> r_s = IpData().trial()->MakeNewIteratesVector();
    if (kkt_residuals_) {
      /* This should be almost zero... */
      r_s->Set_x_NonConst(*IpCq().curr_grad_lag_x()->MakeNewCopy());
      r_s->Set_s_NonConst(*IpCq().curr_grad_lag_s()->MakeNewCopy());
      r_s->Set_y_c_NonConst(*IpCq().curr_c()->MakeNewCopy());
      r_s->Set_y_d_NonConst(*IpCq().curr_d_minus_s()->MakeNewCopy());
      r_s->Set_z_L_NonConst(*IpCq().curr_compl_x_L()->MakeNewCopy());
      r_s->Set_z_U_NonConst(*IpCq().curr_compl_x_U()->MakeNewCopy());
      r_s->Set_v_L_NonConst(*IpCq().curr_compl_s_L()->MakeNewCopy());
      r_s->Set_v_U_NonConst(*IpCq().curr_compl_s_U()->MakeNewCopy());

      r_s->Print(Jnlst(),J_VECTOR,J_USER1,"r_s init");
      delta_u.Print(Jnlst(),J_VECTOR,J_USER1,"delta_u init");
      DBG_PRINT((dbg_verbosity,"r_s init Nrm2=%23.16e\n", r_s->Asum()));

      delta_u_long->Axpy(-1.0, *r_s);
    }

    backsolver_->Solve(&sol, ConstPtr(delta_u_long));

    SmartPtr<IteratesVector> Kr_s;
    if (Do_Boundcheck()) {
      Kr_s = sol.MakeNewIteratesVectorCopy();
    }

    sol.Axpy(1.0, *IpData().trial());

    if (Do_Boundcheck()) {
      DBG_PRINT((dbg_verbosity, "Entering boundcheck"));
      // initialize
      Index new_du_size =0;
      Number* new_du_values;
      std::vector<Index> x_bound_violations_idx;
      std::vector<Number> x_bound_violations_du;
      std::vector<Index> delta_u_sort;
      bool bounds_violated;
      SmartPtr<DenseVectorSpace> delta_u_space = new DenseVectorSpace(0);
      SmartPtr<DenseVector> old_delta_u = new DenseVector(GetRawPtr(delta_u_space));
      SmartPtr<DenseVector> new_delta_u;

      bounds_violated = BoundCheck(sol, x_bound_violations_idx, x_bound_violations_du);
      while (bounds_violated) {
	Driver()->data_A()->Print(Jnlst(),J_VECTOR,J_USER1,"data_A_init");
	Driver()->data_B()->Print(Jnlst(),J_VECTOR,J_USER1,"data_B_init");
	// write new schurdata A
	dynamic_cast<IndexSchurData*>(GetRawPtr(Driver()->data_A_nonconst()))->AddData_List(x_bound_violations_idx, delta_u_sort, new_du_size, 1);
	// write new schurdata B
	dynamic_cast<IndexSchurData*>(GetRawPtr(Driver()->data_B_nonconst()))->AddData_List(x_bound_violations_idx, delta_u_sort, new_du_size, 1);
	Driver()->data_A()->Print(Jnlst(),J_VECTOR,J_USER1,"data_A");
	Driver()->data_B()->Print(Jnlst(),J_VECTOR,J_USER1,"data_B");
	Driver()->SchurBuild();
	Driver()->SchurFactorize();

	old_delta_u->Print(Jnlst(),J_VECTOR,J_USER1,"old_delta_u");
	delta_u_space = NULL; // delete old delta_u space
	delta_u_space = new DenseVectorSpace(new_du_size); // create new delta_u space
	new_delta_u = new DenseVector(GetRawPtr(ConstPtr(delta_u_space)));
	new_du_values = new_delta_u->Values();
	IpBlasDcopy(old_delta_u->Dim(), old_delta_u->Values(), 1, new_du_values, 1);
	for (Index i=0; i<x_bound_violations_idx.size(); ++i) {
	  //	  printf("i=%d, delta_u_sort[i]=%d, x_bound_viol_du[i]=%f\n", i, delta_u_sort[i], x_bound_violations_du[i]);
	  new_du_values[delta_u_sort[i]] = x_bound_violations_du[i];
	}
	SmartPtr<IteratesVector> new_sol = sol.MakeNewIteratesVector();
	new_delta_u->Print(Jnlst(),J_VECTOR,J_USER1,"new_delta_u");

	// solve with new data_B and delta_u
	retval = Driver()->SchurSolve(&sol, ConstPtr(delta_u_long), dynamic_cast<Vector*>(GetRawPtr(new_delta_u)), Kr_s);

	sol.Axpy(1.0, *IpData().trial());

	x_bound_violations_idx.clear();
	x_bound_violations_du.clear();
	delta_u_sort.clear();
	bounds_violated = BoundCheck(sol, x_bound_violations_idx, x_bound_violations_du);
	// copy new vector in old vector ->has to be done becpause otherwise only pointers will be copied and then it makes no sense
	old_delta_u = new_delta_u->MakeNewDenseVector();
	old_delta_u->Copy(*new_delta_u);
      }
    }

    return retval;
  }

  bool StdStepCalculator::BoundCheck(IteratesVector& sol,
				     std::vector<Index>& x_bound_violations_idx,
				     std::vector<Number>& x_bound_violations_du)
  {
    DBG_START_METH("StdStepCalculator::BoundCheck", dbg_verbosity);
    DBG_ASSERT(x_bound_violations_idx.empty());
    DBG_ASSERT(x_bound_violations_du.empty());

    // find bound violations in x vector
    const Number* x_val = dynamic_cast<const DenseVector*>(GetRawPtr(IpData().curr()->x()))->Values();
    const Number* sol_val = dynamic_cast<const DenseVector*>(GetRawPtr(sol.x()))->Values();

    SmartPtr<Vector> x_L_exp = IpData().curr()->x()->MakeNew();
    SmartPtr<Vector> x_U_exp = IpData().curr()->x()->MakeNew();

    SmartPtr<Vector> x_L_comp = IpNLP().x_L()->MakeNew();
    SmartPtr<Vector> x_U_comp = IpNLP().x_U()->MakeNew();

    IpNLP().Px_L()->TransMultVector(1.0, *sol.x(), 0.0, *x_L_comp);
    IpNLP().Px_U()->TransMultVector(1.0, *sol.x(), 0.0, *x_U_comp);

    x_L_comp->Print(Jnlst(),J_VECTOR,J_USER1,"x_L_comp");
    x_U_comp->Print(Jnlst(),J_VECTOR,J_USER1,"x_U_comp");
    //    return false;

    Number* x_L_val = dynamic_cast<DenseVector*>(GetRawPtr(x_L_comp))->Values();
    Number* x_U_val = dynamic_cast<DenseVector*>(GetRawPtr(x_U_comp))->Values();

    const Number* x_L_bound = dynamic_cast<const DenseVector*>(GetRawPtr(IpNLP().x_L()))->Values();
    const Number* x_U_bound = dynamic_cast<const DenseVector*>(GetRawPtr(IpNLP().x_U()))->Values();

    for (Index i=0; i<x_L_comp->Dim(); ++i) {
      x_L_val[i] -= x_L_bound[i];
    }

    for (Index i=0; i<x_U_comp->Dim(); ++i) {
      x_U_val[i] -= x_U_bound[i];
    }

    // project back
    IpNLP().Px_L()->MultVector(1.0, *x_L_comp, 0.0, *x_L_exp);
    IpNLP().Px_U()->MultVector(1.0, *x_U_comp, 0.0, *x_U_exp);

    const Number* x_L_exp_val = dynamic_cast<DenseVector*>(GetRawPtr(x_L_exp))->Values();
    const Number* x_U_exp_val = dynamic_cast<DenseVector*>(GetRawPtr(x_U_exp))->Values();

    for (Index i=0; i<x_L_exp->Dim(); ++i) {
      if (x_L_exp_val[i]<-bound_eps_) {
	x_bound_violations_idx.push_back(i);
	x_bound_violations_du.push_back(-x_L_exp_val[i]+sol_val[i]-x_val[i]); // this is just an awkward way to compute x_bound[i] - x_curr_val[i].
      } else if (-x_U_exp_val[i]<-bound_eps_) {
	x_bound_violations_idx.push_back(i);
	x_bound_violations_du.push_back(-x_U_exp_val[i]+sol_val[i]-x_val[i]);
      }
    }

    // z_L and z_U bound violations -> These are much easier since there is no projecting back and forth
    SmartPtr<const DenseVector> z_L = dynamic_cast<const DenseVector*>(GetRawPtr(sol.z_L()));
    SmartPtr<const DenseVector> z_U = dynamic_cast<const DenseVector*>(GetRawPtr(sol.z_U()));
    z_L->Print(Jnlst(),J_VECTOR,J_USER1,"z_L_boundcheck");
    z_U->Print(Jnlst(),J_VECTOR,J_USER1,"z_U_boundcheck");
    const Number* z_L_val = z_L->Values();
    const Number* z_U_val = z_U->Values();

    SmartPtr<const DenseVector> z_L_trial = dynamic_cast<const DenseVector*>(GetRawPtr(IpData().trial()->z_L()));
    SmartPtr<const DenseVector> z_U_trial = dynamic_cast<const DenseVector*>(GetRawPtr(IpData().trial()->z_U()));
    const Number* z_L_trial_val = z_L_trial->Values();
    const Number* z_U_trial_val = z_U_trial->Values();

    // find absolute index of z_L and z_U in IteratesVector
    Index z_L_ItVec_idx = 0;
    for (Index i=0; i<4; ++i) {
      z_L_ItVec_idx += (sol.GetComp(i))->Dim();
    }
    Index z_U_ItVec_idx = z_L_ItVec_idx + sol.z_L()->Dim();

    for (Index i=0; i<z_L->Dim(); ++i) {
      if (z_L_val[i]<-bound_eps_) {
	x_bound_violations_idx.push_back(i+z_L_ItVec_idx);
	x_bound_violations_du.push_back(-z_L_trial_val[i]);
	//printf("Lower Bound Mult. no. i=%d invalid: delta_u=%f\n", i+z_L_ItVec_idx, z_L_val[i]);
      }
    }

    for (Index i=0; i<z_U->Dim(); ++i) {
      if (z_U_val[i]<-bound_eps_) {
	x_bound_violations_idx.push_back(i+z_U_ItVec_idx);
	x_bound_violations_du.push_back(-z_U_trial_val[i]);
	//printf("Upper Bound Mult. no. i=%d invalid: delta_u=%f\n", i+z_U_ItVec_idx, z_U_val[i]);
      }
    }

    //    if (x_bound_violations_idx.empty() || z_L_bound_violations_idx.empty() || z_U_bound_violations_idx.empty()) {
    if (x_bound_violations_idx.empty()) {
      return false;
    }
    else {
      return true;
    }
  }
}
