// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpUtils.hpp"
#include "IpAlgBuilder.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "AmplTNLP.hpp"
#include "IpTNLPAdapter.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{
  int RunIpoptAlgorithm(const SmartPtr<const Journalist> jnlst, int argv, char**argc)
  {
    jnlst->Printf(J_ERROR, J_MAIN, "\n\n\n*************************************************************\n");
    jnlst->Printf(J_ERROR, J_MAIN, "*** Running Ipopt with AMPL Model ***************************\n");
    jnlst->Printf(J_ERROR, J_MAIN, "*************************************************************\n\n\n");

    // Create the original nlp and the IpoptNLP
    //    NLP* orig_nlp = new __TestNLP();
    SmartPtr<AmplTNLP> ampl_nlp = new AmplTNLP(jnlst, argc);
    SmartPtr<TNLPAdapter> orig_nlp = new TNLPAdapter(GetRawPtr(ampl_nlp));
    SmartPtr<IpoptNLP> ip_nlp = new OrigIpoptNLP(jnlst, GetRawPtr(orig_nlp));

    // Create the IpoptData
    SmartPtr<IpoptData> ip_data = new IpoptData();

    // Create the IpoptCalculators
    SmartPtr<IpoptCalculatedQuantities> ip_cq
    = new IpoptCalculatedQuantities(ip_nlp, ip_data);

    // Get the options
    // TODO: Need method in AmplTNLP that can fill options with the
    //       AMPL user options
    SmartPtr<OptionsList> options = new OptionsList;
    FILE* fp_options = fopen("PARAMS.DAT", "r");
    if (fp_options) {
      // PARAMS.DAT exists, read the content
      options->ReadFromFile(*jnlst, fp_options);
      fclose(fp_options);
      fp_options=NULL;
    }

    // Create the complete Algorithm object
    SmartPtr<IpoptAlgorithm> alg = AlgorithmBuilder(*jnlst, *options, "");

    try {
      // Initialize the algorithm and parse the options
      alg->Initialize(*jnlst, *ip_nlp, *ip_data, *ip_cq, *options, "");
    }
    catch(IpoptException& exc) {
      exc.ReportException(*jnlst);
      return -1;
    }

    if( jnlst->ProduceOutput(J_DETAILED, J_MAIN) ) {
      // Print out the options (including the number of times they were used

      std::string liststr;
      options->PrintList(liststr);
      jnlst->Printf(J_DETAILED, J_MAIN,
                    "\nList of options:\n\n%s", liststr.c_str());
    }

    int retval=-1;

    IpoptAlgorithm::SolverReturn status;
    std::string message("Optimal solution found.");
    try {
      // Run the algorithm
      status = alg->Optimize();
      retval = 0;
      if (status == IpoptAlgorithm::MAXITER_EXCEEDED) {
        message = "Maximal number of iterations exceeded.\n";
        retval = 2;
      }
    }
    catch(IpoptException& exc) {
      exc.ReportException(*jnlst);
      message = exc.Message();
      retval = 1;
    }

    // Write the .sol file
    Index n, m, nnz_jac_g, nnz_h_lag;
    ampl_nlp->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag);
    Number* xsol = new Number[n];
    Number* ysol = new Number[m];
    orig_nlp->ResortX(*ip_data->curr_x(), xsol);
    orig_nlp->ResortG(*ip_data->curr_y_c(), *ip_data->curr_y_d(), ysol);

    message = " \nEXIT: " + message;
    ampl_nlp->write_solution_file(message.c_str(), xsol, ysol);
    delete[] xsol;
    delete[] ysol;

    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "\nNumber of Iterations    = %d\n", ip_data->iter_count());
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Optimal Objective Value = %23.16e\n", ip_cq->curr_f());
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Primal Infeasibility    = %23.16e\n",
                  ip_cq->curr_primal_infeasibility(IpoptCalculatedQuantities::NORM_MAX));
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Dual Infeasibility      = %23.16e\n",
                  ip_cq->curr_dual_infeasibility(IpoptCalculatedQuantities::NORM_MAX));
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Complementarity         = %23.16e\n",
                  ip_cq->curr_complementarity(0., IpoptCalculatedQuantities::NORM_MAX));

    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "\nNumber of objective function evaluations             = %d\n",
                  ip_nlp->f_evals());
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Number of equality constraint evaluations            = %d\n",
                  ip_nlp->c_evals());
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Number of inequality constraint evaluations          = %d\n",
                  ip_nlp->d_evals());
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Number of equality constraint Jacobian evaluations   = %d\n",
                  ip_nlp->jac_c_evals());
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Number of inequality constraint Jacobian evaluations = %d\n",
                  ip_nlp->jac_d_evals());
    jnlst->Printf(J_SUMMARY, J_SOLUTION,
                  "Number of Lagrangian Hessian evaluations             = %d\n",
                  ip_nlp->h_evals());


    jnlst->PrintVector(J_VECTOR, J_SOLUTION, "x", *ip_data->curr_x());
    jnlst->PrintVector(J_VECTOR, J_SOLUTION, "y_c", *ip_data->curr_y_c());
    jnlst->PrintVector(J_VECTOR, J_SOLUTION, "y_d", *ip_data->curr_y_d());
    jnlst->PrintVector(J_VECTOR, J_SOLUTION, "z_L", *ip_data->curr_z_L());
    jnlst->PrintVector(J_VECTOR, J_SOLUTION, "z_U", *ip_data->curr_z_U());
    jnlst->PrintVector(J_VECTOR, J_SOLUTION, "v_L", *ip_data->curr_v_L());
    jnlst->PrintVector(J_VECTOR, J_SOLUTION, "v_U", *ip_data->curr_v_U());

    return retval;
  }

} // namespace Ipopt

int main(int argv, char** argc)
{

  using namespace Ipopt;

  Index retValue = 0;
  Ipopt::SmartPtr<Journalist> jnlst = new Journalist();

# ifdef IP_DEBUG

  DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst));
# endif

  //  Journal* jrnl = jnlst->AddJournal("ConsoleStdOut", "stdout", J_ALL);
  Journal* jrnl = jnlst->AddJournal("ConsoleStdOut", "stdout", J_SUMMARY);
  jrnl->SetPrintLevel(J_DBG, J_NONE);

# ifdef IP_DEBUG
  jrnl = jnlst->AddJournal("Debug", "debug.out", J_DETAILED);
  jrnl->SetPrintLevel(J_DBG, J_ALL);
# endif

  //jrnl = jnlst->AddJournal("All", "all.out", J_MOREVECTOR);
  //jrnl = jnlst->AddJournal("All", "all.out", J_VECTOR);
  //jrnl = jnlst->AddJournal("All", "all.out", J_MOREDETAILED);

  //jrnl = jnlst->AddJournal("LineSearch", "linesearch.out", J_NONE);
  //jrnl->SetPrintLevel(J_LINE_SEARCH, J_ALL);

  retValue = Ipopt::RunIpoptAlgorithm(ConstPtr(jnlst), argv, argc);

  return retValue;
}
