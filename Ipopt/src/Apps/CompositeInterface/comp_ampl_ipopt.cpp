// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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
#include "CompositeNLP.hpp"

extern "C"
{
#include "asl.h"
#include "asl_pfgh.h"
}

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{
  void CreateLinkingEqns(std::vector<SmartPtr<AmplTNLP> >& nlps,
                         std::vector<SmartPtr<AmplSuffixHandler> >& suffix_handlers,
                         std::vector<SmartPtr<VectorSpace> >& linking_eqn_spaces,
                         std::vector<SmartPtr<Matrix> >& Jx_linking_eqns,
                         std::vector<SmartPtr<Matrix> >& Jq_linking_eqns,
                         SmartPtr<VectorSpace>& q_space);

  int RunIpoptAlgorithm(const SmartPtr<const Journalist>& jnlst,
                        int argv, char**argc)
  {
    jnlst->Printf(J_ERROR, J_MAIN, "\n\n\n*************************************************************\n");
    jnlst->Printf(J_ERROR, J_MAIN, "*** Running Ipopt with AMPL Model ***************************\n");
    jnlst->Printf(J_ERROR, J_MAIN, "*************************************************************\n\n\n");

    // For test purposes, every ampl model must have 5 variables and the
    // last variable in all of them is the linked variable
    std::vector<SmartPtr<AmplTNLP> > ampl_nlps;
    std::vector<SmartPtr<NLP> > nlps;
    std::vector<SmartPtr<AmplSuffixHandler> > suffix_handlers;
    // Process the command line options, each one should be a stub file...
    for (Index i=1; i<argv; i++) {
      char** argc_i = new char*[3];
      argc_i[0] = new char[6];
      strcpy(argc_i[0],"bogus");
      Index len = strlen(argc[i]);
      argc_i[1] = new char[len+1];
      strcpy(argc_i[1], argc[i]);
      argc_i[2] = NULL;

      // We use Ampl Suffixes to describe which variables are common
      // If the suffix "common_idx" is non-zero for a particular variable, that variable
      // is linked to the corresponding common variable. i.e. x[10].common_idx = 4
      // adds the linking eqn x[10] - q[4] = 0 (actually x[10] - q[3] :)
      SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
      suffix_handler->AddAvailableSuffix("common_idx", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
      suffix_handlers.push_back(GetRawPtr(suffix_handler));

      SmartPtr<AmplTNLP> ampl_nlp_i = new AmplTNLP(jnlst, argc_i, suffix_handler);
      ampl_nlps.push_back(GetRawPtr(ampl_nlp_i));

      SmartPtr<NLP> nlp_i = new TNLPAdapter(GetRawPtr(ampl_nlp_i));
      nlps.push_back(nlp_i);
    }

    SmartPtr<VectorSpace> q_space;
    std::vector<SmartPtr<VectorSpace> > linking_eqn_spaces;
    std::vector<SmartPtr<Matrix> > Jx_linking_eqns;
    std::vector<SmartPtr<Matrix> > Jq_linking_eqns;
    CreateLinkingEqns(ampl_nlps, suffix_handlers, linking_eqn_spaces, Jx_linking_eqns, Jq_linking_eqns, q_space);

    // Create the composite nlp and the IpoptNLP
    SmartPtr<NLP> nlp = new CompositeNLP(nlps, q_space,
                                         linking_eqn_spaces,
                                         Jx_linking_eqns, Jq_linking_eqns);

    SmartPtr<NLPScalingObject> nlp_scaling = new NoNLPScalingObject();
    SmartPtr<IpoptNLP> ip_nlp = new OrigIpoptNLP(jnlst, nlp, nlp_scaling);

    // Create the IpoptData
    SmartPtr<IpoptData> ip_data = new IpoptData();
    ip_data->Set_iter_count(0);

    // Create the IpoptCalculators
    SmartPtr<IpoptCalculatedQuantities> ip_cq
    = new IpoptCalculatedQuantities(ip_nlp, ip_data);

    // Get the options
    SmartPtr<OptionsList> options = new OptionsList;
    FILE* fp_options = fopen("PARAMS.DAT", "r");
    if (fp_options) {
      // PARAMS.DAT exists, read the content
      options->ReadFromFile(*jnlst, fp_options);
      fclose(fp_options);
      fp_options=NULL;
    }

    // Create the complete Algorithm object
    SmartPtr<IpoptAlgorithm> alg = AlgorithmBuilder::BuildBasicAlgorithm(*jnlst, *options, "");

    // Initialize the algorithm and parse the options
    alg->Initialize(*jnlst, *ip_nlp, *ip_data, *ip_cq, *options, "");

    // Run the algorithm
    SolverReturn status = alg->Optimize();

    int retval=-1;

    if (status == SUCCESS) {
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "\n\nOptimal solution found! \n");
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "Optimal Objective Value = %.16E\n", ip_cq->curr_f());
      jnlst->Printf(J_SUMMARY, J_SOLUTION, "Number of Iterations = %d\n", ip_data->iter_count());
      ip_data->curr()->Print(jnlst, J_VECTOR, J_SOLUTION, "sol");

      retval = 0;
    }
    else {
      jnlst->Printf(J_ERROR, J_MAIN, "Sorry, things failed ?!?!\n");
      retval = 1;
    }

    printf("IN AMPL INTERFACE, THE RESULT NEEDS TO BE STORED IN AMPL STRUCTURE!\n");

    return retval;
  }

  DECLARE_STD_EXCEPTION(AMPL_SUFFIX_ERROR);

  void CreateLinkingEqns(std::vector<SmartPtr<AmplTNLP> >& nlps,
                         std::vector<SmartPtr<AmplSuffixHandler> >& suffix_handlers,
                         std::vector<SmartPtr<VectorSpace> >& linking_eqn_spaces,
                         std::vector<SmartPtr<Matrix> >& Jx_linking_eqns,
                         std::vector<SmartPtr<Matrix> >& Jq_linking_eqns,
                         SmartPtr<VectorSpace>& q_space)
  {
    Index n_nlps = nlps.size();
    DBG_ASSERT(n_nlps == (Index)suffix_handlers.size());

    // first, count the number of required q's
    Index q_dim = 0;
    for (Index nlp_idx = 0; nlp_idx<n_nlps; nlp_idx++) {
      const Index* dp = suffix_handlers[nlp_idx]->GetIntegerSuffixValues("common_idx", AmplSuffixHandler::Variable_Source);
      ASSERT_EXCEPTION(dp, AMPL_SUFFIX_ERROR, "Error in Ampl Suffixes, no linking equations specified");
      SmartPtr<AmplTNLP> ampl_nlp_i = nlps[nlp_idx];
      ASL_pfgh* asl = ampl_nlp_i->AmplSolverObject();
      for (Index i=0; i<n_var; i++) {
        q_dim = Max(q_dim, dp[i]);
      }
    }

    q_space = new DenseVectorSpace(q_dim);

    bool* q_test = new bool[q_dim];

    for (Index nlp_idx = 0; nlp_idx<n_nlps; nlp_idx++) {
      const Index* dp = suffix_handlers[nlp_idx]->GetIntegerSuffixValues("common_idx", AmplSuffixHandler::Variable_Source);
      SmartPtr<AmplTNLP> ampl_nlp_i = nlps[nlp_idx];

      // count the number of linking equations
      ASL_pfgh* asl = ampl_nlp_i->AmplSolverObject();
      Index n_link_eqns = 0;
      for (Index i=0; i<n_var; i++) {
        if (dp[i] != 0) {
          n_link_eqns++;
        }
      }

      // Build the Jx_linking_eqns[i]
      Index* iRows = new Index[n_link_eqns];
      Index* jCols = new Index[n_link_eqns];
      Index curr_eqn_idx = 0;
      for (Index i=0; i<n_var; i++) {
        if (dp[i]) {
          iRows[curr_eqn_idx] = curr_eqn_idx + 1;
          jCols[curr_eqn_idx] = i+1;
          curr_eqn_idx++;
        }
      }

      SmartPtr<GenTMatrixSpace> space = new GenTMatrixSpace(n_link_eqns, n_var, n_link_eqns, iRows, jCols);
      SmartPtr<GenTMatrix> Jx_i = space->MakeNewGenTMatrix();
      Number* values = Jx_i->Values();
      for (Index i=0; i<n_link_eqns; i++) {
        values[i] = 1.0;
      }

      // Build the Jq_linking_eqns
      curr_eqn_idx = 0;
      for (Index i=0; i<n_var; i++) {
        if (dp[i]) {
          iRows[curr_eqn_idx] = curr_eqn_idx+1;
          jCols[curr_eqn_idx] = dp[i];
          q_test[dp[i]-1] = true;
          curr_eqn_idx++;
          ASSERT_EXCEPTION(dp[i] > 0, AMPL_SUFFIX_ERROR, "common_idx must be larger than zero (one based index)");
        }
      }


      space = new GenTMatrixSpace(n_link_eqns, q_dim, n_link_eqns, iRows, jCols);
      SmartPtr<GenTMatrix> Jq_i = space->MakeNewGenTMatrix();
      values = Jq_i->Values();
      for (Index i=0; i<n_link_eqns; i++) {
        values[i] = -1.0;
      }


      linking_eqn_spaces.push_back(new DenseVectorSpace(n_link_eqns));

      Jx_linking_eqns.push_back(GetRawPtr(Jx_i));
      Jq_linking_eqns.push_back(GetRawPtr(Jq_i));

      delete [] iRows;
      delete [] jCols;
    }

    // check the q_test for singularity
    for (Index i=0; i<q_dim; i++) {
      ASSERT_EXCEPTION(q_test[i], AMPL_SUFFIX_ERROR, "One of the common variables q is not linked to any x variables");
    }
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

  SmartPtr<Journal> jrnl =
    jnlst->AddFileJournal("ConsoleStdOut", "stdout", J_SUMMARY);
  jrnl->SetPrintLevel(J_DBG, J_NONE);

  jrnl = jnlst->AddFileJournal("Debug", "debug.out", J_DETAILED);
  jrnl->SetPrintLevel(J_DBG, J_ALL);

  jrnl = jnlst->AddFileJournal("All", "all.out", J_ALL);

  try {
    //***
    // Setup the Journalist
    //***

    retValue = Ipopt::RunIpoptAlgorithm(ConstPtr(jnlst), argv, argc);
  }
  catch(IpoptException& exc) {
    exc.ReportException(*jnlst);
    retValue = 1;
  }
  catch(...) {
    IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
    exc.ReportException(*jnlst);
    retValue = 2;
  }

  return retValue;
}
