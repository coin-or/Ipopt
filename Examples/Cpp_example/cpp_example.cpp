// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpUtils.hpp"
#include "IpJournalist.hpp"
#include "IpSmartPtr.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpIpoptNLP.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptAlg.hpp"
#include "IpAlgBuilder.hpp"

#include "MyNLP.hpp"

using namespace Ipopt;

int main(int argv, char* argc[])
{
  int retValue = 0;
  // Create the journalist for output reporting
  SmartPtr<Journalist> jnlst = new Journalist();

  // Create the output journal for the console
  Journal* jrnl = jnlst->AddJournal("ConsoleStdOut", "stdout", J_SUMMARY);
  jrnl->SetPrintLevel(J_SOLUTION, J_ALL);

  // If you want to log more detailed information to a file, create
  // another journal
  jrnl = jnlst->AddJournal("Detailed", "detailed.out", J_DETAILED);
  jrnl->SetPrintLevel(J_SOLUTION, J_ALL);

  jrnl = jnlst->AddJournal("All", "all.out", J_ALL);


  // start a try catch block for any exceptions
  try
    {
      // create an instance of MyNLP
      SmartPtr<TNLP> my_nlp = new MyNLP();

      // MyNLP is a TNLP (not an NLP), so we need an adapter
      SmartPtr<NLP> nlp = new TNLPAdapter(GetRawPtr(my_nlp));

      // Read the options from IPOPT (TODO, make this easier)
      SmartPtr<OptionsList> options = new OptionsList();
      FILE* fp_options = fopen("PARAMS.DAT", "r");
      if (fp_options) {
	// PARAMS.DAT exists, read the content
	options->ReadFromFile(*jnlst, fp_options);
	fclose(fp_options);
	fp_options=NULL;
      }
      
      // Now create some objects needed by the algorithm
      SmartPtr<IpoptNLP> ip_nlp = new OrigIpoptNLP(ConstPtr(jnlst), nlp);
      SmartPtr<IpoptData> ip_data = new IpoptData();
      SmartPtr<IpoptCalculatedQuantities> ip_cq 
	= new IpoptCalculatedQuantities(ip_nlp, ip_data);

      // Create the algorithm object
      SmartPtr<IpoptAlgorithm> alg = AlgorithmBuilder(*jnlst, *options, "");
      alg->Initialize(*jnlst, *ip_nlp, *ip_data, *ip_cq, *options, "");

      IpoptAlgorithm::SolverReturn status
	= alg->Optimize();

      if (status == IpoptAlgorithm::SUCCESS) {
	jnlst->Printf(J_SUMMARY, J_SOLUTION, "\n\nOptimal solution found! \n");
	jnlst->Printf(J_SUMMARY, J_SOLUTION, "Optimal Objective Value = %.16E\n", ip_cq->curr_f());
	jnlst->Printf(J_SUMMARY, J_SOLUTION, "Number of Iterations = %d\n", ip_data->iter_count());
	jnlst->Printf(J_VECTOR, J_SOLUTION, "\n");
	jnlst->PrintVector(J_VECTOR, J_SOLUTION, "x", *ip_data->curr_x());
	jnlst->Printf(J_VECTOR, J_SOLUTION, "\n");
	jnlst->PrintVector(J_VECTOR, J_SOLUTION, "y_c", *ip_data->curr_y_c());
	jnlst->Printf(J_VECTOR, J_SOLUTION, "\n");
	jnlst->PrintVector(J_VECTOR, J_SOLUTION, "y_d", *ip_data->curr_y_d());
	jnlst->Printf(J_VECTOR, J_SOLUTION, "\n");
	jnlst->PrintVector(J_VECTOR, J_SOLUTION, "z_L", *ip_data->curr_z_L());
	jnlst->Printf(J_VECTOR, J_SOLUTION, "\n");
	jnlst->PrintVector(J_VECTOR, J_SOLUTION, "z_U", *ip_data->curr_z_U());
	jnlst->Printf(J_VECTOR, J_SOLUTION, "\n");
	jnlst->PrintVector(J_VECTOR, J_SOLUTION, "v_L", *ip_data->curr_v_L());
	jnlst->Printf(J_VECTOR, J_SOLUTION, "\n");
	jnlst->PrintVector(J_VECTOR, J_SOLUTION, "v_U", *ip_data->curr_v_U());
	jnlst->Printf(J_VECTOR, J_SOLUTION, "\n");
	
	retValue = 0;
      }
      else {
	jnlst->Printf(J_ERROR, J_MAIN, "Sorry, things failed. Check the detailed output files\n");
	retValue = 1;
      }
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
