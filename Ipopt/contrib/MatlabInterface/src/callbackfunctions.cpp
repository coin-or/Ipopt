// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 18, 2008

#include "callbackfunctions.hpp"
#include "matlabexception.hpp"

// Functions for class CallbackFunctions.
// -----------------------------------------------------------------
CallbackFunctions::CallbackFunctions (const mxArray* ptr) 
  : objfunc(0), gradfunc(0), constraintfunc(0), jacobianfunc(0), 
  jacstrucfunc(0), hessianfunc(0), hesstrucfunc(0), iterfunc(0) {
  const mxArray* p;  // A pointer to a MATLAB array.

  // Check whether we are provided with a structure array.
  if (!mxIsStruct(ptr))
    throw MatlabException("The second input must be a STRUCT");

  // Get the function handle for computing the objective.
  p = mxGetField(ptr,0,"objective");
  if (!p)
    throw MatlabException("You must specify a callback routine for \
computing the value of the objective function");
  if (mxIsEmpty(p) || !isFunctionHandle(p))
    throw MatlabException("You did not provide a valid function handle for \
computing the value of the objective function");
  objfunc = new MatlabFunctionHandle(p);

  // Get the function handle for computing the gradient.
  p = mxGetField(ptr,0,"gradient");
  if (!p)
    throw MatlabException("You must specify a callback routine for \
computing the gradient of the objective");
  if (mxIsEmpty(p) || !isFunctionHandle(p))
    throw MatlabException("You did not provide a valid function handle for \
computing the gradient of the objective");
  gradfunc = new MatlabFunctionHandle(p);  

  // Get the function handle for computing the constraints, if such a
  // function was specified.
  p = mxGetField(ptr,0,"constraints");
  if (p) {
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the response of the constraints");
    constraintfunc = new MatlabFunctionHandle(p);      
  }
  else
    constraintfunc = new MatlabFunctionHandle();

  // Get the function handle for computing the Jacobian. This function
  // is necessary if there are constraints.
  p = mxGetField(ptr,0,"jacobian");
  if (p) {
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the first derivatives (Jacobian) of the constraints");
    jacobianfunc = new MatlabFunctionHandle(p);      
  }
  else {
    if (*constraintfunc)
      throw MatlabException("You must provide a function that returns the \
first derivatives (Jacobian) of the constraints");
    jacobianfunc = new MatlabFunctionHandle();
  }

  // Get the function handle for computing the sparsity structure of
  // the Jacobian. This function is necessary if the Jacobian is being
  // computed.
  p = mxGetField(ptr,0,"jacobianstructure");
  if (p) { 
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the sparsity structure of the Jacobian");
    jacstrucfunc = new MatlabFunctionHandle(p);      
  }
  else {
    if (*jacobianfunc)
      throw MatlabException("You must provide a function that returns the \
sparsity structure of the Jacobian");
    jacstrucfunc = new MatlabFunctionHandle();
  }

  // Get the function handle for computing the Hessian. This function
  // is always optional.
  p = mxGetField(ptr,0,"hessian");
  if (p) {
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the Hessian of the Lagrangian");
    hessianfunc = new MatlabFunctionHandle(p);      
  }
  else 
    hessianfunc = new MatlabFunctionHandle();

  // Get the function handle for computing the sparsity structure of
  // the Hessian of the Lagrangian. This function is necessary if the
  // Hessian is being computed.
  p = mxGetField(ptr,0,"hessianstructure");
  if (p) { 
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for computing the sparsity structure of the Hessian");
    hesstrucfunc = new MatlabFunctionHandle(p);      
  }
  else {
    if (*hessianfunc)
      throw MatlabException("You must provide a function that returns the \
sparsity structure of the Hessian");
    hesstrucfunc = new MatlabFunctionHandle();
  }  

  // Get the iterative callback function handle. This function is
  // always optional.
  p = mxGetField(ptr,0,"iterfunc");
  if (p) { 
    if (mxIsEmpty(p) || !isFunctionHandle(p))
      throw MatlabException("You did not provide a valid function handle \
for the iterative callback");
    iterfunc = new MatlabFunctionHandle(p);      
  }
  else
    iterfunc = new MatlabFunctionHandle();
}

CallbackFunctions::~CallbackFunctions() {
  if (objfunc)        delete objfunc;
  if (gradfunc)       delete gradfunc;
  if (constraintfunc) delete constraintfunc;
  if (jacobianfunc)   delete jacobianfunc;
  if (jacstrucfunc)   delete jacstrucfunc;
  if (hessianfunc)    delete hessianfunc;
  if (hesstrucfunc)   delete hesstrucfunc;
  if (iterfunc)       delete iterfunc;
}

double CallbackFunctions::computeObjective (const Iterate& x) const {
  double         f;  // The return value.
  bool           success;
  const mxArray* inputs[1];
  mxArray*       outputs[1];

  // Call the MATLAB callback function.
  inputs[0] = x;
  success = objfunc->evaluate(1,1,inputs,outputs);
  if (!success)
    throw MatlabException("There was an error when executing the objective \
callback function");

  // Get the output from the MATLAB callback function, which is the
  // value of the objective function at x.
  mxArray* ptr = outputs[0];
  if (!mxIsDouble(ptr) || mxIsComplex(ptr) || mxGetNumberOfElements(ptr) != 1)
    throw MatlabException("The first return value of the objective \
callback function must be a real double scalar");
  if (mxIsSparse(ptr)) {
    // convert sparse objective (unlikely but possible) to full
    mexCallMATLAB(1, &ptr, 1, outputs, "full");
    mxDestroyArray(outputs[0]);
  }
  f = *mxGetPr(ptr);

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);  

  return f;
}

void CallbackFunctions::computeGradient (const Iterate& x, double* g) const {
  bool           success;
  const mxArray* inputs[1];
  mxArray*       outputs[1];

  // Call the MATLAB callback function.
  inputs[0] = x;
  success = gradfunc->evaluate(1,1,inputs,outputs);
  if (!success)
    throw MatlabException("There was an error when executing the \
gradient callback function");

  // Get the output from the MATLAB callback function, which is the
  // value of the gradient of the objective function at x.
  mxArray* ptr = outputs[0];
  if (!mxIsCell(ptr) && (!mxIsDouble(ptr) || mxIsComplex(ptr)))
    throw MatlabException("The gradient callback must return a real double vector");
  if (mxIsSparse(ptr)) {
    // convert sparse gradient to full (simplest method, not fastest)
    mexCallMATLAB(1, &ptr, 1, outputs, "full");
    mxDestroyArray(outputs[0]);
  }
  Iterate grad(ptr);
  if (numvars(x) != numvars(grad))
    throw MatlabException("Invalid gradient passed back from MATLAB \
routine");
  grad.copyto(g);

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);
}

void CallbackFunctions::computeConstraints(const Iterate& x, int m, double* c) const {
  bool           success;
  const mxArray* inputs[1];
  mxArray*       outputs[1];

  // Call the MATLAB callback function.
  inputs[0] = x;
  success = constraintfunc->evaluate(1,1,inputs,outputs);
  if (!success)
    throw MatlabException("There was an error when executing the \
constraints callback function");

  // Get the output from the MATLAB callback function, which is the
  // value of vector-valued constraint function at x.
  mxArray* ptr = outputs[0];
  if ((unsigned) m != mxGetNumberOfElements(ptr))
    throw MatlabException("Invalid constraint values passed back from \
Matlab routine");
  if (!mxIsDouble(ptr) || mxIsComplex(ptr))
    throw MatlabException("The constraint callback must return a real double vector");
  if (mxIsSparse(ptr)) {
    // convert sparse constraint vector (unlikely but possible) to full
    mexCallMATLAB(1, &ptr, 1, outputs, "full");
    mxDestroyArray(outputs[0]);
  }
  copymemory(mxGetPr(ptr),c,m);

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);
}

SparseMatrix* CallbackFunctions::getJacobianStructure (int n, int m) const {
  const mxArray* inputs[1];
  mxArray*       outputs[1];
  bool           success;

  // Call the MATLAB callback function.
  success = jacstrucfunc->evaluate(0,1,inputs,outputs);
  if (!success)
    throw MatlabException("There was an error when getting the structure \
of the Jacobian matrix from MATLAB");

  // Get the output from the MATLAB callback function, which is the
  // sparse matrix specifying the structure of the Jacobian.
  mxArray* ptr = outputs[0];
  if (!mxIsSparse(ptr) || mxIsComplex(ptr))
    throw MatlabException("Jacobian must be a real sparse matrix");
  if ((int) mxGetM(ptr) != m || (int) mxGetN(ptr) != n || 
      !SparseMatrix::inIncOrder(ptr))
    throw MatlabException("Jacobian must be an m x n sparse matrix with \
row indices in increasing order, where m is the number of constraints and \
n is the number of variables");
  if (!mxIsDouble(ptr)) {
    // convert non-double Jacobian structure to double
    mexCallMATLAB(1, &ptr, 1, outputs, "double");
    mxDestroyArray(outputs[0]);
  }
  SparseMatrix* J = new SparseMatrix(ptr);  // The return value.

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);

  return J;
}

SparseMatrix* CallbackFunctions::getHessianStructure (int n) const {
  const mxArray* inputs[1];
  mxArray*       outputs[1];
  bool           success;

  // Call the MATLAB callback function.
  success = hesstrucfunc->evaluate(0,1,inputs,outputs);
  if (!success)
    throw MatlabException("There was an error when getting the structure \
of the Hessian matrix from MATLAB");

  // Get the output from the MATLAB callback function, which is the
  // sparse matrix specifying the structure of the Hessian.
  mxArray* ptr = outputs[0];
  if (!mxIsSparse(ptr) || mxIsComplex(ptr))
    throw MatlabException("Hessian must be a real sparse matrix");
  if ((int) mxGetM(ptr) != n || (int) mxGetN(ptr) != n || 
      !SparseMatrix::isLowerTri(ptr) || !SparseMatrix::inIncOrder(ptr))
    throw MatlabException("Hessian must be an n x n sparse, symmetric and \
lower triangular matrix with row indices in increasing order, where n is \
the number of variables");
  if (!mxIsDouble(ptr)) {
    // convert non-double Hessian structure to double
    mexCallMATLAB(1, &ptr, 1, outputs, "double");
    mxDestroyArray(outputs[0]);
  }
  SparseMatrix* H = new SparseMatrix(ptr);  // The return value.

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);

  return H;
}

void CallbackFunctions::computeJacobian (int m, const Iterate& x, 
					 SparseMatrix& J) const {
  bool           success;
  const mxArray* inputs[1];
  mxArray*       outputs[1];

  // Call the MATLAB callback function.
  inputs[0] = x;
  success = jacobianfunc->evaluate(1,1,inputs,outputs);
  if (!success)
    throw MatlabException("There was an error when executing the \
Jacobian callback function");

  // Get the output from the MATLAB callback function, which is the
  // sparse matrix specifying the value the Jacobian.
  mxArray* ptr = outputs[0];
  if (!mxIsSparse(ptr) || mxIsComplex(ptr))
    throw MatlabException("Jacobian must be a real sparse matrix");
  if ((int) mxGetM(ptr) != m || (int) mxGetN(ptr) != numvars(x) || 
      !SparseMatrix::inIncOrder(ptr))
    throw MatlabException("Jacobian must be an m x n sparse matrix with \
row indices in increasing order, where m is the number of constraints and \
n is the number of variables");
  if (!mxIsDouble(ptr)) {
    // convert non-double Jacobian to double
    mexCallMATLAB(1, &ptr, 1, outputs, "double");
    mxDestroyArray(outputs[0]);
  }
  SparseMatrix Jnew(ptr);
  Jnew.copyto(J);

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);
}

void CallbackFunctions::computeHessian (const Iterate& x, double sigma, int m, 
					const double* lambda, SparseMatrix& H) const {
  bool           success;
  const mxArray* inputs[3];
  mxArray*       outputs[1];

  // Create the input arguments to the MATLAB routine, sigma and lambda.
  mxArray* psigma  = mxCreateDoubleScalar(sigma);
  mxArray* plambda = mxCreateDoubleMatrix(m,1,mxREAL);
  copymemory(lambda,mxGetPr(plambda),m);

  // Call the MATLAB callback function.
  inputs[0] = x;
  inputs[1] = psigma;
  inputs[2] = plambda;
  success = hessianfunc->evaluate(3,1,inputs,outputs);
  if (!success)
    throw MatlabException("There was an error when executing the Hessian \
callback function");

  // Get the output from the MATLAB callback function, which is the
  // sparse matrix specifying the value the Hessian.
  mxArray* ptr = outputs[0];
  if (!mxIsSparse(ptr) || mxIsComplex(ptr))
    throw MatlabException("Hessian must be a real sparse matrix");
  if ((int) mxGetM(ptr) != numvars(x) || (int) mxGetN(ptr) != numvars(x) || 
      !SparseMatrix::isLowerTri(ptr) || !SparseMatrix::inIncOrder(ptr))
    throw MatlabException("Hessian must be an n x n sparse, symmetric and \
lower triangular matrix with row indices in increasing order, where n is \
the number of variables");
  if (!mxIsDouble(ptr)) {
    // convert non-double Hessian to double
    mexCallMATLAB(1, &ptr, 1, outputs, "double");
    mxDestroyArray(outputs[0]);
  }
  SparseMatrix Hnew(ptr);
  Hnew.copyto(H);

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);
  mxDestroyArray(psigma);
  mxDestroyArray(plambda);
}

bool CallbackFunctions::iterCallback (int t, double f, 
				      double inf_pr, double inf_du, 
				      double mu, double d_norm,
				      double regularization_size,
				      double alpha_du, double alpha_pr,
				      int ls_trials, const Ipopt::IpoptData* ip_data, 
				      Ipopt::IpoptCalculatedQuantities* ip_cq,
				      int n) const {
  bool           success;
  const mxArray* inputs[3];
  mxArray*       outputs[1];

  // Create the input arguments to the MATLAB routine.
  mxArray* pt = mxCreateDoubleScalar(t);
  mxArray* pf = mxCreateDoubleScalar(f);

  // Create structure to hold extra IPOPT variables
  const char* varfields[9];
  varfields[0] = "x";
  varfields[1] = "inf_pr";
  varfields[2] = "inf_du";
  varfields[3] = "mu";
  varfields[4] = "d_norm";
  varfields[5] = "regularization_size";
  varfields[6] = "alpha_du";
  varfields[7] = "alpha_pr";
  varfields[8] = "ls_trials"; 
  mxArray *varStruct = mxCreateStructMatrix(1,1,9,varfields);
  mxSetField(varStruct,0,"inf_pr",mxCreateDoubleScalar(inf_pr));
  mxSetField(varStruct,0,"inf_du",mxCreateDoubleScalar(inf_du));
  mxSetField(varStruct,0,"mu",mxCreateDoubleScalar(mu));
  mxSetField(varStruct,0,"d_norm",mxCreateDoubleScalar(d_norm));
  mxSetField(varStruct,0,"regularization_size",mxCreateDoubleScalar(regularization_size));  
  mxSetField(varStruct,0,"alpha_du",mxCreateDoubleScalar(alpha_du));
  mxSetField(varStruct,0,"alpha_pr",mxCreateDoubleScalar(alpha_pr));
  mxSetField(varStruct,0,"ls_trials",mxCreateDoubleScalar(ls_trials));
  
  //The following code translates IPOPT's NLP to the Original NLP, so we can extract x
  //Original code by Steven Dirske, Stefan Vigerske [GAMS]
  Ipopt::TNLPAdapter* tnlp_adapter = NULL;
  if(ip_cq != NULL) {
      Ipopt::OrigIpoptNLP* orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
      if(orignlp != NULL) 
          tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
  }
  //If we successfully converted the NLP, extract x [ResortX auto sorts and fills in fixed vars]
  if(tnlp_adapter != NULL && ip_data != NULL && IsValid(ip_data->curr())) {
      mxSetField(varStruct,0,"x",mxCreateDoubleMatrix(n,1,mxREAL)); //ip_data->curr()->x()->Dim()+1
      tnlp_adapter->ResortX(*ip_data->curr()->x(),mxGetPr(mxGetField(varStruct,0,"x")));
  }
  
  // Call the MATLAB callback function.
  inputs[0] = pt;
  inputs[1] = pf;
  inputs[2] = varStruct;
  success = iterfunc->evaluate(3,1,inputs,outputs);
  if (!success)
    throw MatlabException("There was an error when executing the iterative \
callback function");

  // Get the output from the MATLAB callback function, which is a
  // boolean value telling whether or not IPOPT should continue.
  mxArray* ptr = outputs[0];
  if (!mxIsLogicalScalar(ptr))
    throw MatlabException("The return value for the iterative callback must \
either be TRUE or FALSE");
  bool b = mxIsLogicalScalarTrue(ptr);

  // Free the dynamically allocated memory.
  mxDestroyArray(ptr);
  mxDestroyArray(pt);
  mxDestroyArray(pf);
  mxDestroyArray(varStruct);

  return b;
}
