// Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         September 15, 2008

#include "ipoptoptions.hpp"
#include "matlabexception.hpp"

using Ipopt::IpoptApplication;
using Ipopt::SmartPtr;
using Ipopt::IsValid;
using Ipopt::RegisteredOption;
using Ipopt::SmartPtr;
using Ipopt::Snprintf;

// need 

// Function definitions for class IpoptOptions.
// -----------------------------------------------------------------
IpoptOptions::IpoptOptions (IpoptApplication& app, const mxArray* ptr) 
  : app(app) {

  if (ptr != NULL) { // avoid segfault if no ipopt field in options
    // Check to make sure the MATLAB array is a structure array.
    if (!mxIsStruct(ptr))
      throw MatlabException("The OPTIONS input must be a structure array; \
type HELP STRUCT in the MATLAB console for more information");

    // Each field in the structure array should correspond to an option
    // in IPOPT. Repeat for each field.
    int n = mxGetNumberOfFields(ptr);
    for (int i = 0; i < n; i++) {
      const char* label = mxGetFieldNameByNumber(ptr,i);
      mxArray*    p     = mxGetFieldByNumber(ptr,0,i);
      setOption(label,p);
    }
  }
}

bool IpoptOptions::useQuasiNewton() const {
  bool        b;  // The return value.
  std::string value;

  app.Options()->GetStringValue("hessian_approximation",value,"");
  b = !value.compare("limited-memory");
  return b;
}

bool IpoptOptions::useDerivChecker() const {
  bool        b;  // The return value.
  std::string value;

  app.Options()->GetStringValue("derivative_test",value,"");  
  b = value.compare("none");
  return b;
}

bool IpoptOptions::userScaling() const {
  bool        b;  // The return value.
  std::string value;

  app.Options()->GetStringValue("nlp_scaling_method",value,"");  
  b = !value.compare("user-scaling");
  return b;
}

int IpoptOptions::printLevel() const {
  int value;  // The return value.
  app.Options()->GetIntegerValue("print_level",value,"");
  return value;
}

double IpoptOptions::getPosInfty() const {
  double value;  // The return value.
  app.Options()->GetNumericValue("nlp_upper_bound_inf",value,"");
  return value;
}

double IpoptOptions::getNegInfty() const {
  double value;  // The return value.
  app.Options()->GetNumericValue("nlp_lower_bound_inf",value,"");
  return value;
}

void IpoptOptions::setOption (const char* label, const mxArray* ptr) {

  // Check to make sure we have a valid option.
  SmartPtr<const RegisteredOption> option = app.RegOptions()->GetOption(label);
  if (!IsValid(option)) {
    char buf[256];
    Snprintf(buf, 255, "You have specified a nonexistent IPOPT option (\"%s\")", label);
    throw MatlabException(buf);
  }

  Ipopt::RegisteredOptionType type = option->Type();
  if (type == Ipopt::OT_String)
    setStringOption(label,ptr);
  else if (type == Ipopt::OT_Integer)
    setIntegerOption(label,ptr);
  else
    setNumberOption(label,ptr);
}

void IpoptOptions::setStringOption (const char* label, const mxArray* ptr) {

  // Check whether the option value is a string.
  if (!mxIsChar(ptr)) {
    char buf[256];
    Snprintf(buf, 255, "IPOPT option value for option \"%s\" should be a string", label);
    throw MatlabException(buf);
  }

  // Get the option value.
  char* value = mxArrayToString(ptr);

  // Set the option.
  bool success = app.Options()->SetStringValue(label,value);
  if (!success) {
    char buf[256];
    Snprintf(buf, 255, "Invalid value for IPOPT option \"%s\"", label);
    throw MatlabException(buf);
  }

  // Free the dynamically allocated memory.
  mxFree(value);
}

void IpoptOptions::setIntegerOption (const char* label, const mxArray* ptr) {
  
  // Check whether the option value is a number.
  if (!mxIsDouble(ptr)) {
    char buf[256];
    Snprintf(buf, 255, "IPOPT option value for option \"%s\" should be an integer", label);
    throw MatlabException(buf);
  }
  
  // Set either the integer option.
  double value   = mxGetScalar(ptr);
  bool   success = app.Options()->SetIntegerValue(label,(int) value);
  if (!success) {
    char buf[256];
    Snprintf(buf, 255, "Invalid value for integer IPOPT option \"%s\"", label);
    throw MatlabException(buf);
  }
}

void IpoptOptions::setNumberOption (const char* label, const mxArray* ptr) {
  
  // Check whether the option value is a number.
  if (!mxIsDouble(ptr)) {
    char buf[256];
    Snprintf(buf, 255, "IPOPT option value for option \"%s\" should be a number", label);
    throw MatlabException(buf);
  }
  
  // Set either the numeric option.
  double value   = mxGetScalar(ptr);
  bool   success = app.Options()->SetNumericValue(label,value);
  if (!success) {
    char buf[256];
    Snprintf(buf, 255, "Invalid value for numeric IPOPT option \"%s\"", label);
    throw MatlabException(buf);
  }
}
