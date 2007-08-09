// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_MATLABEXCEPTION
#define INCLUDE_MATLABEXCEPTION

#include <exception>

// Class MatlabException
// -----------------------------------------------------------------
// It is assumed that the argument passed to the constructor persists
// as long as the MatlabException object is in scope. Usually, this
// means that it should persist for the duration of the entire
// program. This is always the case if the input "message" is a literal.
class MatlabException : public std::exception {
public:
  MatlabException (const char* message) throw();
  ~MatlabException()                    throw() { };

  // The copy constructor makes a shallow copy.
  MatlabException (const MatlabException& source) throw();

  // The copy assignment operator makes a shallow copy as well.
  MatlabException& operator= (const MatlabException& source);
    
  // Return the message string.
  virtual const char* what () const throw() { return message; };
    
  private:
    const char* message;  // The error message.
};

#endif
