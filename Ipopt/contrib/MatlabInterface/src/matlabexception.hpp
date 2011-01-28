// Copyright (C) 2007, 2009 Peter Carbonetto. All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_MATLABEXCEPTION
#define INCLUDE_MATLABEXCEPTION

#include <exception>

#define ME_BUFLEN 512

// Class MatlabException
// -----------------------------------------------------------------
// It is assumed that the argument passed to the constructor persists
// as long as the MatlabException object is in scope. Usually, this
// means that it should persist for the duration of the entire
// program. This is always the case if the input "message" is a literal.
//
// AW: Since I would like to include more detailed information (e.g.,
// which option is the one that is unknown), I changed this so that
// this object keeps a copy of the exception
class MatlabException : public std::exception {
public:
  MatlabException (const char* message) throw();
  ~MatlabException()                    throw() { };

  // The copy constructor makes a copy.
  MatlabException (const MatlabException& source) throw();

  // The copy assignment operator makes a copy as well.
  MatlabException& operator= (const MatlabException& source);
    
  // Return the message string.
  virtual const char* what () const throw() { return message; };
    
private:
  char message[ME_BUFLEN];  // The error message.
};

#endif
