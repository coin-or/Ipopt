// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef INCLUDE_ARRAY
#define INCLUDE_ARRAY

#include "matlabexception.h"
#include <string.h>

// Function definitions.
// ---------------------------------------------------------------
template <class Type> void copymemory (const Type* source, Type* dest, 
				       int length) {
  memcpy(dest,source,sizeof(Type)*length);
}

// Class Array
// ---------------------------------------------------------------
// An Array object stores a collection of ordered elements. The key
// aspect of Array objects is that they do not necessarily achieve
// encapsulation; they do not necessarily retain independent
// ownership of their data. That is, the data could be modified
// externally by another agent. This behaviour is determined by the
// data member "owner". If "owner" is false, the array does not take
// care of allocation and deallocation of the data in memory.
//
// Copy assignment behaves different than usual---it copies the
// data, but requires that the destination already have the proper
// resources allocated in memory. 
//
// The copy constructor performs a shallow copy.
template <class Type> class Array { 
  public:
  
  // This constructor allocates memory for an array of elements of
  // the specificed length. 
  explicit Array (int length);

  // This constructor is set to point to an already existing
  // array. As such, it does not gain ownership of the elements.
  Array (Type* data, int length);

  // The copy constructor performs a shallow copy, so that both
  // arrays share the same copy of the data.
  Array (const Array<Type>& source);

  // The destructor.
  ~Array();

  // Set all the elements of the array to the specified value.
  void setvalue (const Type& value);

  // Copy the elements to from the source array. Both the source and
  // destination array must have the same length.
  void         inject    (const Type* source);
  void         inject    (const Array<Type>& source);
  Array<Type>& operator= (const Type* source);
  Array<Type>& operator= (const Array<Type>& source);

  // Copy the elements to the location in memory pointed to by
  // "dest". It is assumed that sufficient memory is allocated for
  // the destination.
  void copyto (Type* dest) const;

  // Get the number of elements in the array.
  int length() const { return n; };

  // Elements of an array A can be accessed the subscript operator;
  // e.g. A[i]. For the sake of efficiency, no checks are made to
  // ensure that the subscripting is legal. It is up to the user to
  // ensure that i is non-negative and less than A.length(). Note
  // that subscripts, as per C++ convention, start at 0.
  Type&       operator[] (int i);
  const Type& operator[] (int i) const;

  // Returns true if the two arrays have the same dimensions
  // (i.e. the same lengths).
  bool operator== (const Array<Type>& a) const { return n == a.n; };
  bool operator!= (const Array<Type>& a) const { return !(*this == a); };

protected:
  Type* elems;
  int   n;      // The number of elements in the array.
  bool  owner;  // Whether the object has ownership of the data.
};

// Function definitions for class Array.
// -----------------------------------------------------------------
template <class Type> Array<Type>::Array (int length) {
  n     = length;
  owner = true;
  elems = 0;
  elems = new Type[length];
}  

template <class Type> Array<Type>::Array (Type* data, int length) {
  n     = length;
  owner = false;
  elems = data;
}

template <class Type> Array<Type>::Array (const Array<Type>& source) {
  n     = source.n;
  owner = false;
  elems = source.elems;
}

template <class Type> Array<Type>::~Array () { 
  if (owner && elems) 
    delete[] elems;
}

template <class Type> void Array<Type>::setvalue (const Type& value) { 
  for (int i = 0; i < n; i++)
    elems[i] = value;
}

template <class Type> void Array<Type>::inject (const Type* source) {

  // Check for self-assignment. If there is self-assignment, there
  // is no need to copy because the objects share the same data!
  if (elems != source)
    copymemory<Type>(source,elems,n);
}

template <class Type> void Array<Type>::inject (const Array<Type>& source) {
  if (n != source.n)
    throw MatlabException("Unable to perform copy; arrays have \
different lengths");    
  inject(source.elems);
}

template <class Type> Array<Type>& Array<Type>::operator= 
(const Type* source) {
  inject(source);
  return *this;
}

template <class Type> Array<Type>& Array<Type>::operator= 
(const Array<Type>& source) {
  inject(source);
  return *this;
}

template <class Type> void Array<Type>::copyto (Type* dest) const {
  Array destarray(dest,n);
  destarray.inject(*this);
}

template <class Type> Type& Array<Type>::operator[] (int i) {
  return elems[i]; 
}

template <class Type> const Type& Array<Type>::operator[] (int i) const {
  return elems[i]; 
}

#endif
