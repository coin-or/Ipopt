// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-16

#ifndef __AS_MEASUREMENT_HPP__
#define __AS_MEASUREMENT_HPP__

#include "IpReferenced.hpp"
#include "IpDenseVector.hpp"
#include "IpIteratesVector.hpp"

namespace Ipopt
{

class SIPOPTLIB_EXPORT Measurement: public ReferencedObject
{
   /** This class provides an abstraction for the measurements of the states coming in
    * and the solutions of the controller.
    *
    * It basically acts as the "plant" of the controller.
    */
public:

   Measurement()
   { }

   virtual ~Measurement()
   { }

   /** This function returns a std::vector holding the indices in IteratesVector of the
    *  equations that are to be "slacked" to free the initial values for sIPOPT.
    *
    *  This std::vector is used in the construction of the A-SchurData for the Schur Decomposition.
    */
   virtual std::vector<Index> GetInitialEqConstraints() = 0;

   /** This function returns delta_u.
    *
    * It should use the values of IpData().trial()->x()
    */
   virtual SmartPtr<DenseVector> GetMeasurement(
      Index measurement_number
   ) = 0;

   /** This function does whatever the measurement machine does with the solution of the SensAlgorithm */
   virtual void SetSolution(
      Index                    measurement_number,
      SmartPtr<IteratesVector> sol
   ) = 0;
};

}

#endif
