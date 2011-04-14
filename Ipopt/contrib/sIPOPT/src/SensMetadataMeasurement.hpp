// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-07-10

#ifndef __AS_METADATAMEASUREMENT_HPP__
#define __AS_METADATAMEASUREMENT_HPP__

#include "SensMeasurement.hpp"
#include "SensSuffixHandler.hpp"
#include "IpAlgStrategy.hpp"


namespace Ipopt
{

  class MetadataMeasurement : public Measurement, public SuffixHandler, public AlgorithmStrategyObject
  {
  public:
    MetadataMeasurement();

    virtual ~MetadataMeasurement();

    /* AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /* measurement methods */
    virtual std::vector<Index> GetInitialEqConstraints();

    virtual SmartPtr<DenseVector> GetMeasurement(Index measurement_number);

    virtual void SetSolution(Index measurement_number, SmartPtr<IteratesVector> sol);

    /** suffix handler methods */

    virtual std::vector<Index> GetIntegerSuffix(std::string suffix_string);

  private:

    /** Number of sens_indices */
    Index n_idx_;

    /** owner space of x */
    SmartPtr<const DenseVectorSpace> x_owner_space_;
    /** owner space of s */
    SmartPtr<const DenseVectorSpace> s_owner_space_;
    /** owner space of y_c */
    SmartPtr<const DenseVectorSpace> y_c_owner_space_;
    /** owner space of y_d */
    SmartPtr<const DenseVectorSpace> y_d_owner_space_;
    /** owner space of z_L */
    SmartPtr<const DenseVectorSpace> z_L_owner_space_;
    /** owner space of z_U */
    SmartPtr<const DenseVectorSpace> z_U_owner_space_;

  };

}

#endif
