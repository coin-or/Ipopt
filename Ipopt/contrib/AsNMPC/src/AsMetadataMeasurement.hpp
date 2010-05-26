// Copyright 2009 Hans Pirnay
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Date   : 2009-07-10

#ifndef __AS_METADATAMEASUREMENT_HPP__
#define __AS_METADATAMEASUREMENT_HPP__

#include "AsMeasurement.hpp"
#include "AsSuffixHandler.hpp"
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

    virtual std::vector<Index> GetNmpcState(Index i);

    virtual SmartPtr<DenseVector> GetMeasurement(Index measurement_number);

    virtual void SetSolution(Index measurement_number, SmartPtr<IteratesVector> sol);    

    /* suffix handler methods */
    
    virtual std::vector<Index> GetIntegerSuffix(std::string suffix_string);

  private:

    /* Number of nmpc_indices */
    Index n_idx_;
    
    std::string select_step_;
    /* owner space of x */
    SmartPtr<const DenseVectorSpace> x_owner_space_;
    /* owner space of y_c */
    SmartPtr<const DenseVectorSpace> c_owner_space_; 
  };

}

#endif
