// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpNLPScaling.hpp,v 1.3 2004/11/19 00:55:18 andreasw Exp $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPGRADIENTCALING_HPP__
#define __IPGRADIENTSCALING_HPP__

#include "IpNLPScaling.hpp"
#include "IpNLP.hpp"
#include "IpIpoptType.hpp"

namespace Ipopt
{

  DeclareIpoptType(GradientScaling);

  /** This class does problem scaling by setting the 
   *  scaling parameters based on the maximum of the
   *  gradient at the user provided initial point.
   */
  class GradientScaling : public StandardScalingBase
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    GradientScaling(const SmartPtr<NLP>& nlp)
        :
      StandardScalingBase(),
      nlp_(nlp)
    {}

    /** Default destructor */
    virtual ~GradientScaling()
    {}
    //@}

    /** Initialize the object from the options */
    bool Initialize(const Journalist& jnlst,
		    const OptionsList& options,
		    const std::string& prefix);

    /** Methods for IpoptType */
    //@{
    /** Register the options for this class */
    static void RegisterOptions(SmartPtr<RegisteredOptions>& roptions);
    //@}

  protected:
    virtual void DetermineScalingParametersImpl(
	const SmartPtr<const VectorSpace> x_space,
	const SmartPtr<const VectorSpace> c_space,
	const SmartPtr<const VectorSpace> d_space,
	const SmartPtr<const MatrixSpace> jac_c_space,
	const SmartPtr<const MatrixSpace> jac_d_space,
	const SmartPtr<const SymMatrixSpace> h_space,
	Number& df, Vector& dx, 
	Vector& dc, Vector& dd);

  private:

    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{

    /** Copy Constructor */
    GradientScaling(const GradientScaling&);

    /** Overloaded Equals Operator */
    void operator=(const GradientScaling&);
    //@}

    /** pointer to the NLP to get scaling parameters */
    SmartPtr<NLP> nlp_;

    /** maximum allowed gradient before scaling is performed */
    Number scaling_max_gradient_;
  };
} // namespace Ipopt
#endif
