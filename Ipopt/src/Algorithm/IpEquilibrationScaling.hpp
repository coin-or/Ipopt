// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2007-05-21

#ifndef __IPEQUILIBRATIONSCALING_HPP__
#define __IPEQUILIBRATIONSCALING_HPP__

#include "IpNLPScaling.hpp"
#include "IpNLP.hpp"

namespace Ipopt
{
  /** This class does problem scaling by setting the
   *  scaling parameters based on the maximum of the
   *  gradient at the user provided initial point.
   */
  class EquilibrationScaling : public StandardScalingBase
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    EquilibrationScaling(const SmartPtr<NLP>& nlp)
        :
        StandardScalingBase(),
        nlp_(nlp)
    {}

    /** Default destructor */
    virtual ~EquilibrationScaling()
    {}
    //@}

    /** Methods for IpoptType */
    //@{
    /** Register the options for this class */
    static void RegisterOptions(const SmartPtr<RegisteredOptions>& roptions);
    //@}

  protected:
    /** Initialize the object from the options */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    virtual void DetermineScalingParametersImpl(
      const SmartPtr<const VectorSpace> x_space,
      const SmartPtr<const VectorSpace> c_space,
      const SmartPtr<const VectorSpace> d_space,
      const SmartPtr<const MatrixSpace> jac_c_space,
      const SmartPtr<const MatrixSpace> jac_d_space,
      const SmartPtr<const SymMatrixSpace> h_space,
      const Matrix& Px_L, const Vector& x_L,
      const Matrix& Px_U, const Vector& x_U,
      Number& df,
      SmartPtr<Vector>& dx,
      SmartPtr<Vector>& dc,
      SmartPtr<Vector>& dd);

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
    EquilibrationScaling(const EquilibrationScaling&);

    /** Overloaded Equals Operator */
    void operator=(const EquilibrationScaling&);
    //@}

    /** pointer to the NLP to get scaling parameters */
    SmartPtr<NLP> nlp_;

    /** maximal radius for the random perturbation of the initial
     *  point. */
    Number point_perturbation_radius_;
  };

  /** This class is a simple object for generating randomly perturbed
   *  points that are withing the NLP bounds.  The
   *  ramdon_perturb_radius gives the upper bound of the
   *  perturbation. */
  class PointPerturber : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    PointPerturber(const Vector& reference_point,
                   Number random_pert_radius,
                   const Matrix& Px_L, const Vector& x_L,
                   const Matrix& Px_U, const Vector& x_U);

    /** Default destructor */
    virtual ~PointPerturber()
    {}
    //@}

    /** Return a new perturbed point */
    SmartPtr<Vector> MakeNewPerturbedPoint() const;

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
    PointPerturber(const PointPerturber&);

    /** Overloaded Equals Operator */
    void operator=(const PointPerturber&);
    //@}

    /** pointer to the midpoint of the perturbation */
    SmartPtr<Vector> ref_point_;

    /** pointer to the perturbation vector */
    SmartPtr<Vector> pert_dir_;
  };

} // namespace Ipopt
#endif
