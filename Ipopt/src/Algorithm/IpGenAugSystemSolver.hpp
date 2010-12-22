// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter     IBM    2007-03-01

#ifndef __IP_STDAUGSYSTEMSOLVER_HPP__
#define __IP_STDAUGSYSTEMSOLVER_HPP__

#include "IpAugSystemSolver.hpp"
#include "IpGenKKTSolverInterface.hpp"

namespace Ipopt
{
  /** Solver for the augmented system using GenKKTSolverInterfaces.
   *
   *  This takes any Vector values out and provides Number*'s, but
   *  Matrices are provided as given from the NLP.
   */
  class GenAugSystemSolver : public AugSystemSolver
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor using only a linear solver object */
    GenAugSystemSolver(GenKKTSolverInterface& SolverInterface);

    /** Default destructor */
    virtual ~GenAugSystemSolver();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    /** Set up the augmented system and solve it for a set of given
     *  right hand side - implementation for GenTMatrices and
     *  SymTMatrices.
     */
    virtual ESymSolverStatus MultiSolve(
      const SymMatrix* W,
      double W_factor,
      const Vector* D_x,
      double delta_x,
      const Vector* D_s,
      double delta_s,
      const Matrix* J_c,
      const Vector* D_c,
      double delta_c,
      const Matrix* J_d,
      const Vector* D_d,
      double delta_d,
      std::vector<SmartPtr<const Vector> >& rhs_xV,
      std::vector<SmartPtr<const Vector> >& rhs_sV,
      std::vector<SmartPtr<const Vector> >& rhs_cV,
      std::vector<SmartPtr<const Vector> >& rhs_dV,
      std::vector<SmartPtr<Vector> >& sol_xV,
      std::vector<SmartPtr<Vector> >& sol_sV,
      std::vector<SmartPtr<Vector> >& sol_cV,
      std::vector<SmartPtr<Vector> >& sol_dV,
      bool check_NegEVals,
      Index numberOfNegEVals);

    /** Number of negative eigenvalues detected during last
     * solve.  Returns the number of negative eigenvalues of
     * the most recent factorized matrix.  This must not be called if
     * the linear solver does not compute this quantities (see
     * ProvidesInertia).
     */
    virtual Index NumberOfNegEVals() const;

    /** Query whether inertia is computed by linear solver.
     * Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const;

    /** Request to increase quality of solution for next solve.  Ask
     *  underlying linear solver to increase quality of solution for
     *  the next solve (e.g. increase pivot tolerance).  Returns
     *  false, if this is not possible (e.g. maximal pivot tolerance
     *  already used.)
     */
    virtual bool IncreaseQuality();

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default constructor. */
    GenAugSystemSolver();
    /** Copy Constructor */
    GenAugSystemSolver(const GenAugSystemSolver&);

    /** Overloaded Equals Operator */
    void operator=(const GenAugSystemSolver&);
    //@}

    /** Check the internal tags and decide if the passed variables are
     *  different from what is in the augmented_system_ */
    bool AugmentedSystemChanged(const SymMatrix* W,
                                double W_factor,
                                const Vector* D_x,
                                double delta_x,
                                const Vector* D_s,
                                double delta_s,
                                const Matrix& J_c,
                                const Vector* D_c,
                                double delta_c,
                                const Matrix& J_d,
                                const Vector* D_d,
                                double delta_d);

    void UpdateTags(const SymMatrix* W,
                    double W_factor,
                    const Vector* D_x,
                    double delta_x,
                    const Vector* D_s,
                    double delta_s,
                    const Matrix& J_c,
                    const Vector* D_c,
                    double delta_c,
                    const Matrix& J_d,
                    const Vector* D_d,
                    double delta_d);

    /** The linear solver object that is to be used to solve the
     *  linear systems.
     */
    SmartPtr<GenKKTSolverInterface> solver_interface_;

    /**@name Tags and values to track in order to decide whether the
       matrix has to be updated compared to the most recent call of
       the Set method.
     */
    //@{
    /** Tag for W matrix.  If W has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag w_tag_;
    /** Most recent value of W_factor */
    double w_factor_;
    /** Tag for D_x vector, representing the diagonal matrix D_x.  If
     *  D_x has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_x_tag_;
    /** Most recent value of delta_x from Set method */
    double delta_x_;
    /** Tag for D_s vector, representing the diagonal matrix D_s.  If
     *  D_s has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_s_tag_;
    /** Most recent value of delta_s from Set method */
    double delta_s_;
    /** Tag for J_c matrix.  If J_c has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag j_c_tag_;
    /** Tag for D_c vector, representing the diagonal matrix D_c.  If
     *  D_c has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_c_tag_;
    /** Most recent value of delta_c from Set method */
    double delta_c_;
    /** Tag for J_d matrix.  If J_d has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag j_d_tag_;
    /** Tag for D_d vector, representing the diagonal matrix D_d.  If
     *  D_d has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_d_tag_;
    /** Most recent value of delta_d from Set method */
    double delta_d_;
    //@}

    /** @name Space for storing the diagonal matrices.  If the matrix
     *  hasn't changed, we can use it from the last call. */
    //@{
    Number* dx_vals_copy_;
    Number* ds_vals_copy_;
    Number* dc_vals_copy_;
    Number* dd_vals_copy_;
    //@}

    /** @name Algorithmic parameters */
    //@{
    /** Flag indicating whether the TNLP with identical structure has
     *  already been solved before. */
    bool warm_start_same_structure_;
    //@}
  };

} // namespace Ipopt

#endif
