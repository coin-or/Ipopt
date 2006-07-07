// Copyright (C) 2006 Damien Hocking, KBC Advanced Technologies
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors: Damien Hocking                 KBC    2006-03-20
//        (included his original contribution into Ipopt package on 2006-03-25)
//          Andreas Waechter               IBM    2006-03-25
//           (minor changes and corrections)


#ifndef __IPMUMPSSOLVERINTERFACE_HPP__
#define __IPMUMPSSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

extern "C"
{
#include "dmumps_c.h"
}

namespace Ipopt
{

  /** Interface to the linear solver Mumps, derived from
   *  SparseSymLinearSolverInterface.  For details, see description of
   *  SparseSymLinearSolverInterface base class.
   */
  class MumpsSolverInterface: public SparseSymLinearSolverInterface
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor */
    MumpsSolverInterface();

    /** Destructor */
    virtual ~MumpsSolverInterface();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);


    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Method for initializing internal stuctures. */
    virtual ESymSolverStatus InitializeStructure(Index dim, Index nonzeros, const Index *ia, const Index *ja);

    /** Method returing an internal array into which the nonzero
     *  elements are to be stored. */
    virtual double* GetValuesArrayPtr();

    /** Solve operation for multiple right hand sides. */
    virtual ESymSolverStatus MultiSolve(bool new_matrix,
                                        const Index* ia,
                                        const Index* ja,
                                        Index nrhs,
                                        double* rhs_vals,
                                        bool check_NegEVals,
                                        Index numberOfNegEVals);

    /** Number of negative eigenvalues detected during last
     *  factorization.
     */
    virtual Index NumberOfNegEVals() const;
    //@}

    //* @name Options of Linear solver */
    //@{
    /** Request to increase quality of solution for next solve.
     */
    virtual bool IncreaseQuality();

    /** Query whether inertia is computed by linear solver.
     *  Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const
    {
      return true;
    }
    /** Query of requested matrix type that the linear solver
     *  understands.
     */
    EMatrixFormat MatrixFormat() const
    {
      return Triplet_Format;
    }
    //@}

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

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
    MumpsSolverInterface(const MumpsSolverInterface&);

    /** Overloaded Equals Operator */
    void operator=(const MumpsSolverInterface&);
    //@}

    /** @name Information about the matrix */
    //@{
    /** @name Information about the matrix */
    //@{
    /** Number of rows and columns of the matrix */
    Index n;

    /** Number of nonzeros of the matrix in triplet representation. */
    Index nz;

    /** Array for storing the values of the matrix. */
    double* a;

    /** Array for storing the row indices of the matrix */
    int* irn_;
    /** Array for storing the column indices of the matrix */
    int* jcn_;
    //@}
    /** @name Information about most recent factorization/solve */
    //@{
    /** Number of negative eigenvalues */
    Index negevals;
    //@}

    /** @name Solver specific options */
    //@{
    //@}

    /** @name Initialization flags */
    //@{
    //@}

    /** @name Solver specific information */
    //@{
    /**@name Some counters for debugging */
    //@{
    //@}

    /** @name Internal functions */
    //@{
    /** Call Mumps to do the analysis phase.
     */
    /** Call Mumps to factorize the Matrix.
     */
    ESymSolverStatus Factorization(const Index* ia,
                                   const Index* ja,
                                   bool check_NegEVals,
                                   Index numberOfNegEVals);

    /** Call Mumps to do the Solve.
     */
    ESymSolverStatus Solve(const Index* ia,
                           const Index* ja,
                           Index nrhs,
                           double *rhs_vals);
    //@}
    //MUMPS data structure
    DMUMPS_STRUC_C mumps_data;

  };

} // namespace Ipopt
#endif
