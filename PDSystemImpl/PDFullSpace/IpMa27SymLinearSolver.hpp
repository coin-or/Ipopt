// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPMA27SYMLINSOLVER_HPP__
#define __IPMA27SYMLINSOLVER_HPP__

#include "IpSymLinearSolver.hpp"
#include "IpSymMatrix.hpp"
#include <vector>

namespace Ipopt
{

  /** Implementation of a SymLinearSolver Class, using the Harwell
   * Routines MA27.
   */
  class Ma27SymLinearSolver: public SymLinearSolver
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor */
    Ma27SymLinearSolver();

    /** Destructor */
    virtual ~Ma27SymLinearSolver();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Solve operation for multiple right hand sides.  For details
     * see the description in the base class SymLinearSolver.
     */
    virtual ESolveStatus MultiSolve(const SymMatrix &A,
                                    std::vector<const Vector*>& rhsV,
                                    std::vector<Vector*>& solV,
                                    bool check_NegEVals,
                                    Index numberOfNegEVals);

    /** Number of negative eigenvalues detected during last
     * factorization.  Returns the number of negative eigenvalues of
     * the most recent factorized matrix.
     */
    virtual Index NumberOfNegEVals() const;
    //@}

    //* @name Options of Linear solver */
    //@{
    /** Request to increase quality of solution for next solve.
     * Ask linear solver to increase quality of solution for the next
     * solve (e.g. increase pivot tolerance).  Returns false, if this
     * is not possible (e.g. maximal pivot tolerance already used.)
     */
    virtual bool IncreaseQuality();

    /** Query whether inertia is computed by linear solver.
     * Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const
    {
      return true;
    }
    //@}

    /** @name Setting specific options */
    //@{
    /** Set the default pivot tolerance */
    void SetPivTol(Number PivTol)
    {
      pivtol_ = PivTol;
    }

    /** Set the maximal pivot tolerance */
    void SetPivTolMax(Number PivTolMax)
    {
      pivtolmax_ = PivTolMax;
    }

    /** Set factor for estimating initial value of liw */
    void SetLiwInitFactor(Number LiwInitFactor)
    {
      liw_init_factor_ = LiwInitFactor;
    }

    /** Set factor for estimating initial value of la */
    void SetLaInitFactor(Number LaInitFactor)
    {
      la_init_factor_ = LaInitFactor;
    }
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
    Ma27SymLinearSolver(const Ma27SymLinearSolver&);

    /** Overloaded Equals Operator */
    void operator=(const Ma27SymLinearSolver&);
    //@}

    /** @name Information about the matrix */
    //@{
    /** Tag for the incoming matrix */
    TaggedObject::Tag atag_;

    /** Number of rows and columns of the matrix */
    Index dim_;

    /** Number of nonzeros of the matrix */
    Index nonzeros_;
    //@}

    /** @name Information about most recent factorization/solve */
    //@{
    /** Number of negative eigenvalues */
    Index negevals_;
    //@}

    /** @name Initialization flags */
    //@{
    /** Flag indicating if internal data is initialized.
     *  For initialization, this object needs to have seen a matrix */
    bool initialized_;
    /** Flag indicating if matrix has been factorizated.
     *  Is set to true, if data from a factorization is available */
    bool factorized_;
    //@}

    /** @name Solver specific data/options */
    //@{
    /** Pivol tolerance */
    Number pivtol_;

    /** Maximal pivot tolerance */
    Number pivtolmax_;

    /** Factor for estimating initial value of liw */
    Number liw_init_factor_;
    /** Factor for estimating initial value of la */
    Number la_init_factor_;
    /** Factor for increaseing memory */
    Number meminc_factor_;
    //@}

    /** @name Data for the linear solver.
     * Storing factorization and other solver specific data structure.
     */
    //@{
    /** integer control values */
    ipfint icntl_[30];
    /** real control values */
    double cntl_[5];

    /** row indices of matrix.
     * We keep a copy of this in case type ipfint is not same as type Index
     */
    ipfint* airn_;
    /** column indices of matrix.
     * We keep a copy of this in case type ipfint is not same as type Index
     */
    ipfint* ajcn_;

    /** length of integer work space */
    ipfint liw_;
    /** integer work space */
    ipfint* iw_;

    /** MA28's IKEEP */
    ipfint* ikeep_;
    /** MA27's NSTEPS */
    ipfint nsteps_;
    /** MA27's MAXFRT */
    ipfint maxfrt_;

    /** length LA of A */
    ipfint la_;
    /** factor A of matrix */
    double* a_;

    /** flag indicating that la should be increased before next factorization
     */
    bool la_increase_;
    /** flag indicating that liw should be increased before next factorization
     */
    bool liw_increase_;
    //@}

    /** @name Internal functions */
    //@{
    /** Initialize nonzero structure.
     *  Set dim_ and nonzeros_, and copy the nonzero structure of symT_A
     *  into airn_ and ajcn_
     */
    void InitializeStructure(const SymMatrix& symT_A);

    /** Call MA27AD and reserve memory for MA27 data.
     *  Reserve memory for iw_ and ikeep_, call MA27AD to perform
     *  symbolic manipulations, and reserve all the remaining data memory
     */
    ESolveStatus SymbolicFactorization();

    /** Call MA27BD to factorize the Matrix.
     *  It is assumed that the first nonzeros_ element of a_ contain the values
     *  of the matrix to be factorized.
     */
    ESolveStatus Factorization(const SymMatrix& A,
                               bool check_NegEVals,
                               Index numberOfNegEVals);

    /** Call MA27CD to do the backsolve.
     */
    ESolveStatus Backsolve(std::vector<const Vector*>& rhsV,
                           std::vector<Vector*>& solV);
    //@}
  };

} // namespace Ipopt
#endif
