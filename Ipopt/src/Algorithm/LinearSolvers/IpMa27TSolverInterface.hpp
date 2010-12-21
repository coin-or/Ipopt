// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17

#ifndef __IPMA27TSOLVERINTERFACE_HPP__
#define __IPMA27TSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

namespace Ipopt
{
  /** Interface to the symmetric linear solver MA27, derived from
   *  SparseSymLinearSolverInterface.
   */
  class Ma27TSolverInterface: public SparseSymLinearSolverInterface
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor */
    Ma27TSolverInterface();

    /** Destructor */
    virtual ~Ma27TSolverInterface();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);


    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Method for initializing internal stuctures.  Here, ndim gives
     *  the number of rows and columns of the matrix, nonzeros give
     *  the number of nonzero elements, and airn and acjn give the
     *  positions of the nonzero elements.
     */
    virtual ESymSolverStatus InitializeStructure(Index dim, Index nonzeros,
        const Index *airn,
        const Index *ajcn);

    /** Method returing an internal array into which the nonzero
     *  elements (in the same order as airn and ajcn) are to be stored
     *  by the calling routine before a call to MultiSolve with a
     *  new_matrix=true.  The returned array must have space for at least
     *  nonzero elements. */
    virtual double* GetValuesArrayPtr();

    /** Solve operation for multiple right hand sides.  Overloaded
     *  from SparseSymLinearSolverInterface.
     */
    virtual ESymSolverStatus MultiSolve(bool new_matrix,
                                        const Index* airn,
                                        const Index* ajcn,
                                        Index nrhs,
                                        double* rhs_vals,
                                        bool check_NegEVals,
                                        Index numberOfNegEVals);

    /** Number of negative eigenvalues detected during last
     *  factorization.  Returns the number of negative eigenvalues of
     *  the most recent factorized matrix.  This must not be called if
     *  the linear solver does not compute this quantities (see
     *  ProvidesInertia).
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
    Ma27TSolverInterface(const Ma27TSolverInterface&);

    /** Overloaded Equals Operator */
    void operator=(const Ma27TSolverInterface&);
    //@}

    /** @name Information about the matrix */
    //@{
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
    /** Flag indicating if the matrix has to be refactorized because
     *  the pivot tolerance has been changed. */
    bool pivtol_changed_;
    /** Flag that is true if we just requested the values of the
     *  matrix again (SYMSOLVER_CALL_AGAIN) and have to factorize
     *  again. */
    bool refactorize_;
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
    /** Flag indicating whether the TNLP with identical structure has
     *  already been solved before. */
    bool warm_start_same_structure_;
    /** Flag indicating if the interia is always assumed to be
     *  correct. */
    bool skip_inertia_check_;
    /** Flag indicating if MA27 should continue if a singular matrix
    is detected, but right hands sides are still accepted. */
    bool ignore_singularity_;
    //@}

    /** @name Data for the linear solver.
     * Storing factorization and other solver specific data structure.
     */
    //@{
    /** integer control values */
    ipfint icntl_[30];
    /** real control values */
    double cntl_[5];

    /** length of integer work space */
    ipfint liw_;
    /** integer work space */
    ipfint* iw_;

    /** MA27's IKEEP */
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
    /** Call MA27AD and reserve memory for MA27 data.
     *  Reserve memory for iw_ and ikeep_, call MA27AD to perform
     *  symbolic manipulations, and reserve all the remaining data memory
     */
    ESymSolverStatus SymbolicFactorization(const Index* airn,
                                           const Index* ajcn);

    /** Call MA27BD to factorize the Matrix.
     *  It is assumed that the first nonzeros_ element of a_ contain the values
     *  of the matrix to be factorized.
     */
    ESymSolverStatus Factorization(const Index* airn,
                                   const Index* ajcn,
                                   bool check_NegEVals,
                                   Index numberOfNegEVals);

    /** Call MA27CD to do the backsolve.
     */
    ESymSolverStatus Backsolve(Index nrhs,
                               double *rhs_vals);
    //@}
  };

} // namespace Ipopt
#endif
