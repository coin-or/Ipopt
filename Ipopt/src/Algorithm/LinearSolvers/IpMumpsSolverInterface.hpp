// Copyright (C) 2006, 2007 Damien Hocking, KBC Advanced Technologies
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Damien Hocking                 KBC    2006-03-20
//        (included his original contribution into Ipopt package on 2006-03-25)
//          Andreas Waechter               IBM    2006-03-25
//           (minor changes and corrections)
//          Scott Turnberg                 CMU    2006-05-12
//           (major revision)
//           (incorporated by AW on 2006-11-11 into Ipopt package)


#ifndef __IPMUMPSSOLVERINTERFACE_HPP__
#define __IPMUMPSSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

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

    /** Query whether the indices of linearly dependent rows/columns
     *  can be determined by this linear solver. */
    virtual bool ProvidesDegeneracyDetection() const;

    /** This method determines the list of row indices of the linearly
     *  dependent rows. */
    virtual ESymSolverStatus DetermineDependentRows(const Index* ia,
        const Index* ja,
        std::list<Index>& c_deps);

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
    /** Primary MUMP data structure */
    void* mumps_ptr_;
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
    /** Counter on number of alive Mumps interface objects, if we have called MPI_Initialize.
     *
     * When the last object is destroyed, we will call MPI_Finalize.
     */
    static int instancecount_mpi;
    //@}

    /** @name Solver specific data/options */
    //@{
    /** Pivol tolerance */
    Number pivtol_;

    /** Maximal pivot tolerance */
    Number pivtolmax_;

    /** Percent increase in memory */
    Index mem_percent_;

    /** Permution and scaling method in MUMPS */
    Index mumps_permuting_scaling_;

    /** Pivot order in MUMPS. */
    Index mumps_pivot_order_;

    /** Scaling in MUMPS */
    Index mumps_scaling_;

    /** Threshold in MUMPS to stay that a constraint is linearly
     *  dependent */
    Number mumps_dep_tol_;

    /** Flag indicating whether the TNLP with identical structure has
     *  already been solved before. */
    bool warm_start_same_structure_;
    //@}

    /** Flag indicating if symbolic factorization has already been
     *  called */
    bool have_symbolic_factorization_;

    /** @name Internal functions */
    //@{
    /** Call MUMPS (job=1) to perform symbolic manipulations, and reserve
     *  memory.
     */
    ESymSolverStatus SymbolicFactorization();

    /** Call MUMPS (job=2) to factorize the Matrix.
     *  It is assumed that the first nonzeros_ element of a_ contain the values
     *  of the matrix to be factorized.
     */
    ESymSolverStatus Factorization(bool check_NegEVals,
                                   Index numberOfNegEVals);

    /** Call MUMPS (job=3) to do the solve.
     */
    ESymSolverStatus Solve(Index nrhs, double *rhs_vals);
    //@}
  };

} // namespace Ipopt
#endif
