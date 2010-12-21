// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-18
//            based on IpPardisoSolverInterface.hpp rev 1119


#ifndef __IPITERATIVEPARDISOSOLVERINTERFACE_HPP__
#define __IPITERATIVEPARDISOSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpInexactCq.hpp"
#include "IpIterativeSolverTerminationTester.hpp"

namespace Ipopt
{

  /** Interface to the linear solver Pardiso, derived from
   *  SparseSymLinearSolverInterface.  For details, see description of
   *  SparseSymLinearSolverInterface base class.
   */
  class IterativePardisoSolverInterface: public SparseSymLinearSolverInterface
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor */
    IterativePardisoSolverInterface(IterativeSolverTerminationTester& normal_tester,
                                    IterativeSolverTerminationTester& pd_tester);

    /** Destructor */
    virtual ~IterativePardisoSolverInterface();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);


    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Method for initializing internal stuctures. */
    virtual ESymSolverStatus InitializeStructure(Index dim, Index nonzeros,
        const Index *ia,
        const Index *ja);

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
      return CSR_Format_1_Offset;
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
    /** Default Constructor */
    IterativePardisoSolverInterface();

    /** Copy Constructor */
    IterativePardisoSolverInterface(const IterativePardisoSolverInterface&);

    /** Overloaded Equals Operator */
    void operator=(const IterativePardisoSolverInterface&);
    //@}

    /** @name Information about the matrix */
    //@{
    /** Number of rows and columns of the matrix */
    Index dim_;

    /** Number of nonzeros of the matrix in triplet representation. */
    Index nonzeros_;

    /** Array for storing the values of the matrix. */
    double* a_;
    //@}

    /** @name Information about most recent factorization/solve */
    //@{
    /** Number of negative eigenvalues */
    Index negevals_;
    //@}

    /** @name Solver specific options */
    //@{
    /** Type for mathcing strategies */
    enum PardisoMatchingStrategy
    {
      COMPLETE,
      COMPLETE2x2,
      CONSTRAINT
    };
    /** Option that controls the matching strategy. */
    PardisoMatchingStrategy match_strat_;
    /** Flag indicating if symbolic factorization has already been
     *  performed. */
    bool have_symbolic_factorization_;
    /** Flag indicating whether the symbolic factorization should only
     *  be done after perturbed elements, if the inertia was wrong */
    bool pardiso_redo_symbolic_fact_only_if_inertia_wrong_;
    /** Flag indicating whether repeated perturbed elements even after
     *  a new symbolic factorization should be interpreted as a
     *  singular matrix */
    bool pardiso_repeated_perturbation_means_singular_;
    /** Flag indicating if the interia is always assumed to be
      *  correct. */
    bool skip_inertia_check_;
    /** Maximal number of decreases of drop tolerance during one solve. */
    Index pardiso_max_droptol_corrections_;
    //@}

    /** Options for the preconditioner */
    //@{
    Index pardiso_max_iter_;
    Number pardiso_iter_relative_tol_;
    Index pardiso_iter_coarse_size_;
    Index pardiso_iter_max_levels_;
    Number pardiso_iter_dropping_factor_;
    Number pardiso_iter_dropping_schur_;
    Index pardiso_iter_max_row_fill_;
    Number pardiso_iter_inverse_norm_factor_;

    Index normal_pardiso_max_iter_;
    Number normal_pardiso_iter_relative_tol_;
    Index normal_pardiso_iter_coarse_size_;
    Index normal_pardiso_iter_max_levels_;
    Number normal_pardiso_iter_dropping_factor_;
    Number normal_pardiso_iter_dropping_schur_;
    Index normal_pardiso_iter_max_row_fill_;
    Number normal_pardiso_iter_inverse_norm_factor_;
    //@}

    /** Decrease factor for dropping tolerances */
    Number decr_factor_;

    /** Actualy used dropping tolerances */
    //@{
    Number pardiso_iter_dropping_factor_used_;
    Number pardiso_iter_dropping_schur_used_;
    Number normal_pardiso_iter_dropping_factor_used_;
    Number normal_pardiso_iter_dropping_schur_used_;
    //@}

    /** @name Initialization flags */
    //@{
    /** Flag indicating if internal data is initialized.
     *  For initialization, this object needs to have seen a matrix */
    bool initialized_;
    //@}

    /** @name Solver specific information */
    //@{
    /** Internal data address pointers. */
    void** PT_;
    /** Maximal number of factors with identical nonzero
     *  structure. Here, we only store one factorization. Is always 1.*/
    ipfint MAXFCT_;
    /** Actual matrix for the solution phase. Is always 1.*/
    ipfint MNUM_;
    /** Matrix type; real and symmetric indefinite.  Is always -2.*/
    ipfint MTYPE_;
    /** Parameter and info array for Pardiso. */
    ipfint* IPARM_;
    /** Parameter and info array for Pardiso. */
    double* DPARM_;
    /** Message level. */
    ipfint MSGLVL_;
    //@}

    /**@name Some counters for debugging */
    //@{
    Index debug_last_iter_;
    Index debug_cnt_;
    //@}

    /** @name Internal functions */
    //@{
    /** Call Pardiso to do the analysis phase.
     */
    ESymSolverStatus SymbolicFactorization(const Index* ia,
                                           const Index* ja);

    /** Call Pardiso to factorize the Matrix.
     */
    ESymSolverStatus Factorization(const Index* ia,
                                   const Index* ja,
                                   bool check_NegEVals,
                                   Index numberOfNegEVals);

    /** Call Pardiso to do the Solve.
     */
    ESymSolverStatus Solve(const Index* ia,
                           const Index* ja,
                           Index nrhs,
                           double *rhs_vals);
    //@}

    /** Method to easily access Inexact data */
    InexactData& InexData()
    {
      InexactData& inexact_data =
        static_cast<InexactData&>(IpData().AdditionalData());
      DBG_ASSERT(dynamic_cast<InexactData*>(&IpData().AdditionalData()));
      return inexact_data;
    }

    /** Method to easily access Inexact calculated quantities */
    InexactCq& InexCq()
    {
      InexactCq& inexact_cq =
        static_cast<InexactCq&>(IpCq().AdditionalCq());
      DBG_ASSERT(dynamic_cast<InexactCq*>(&IpCq().AdditionalCq()));
      return inexact_cq;
    }

    /** Termination tester for normal step computation */
    SmartPtr<IterativeSolverTerminationTester> normal_tester_;

    /** Termination tester for primal-dual step computation */
    SmartPtr<IterativeSolverTerminationTester> pd_tester_;

  };

} // namespace Ipopt
#endif
