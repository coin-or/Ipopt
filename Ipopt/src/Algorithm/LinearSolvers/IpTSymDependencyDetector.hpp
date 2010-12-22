// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2007-04-18

#ifndef __IPTSYMDEPENDENCYDETECTOR_HPP__
#define __IPTSYMDEPENDENCYDETECTOR_HPP__

#include "IpTDependencyDetector.hpp"
#include "IpTSymLinearSolver.hpp"

namespace Ipopt
{

  /** Base class for all derived algorithms for detecting linearly
   *  dependent rows in the constraint Jacobian. */
  class TSymDependencyDetector: public TDependencyDetector
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    TSymDependencyDetector(TSymLinearSolver& tsym_linear_solver);

    virtual ~TSymDependencyDetector()
    {}
    //@}

    /** Has to be called to initialize and reset these objects. */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Method determining the number of linearly dependent rows in
     *  the matrix and the indices of those rows.  We assume that the
     *  matrix is available in "Triplet" format (MA28 format), and
     *  that the arrays given to this method can be modified
     *  internally, i.e., they are not used by the calling program
     *  anymore after this call.  This method returns false if there
     *  was a problem with the underlying linear solver.
     */
    virtual bool DetermineDependentRows(Index n_rows, Index n_cols,
                                        Index n_jac_nz,
                                        Number* jac_c_vals,
                                        Index* jac_c_iRow,
                                        Index* jac_c_jCol,
                                        std::list<Index>& c_deps);

    /** This must be called to make the options for this class
     *  known */
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);

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
    TSymDependencyDetector();

    /** Copy Constructor */
    TSymDependencyDetector(const TSymDependencyDetector&);

    /** Overloaded Equals Operator */
    void operator=(const TSymDependencyDetector&);
    //@}

    SmartPtr<const Journalist> jnlst_;

    SmartPtr<TSymLinearSolver> tsym_linear_solver_;
  };

} // namespace Ipopt

#endif
