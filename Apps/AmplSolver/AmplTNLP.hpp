// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPAMPLTNLP_HPP__
#define __IPAMPLTNLP_HPP__

#include "IpUtils.hpp"
#include "IpTNLP.hpp"
#include "IpJournalist.hpp"
#include "IpException.hpp"
#include "IpSmartPtr.hpp"

/* non Ipopt forward declaration */
struct ASL_pfgh;
struct SufDecl;
struct SufDesc;

namespace Ipopt
{
  /* forward declarations */
  class AmplSuffixHandler;

  /** Ampl Interface.
   *  Ampl Interface, implemented as a TNLP.
   */
  class AmplTNLP : public TNLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    AmplTNLP(const SmartPtr<Journalist>& jnlst, char**& argv, SmartPtr<AmplSuffixHandler> suffix_handler = NULL, bool allow_discrete = false);

    /** Default destructor */
    virtual ~AmplTNLP();
    //@}

    /** Exceptions */
    DECLARE_STD_EXCEPTION(NONPOSITIVE_SCALING_FACTOR);

    /**@name methods to gather information about the NLP. These
    * methods are overloaded from TNLP. See TNLP for their more
    * detailed documentation. */
    //@{
    /** returns dimensions of the nlp. Overloaded from TNLP */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag);

    /** returns bounds of the nlp. Overloaded from TNLP */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);

    /** provides a starting point for the nlp variables. Overloaded
    from TNLP */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda, Number* lambda);

    /** evaluates the objective value for the nlp. Overloaded from TNLP */
    virtual bool eval_f(Index n, const Number* x, bool new_x,
                        Number& obj_value);

    /** evaluates the gradient of the objective for the
    nlp. Overloaded from TNLP */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
                             Number* grad_f);

    /** evaluates the constraint residuals for the nlp. Overloaded from TNLP */
    virtual bool eval_g(Index n, const Number* x, bool new_x,
                        Index m, Number* g);

    /** specifies the jacobian structure (if values is NULL) and
     *  evaluates the jacobian values (if values is not NULL) for the
     *  nlp. Overloaded from TNLP */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow,
                            Index *jCol, Number* values);

    /** specifies the structure of the hessian of the lagrangian (if
     *  values is NULL) and evaluates the values (if values is not
     *  NULL). Overloaded from TNLP */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values);

    virtual void get_scaling_parameters(Number& obj_scaling,
                                        Index n, Number* x_scaling,
                                        Index m, Number* g_scaling);
    //@}

    /** @name Solution Methods */
    //@{
    virtual void finalize_solution(ApplicationReturnStatus status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value);
    //@}

    /**@name Ampl specific methods */
    //@{
    /** Return the ampl solver object (ASL*) */
    ASL_pfgh* AmplSolverObject()
    {
      return asl_;
    }

    /** Write the solution file.  This is a wrapper for AMPL's
     *  write_sol.  TODO Maybe this should be at a different place, or
     *  collect the numbers itself? */
    void write_solution_file(const std::string& message) const;

    /** ampl orders the variables like (continuous, binary, integer).
     *  This method gives the number of binary and integer variables.
     *  For details, see Tables 3 and 4 in "Hooking Your Solver to
     *  AMPL"
     */
    void get_discrete_info(Index& nlvb_,
                           Index& nlvbi_,
                           Index& nlvc_,
                           Index& nlvci_,
                           Index& nlvo_,
                           Index& nlvoi_,
                           Index& nbv_,
                           Index& niv_) const;
    //@}

  private:
    /** Journlist */
    SmartPtr<const Journalist> jnlst_;

    /** pointer to the main ASL structure */
    ASL_pfgh* asl_;

    /** Sign of the objective fn (1 for min, -1 for max) */
    double obj_sign_;

    /**@name Problem Size Data*/
    //@{
    Index nz_h_full_; // number of nonzeros in the full_x hessian
    /* the rest of the problem size data is available easily through the ampl variables */
    //@}

    /**@name Internal copies of data */
    //@{
    /** A non-const copy of x - this is kept up-to-date in apply_new_x */
    Number* non_const_x_;

    /** Solution Vectors */
    Number* x_sol_;
    Number* z_L_sol_;
    Number* z_U_sol_;
    Number* g_sol_;
    Number* lambda_sol_;
    Number obj_sol_;
    //@}

    /**@name Flags to track internal state */
    //@{
    /** true when the objective value has been calculated with the
     *  current x, set to false in apply_new_x, and set to true in
     *  internal_objval */
    bool objval_called_with_current_x_;
    /** true when the constraint values have been calculated with the
     *  current x, set to false in apply_new_x, and set to true in
     *  internal_conval */
    bool conval_called_with_current_x_;
    //@}


    /** Suffix Handler */
    SmartPtr<AmplSuffixHandler> suffix_handler_;

    /** Make the objective call to ampl */
    bool internal_objval(Number& obj_val);

    /** Make the constraint call to ampl*/
    bool internal_conval(Index m, Number* g=NULL);

    /** Internal function to update the internal and ampl state if the
    x value changes */
    void apply_new_x(bool new_x, Index n, const Number* x);


    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    AmplTNLP();

    /** Copy Constructor */
    AmplTNLP(const AmplTNLP&);

    /** Overloaded Equals Operator */
    void operator=(const AmplTNLP&);
    //@}

  };


  class AmplSuffixHandler : public ReferencedObject
  {
  public:
    AmplSuffixHandler();

    ~AmplSuffixHandler();

    enum Suffix_Type
    {
      Index_Type,
      Number_Type
    };

    enum Suffix_Source
    {
      Variable_Source,
      Constraint_Source,
      Objective_Source,
      Problem_Source
    };

    void AddAvailableSuffix(std::string suffix_string, Suffix_Source source, Suffix_Type type)
    {
      suffix_ids_.push_back(suffix_string);
      suffix_types_.push_back(type);
      suffix_sources_.push_back(source);
      //      suffix_values_.push_back();
    }

    const Index* GetIntegerSuffixValues(std::string suffix_string, Suffix_Source source) const;

    const Number* GetNumberSuffixValues(std::string suffix_string, Suffix_Source source) const;

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
    //AmplSuffixHandler();

    /** Copy Constructor */
    AmplSuffixHandler(const AmplSuffixHandler&);

    /** Overloaded Equals Operator */
    void operator=(const AmplSuffixHandler&);
    //@}

    mutable ASL_pfgh* asl_;

    SufDecl* suftab_;

    std::vector<std::string> suffix_ids_;
    std::vector<Suffix_Type> suffix_types_;
    std::vector<Suffix_Source> suffix_sources_;

    /** Method called by AmplTNLP to prepare the asl for the suffixes */
    void PrepareAmplForSuffixes(ASL_pfgh* asl);

    /** Method called by AmplTNLP to retrieve the suffixes from asl */
    //    void RetrieveSuffixesFromAmpl(ASL_pfgh* asl);

    friend class AmplTNLP;
  };


} // namespace Ipopt

#endif
