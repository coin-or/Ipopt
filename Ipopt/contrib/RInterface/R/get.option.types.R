# Copyright (C) 2010 Jelmer Ypma. All Rights Reserved.
# This code is published under the Eclipse Public License.
#
# File:   get.option.types.R
# Author: Jelmer Ypma
# Date:   18 April 2010
#
# This function converts a list with ipopt options into 
# three sub-lists, where the options are sorted into 
# the different value types (integer, numeric, string).
#
# Input: list of ipopt options and their values
# Output: list containing three sub-lists by type with ipopt options and their values

get.option.types <- function(opts) {

	# define types of ipopt options
    ipopt.option.types <- list(

        # Output
              "print_level"="integer",
              "print_user_options"="string",
              "print_options_documentation"="string",
              "output_file"="string",
              "file_print_level"="integer",
              "option_file_name"="string",

        # Termination
              "tol"="numeric",
              "max_iter"="integer",
              "max_cpu_time"="numeric",
              "dual_inf_tol"="numeric",
              "constr_viol_tol"="numeric",
              "compl_inf_tol"="numeric",
              "acceptable_tol"="numeric",
              "acceptable_iter"="integer",
              "acceptable_constr_viol_tol"="numeric",
              "acceptable_dual_inf_tol"="numeric",
              "acceptable_compl_inf_tol"="numeric",
              "acceptable_obj_change_tol"="numeric",
              "diverging_iterates_tol"="numeric",

        # NLP Scaling
              "obj_scaling_factor"="numeric",
              "nlp_scaling_method"="string",
              "nlp_scaling_max_gradient"="numeric",

        # NLP
              "bound_relax_factor"="numeric",
              "honor_original_bounds"="string",
              "check_derivatives_for_naninf"="string",
              "nlp_lower_bound_inf"="numeric",
              "nlp_upper_bound_inf"="numeric",
              "fixed_variable_treatment"="string",
              "jac_c_constant"="string",
              "jac_d_constant"="string",
              "hessian_constant"="string",

        # Initialization
              "bound_frac"="numeric",
              "bound_push"="numeric",
              "slack_bound_frac"="numeric",
              "slack_bound_push"="numeric",
              "bound_mult_init_val"="numeric",
              "constr_mult_init_max"="numeric",
              "bound_mult_init_method"="string",

        # Barrier Parameter
              "mehrotra_algorithm"="string",
              "mu_strategy"="string",
              "mu_oracle"="string",
              "quality_function_max_section_steps"="integer",
              "fixed_mu_oracle"="string",
              "mu_init"="numeric",
              "mu_max_fact"="numeric",
              "mu_max"="numeric",
              "mu_min"="numeric",
              "barrier_tol_factor"="numeric",
              "mu_linear_decrease_factor"="numeric",
              "mu_superlinear_decrease_power"="numeric",

        # Multiplier Updates
              "alpha_for_y"="string",
              "alpha_for_y_tol"="numeric",
              "recalc_y"="string",
              "recalc_y_feas_tol"="numeric",

        # Line Search
              "max_soc"="integer",
              "watchdog_shortened_iter_trigger"="integer",
              "watchdog_trial_iter_max"="integer",
              "accept_every_trial_step"="string",
              "corrector_type"="string",

        # Warm Start
              "warm_start_init_point"="string",
              "warm_start_bound_push"="numeric",
              "warm_start_bound_frac"="numeric",
              "warm_start_slack_bound_frac"="numeric",
              "warm_start_slack_bound_push"="numeric",
              "warm_start_mult_bound_push"="numeric",
              "warm_start_mult_init_max"="numeric",

        # Restoration Phase
              "expect_infeasible_problem"="string",
              "expect_infeasible_problem_ctol"="numeric",
              "expect_infeasible_problem_ytol"="numeric",
              "start_with_resto"="string",
              "soft_resto_pderror_reduction_factor"="numeric",
              "required_infeasibility_reduction"="numeric",
              "bound_mult_reset_threshold"="numeric",
              "constr_mult_reset_threshold"="numeric",
              "evaluate_orig_obj_at_resto_trial"="string",

        # Linear Solver
              "linear_solver"="string",
              "linear_system_scaling"="string",
              "linear_scaling_on_demand"="string",
              "max_refinement_steps"="integer",
              "min_refinement_steps"="integer",

        # Hessian Perturbation
              "max_hessian_perturbation"="numeric",
              "min_hessian_perturbation"="numeric",
              "first_hessian_perturbation"="numeric",
              "perturb_inc_fact_first"="numeric",
              "perturb_inc_fact"="numeric",
              "perturb_dec_fact"="numeric",
              "jacobian_regularization_value"="numeric",

        # Quasi-Newton
              "hessian_approximation"="string",
              "limited_memory_max_history"="integer",
              "limited_memory_max_skipping"="integer",

        # Derivative Test
              "derivative_test"="string",
              "derivative_test_perturbation"="numeric",
              "derivative_test_tol"="numeric",
              "derivative_test_print_all"="string",
              "point_perturbation_radius"="numeric",

        # MA27 Linear Solver
              "ma27_pivtol"="numeric",
              "ma27_pivtolmax"="numeric",
              "ma27_liw_init_factor"="numeric",
              "ma27_la_init_factor"="numeric",
              "ma27_meminc_factor"="numeric",

        # MA57 Linear Solver
              "ma57_pivtol"="numeric",
              "ma57_pivtolmax"="numeric",
              "ma57_pre_alloc"="numeric",
              "ma57_pivot_order"="integer",

        # MUMPS Linear Solver
              "mumps_pivtol"="numeric",
              "mumps_pivtolmax"="numeric",
              "mumps_mem_percent"="integer",
              "mumps_permuting_scaling"="integer",
              "mumps_pivot_order"="integer",
              "mumps_scaling"="integer",

        # Pardis"Linear Solver
              "pardiso_msglvl"="integer",
              "pardiso_matching_strategy"="string",
              "pardiso_out_of_core_power"="integer",

        # WSMP Linear Solver
              "wsmp_num_threads"="integer",
              "wsmp_ordering_option"="integer",
              "wsmp_pivtol"="numeric",
              "wsmp_pivtolmax"="numeric",
              "wsmp_scaling"="integer",
              "wsmp_singularity_threshold"="numeric"
    )
	
	 
	
	# initialize list with options sorted by type
	converted.opts <- list( "integer"=list(), "string"=list(), "numeric"=list() )
	
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

	# check if we have at least 1 element in the list, otherwise the 
    # loop runs from 1 to down 0 and we get errors
    if ( length( opts ) > 0 ) {

        # loop over all options and give them the correct type
        for ( i in 1:length( opts ) ) {
            tmp.type <- ipopt.option.types[[match( names(opts)[i], names(ipopt.option.types) )]]
            if ( is.null( tmp.type ) ) {
                # determine type
                if ( is.character(opts[[i]]) ) {
                    tmp.type <- "string"
                } else if ( is.wholenumber(opts[[i]]) ) {
                    tmp.type <- "integer"
                } else {
                    tmp.type <- "numeric"
                }
                cat( paste( "Warning: ", names(opts)[i], " is not a recognized option, we try to pass it to Ipopt as ", tmp.type, "\n" ) )
            }
            
            if ( tmp.type=="string" ) {
                converted.opts$string[[ names(opts)[i] ]] <- as.character(opts[[i]])
            } else if ( tmp.type=="integer" ) {
                converted.opts$integer[[ names(opts)[i] ]] <- as.integer(opts[[i]])
            } else if ( tmp.type=="numeric" ) {
                converted.opts$numeric[[ names(opts)[i] ]] <- as.numeric(opts[[i]])
            } else {
                stop(paste("Type of option ", names(opts)[i], " not recognized"))
            }
        }
    }
	
	return ( converted.opts )
}
