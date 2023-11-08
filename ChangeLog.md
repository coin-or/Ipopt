# ChangeLog

Here we list changes of Ipopt since the release of version 2.2.0.
More detailed information about incremental changes can be found in the
[commit history](https://github.com/coin-or/Ipopt/commits).
[TOC]

## 3.14

### 3.14.13 (2023-11-08)

- Reduced priority for making Spral the default value for option linear_solver [#677].
- Adapted to change of integer types in Spral interface. Minimal required Spral version is now v2023.03.29.
- Fixed that return code from `TNLP::eval_jac_g()` was ignored at first call [#697, by
  Christoph Hansknecht].
- Print additional messages when reallocation of MA27 working space failed [#671, by
  Daniel Oliveira].
- Added option `file_append` to specify whether to append to `output_file` instead of
  truncating it. Added default argument `file_append` to `Journalist::AddFileJournal()`,
  added default argument `fappend` to `FileJournal::Open()`, and added default argument
  `file_append` to `IpoptApplication::OpenOutputFile()`. [#720]

### 3.14.12 (2023-04-05)

- Fix that a source file was installed and install more header files.
  [#641, #642, by Joao Sousa Pinto]
- Fixed crash of GetIpoptCurrentIterate() and GetIpoptCurrentViolations() in
  C interface when called before or after IpoptSolve(). [#644, #645, by Robbybp]
- Updated HSL_MA97 header file to the one from HSL MA97 2.8.0 [#646, by Jari Fowkes].
- Fixed crash when trying to solve problem without variables and constraints. [#648]
- Added optional argument to `AlgorithmBuilder` constructor to provide name of custom solver. [#618]
- Changed handling of dual solution for square problems: When solving a problem
  with as many equations as variables, Ipopt used to ignore the violation of
  dual feasibility and complementarity in the convergence check and computed
  a final dual solution via a least-square estimate. If this failed, Ipopt would
  claim a square problem to be solved to optimality without providing a solution
  that satisfies optimality conditions. With this version, the behavior has
  been changed so that dual feasibility is no longer ignored by the convergence
  check, the least-square estimate is only computed if optimality is not proven
  already, and the normal Ipopt algorithm continues if the least-square estimate
  does not satisfy optimality conditions.
- Updated HSL_MC86 header file to the one from HSL MC68 3.3.3 [#653, by Jari Fowkes].

### 3.14.11 (2023-02-07)

- Added `IpoptData::TimingStats() const` [#611]
- Assume DLL library extension in linear solver library loader on Windows
  also when building with other compiler than MSVC/Intel [#628].
- Updated buildsystem files after upgrading to most recent versions of autotools.
- Install some additional header files [#637].

### 3.14.10 (2022-10-11)

- Added option `grad_f_constant` to specify that objective function is linear.
  If set, the gradient of the objective will be requested by Ipopt only once. [#597]
- Added `OrigIpoptNLP::orig_d_L()` and `OrigIpoptNLP::orig_d_U()` to get
  original constraint sides (before relaxation due to bound_relax_factor > 0).
- `TNLP::get_curr_violations()` now returns the constraint violation and
  complementarity w.r.t. the original (non-relaxed) constraint sides. [#603]

### 3.14.9 (2022-07-21)

- Fixed mapping of meta data for variable bounds, e.g., variable names,
  from TNLP to Ipopts internal NLP [#590].

### 3.14.8 (2022-07-13)

- Added options ma27_print_level, ma57_print_level, and mumps_print_level
  to enable output from these linear solvers.

### 3.14.7 (2022-06-24)

- Fixed that ComputeSensitivityMatrix() of sIpopt assumed that there are
  no more than two parameters [#578, by Andrea Vescovini].
- For completeness, added option `gradient_approximation` to enable approximation
  of gradient of objective function by finite differences. Do not use. [#573]
- Added function `IPSETPROBLEMSCALING` to Fortran interface to set problem
  scaling [#577, by Steven R. Hall]

### 3.14.6 (2022-05-02)

- Fixed mapping of meta data for inequalities, e.g., constraint names,
  from TNLP to Ipopts internal NLP [#570, by Daniel Oliveira].
- Fixed that MC68 ordering time was not accounted in symbolic factorization
  time of HSL MA86 [#571].
- Include more header files in IpIpoptCalculatedQuantities.hpp for setups
  where forward declarations are not sufficients [#572].

### 3.14.5 (2022-02-09)

- Tried to fix recognition of JNI headers on macOS >= 11 [#516].
- Fixed that only primal variable values where passed to `finalize_solution()`
  when a timelimit was reached [#552].

### 3.14.4 (2021-09-20)

- Skip build of Java interface if either java or jar is not found [#510].
  Only give warning if javac and jar are found, but no java or javadoc.
- Fixed that `--with-lapack-lflags` was ignored if `--with-lapack` was not
  specified explicitly [#512,#515].

### 3.14.3 (2021-09-03)

- Fixed timing for iterate initialization if initialization fails due to
  an evaluation error.
- Fixed possible integer overflow when reserving space for indices of Jacobian
  belonging to fixed variables (introduced with 3.14.0) and reduced memory
  usage for indices of Jacobian belonging to fixed variables.

### 3.14.2 (2021-07-21)

- Added `OptionsList::UnsetValue()` to remove an option setting.
- Added missing translation of some Ipopt status codes into AMPL solve result codes.
- If using the MPI-parallel version of MUMPS: Moved calls to
  `MPI_Init()`/`MPI_Finalize()` in MUMPS interface into global constructor/destructor
  of Ipopt library (if building with GCC/clang). Use configure flag
  `--disable-mpiinit` to disable. [#500]

### 3.14.1 (2021-06-25)

- Fixed invalidation of cached Hessians when reoptimizing with same structure.
- Added `OptionsList::SetBoolValue()` and `OptionsList::SetBoolValueIfUnset()`. [#492]
- Skip check for and link against libdl if linear solver loader is disabled.
- Fixed missing initialization of `BacktrackingLineSearch::in_watchdog_`.
- Fixed a problem with the current solution not being reset when
  initialization of a NLP fails in reoptimization.
- Fixed that C++11 capability was not correctly identified with MS/Intel compilers.

### 3.14.0 (2021-06-15)

#### Data Types

- Due to a contribution by Mitchell Clement [#428], it is now possible
  to build Ipopt in a variant that uses single-precision floating
  point arithmetic instead of the default double-precision.
  This can be enabled by specifying the configure flag `--with-precision=single`.
  Doing so has a number of consequences on Ipopt interfaces and Ipopt dependencies.
  See the Ipopt installation instructions for more details.
- The name of all functions in `IpBlas.hpp` and `IpLapack.hpp` has been
  changed to drop the precision-specifier "D" in the name.
  Wrapper using the old names are available, but their use is deprecated.
- It is now possible to build Ipopt in a variant that uses 64-bit
  integers for its `Index` type. This can be enabled by specifying
  the configure flag `--with-intsize=64`. Doing so has a number of
  consequences on Ipopt interfaces and Ipopt dependencies. See
  the Ipopt installation instructions for more details. [#259]
- Added new header `IpTypes.h` that defines Ipopt types for integer
  and real numbers for C users.
- Deprecated almost never used type `Ipopt::Int`, use `int` instead.
- Deprecated `Number`, `Index`, and `Int` (on global namespace) in
  `IpStdCInterface.h`, use `ipnumber`, `ipindex`, and `int` instead, respectively.
- Deprecated `Bool` in `IpStdCInterface.h` and replaced it by `bool`
  from `stdbool.h`. Note, that `Bool` was defined to be `int`, but `bool`
  is likely a shorter type, e.g., `signed char`. Deprecated `TRUE` and
  `FALSE`, use `true` and `false` instead.
- Deprecated `IPOPT_EXPORT` macro and introduced `IPOPT_CALLCONV`.
- Deprecated `IPOPT_FORTRAN_INTEGER_TYPE` and `ipfint`. These were
  always assumed to be `int`. Use `ipindex` or `ipopt::Index` instead.

#### Linear Solver Interfaces

- Due to a contribution by Byron Tasseff [#446], it is now possible to use
  the linear solver [SPRAL](https://github.com/ralna/spral) (Sparse Parallel
  Robust Algorithms Library) with Ipopt.
  SPRAL is open-source and can, optionally, make use of NVIDIA GPUs.
  If Ipopt has been build with SPRAL, then option `linear_solver` can be
  set to `spral` to enable use of SPRAL.
  See the installation instructions on how to build the Ipopt/SPRAL interface
  and the options documentation for new options that are available for
  the Ipopt/SPRAL interface.
- Added `IpoptLinearSolvers.h` with function `IpoptGetAvailableLinearSolvers()`
  to retrieve information which linear solvers are available for Ipopt
  (linked in or loaded at runtime). Options `linear_solver` and
  `linear_system_scaling` can now only be set to values which corresponding
  code is available (linked in or for load at runtime).
- Revised and streamlined implementation of feature that loads libraries
  with HSL or Pardiso routines at runtime. Added options `hsllib` and `pardisolib`
  to specify name of of HSL and Pardiso libraries, respectively.
- Separated interfaces for Pardiso from [pardiso-project.org](https://pardiso-project.org)
  and Pardiso from Intel MKL and allow to have both Pardiso versions available
  with the same Ipopt libraries.
  - Pardiso from pardiso-project.org can be selected via `linear_solver=pardiso`.
    configure option <b>`--with-pardiso` should now specify only a Pardiso
    library to be loaded at runtime</b> (the value for `--with-pardiso` decides
    the default value for option `pardisolib`).
    To avoid conflicts with Pardiso from MKL, it is **no longer possible to
    link against Pardiso from pardiso-project.org**.
  - Pardiso from MKL can be selected via <b>`linear_solver=pardisomkl`</b>.
  - Options that influence Pardiso from MKL are named <b>`pardisomkl_msglvl`, `pardisomkl_order`, etc.</b>
  - configure option `--disable-pardisomkl` can be used to disable the check
    for MKL Pardiso. [#454]
- Fixed that return code (`info%flag`) of `ma97_solve()` was not checked
  in `Ma97SolverInterface::MultiSolve()`.
- Calls into MUMPS are now protected by a mutex if compiling for C++11 or higher.
  This prevents Ipopt from calling MUMPS concurrently from several threads.
- An insufficient memory return status is now also given if the memory
  required for a working space array of a linear solver exceeds the maximal
  value for the working space length variable, e.g., if MA27 requires a
  working space array of length higher than 2^31.
- Avoid floating point overflow when computing memory increase in interfaces
  to MA27 and MA57.

#### Algorithm

- Changed treatment of NLPs with all variables fixed and `fixed_variable_treatment`
  set to `make_parameter`:
  Ipopt will still terminate early, but initialize its data structures
  first and print a solve summary (49a0d3a56).
  Changed treatment of NLPs with inconsistent variable bounds or constraint
  sides: These will now result in an invalid problem definition error (-11)
  (5cdb2624bb).
  As a consequence, solve statistics should now always be available if the
  return status from `(Re)Optimize()` is larger than -10.
- Bound multipliers are now computed for fixed variables if
  `fixed_variable_treatment` is set to `make_parameter` (the default).
  This can trigger a reevaluating of the gradient of the objective
  function or the Jacobian of the constraint functions. If this
  is not desired, then option `fixed_variable_treatment` can be set
  to the new value `make_parameter_nodual`. [#308]
- Changed default for `honor_original_bounds` to no.
  Variable bounds (see option `bound_relax_factor`) are now relaxed
  by at most the value of `constr_viol_tol`.
  The solve summary now reports the violation of original bounds
  and `SolveStatistics::(Scaled)Infeasibilities()` has been extended
  to report the violation of original variable bounds. [#312]
- If Ipopt hits a time or iteration limit during watchdog phase,
  the iterate from before the watchdog phase is now restored and
  passed to `finalize_solution`. Note that this does not apply if
  Ipopt is stopped due to a user-interrupt (`intermediate_callback`).
  [#289]
- When `BacktrackingLinesearch` could not find a trial point that
  provided sufficient progress, it may resort to call the feasibility
  restoration phase. If, however, the current (scaled) constraint
  violation was below tol/100, this would result in the infamous
  "Restoration phase is called at almost feasible point..." abort
  with status code `Restoration_Failure`. Now, this message has been
  changed to "Linesearch failed, but no restoration phase or other
  fall back is available." and the status code to `Error_In_Step_Computation`
  to better reflect that the linesearch failed and not the restoration
  phase. Further, the unscaled constraint violation tolerance now
  needs to be below `constr_viol_tol/10` as well in order to trigger
  this abort.
- When a square problem is solved and the restoration phase only succeeded to
  find a point that is feasible w.r.t. constr_viol_tol, but not w.r.t. tol,
  then status Feasible_Point_Found is returned now.

#### Ipopt interfaces

- Added `TNLP::get_curr_iterate()` and `TNLP::get_curr_violations()`
  to request the current iterate (primal and dual variable values)
  and primal and dual infeasibility w.r.t. the TNLP. The methods
  are meant to be called during `intermediate_callback` to inspect
  the current iterate. The C, Fortran, and Java interfaces have
  been extended by corresponding functions, too. The hs071 examples
  have been extended to show use of the new functions. Added test
  `getcurr` to test new functions. [#382, #451]
- Extended the Java interface by the possibility to specify an `intermediate_callback`.
- Added flag to `TNLPAdapter::ResortX()` to specify how to handle fixed variables.
  Added flag to `TNLPAdapter::ResortG()` to specify whether to correct by
  right-hand-side of equality constraints.
  Added `TNLPAdapter::ResortBoundMultipliers()` to generate correct duals
  for fixed variables if `fixed_variable_treatment` is `make_constraint`.
  Added `TNLPAdapter::GetFullDimensions()`, `TNLPAdapter::GetFixedVariables()`,
  `TNLPAdapter::GetPermutationMatrices()`, and `TNLPAdapter::GetC_Rhs()`
  to retrieve more information of the transformation in a `TNLPAdapter`.
- Added `TNLPAdapter::ResortBounds()` and deprecated `TNLPAdapter::ResortBnds()`.
- `AmplTNLP` constructor and `AmplTNLP::get_options()` now expect to
  receive a pointer to a `RegisteredOptions` object as well.
  Previous versions of these methods are still available but deprecated.
  Using a NULL pointer for the `RegisteredOptions` argument is
  also deprecated.
- `std::overflow_error` exceptions are now caught by Ipopt even if rethrowing
  of non-Ipopt exceptions is enabled.

#### Timing

- Added option `max_wall_time` to specify a wallclock time limit.
  Added `SolverReturn` code `WALLTIME_EXCEEDED` and `ApplicationReturnStatus`
  code `Maximum_WallTime_Exceeded`.
- Changed default for `max_cpu_time` to 1e20 to indicate no CPU timelimit.
- Detailed timing statistics are no longer collected by default.
  To reenable, set the new option `timing_statistics` to yes or set the
  option `print_timing_statistics` to yes [#299].
- Removed `IpoptData::ResetCpuStartTime()`. Deprecated `IpoptData::cpu_time_start()`,
  use `IpoptData::TimingStats()::OverallAlgorithm()::StartCpuTime()` instead.
- Deprecated `SolveStatistics::TotalCPUTime()`, use `SolveStatistics::TotalCpuTime()`
  instead.

#### Option handling

- All Ipopt options are now also available via the AMPL interface,
  that is, can be set in AMPL via an option `ipopt_options "..."` statement.
- Added flag "advanced" for `RegisteredOption` to store whether an
  option is rather meant for expert users.
- Added class `RegisteredCategory` to store information on category
  of registered option. Next to the name, this is a priority, which
  determines the order in which categories are printed. As a
  consequence, methods `RegisteredOption::RegisteringCategory()` and
  `RegisteredOptions::RegisteringCategory()` now return a
  `RegisteredCategory` object instead of a string. Further,
  `RegisteredOption::SetRegisteringCategory()` has been removed.
- Deprecated existing `RegisteredOptions::OutputXyz()` methods and
  added new variant of `RegisteredOptions::OutputOptionDocumentation()`.
  Added parameter `print_advanced_options` to indicate whether to
  print documentation for advanced options, too.
- Added previously undocumented advanced options to options reference.

#### Miscellaneous

- Updated buildsystem, including improved checks for dependencies,
  use of current autotools, and skipping the build of intermediate
  non-distributed libraries.
- Fixed missing initialization of `IpoptApplication::smart_jnlst` in
  second `IpoptApplication` constructor, if `ipopt verbosity > 0`.
- Add GCC format attribute to `Journalist::Printf` functions to enable
  `printf`-formatter check.
- Fixed `DenseVector::SumLogsImpl()` such that it returns 0.0 for
  a vector of dimension 0. Returned `nan` for homogeneous 0-vector
  of dimension 0 before, which may have caused the restoration phase to
  fail for problems with only equality or only inequality constraints.
  Also other `DenseVector` methods now skip calculations when dimension
  is 0 to avoid (probably harmless) divisions by zero.
- Fixed a problem where moving slack away from 0 did not succeed
  when mu was very small. [#212]
- Fixed a problem where moving slacks away from 0 resulted in `nan`
  if multipliers were zero. Added `Vector::ElementWiseSelect()`.
- Various tiny bugfixes and improvements in performance and code style
  by following suggestions of `cppcheck`.
- Added documentation on some available C preprocessor flags for expert users.
- Fixed static build of sIpopt without GCC. Fixed that installed sIpopt
  headers were not usable (`SIPOPTLIB_EXPORT` not defined).
- Fixed wrong gradient of objective function and Lagrangian Hessian in
  restoration problem [#478, by Nai-Yuan Chiang].
- If Ipopt is compiled for checklevel 2 or higher and the GLIBC extension `feenableexcept()`
  is available, then floating-pointing exceptions divbyzero, overflow, and
  invalid are raised while `IpoptAlgorithm::Optimize()` is running.
- Fixed that norm on unscaled complementarity or scaled complementarity tolerance
  were negative when maximizing (by using a negative scaling factor for the
  objective).
- Changed formula for relative error in derivative checker. The absolute
  error is now scaled up if the approximate derivative value is between
  derivative_test_tol and 1. [#487].
- The second-order derivative checker now uses values for obj_factor and
  lambda that are different from 1.


## 3.13

### 3.13.4 (2021-02-24)

- Fixed a linking issue for `ipopt_sens` [#418]
- Fixed Makefile for Java example regarding location of jar file
- Fixed build of R interface if using `-fvisibility=hidden`.

### 3.13.3 (2020-10-16)

- Members of `AmplTNLP` class are now protected instead of private.
- Updated Eclipse Public License from 1.0 to 2.0.
- Fixed dangling pointer problems with Journalist used for debugging
  (`--with-ipopt-verbosity > 0`) when more than one `IpoptApplication`
  is used. [#393, thanks to Brad Bell]
- Fixed build problem when using HSL library that does not include
  MA27, MA57, or MC19. [#395]
- Added example `recursive_nlp` that uses Ipopt to solves an
  optimization problem for the evaluation of the objective function.
  [contributed by Brad Bell]
- Fixed build of linear-solver loader on Windows [#408]

### 3.13.2 (2020-04-30)

- The C-preprocessor defines `COIN_IPOPT_CHECKLEVEL`,
  `COIN_IPOPT_VERBOSITY`, and `FORTRAN_INTEGER_TYPE`, which are defined
  by `IpoptConfig.h`, have been renamed to `IPOPT_CHECKLEVEL`,
  `IPOPT_VERBOSITY`, and `IPOPT_FORTRAN_INTEGER_TYPE`, respectively.
  They are still available under their previous name, but these
  will be removed in Ipopt 3.14.
- Changed dependencies as used by coinbrew to use new versions (2.1)
  of ThirdParty/HSL and ThirdParty/MUMPS and dropped ThirdParty/Metis.
  The new versions of the HSL and MUMPS build scripts now look
  for a Metis library in the system and should work with both
  Metis 4 and Metis 5.
- Changed location where Java interface jar gets installed from
  `$libdir` to `$datadir/java/`.
- minor fixes to buildsystem

### 3.13.1 (2020-03-11)

- Added asserts that check whether sparsity pattern of Jacobian
  and Hessian as returned by TNLP are within range w.r.t. number
  of variables and constraints. [#350]
- `TNLPAdapter::ResortBnds` now initializes complete output arrays
  with 0.0 before filling in values for non-fixed variables. Use
  new argument `clearorig` to turn this off. [#352]
- bring back configure variables `ADD_{C,CXX,F}FLAGS`
- added configure option `--enable-relocatable` to make prefix in
  pkg-config files relative to pcfiledir (assuming that `--libdir`
  hasn't been set)
- bring back `configall_system.h` for build without config header
- minor fixes to buildsystem

### 3.13.0 (2019-10-19)

This major release comes with a larger renovation of the build
system and a changed directory structure (eliminated top directory),
which is the result of a long and still on-going effort to use
recent autotools versions for various COIN-OR projects, reduce
future maintenance efforts, and adapting behaviors of standard
autotools-based projects.
As a consequence, a monolithic build of Ipopt, which builds Ipopt
with all its dependencies in one run of configure and make is no
longer possible. Dependencies should now be build and installed
before building Ipopt.
Additionally, support for some outdated versions of dependencies
and unmaintained components of Ipopt has been dropped and some
improvements that may require changes on the users side have been
applied.

A more detailed, probably incomplete, list of changes follows:
- Removed git submodules. Dependencies (HSL, Mumps, ASL, etc) now
  need to be build and installed in advance, either manually or
  by using coinbrew.
- Dropped support for HSL < 2013.
- Dropped support for MA28 in the linear solver loader.
- Dropped support for Pardiso < 4.0 from pardiso-project.org.
- Added support for Mumps 5.2.x, though initial experiments on
  CUTEst indicated that, on average, performance is worse than
  when using Mumps 4.10.0.
- Dropped CUTEr interface, the successor CUTEst includes an
  interface to Ipopt.
- Dropped Matlab interface as it is unmaintained and it was
  reported that it stopped functioning.
  Use https://github.com/ebertolazzi/mexIPOPT instead.
- Dropped MSVS project files as unmaintained and not functioning
  with current Ipopt anymore.
- Integrated Java interface into the main Ipopt library, that is,
  it is handled equivalently to the C and Fortran interfaces:
  - The source moved into `src/Interfaces`.
  - The JNI functions are now included in the main Ipopt library,
    thus an extra jipopt library is no longer build or necessary.
  - The Java class and `org.coinor.ipopt.jar` package are build and
    installed as part of the main Ipopt build.
  - The examples moved into `examples/*_java`.
  - A Java interface test is executed by `make test`.
  - To build javadoc, run `make javadoc` in the main build directory.
  - The configure flag `--disable-java` can be used to disable the
    check for Java and build of the Java interface.
  - `DLLPATH` and `DLLNAME` have been removed from the Ipopt class and
    constructors that works without arguments and with only one
    argument (specifying the Ipopt library namestem) have been added.
  - Method `Ipopt::finalize` has been marked as deprecated and will
    be removed in some future Ipopt version. Users must call
    `dispose()` explicitly.
- Integrated sIpopt into the main Ipopt build, that is, it is now
  build together with Ipopt, but installed as separate library
  and executable. Use `--disable-sipopt` to disable building sIpopt.
- `IPOPT_THREAD_LOCAL` now uses C++11's `thread_local` keyword if C++11
  is available.
- When using a GCC-compatible compiler, Ipopt and sIpopt interface
  functions are now declared with `visibility(default)`-attribute,
  thus building Ipopt with `-fvisibility=hidden` still produces a
  usable library.
- When using a MSVC-compatible compiler, Ipopt and sIpopt interface
  functions are now declared with `dllimport`-attribute, so that an
  Ipopt C++ DLL can be used.
- Under Windows/Msys2, DLLs are now build by default.
- Cygwin and MSys1 are not supported.
- pkg-config is now mandatory to use dependencies like ASL or HSL.
  On Windows, make sure to use a pkg-config version that produces
  Unix-style paths.
- Script "`compile`" is now used to wrap around calls of cl/icl/ifort
  and translate GCC-style compiler flags to MSVC style.
- "Addlibs" files have been removed, pkg-config should be used instead.
- Header files are now installed in the better named
  `$prefix/include/coin-or` instead of `$prefix/include/coin`.
- The default for `--prefix` is no longer the build directory, but
  the autotools-default, probably `/usr/local`.
- The check for a Fortran compiler can be disabled via `--disable-f77`
  and Ipopt can easier be build without a Fortran compiler.
- Lapack is no longer optional, but required. The separate check
  for Blas and the `--with-blas` flags have been removed.
- `--enable-debug` does not imply `--disable-shared` anymore.
- Removed `--enable-debug-ipopt`, use `--enable-debug` instead.
- Removed configure variables `{ADD,OPT,DBG}_{C,CXX,F77}FLAGS`.
  Use `{C,CXX,F77}FLAGS` instead.
- Silent build output is now enabled by default, use configure
  flag `--disable-silent-rules` or call make with `V=1` to disable.
- Also for static builds, PIC objects are now generated by default,
  use `--without-pic` to disable.
- The `--with-*-incdir` and `--with-*-lib` configure flags have been
  replaced by corresponding `--with-*-cflags` and `--with-*-lflags`
  flags. Note that the include directories need to be specified
  via `-I<dir>` in `--with-*-cflags`.
- Fixed handling of `ma77_default_control` in `LSL_setMA77()`.
- Fixed calculation of quality function when setting option
  `quality_function_centrality` to `reciprocal`.
- Fixed compiler warnings, in particular when using `-Wunused-parameter`.
- Changed default for `ma97_print_level` to -1. This avoids messages
  about numerical singular systems written to stdout by default.

## 3.12

### 3.12.13 (2019-04-08)

- fixed Pardiso settings when using Pardiso from Pardiso project
  website (by Olaf Schenk): the new settings should provide much
  better performance; the default for option `pardiso_order` changed
  from `five` to `metis`.
- changed distinction of MKL and Basel Pardiso in configure: to
  use MKL Pardiso, only specify MKL for Blas; to use Basel Pardiso,
  use `--with-pardiso`

### 3.12.12 (2018-11-17)

- allow for `--without-matlab-home` to disable check for Matlab [r2748]
- add `dppsv` to `v8-ifort` [r2746]
- disable error in `LibraryHandler.c` if `snprintf` detection failed [r2751]

### 3.12.11 (2018-09-16)

- fill MUMPS struct with zeros when allocating in MUMPS interface [r2724]
- minor fix in build-system of ThirdParty/ASL

### 3.12.10 (2018-06-02)

- fixed setting for parallel solve when using MKL Pardiso
  (by t1393988511) [r2711]: parallel solve was disabled (which
  is not the default); note, that the setting for parallel
  factorization was not affected
- fixed invalid read in AMPL interface for problems without
  objective function [r2715, #305]
- updated ThirdParty/ASL to retrieve updated ASL (20180528) [#305]
- name JIpopt library `libjipopt.dylib` on Mac OS X [r2718, #275]

### 3.12.9 (2018-01-15)

- fixed memory leak in MA86 interface (by mhahn) [r2700,#283]
- fixed handling of time limit when reoptimizing: CPU time spend
  was accumulated when reoptimizing, while it should have been
  reset for each solve (by paul-scott) [r2702,r2703]
- fixed sign in Jacobian finite-difference approximation when point
  was close to variable upper bounds (by Enrico Bertolazzi) [r2704]

### 3.12.8 (2017-06-12)

- add define for `FORTRAN_INTEGER_TYPE` to `config_ipopt_default.h`
- `IpoptApplication::RethrowNonIpoptException()` now returns whether
  non-ipopt exceptions were rethrown before the method was called.

### 3.12.7 (2017-02-25)

- removed compiler flag `-pedantic-errors` to avoid problems with some
  configure tests when using recent GCC versions
- fixed rare bug in handling variable/constraint names in `AmplTNLP`
  (by G. Hackebeil) [r2673]
- the `get.Mumps` script in ThirdParty/Mumps now renames `libseq/mpi.h`
  to `libseq/mumps_mpi.h` to avoid conflicts when building in a MPI
  environment (by T. Ralphs); note that if updating an existing
  checkout/download of Ipopt, you may have to rerun get.Mumps

### 3.12.6 (2016-07-20)

- better support for custom algorithm development [r2659] (by Gabriel Hackebeil):
  &quot;Reorganization of the `AlgorithmBuilder` class to allow easier
  customization of the Ipopt algorithm. In particular, we wanted to
  make use of the code that creates the `SymLinearSolver` object to
  implement our own `SymLinearSolver` without copy-pasting everything
  in AlgorithmBuilder.
  `AlgorithmBuilder::BuildBasicAlgorithm` now consists of 8 method calls
  that build the core components passed into the arguments of the
  `IpoptAlgorithm` class. These calls are ordered based on any dependencies
  they might have. In addition, all code for creating the `PDSystemSolver`,
  `AugSystemSolver`, and `SymLinearSolver` has been moved into separate factory
  methods.
  Also, included is a change to install a few more header files with Ipopt.
  Some of these are required to subclass AlgorithmBuilder, and the others
  are simply some matrix types that we require.&quot;
- extend build system to work without Fortran compiler [r2660,r2661]:
  If no Fortran compiler is available (`F77=unavailable`), then
  the build system checks for functions in Blas, Lapack, and
  Pardiso via C linkage. This seems to work when using the Intel MKL,
  thus allowing to build Ipopt with C/C++ compilers and MKL only.
  The linear solver loader and the CuteR interface are disabled when
  no Fortran compiler is available. A user may have to adjust the
  definition of `F77_FUNC` in `Ipopt/src/Common/IpoptConfig.h`.

### 3.12.5 (2016-04-30)

- changed `fptr` from `long` to `void*`: the Fortran side needs to
  make sure that it uses a big enough integer type to store a
  C pointer, thus `void*` can be used on the C side [r2599]
- added additional second-order-correction method, which can be
  selected by setting the new option `soc_method` to 1 (by Wei Wan)
  [r2606, r2607]
- added parameter `allow_clobber` with default value false to
  `IpoptApplication::Initialize()` and `OptionsList::ReadFromStream()`

### 3.12.4 (2015-08-09)

- option to use regularized Hessian when doing a curvature test
  without inertia information (`neg_curv_test_tol` > 0), new
  option `neg_curv_test_reg` to switch back to original behavior
  (by N.-Y. Chiang and V. Zavala Tejeda) [r2579]
- sIpopt: Added access to sensitivity directional derivative
  vector (`ds/dp*(p-p0)` Eq. 14 sIpopt implementation paper). Also,
  added an option to compute the sensitivity matrix and provide
  access to it. Finally, added an example that shows how to
  access the new information. (by R. Lopez-Negrete)
- use workaround for failing check for random number generator
  with any gcc 4.8.x, x >= 2

### 3.12.3 and 3.11.11 (2015-04-15)

- fixed bug in MA97 interface that lead to conversion issues
  (by J. Hogg) [r2566, #260]

### 3.12.2 (2015-04-04)

- revised integration of doxygen-generated documentation into build system
  (by T. Ralphs)

### 3.12.1 (2015-02-13)

- fixes to build system for dependency linking and library versioning
- Ipopt will now report an NLP with inconsistent variable bounds
  or inconsistent constraints sides as infeasible instead of
  throwing an invalid TNLP exception (by T. Kelman) [r2548]

### 3.12.0 (2015-01-23)

- Library dependencies are now recorded in shared library builds,
  which is intended to simplify linking against the Ipopt library.
  However, the pkg-config and `ipopt_addlibs` files do not reflect
  this change yet (it is rather experimental, imho). To restore
  the previous behavior, use `--disable-dependency-linking` as
  configure option.
- If linking against Intel MKL for Blas/lapack, use of Pardiso
  from MKL is now automatically enabled. Note, that this will
  change the default solver on Ipopt builds without any of the
  linear solvers MA27, MA57, MA97, and MA86 (these take preference
  over Pardiso). [#216]
- dropped support for old HSL sources (<2013; ThirdParty/HSLold)
- updated ASL sources, now downloaded from AMPL-MP (github)
- some internal changes to data structures (improved use of compound
  component spaces) and addition of `IpLapackDppsv` (by Gabe Hackebeil)

## 3.11

### 3.11.10 (2015-01-18)

- fix a memory allocation in Java interface in cases where `jint`
  has a different size than `int` [r2513]
- the buildsystem now tries the `Accelerate` framework instead of
  `vecLib` for finding Blas/Lapack on MacOS X

### 3.11.9 (2014-08-16)

- fix compilation issue of Java interface on systems where `Index`
  and `jint` have different size [r2498, #241]
- work around failing check for random number generator with gcc
  4.8.3 [r2495, r2496]
- readded `IpTaggedObject.cpp` to list of sources to compile in
  MSVS `v8-ifort` project file [r2492]
- work around missing support for thread-local storage with gcc < 4.5
  on MacOS X [r2491, #243]
- fix call to MKL Pardiso init function [r2489]
- speed up Triplet to CSR converter [r2487, #234]
- fixed a bug in equilibration scaling where average values were
  computed incorrectly (by V. Zverovich) [r2483]

### 3.11.8 (2014-04-08)

- fixed a bug, introduced with Ipopt 3.11.0, where the tag in the
  Ipopt's caching mechanism was not unique over time, which lead
  to failures of Ipopt that were difficult to debug or recognize
  (e.g., Ipopt may have stopped with an restoration failure for
  instances that solved fine with Ipopt 3.10) [r2472, r2473]
  I'm very thankful to Gabriel Hackebeil and Kurt Majewski for
  their debugging effort on this issue.
- building Mumps with pthreads is now disabled by default [#229]
- fixed setting of `LD` on Windows (now set to link only iff using
  MS/Intel compilers) [#230]
- fixed download link for Gnumex [r2471]
- for some messages about too-few-degrees-of-freedom and restoration
  failure, the message level changed from error to strong-warning
  [r2460, r2469]
- revised calls to `MPI_Init` and `MPI_Finalize` in MUMPS interface [r2467]
  (`MPI_Init` is now called only if function `MPI_Initialized` is available
  and MPI has not been initialized already; `MPI_Finalize` is only called
  if Ipopt also called `MPI_Init`; ...)

### 3.11.7 (2013-12-18)

- adapted PARDISO parameters when using MKL PARDISO to be close
  to using Basel PARDISO
- added options `pardiso_max_iterative_refinement_steps` and
  `pardiso_order`; the former defaults to 1 in case of MKL PARDISO,
  which may help on instances that otherwise fail due to numerical issues
- removed duplicate code in `IpQualityFunctionMuOracle.cpp` [#225, r2445]
- fixed bug in triplet to csr converter [#226, r2446]
- minor changes in buildsystem

### 3.11.6 (2013-11-16)

- updates to Matlab Interface build system (by T. Kelman)
- fix to updates of R Interface [r2416, #223]
- fixed `SHAREDLIBEXT` in `v8-ifort`'s `config.h` [r2426, #224]
- minor fixes to the buildsystem

### 3.11.5 (2013-10-26)

- added method `IpoptApplication::RethrowNonIpoptException()` to enable
  rethrowing of non-ipopt and non-`bad_alloc` exceptions catched in
  the `*Optimize()` and `Initialization()` methods; default is still to
  return with `NonIpopt_Exception_Thrown` status
- minor fixes to the buildsystem [#215, #222]

### 3.11.4 (2013-09-12)

- hopefully fixed non-working linear solver loader in DLLs build with
  MSVS/`v8-ifort` project files [r2365]
- allow MC19 to be loaded via linear solver loader (by J. Currie) [r2366]
- fixed new point flag when running dependency detector [r2368]
- experimental: adapt Pardiso interface to work with MKL Pardiso
  (by J. Currie, T. Kelman) [r2369, #216]:
  - in a few tests it has been found that Pardiso from Intel MKL nowadays
    seems to work fine with Ipopt
  - to use Intel MKL with Ipopt 3.11, one has to specify the MKL libs via
    `--with-pardiso` and add `-DHAVE_PARDISO_MKL -DHAVE_PARDISO_PARALLEL`
    to the compiler flags
  - note that this is still an experimental feature (and thus not enabled
    by default)
- updated Ipopt/R interface to version 0.8.4 [r2373]
  - additional variables have been included in the object returned from `ipoptr`:
    - `z_L`: final values for the lower bound multipliers
    - `z_U`: final values for the upper bound multipliers
    - `constraints`: final values for the constraints
    - `lambda`: final values for the Lagrange multipliers
  - removed `ipoptr_environment` as argument in `ipoptr` (see also r2372)
- fixed bug in penalty term for centrality in quality function (not used by
  default) [#219, r2374]
- minor bugfixes in AMPL interface, debug print statements, and compound matrix
  (by G. Hackebeil) [#218, r2371, r2377, r2378, r2379]
- download scripts for ASL, Blas, and Lapack now first try to download tarball
  copies from the COIN-OR server

### 3.11.3 (2013-08-08)

- `get.*` scripts for ThirdParty/{ASL,Blas,Lapack} now work around broken
  ftp access to www.netlib.org.

### 3.11.2 (2013-07-01)

- changed default for option `option_file_name` to `ipopt.opt`; specifying an
  empty string for this option now disables reading of an option file [r2339]
- missing initial values are now set to 0.0, projected onto variable bounds,
  in AMPL interface [r2340, #205]
- fixed missing variable initialization in MA97 interface [r2341, #206]

### 3.11.1 (2013-06-14)

- the setup for the `v8-ifort` MSVS project changed to use dynamic runtime
  DLLs instead of static linking, which caused crashes in debug mode
  (by M. Roelofs) [r2301]
- fixed memory leaks in Java Interface (by javier) [#200, r2312]
- updates and fixes to MA77 and MA87 interfaces, adding support of
  HSL 2013 codes (by J. Hogg);
  HSL 2012 still supported when compiled with Ipopt, but the linear solver
  loader to dynamically load a HSL library at runtime now assumes HSL 2013
- added option `ma97_solve_blas3` (by J. Hogg) [r2329]
- changed default for option `ma27_meminc_factor` from 10.0 to 2.0 [r2330]
- fixed bug in `ipopt_auxdata` of MATLAB Interface related to `iterfunc` [r2325]

### 3.11.0 (2013-05-07)

#### Miscellaneous

- update and extension of Ipopt documentation
- updated build of doxygen-generated documentation to comply with other
  COIN-OR projects
- localized global variables in `TaggedObject` and `RegisteredOption`,
  so that Ipopt should now be threadsafe as long as Ipopt objects
  (esp. `SmartPtr`'s) are not shared between threads and a threadsafe
  linear solver is used (e.g., MA27) [#167]
- no more need for whitespace character at end of options file
- added options `print_frequency_iter` and `print_frequency_time` to regulate
  which iteration summary lines should be printed [#161]
- function evaluation timings are now available in `OrigIpoptNLP` [#86]
- some fixes to uncommon issues with the Ipopt `SmartPtr` [#162]

#### Linear Solver Interfaces

- new build system for Harwell codes (ThirdParty/HSL), which requires
  the coin-hsl archives from http://www.hsl.rl.ac.uk/ipopt/;
  previously downloaded HSL sources can still be used by placing them
  into ThirdParty/HSLold, but this option will be removed in a future
  Ipopt version
- new interfaces for Harwell codes HSL_MA77, HSL_MA86, and HSL_MA97;
  see http://www.hsl.rl.ac.uk/ipopt/ about help on when to use which solver;
  especially MA57 and HSL_MA97 should be considered as replacement for MA27;
  however, MA27 is still the default for now
- changed default of `ma57_automatic_scaling` to `no` (faster in general,
  but for higher reliability, you may want to set this option to yes)

#### Ipopt Interfaces

- major improvements for building the MATLAB interface (see documentation)
- MATLAB interface now returns number of function evaluations, too
- the MA57 interface can now be used with the MA57 library that comes with
  MATLAB (configure option `--enable-matlab-ma57`; cannot use Metis)
- `auxdata` is now handled by a wrapper function `ipopt_auxdata.m` instead
  of internally within the MEX code (replace Matlab call to `ipopt` with
  `ipopt_auxdata` if using auxiliary data in any callbacks) [r2282]
- exposed more intermediate Ipopt information to `iterfunc` callback [r2283]

- fixes to JIpopt buildsystem (now may work on windows and uses libtool)
- JIpopt now reads options file `ipopt.opt` by default, if present
- removed predefined `KEY_` strings in JIpopt
- renamed API functions that retrieve solution values in JIpopt

- simplified installation of R interface

## 3.10

### 3.10.4 (2013-05-05)

- fixed sign of dual values in AMPL solution again (with help of Gabe)
  [r2169, r2170, r2184, #183]
- fixed bugs with reoptimizing a TNLP with all variables fixed [r2172, r2173, #179]
- fixed more issues with sparse data structures and non-double numbers
  in Matlab interface (by T. Kelman) [r2191]
- added missing `final int` for Ipopt return code `Maximum_CpuTime_Exceeded`
  in Java interface [r2216]
- fixed bug when trying to warmstart Ipopt in Java interface [r2253]
- fixed wrong use of `SmartPtr`'s in Java interface [r2255, r2263]
- fixed bug in returning final solution in Java interface [r2258]
- included patch in ThirdParty/Mumps to work around bugs in Mumps
  matrix ordering routines AMF and QAMD (now give preference to AMD and METIS)

### 3.10.3 (2012-11-19)

- minor fixes in MA86 interface (by Jonathan Hogg) [r2069, r2086]
- fix in `IpTripletToCSRConverter` for CSR forms with
  extra entries (by Jonathan Hogg) [r2087]
- workaround problems with Mac OS X Lion's blas library
  (by Frederic Hetch) [r2102, #181]
- the C interface now catches also Ipopt exceptions thrown
  within the `OptimizeTNLP` call and returns `Ipopt::Unrecoverable_Exception`
  as status [r2094, #144]
- fixed segmentation fault in adaptive barrier parameter update rule
  when using the mehrotra option on unconstrained problems [r2114, #114]
- fixed segmentation fault in case no iterate is available in case of
  catastrophic failure in restoration phase [r2115]
- fixed default for `mumps_dep_tol` to work with current Mumps versions [r2116]
- fixed sign of dual values in AMPL solution [r2123, #183]
- fixed issue with sparse gradients in Matlab interface
  (by T. Kelman) [r2133, #187]
- sIPOPT (by H. Pirnay):
  - starting values in C++ version of parametric example now
    match AMPL version [r2098]
  - numerical values in parametric example now match publication [r2099]
  - added options `n_sens_steps` and `sens_boundcheck` as AMPL options [r2099]
  - any non-zero suffix value is now accepted for `sens_init_constr` [r2100]
  - fix typo in AMPL interface (by Weifeng Chen) [r2110]
  - fix bug when compiling without AMPL interface [r2113]
- build system:
  - updated instruction on using nowadays HSL sources (by T. Kelman)
  - fixed issue with libdir being `<prefix>/lib64`
- other minor fixes

### 3.10.2 (2012-02-12)

- updates to HSL interface (by Jonathan Hogg):
  - MC68 can now be loaded dynamically, too
  - MA86 exploits built-in scaling
  - MA86 can choose best ordering from AMD and Metis
  - fix for return code of MA86 for singular matrices
- corrected computation of Upsilon (norm of step SQUARED)
- updates to MSVS `v8-ifort` project files and addition of vc10
  project files (using vc8-generated `IpoptFSS.dll`) (by Marcel Roelofs)
- minor bugfixes, include updates in BuildTools

### 3.10.1 (2011-09-20)

- include updates in BuildTools, including new ThirdParty/Metis
   (fix for URL to download Metis 4.0.3 release)
- linear solver loader prints error code when failing to load
  library under Windows
- message on failure when reallocating memory in Mumps now includes
  size of memory that was tried to be allocated
- added missing include of `cstdio/stdio.h` in `IpJournalist.hpp`
- minor fixes to build system

### 3.10.0 (2011-06-20)

- move to new COIN-OR configuration and installation convention
- primal infeasibility output is now true infeasibility in original
  problem formulation

## 3.9

### 3.9.3 (2011-04-07)

- include updates in BuildTools, including new ThirdParty/Metis
   (required to work with current metis release)

### 3.9.2 (2010-12-22)

- converted from Common Public License to Eclipse Public License
- some bugfixes from BuildTools

### 3.9.1 (2010-11-26)

- improved Hessian update for restoration phase
- added intermediate callback feature to C and Fortran interface

### 3.9.0 (2010-11-05)

- switching to new BuildTools system
- added R interface (contributed by Jelmer Ypma)
- updates in AsNMPC (by Hans Pirnay)

## 3.8

### 3.8.3 (2010-06-29)

- restated `SolveStatistics::TotalCPUTime` method for backward
  compatibility

### 3.8.2 (2010-06-16)

- uses MUMPS version 4.9 and Lapack version 3.2.1
- added AsNMPC contribution made by Hans Pirnay
- enhanced MA57 options
- several bug fixes and minor additions

### 3.8.1 (2009-10-30)

- Bugfix in NLP function evaluation timing measurement.  The
  time for the objective function gradient had been forgotten.

### 3.8.0 (2009-10-30)

- Added MSVC solution with Intel Fortran compiler to generate DLLs
  (contributed by Marcel Roelofs).  To make this work, a lot of methods
  in exported headers have been made virtual
- changed default convergence tolerance in restoration phase (now same
  as regular tolerance)
- output is flushed after each iteration

## 3.7

### 3.7.1 (2009-10-06)

- bugfix for square problems
- correct timing information (obj gradient was forgotten)
- flush output buffer after each iteration
- first code for iterative WSMP version (experimental and undocumented)

### 3.7.0 (2009-07-16)

- a number of fixes (including those from 2009 COIN-OR Bug Squashing Party)
- changes in some exposed header files to provide access to internal
  data structures for specific applications

## 3.6

### 3.6.1 (2009-05-01)

- minor corrections in tutorial files

### 3.6.0 (2009-04-29)

- new Matlab interface
- added new option to limit cpu time: `max_cpu_time`
- added ThirdParty directory for Metis to be used with MUMPS or MA57
- updated CUTEr Makefile to make it work with CUTEr2
- added files for a tutorial (including coding exercise)

## 3.5

### 3.5.5 (2009-01-13)

- minor fixes regarding compilation
- undocumented version of inexact method

### 3.5.4 (2008-09-29)

- changed to MUMPS version 4.8.3 in externals (Mumps developers
  seem to have removed 4.8.1).

### 3.5.3 (2008-09-19)

- changed back to MUMPS version 4.8.1 since there seem to be issues
  on Windows

### 3.5.2 (2008-09-18)

- changed to latest version of MUMPS (4.8.2)
- some bugfixes (linear algebra objects, automatic problem scaling)
- made sure the subversion revision number is correct in all files
- allowed general additional data and cq in `IpData` and `IpCq`

### 3.5.1 (2008-08-26)

- changed to latest version of MUMPS (4.8.1)

### 3.5.0 (2008-08-25)

- added `ComputeRowAMax` and `ComputeColAMax` methods to Matrix base class
- changed externals for MUMPS to stable/1.1
- call `finalize_solution` in more failure cases
  (this means that AMPL writes .sol file in more situations)
- added `IpTNLPReducer` as simple way to exclude constraints from problem
- several fixes, also from COIN-OR Bug Squashing Party 2008

## 3.4

### 3.4.2 (2008-07-18)

- some bug fixes
- added wallclock time routine
- penalty function version does no longer crash if it
  wants to go to restoration phase (not that this really helps
  convergence though)

### 3.4.1 (2008-05-30)

- some bug fixes
- deleted `v9` MSVC files again (since `v8` works fine for `v9`)
- print Ipopt version in default print level
- added option that allows to change name of options file
  (`option_file_name`)

### 3.4.0 (2008-04-25)

- added support to dynamically load HSL or Pardiso:
  If Ipopt has been compiled without some HSL or Pardiso solver,
  it can now load those solvers from a shared library at runtime
  without recompilation.  This will make distribution of binaries
  easier.  Does not work on all platforms yet.
- several bugfixes
- ensured compilation of MSVS project files (`v8` and `v9`)
- new special return code for square problems
  (`Feasible_Point_Found` returned if dual inf not small)
- new initialization option for bound multipliers
  (see option `bound_mult_init_method`)
- added simple penalty function line search option
  (`line_search_method=penalty`) - not guaranteed to converge, see
  Ipopt implementation paper (in MathProg)
- some very basic method to approximate constraint Jacobian by
  finite differences (not efficient, but will hopefully be extended)

## 3.3

### 3.3.5 (2008-02-28)

- corrected links for Ipopt mailing list
- added missing `Makefile.in` for Matlab interface
- the `addlibs*` files are now installed in `share/doc/coin/Ipopt`
  instead of lib
- updates in Matlab interface
- bugfix for ticket #56

### 3.3.4 (2007-12-27)

- headers are now installed in `include/coin` (no longer in `include/ipopt`)
- default for `dual_inf_tol` is now 1 (instead of 1e-4)
- In matlab interface, here the text from Peter Carbonetto:
  There have been several significant changes made to the MATLAB interface
  since the last release. The most important two changes are:
  1. Following the "warm start" feature of IPOPT, you may pass in initial
     estimates for the Lagrange multipliers.
  2. Callback routines for computing the objective, gradient (etc.) are now
     specified using function handles instead of strings.
     (Thanks to Marcus Brubaker at the University of Toronto for the initial suggestion.)

### 3.3.3 (2007-09-25)

- minor changes, bug fixes

### 3.3.1 (2007-06-20)

Synchronized with all changes in trunk; probably more than to be
remembered.  In the following a few:
- support for Mumps linear solver (contributed by Damian Hocking)
- `--print-options` flag for ipopt ASL solver executable to see all
  Ipopt options (available through `ipopt.opt` file)
- added Matlab interface (contributed by Peter Carbonetto)
- added support for `f2c` compiler to compiler Fortran code with
  MSVC++ compiler
- new MSVS support (now within MSVS project and also with `f2c`)
- a number of small changes/bug fixes/improvements
- small change in interface (e.g., `FinalizeSolution` method)

## 3.2

### 3.2.4 (2007-04-24)

- updated download script for Blas to fit netlib's recent changes
- using a more recent version of BuildTools

### 3.2.3 (2006-11-29)

- updated download script for Lapack to fit to netlib's recent changes

### 3.2 r795 (2006-10-11)

- Bugfix in L-BFGS update
- fix in configure with detection of `sizeof(int*)` on Cygwin

### 3.2.1 (2006-07-14) - dev release number 764

- Bugfix in least square multiplier estimate.
  It mainly showed up in LBFGS with restoration phase as seg fault

### 3.2.0 (2006-07-07) - dev release number 757

- changed installation procedure and directory structure to
  conform with new COIN-OR convention

## 3.1

### 3.1.0 (2006-04-08) - dev release number 714

Several bug-fixes, improvements and additions. In particular:
- new quasi-Newton approximation using L-BFGS
- interfaces to linear solver MA57, WSMP, Pardiso
  (MUMPS and TAUCS not yet completed)
- derivative checker
- unit test
- configure supports compilation under Cygwin with native Windows compilers
- ScalableExample
- user call-back method in TNLP

## 3.0

### 3.0.1 (2005-12-04)

- Several corrections to Windows files
- Fix termination if number of iterations is exceeded in restoration phase

### 3.0.0 (2005-08-26) - dev release number 510

- First official release of the new C++ implementation of Ipopt.

## 2.2

### no new release (2005-08-19)

- corrected detection of BLAS libraries for SUN (make sure the example Makefiles work)
- upgrade LICENSE file to CPL version 1.0 as retrieved from www.opensource.org

### 2.2.1e (2005-05-30)

- fixed sign of multipliers returned to AMPL
  (bug reported by Rhoda Baker and Karsten Theissen)
- switched to automake 1.9.5

### no new release (2005-01-07)

- bugfix for the limited memory BFGS in case of square problems
  (bug reported by Wanhe Zhang and Ned Nedialkov)

### 2.2.1d (2004-10-05)

- Added `outlev` as an option to the AMPL solver as a synonym for `iprint`
- For `iprint` = 0, the output lines per iteration are now suppressed
- corrected two bugs in `configure` script (test for size of `long` etc
  before Fortran libraries as added to `LIBS`; prevent cycling in
  `make -j 1` test)
- internally renamed subroutine `ERROR` to `OPTERROR` (`ERROR` had a
  name clash for some Fortran compiler)
- avoid uninitialized variable in `update_b_lm.f`
- minor correction in computation of residual in `get_step_full.F`
- minor change for slack correction in `filter.F`

### 2.2.1c (2004-07-20)

- corrected bug leading to very small `QTAU` in rare circumstances

### 2.2.1b (2004-05-21)

- Make `DFILLINFACT` option available through AMPL interface
- Now, later increase of memory requirement for Harwell solvers is also
  based on `DFILLINFACT`, instead of using a fixed values of 10.

### 2.2.1a (2004-05-13)

- fix in `IPOPT/ipopt/mainloop.F`:
  The multipliers were not scaled back for low printlevel.

### no new release (2004-04-28)

- fix in `IPOPT/AMPL_interface/ipoptAMPL.c`:
  Now the mulitpliers for the constraints are passed back to AMPL.
- Added download scripts using `wget` to get ASL, BLAS, and LAPACK more easily.
  Thanks to Frank Wuebbeling for the hint.

### 2.2.1 (2004-04-25)

- AMPL solver executable is now called `ipopt` (instead of `ipoptAMPL`).
  This fixed also problem with assigning IPOPT options from within
  AMPL. (reported by Karsten Theissen)
- default value for number of iterations is now 10000 (instead of 1000)
- new option: `IMAXCPUSEC` to be set to the maximum number of CPU seconds
  after which the algorithm should stop.  The check is performed only at
  certain points in the algorithm, so that the executable might run longer
  than specified.
  The algorithm also stops (as before) when the file "STOP" is detected in
  the current directory.  Finally, a call to `USER_REQUESTED_STOP` has
  been added if the preprocessor macro `USE_USER_REQUESTED_STOP` has been
  defined. If this `LOGICAL` function returns `.TRUE.`, the algorithm also
  stops.  This feature was requested by David Ternet.
- if IPOPT is run several time in a row, the counting of function
  evaluations is restarted after every new call of IPOPT.
- in `get_step_full`, now check if the number of negative eigenvalues is LESS
  than then number constraints.  If so, increase the pivot tolerance,
  and if that doesn't help, pretend that the system is singular.
  (this fixed a problem reported by Hans Mittelmann)
- decrease default values for `DPIVTOLMAX`
- suppress superflous leading zeros in iteration output.  For long output
  (`ioutput=1`) include CPU time
- a few changes regarding the inertia correction (`get_step_full.F`).
  This decreases CPU time significantly in a few cases.
- corrected problem in C-files related to names for `struct`'s (some GNU
  compiler complained)
- changed default options for GNU compilers (now `-O3 -funroll-loops`,
  no longer `-O2` and `-mieee-fp`)
- no reference to MC35 in `resto_tron` if `HAVE_MC35` is not defined
- `#error` preprocessor directive removed from `*.F` files, since not all
  compilers understand it (reported by Hans Mittelmann)
- switch to automake 1.8.3

### 2.2.0 (2004-03-10)

Many things have changed since the last official release.
Here a few highlights:
- easier installation procedure with autoconf
- algorithm made more robust and efficient
- new restoration phase for filter method (TRON no longer needed for
  full-space option anymore)
- C-interface
