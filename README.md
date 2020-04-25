Ipopt
=====

Introduction
------------

Ipopt (Interior Point OPTimizer, pronounced eye-pea-Opt) is a software package for large-scale [nonlinear optimization](http://wiki.mcs.anl.gov/NEOS/index.php/Nonlinear_Programming_FAQ).
It is designed to find (local) solutions of mathematical optimization problems of the form

```
   min     f(x)
  x ∈ Rⁿ

s.t.       g_L ≤ g(x) ≤ g_U
           x_L ≤  x   ≤ x_U
```
where ```f(x): Rⁿ --> R``` is the objective function, and ```g(x): Rⁿ --> Rᵐ```
are the constraint functions.  The vectors `g_L` and `g_U` denote the lower and upper bounds on the constraints, and the vectors `x_L` and `x_U` are the bounds on the variables `x`.
The functions `f(x)` and `g(x)` can be nonlinear and nonconvex, but should be twice continuously differentiable.
Note that equality constraints can be formulated in the above formulation by setting the corresponding components of `g_L` and `g_U` to the same value.

Ipopt is part of the [COIN-OR Initiative](http://www.coin-or.org).
The Ipopt project webpage is <https://github.com/coin-or/Ipopt>.


Background
----------

Ipopt is written in C++ and is released as open source code under the [Eclipse Public License (EPL)](LICENSE).
The code has been written by [Andreas Wächter](http://www.mccormick.northwestern.edu/directory/profiles/Andreas-Waechter.html) and [Carl Laird](http://allthingsoptimal.com/biography/).
The COIN-OR project managers for Ipopt are [Andreas Wächter](http://users.iems.northwestern.edu/~andreasw) und [Stefan Vigerske](https://www.gams.com/~stefan).
For a list of **all contributors**, see the [AUTHORS file](AUTHORS).

The C++ version has first been [released on Aug 26, 2005](http://list.coin-or.org/pipermail/ipopt/2005-August/000331.html) as version 3.0.0.
The previously released [pre-3.0 Fortran version](https://github.com/coin-or/Ipopt/tree/stable/2.3) is no longer maintained.


The Ipopt distribution can be used to generate a library that can be linked to one's own C++, C, Fortran, or Java code, as well as a solver executable for the [AMPL](http://www.ampl.com) modeling environment.
The package includes an interface to the [R](http://www.r-project.org/) programming environment.
IPOPT can be used on Linux/UNIX, Mac OS X, and Windows platforms.

As open source software, the source code for Ipopt is provided without charge.
You are free to use it, also for commercial purposes.
You are also free to modify the source code (with the restriction that you need to make your changes public if you decide to distribute your version in any way, e.g. as an executable); for details see the EPL license.
And we are certainly very keen on feedback from users, including contributions!

In order to compile Ipopt, certain third party code is required (such as some linear algebra routines).
Those are available under different conditions/licenses.

If you want to learn more about Ipopt, you can find references in the [bibliography of the documentation](https://coin-or.github.io/Ipopt/citelist.html) and this ["Papers about Ipopt" page](https://github.com/coin-or/Ipopt/wiki/IpoptPapers).

For information on projects that use Ipopt, refer to the [Success Stories page](https://github.com/coin-or/Ipopt/wiki/SuccessStories).


Getting Started
---------------

On sufficiently prepared systems, a quick way to build and install Ipopt
is to get the coinbrew script from https://coin-or.github.io/coinbrew/
and running

    /path/to/coinbrew fetch Ipopt --no-prompt
    /path/to/coinbrew build Ipopt --prefix=/dir/to/install --test --no-prompt --verbosity=3
    /path/to/coinbrew install Ipopt --no-prompt

The coinbrew script will take care of building and installing the
dependencies ASL and Mumps before building Ipopt.

More details on using coinbrew can be found at the instructions on
[Getting Started with the COIN-OR Optimization Suite](https://coin-or.github.io/user_introduction).

If using coinbrew is not sufficient, then the
[installation instructions in the Ipopt documentation](https://coin-or.github.io/Ipopt/INSTALL.html)
should be studied.

Some precompiled binaries of Ipopt are also available:

- **[JuliaOpt provides Ipopt binaries](https://github.com/JuliaOpt/IpoptBuilder/releases)**
- **[AMPL provides binaries](http://ampl.com/products/solvers/open-source/#ipopt)** for using Ipopt through AMPL
- **[Pardiso project provides binaries](https://pardiso-project.org/index.html#binaries)** for using Ipopt with Pardiso through Matlab


Getting Help
------------

 * **[Ipopt Documentation](https://coin-or.github.io/Ipopt/)** with installation instructions, options reference, and more
 * **[Issue tracking system](https://github.com/coin-or/Ipopt/issues/)**: If you believe you found a **bug** in the code, please use the issue tracking system.
   Please include as much information as possible, and if possible some (ideally simple) example code so that we can reproduce the error.
 * **[Mailing list](http://list.coin-or.org/mailman/listinfo/ipopt)**: subscribe to get notifications about updates and to post questions and comments regarding Ipopt
 * **[Mailing list archive](http://list.coin-or.org/pipermail/ipopt/)**
 * **[Ipopt Wiki](https://github.com/coin-or/Ipopt/wiki)** with hints and tricks
 * [short Ipopt tutorial](http://drops.dagstuhl.de/volltexte/2009/2089/pdf/09061.WaechterAndreas.Paper.2089.pdf)

Please Cite Us
--------------

We provide this program in the hope that it may be useful to others, and we would very much like to hear about your experience with it.
If you found it helpful and are using it within our software, we encourage you to add your feedback to the [Success Stories page](https://github.com/coin-or/Ipopt/wiki/SuccessStories).

Since a lot of time and effort has gone into Ipopt's development, **please cite the following publication if you are using Ipopt for your own research**:

* A. Wächter and L. T. Biegler, **[On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm for Large-Scale Nonlinear Programming](http://dx.doi.org/10.1007/s10107-004-0559-y)**, _Mathematical Programming_ 106(1), pp. 25-57, 2006
  ([preprint](http://www.optimization-online.org/DB_HTML/2004/03/836.html))
