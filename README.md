# IPOPT 2.3.x

Ipopt is a a software package for large-scale nonlinear optimization.

**NOTE: This is the pre-3.0 (Fortran) version of Ipopt.
This version of Ipopt is no longer maintained and cannot be build
anymore out-of-the-box. The project page for the new (C++) Ipopt 3.x version
is <https://github.com/coin-or/Ipopt>.**

## What is IPOPT?

IPOPT is an open source software package for large-scale [nonlinear
optimization](http://en.wikipedia.org/wiki/Nonlinear_programming) (NLP).
It can be used to solve optimization problems of the form

```
   min     f(x)
  x ∈ Rⁿ

s.t.       c(x) = 0
           x_L ≤  x   ≤ x_U
```

where `x` are the optimization variables (possibly with lower and upper
bounds `x_L` and `x_U`), `f(x)` is the objective function, and `c(x)`
are constraint functions. The functions `f(x)` and `c(x)`
can be nonlinear. (Note that nonlinear inequality constraints can be
formulated in the above statement using slack variables). IPOPT aims to
find a local solution of such a problem.

The IPOPT distribution can be used to generate a library that can be
linked to one's own FORTRAN or C code, as well as a solver executable
for the [AMPL](http://www.ampl.com/) modeling language. It also includes
an interface to CUTEr. IPOPT can be used on Linux/UNIX platforms and
Windows.

IPOPT is also available as NLP solver at the
[NEOS Server](http://www.neos-server.org) at Argonne National Laboratories.
There you can submit and solve your
[AMPL](http://www.neos-server.org/neos/solvers/nco:Ipopt/AMPL.html) and
[GAMS](http://www.neos-server.org/neos/solvers/nco:Ipopt/GAMS.html)
models online.

And to say it once technically: IPOPT implements an interior point line
search filter method.

## How can I get it?

IPOPT is written in Fortran 77 (and a little bit of C) and is released
as open source code under the
[Common Public License (CPL)](http://opensource.org/licenses/cpl1.0.php).
It is available from the [COIN-OR](http://www.coin-or.org) repository.
You can obtain the source code from the branch stable/3.2 of the
[Ipopt repository at GitHub](https://github.com/coin-or/Ipopt).

As [open source](http://www.opensource.org/) software, the source code
for IPOPT is provided without charge. You are free to use it, also for
commercial purposes. You are also free to modify the source code (with
the restriction that you need to make your changes public if you decide
to distribute your version in any way, e.g. as an executable); for
details see the LICENSE file. And we are certainly very keen on
[feedback](http://list.coin-or.org/mailman/listinfo/ipopt) from users,
including contributions\!

In order to compile IPOPT, certain third party code is required (such as
some linear algebra routines, or the AMPL Solver Library). Those are
available under different conditions/licenses.

We provide this program in the hope that it may be useful to others, and
we would very much like to hear about your experience with it. If you
found it helpful and are using it within our software, we would like to
include you in our user list below.

Since a lot of time and effort has gone into IPOPT's development,
<span style="font-weight: bold;">please cite the following publication
if you are using IPOPT for your own research</span>:

A. Wächter and L. T. Biegler,
**On the Implementation of a Primal-Dual Interior Point Filter Line
Search Algorithm for Large-Scale Nonlinear Programming**,
[Mathematical Programming 106(1), pp. 25-57, 2006](http://dx.doi.org/10.1007/s10107-004-0559-y)
([preprint](http://www.optimization-online.org/DB_HTML/2004/03/836.html))

## How can I learn more about it?

  - Some more general information is available in the README file.
    Installation instructions are provided in the INSTALL file (or its
    short version QUICKINSTALL). A detailed documentation of the
    interface and algorithmic options is given in the README.IPOPT file.
    Also, some remarks on using IPOPT on Windows can be found in the
    README.Windows file. Finally, a table with the most important
    options is given below.
  - If you are interested in the mathematical details of the algorithm,
    you may want to have a look at the above mentioned paper. Additional
    aspects of the method are described in
    [this PhD thesis](http://researcher.watson.ibm.com/researcher/files/us-andreasw/thesis.pdf).
  - There also exists a [mailing list](http://list.coin-or.org/mailman/listinfo/ipopt)
    for the IPOPT  project. Here we post announcements, and users can submit their
    comments, bug reports, questions etc.
    (<span style="text-decoration: underline;">Note</span>: You need to
    subscribe to the mailing list before you can submit a message - this
    way we avoid that subscribers are overloaded by spam messages...)


## Who is involved in IPOPT?

The main author and project leader is
[Andreas Wächter](http://users.iems.northwestern.edu/~andreasw).
The other original authors are
[Lorenz T. Biegler](http://www.cheme.cmu.edu/who/faculty/biegler.html),
Yi-Dong Lang, and Arvind Raghunathan.
A new C++ re-implementation is being developed by
[Andreas Wächter](http://users.iems.northwestern.edu/~andreasw) and
[Carl Laird](http://allthingsoptimal.com/biography/).

Further contributors are:

  - Dominique Orban (update of CUTEr interface)
  - Krik Abbott (C-interface)

We are also very greatful for helpful feedback and comments from: Juan
Arrieta, Hande Benson, Gerd Bürger, Andrew Conn, Antonio Flores,
Sebastien Gros, Tobias Jockenhövel, Dieter Kraft, Carl Laird, Petros
Mamales, Hans Mittelmann, Ned Nedialkov, Jorge Nocedal, Jorge Paloschi,
Richard Waltz.

## Who is using IPOPT?

Yes, this is what we would like to know, too\! :) So, if you are using
it, please let us know\!

Currently, IPOPT is/has been used in the following projects/products:

  - other [COIN-OR](http://www.coin-or.org) projects:
    [NLPAPI](https://projects.coin-or.org/NLPAPI/) and
    [DFO](https://projects.coin-or.org/Dfo)
  - several research projects in
    [Larry Biegler's group at CMU](http://www.cheme.cmu.edu/people/faculty/index.htm),
    including one with Ipopt as optimization engine in the process simulator
    [ROMeo<sup>TM</sup>](http://www.simsci-esscor.com/us/eng/simsciProducts/productlist/romeo/default.htm)
  - in the IBM-internal circuit tuning tool EinsTuner
  - in [ABB](http://www.abb.com)'s Pulp and Paper Online Production
    Optimizer
  - Tobias Jockenhövel's [OptControlCentre](http://www.optcontrolcentre.com/)
  - by [Ned Nedialkov](http://www.cas.mcmaster.ca/~nedialk) in the
    context of solving Differential-Algebraic Equations
  - at the Munich University of Applied Sciences (for the course on
    "Modeling, Simulation and Optimization" by
    [Dieter Kraft](http://dkraft.userweb.mwn.de/))
  - by Juan Arrieta, in the computation of optimal trajectories for
    spacecraft and related astronautics problems
  - *who/what did we forget? Let us know\!*

## Contact

If you have any questions or comments please send a message to the
[mailing list](http://list.coin-or.org/mailman/listinfo/ipopt). (Note:
You need to subscribe to the mailing list before you can post a
message.) However, note that the Fortran version is no longer maintained
and compilation or execution problems are not fixed.

## AMPL Options for IPOPT

You can set options for IPOPT from AMPL using the
```
option ipopt_options 'option1=value1 option2=value2 ...';
```
command in your AMPL model. For example, if you want to set the maximal
number of iterations to 1000, and want the problem to be solved to an
accuracy of 1e-6, you would add
```
option ipopt_options 'imaxiter=1000 dtol=1e-6';
```

(Note: Even if you rename your AMPL solver executable, you still need to
use the `ipopt_options` keyword.)

The following table describes some of the available options. If you are
looking for more options, please consult the [README.IPOPT](doc/README.IPOPT)
file.


| Parameter | Default Value | Description |
| :-------- | ------------: | :---------- |
| imaxiter  | 10000         | Maximum number of iterations |
| imaxcpusec| 999999999     | Upper limit on computation time (CPU seconds) |
| dtol      | 1e-8          | Overall desired (scaled) error tolerance |
| dcmaxtol  | 1e300         | Absolute tolerance for (unscaled) primal infeasibility |
| dinfmaxtol| 1e300         | Absolute tolerance for (unscaled) dual infeasibility |
| iprint    | 0             | Amount of output. The larger, the more detailed output (-1: minimal output) |
| ioutput   | 0             | If set to 1, a longer information line is printed per iteration |
| iscale    | 2             | Options for automatic scaling of the problem statement (0: no scaling) |
| dfscale   | 1.0           | Scaling factor for the objective function |
| dmu0      | 0.1           | Initial value of the barrier parameter |
| dbndfrac  | 0.01          | Relative distance of starting point from closest bound |
| dbndpush  | 0.01          | Absolute distance of starting point from closest bound |
| dmovebounds | 1e-8        | Initial relative perturbation of the bounds |
| dpivtol   | 1e-8          | Pivot tolerance for linear solver |
| dfillinfact | 5.0         | Factor for estimating memory requirement for linear solver (has to be >=1). If IPOPT runs out of memory, a smaller value (such as 1, 1.2, or 2) might correct the problem. |

**<span class="underline">COPYRIGHT:</span>**

Copyright (C) 2002, 2004, 2005 Carnegie Mellon University, IBM, and
others.
All Rights Reserved.
This code is published under the Common Public License.
