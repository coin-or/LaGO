# LaGO

## **Note**: The development of LaGO has ceased.

## Introduction

LaGO (*La*grangian *G*lobal *O*ptimizer) is a software-package for the global optimization of nonconvex mixed-integer nonlinear programs (MINLP).
It is written in C++ and is released as open source code under the [Common Public Licence (CPL)](http://www.opensource.org/licenses/cpl.php).
The code has been written by Ivo Nowak and [Stefan Vigerske](http://www.math.hu-berlin.de/~stefan) (then Humboldt-University Berlin), who is the COIN-OR project leader for LaGO.

LaGO is designed to find global solutions of mathematical optimization programs of the form
```
    min   f(x)
    s.t.         g(x) <=  0
                 h(x)  =  0
          x_L <=  x   <=  x_U
                  x_i integer   i\in B
```
where ` f: R^n --> R ` is the objective function and ` g: R^n --> R^m ` and ` h: R^n --> R^k ` are the constraint functions.
` x_L ` and ` x_U ` denote lower and upper bounds on the variables, and ` B ` the indices of the discrete variables.
The functions ` f(x), g(x),` and ` h(x) ` can be nonlinear and nonconvex, but have to be twice continuously differentiable.

LaGO can work with black-box formulations of the functions ` f(x), g(x),` and ` h(x).`
Only methods for the evaluation of function values, gradients, and Hessian-vector-multiplications, and information about sparsity of the functions have to be provided.

The bounds `x_L` and `x_U` should be finite.
LaGO has some methods to deal with unbounded variables, but they are likely to fail on more difficult optimization problems.


## Download / Installation

You can obtain the LaGO code from [GitHub](https://github.com/coin-or/LaGO).
The LaGO distribution can be used to generate solver executables for the [AMPL](http://www.ampl.com) or [GAMS](http://www.gams.com) modeling environments.

LaGO has only been used on Linux 32- and 64-bit systems with Intel compatible processors so far.
The installation is likely to fail on other platforms since it has never been tested.


### Latest stable version

The latest stable release is version 0.3. You can obtain it via Git by
```
  git clone -b stable/0.3 https://github.com/coin-or/LaGO.git
```

The build system can be used as documented on the [BuildTools webpage](https://github.com/coin-or-tools/BuildTools).

LaGO uses the following COIN-OR packages:

 * [BuildTools](https://github.com/coin-or-tools/BuildTools)
 * [CoinUtils](https://github.com/coin-or/CoinUtils)
 * [Clp](https://github.com/coin-or/Clp)
 * [Cgl](https://github.com/coin-or/Cgl)
 * [Osi](https://github.com/coin-or/Osi)
 * [Ipopt](https://github.com/coin-or/Ipopt)
 * ThirdParty/HSL or ThirdParty/Mumps, ThirdParty/Blas, ThirdParty/Lapack (required by Ipopt)
 * ThirdParty/ASL and ThirdParty/GAMSIO (one or the other)


### Third Party packages

Further, LaGO depends on the following third party packages:

 * [GAMS](http://www.gams.com) I/O libraries (alternative to ASL)
 * NIST Template Numerical Toolkit: [TNT](http://math.nist.gov/tnt)
 * Serial Graph Partitioning and Fill-reducing Matrix Ordering: [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
 * Random Number Generation: [RANLIB.C](http://www.netlib.org/random)
 * Interval Arithmetic Library: [FILIB++](http://www.math.uni-wuppertal.de/wrswt/software/filib.html) (optional)
 * ILOG LP Solver: [CPLEX](http://www.ilog.com/products/cplex) (optional)

For the first five packages and the ASL, scripts are provided which download and install these packages.
These scripts have to be executed before `configure` is called.
Alternatively, you can give the location of the compiled packages as parameters in your `configure` call.


#### GAMS interface

If you want to use LaGO with GAMS interface, you can download the required GAMS I/O libraries using the `get.LX3` (32-bit Linux) or `get.LEI` (64-bit Linux) script in ThirdParty/GAMS_IO.

Furthermore, LaGO requires some libraries from an installed GAMS distribution for linking.
Hence, the `configure` script searches for a GAMS system in your `PATH` environment variable.
Alternatively, you can give the path to a GAMS system as a parameter to the `configure` script.
A GAMS system can be downloaded [here](http://www.gams.com/download).

To make GAMS known of LaGO, you copy the file `bin/gmsla_.zip` into your GAMS distribution.
To inform GAMS about the new solver, you now need to call the program `gamsinst` of your GAMS distribution.
From now on you can use LaGO by giving the argument `MINLP=LAGO` when solving a MINLP via GAMS.


#### AMPL interface

If the `configure` script does not find GAMS I/O libraries or a GAMS system, it will check whether the AMPL Solver Library is installed.
If so, a LaGO binary with AMPL interface will be generated.
You use this binary by first generating a `.nl` file from you AMPL model and then parsing the name of the .nl file as argument to LaGO.


#### FILIB++

If you make the FILIB++ available to LaGO, features based on interval arithmetic are activated.
Note, that if you use the FILIB++, also the interface to the optimization problem need to provide methods to evaluate functions over an interval.
This is only the case for the GAMS interface.
Thus, the interval arithmetic features in LaGO will not work when you use the AMPL interface (exception: If your MINLP is purely linear and quadratic (MIQQP), then it will work, since quadratic and linear functions are handled separately by LaGO.).


#### CPLEX

If you want to use CPLEX, you have to provide the linker flag and the path of the include file `cplex.h` to the configure script (`--with-cplex-lib` and `-with-cplex-incdir`), e.g.,
```
--with-cplex-incdir=<path-to-cplex.h> --with-cplex-lib="-L<path-to-cplex-library> -lcplex -lpthread"
```
If CPLEX is not available for LaGO, COIN/Clp is used.


## Usage / Documentation

LaGO can be used with either AMPL or GAMS, depending on your configuration.

If you use GAMS, you can install LaGO as MINLP solver in your GAMS system, see the documentation above.

If you use AMPL, you have to pass the name of a `.nl` file as argument to LaGO. Optional you can pass the name of parameter file as second argument.

Due to the early state of the LaGO project, documentation is available only in a very limited form.
The directory doc contains a makefile which calls flex to generate a program that generates a parameter documentation.
The easiest is to call make in the doc directory. Then point your webbrower to one of the files `doc/paramdoc.html`, `doc/shortparamdoc.html`, `doc/veryshortparamdoc.html`.
These files differ in the amount of LaGO parameters that are listed.
  * `doc/veryshortparamdoc.html` documents those parameters which might be most interesting for you when you are starting to use LaGO.
  * `doc/shortparamdoc.html` documents those parameters that can also be useful to influence the behaviour of LaGO.
  * `doc/paramdoc.html` documents all parameters of LaGO.
    This documentation is more directed to developers, since many features that you can switch on there are not tested and are in a very early state.

The parameters in a LaGO optionfile have the format `parametername : parametervalue`.
An exemplary `lago.opt` has the form
```
  MinlpBCP max iter: 1000
  BCP subdiv type: Violation
  #nlp solver for local minimization
  GAMS LocOpt solver  : snopt
```

A doxygen documentation of the classes in LaGO can be generated by calling `make doxydoc` in the main directory of the LaGO package (after configure has been called).

The paper [LaGO - a (heuristic) Branch and Cut algorithm for nonconvex MINLPs](http://www.mathematik.hu-berlin.de/publ/pre/2006/P-06-24.ps) and the slides
[LaGO - Branch and Cut for nonconvex MINLPs](http://www.math.hu-berlin.de/~eopt/papers/ago07.pdf) illustrate the Branch and Cut algorithm implemented in LaGO.

We have also a [mailing list](http://list.coin-or.org/mailman/listinfo/lago) where you can subscribe and post questions and comments regarding LaGO.
If you believe you found a bug in the code, then you can keep it.


## Tests

The LaGO distribution contains a directory `test`.
After installation, you find in this directory the script `run_tests`.
This script runs LaGO with several models that are taken from the [GAMS MINLPLib](http://www.gamsworld.org/minlp/minlplib.htm).
All these tests should work.
If you experience that some test is not working, you may want to try an actively maintained MINLP solver instead.


## Project Links

* [COIN-OR Initiative](http://www.coin-or.org)

* [mailing list](http://list.coin-or.org/mailman/listinfo/lago)
