# OFF

~~~
OFF, Open source Finite volumes Fluid dynamics code
~~~

_OFF_ is a CFD code designed to be accurate, efficient and modular for solving, numerically, the Navier-Stokes equations of fluid dynamics by means of Finite Volume technique. It is written in standard (compliant) Fortran 2003 by means of OOP paradigm. _OFF_ can be executed on parallel CPU-based architecture (shared memory multi-cores workstation, distributed memory cluster and hybrid distributed memory cluster based on shared memory nodes).

### Obtaining OFF

_OFF_ can be obtained only from github repository. You can use `Download ZIP` button provided from github or you can locally clone the github repository:

~~~
git clone https://github.com/szaghi/OFF
~~~

Two branches are generally present: `master` branch that holds the stable version and `testing` one for developing purposes. In general, the git `tags` are used for referencing versions updates.

### Copyrights

_OFF_ is an open source project, it is distributed under the [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html). Anyone is interest to use, to develop or to contribute to _OFF_ is welcome. Take a look at the [contributing guidelines](CONTRIBUTING.md) for starting to contribute to the project.

## Main features

### API

* Clean and lightweight programming style:
    - explicit declaration imposed (pervasive `implicit none` usage);
    - Fortran free form syntax;
    - clear and comprehensive comments;
* easy to maintain and to extend:
    - extensive use of OOP paradigm;
* high quality API documentation:
    - doxygen-based high quality html pages;
* simple development environment:
    - easy source files maintenance by means of Git, one of the best distributed versioning system;
    - collaborative framework by means of GitHub repository.

### Mathematical and numerical models

* Finite Volume, Godunov-like scheme based on Euler conservation Laws written in fully conservative formulation:
    - the extension to viscous Navier-Stokes equations is under developing;
* Underling Riemann Problem solver for convective fluxes:
    - Approximate Riemann solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm;
    - Approximate Riemann solver based on Primitive Variables Linearization algorithm;
    - Approximate Riemann solver based on Two Rarefactions algorithm;
    - Approximate Riemann solver based on Two Shocks algorithm;
    - Approximate Riemann solver based on Adaptive (non iterative) PVL-TR-TS algorithm;
    - Approximate Riemann solver based on Adaptive (non iterative) LF-TR algorithm;
    - Approximate Riemann solver based on HLLC algorithm;
    - Approximate Riemann solver based on Roe linearization.
    - Exact Riemann solver based on iterative solution of u-function;
* Multi-Species fluids models:
    - Partial Densities species conservation (Standard Thermodynamic Model);
    - New multi-dimensional conservation models of Favini, B. et al (under developing);
* Multi-Phases fluids models:
    - Fully-coupled Lagrangian particles transport model (under developing);
* Space numerical integration models:
    - 1-st order piece-wise constant reconstruction;
    - 2-nd order TVD linear-wise reconstruction;
    - 3-rd,5-th,7-th orders WENO non-linear reconstruction;
* Time approximation models:
    - 1-st order forward Euler integration;
    - 2-nd,3-rd,4-th orders Strong-Stability-Preserving explicit Runge-Kutta integration;
* Local pseudo-time convergence acceleration for steady simulations;
* Multi-grid time convergence acceleration:
    - Multi-grid model has been already developed, but it is affected by some not still recognized bugs. Testing and bugs fixing
      are in progress.
* Underling numerical grid models:
    - 3D, general curvilinear, body-fitted, structured multi-blocks mesh;
    - Adaptive Mesh Refinement, AMR model (under developing);
    - Blocks overlapping, overset (Chimera) model (to be developed in future);
* Computational parallelism ability:
    - Domain decomposition by means of Message Passing Interface (MPI) paradigm providing the ability to use distributed-memory
      cluster facilities;
    - Fine, local parallelism by means of OpenMP paradigm providing the ability to use shared-memory cluster facilities;
    - Fine, local parallelism by means of GPU programming (e.g. CUDA framework) providing the ability to use GPUs cluster
      facilities (to be developed in future).

## Documentation

Detailed documentation can be found on the [GitHub pages](http://szaghi.github.com/OFF/index.html) of the project.
