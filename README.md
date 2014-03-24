# OFF

OFF, Open source Finite volumes Fluid dynamics code [see documentation](http://szaghi.github.com/OFF/index.html).

It is written in in standard (compliant) Fortran 2003 with highly modularity as design target.

The aim of _OFF_ is to solve, numerically, the Navier-Stokes equations of fluid dynamics by means of Finite Volume technique.

The main features of _OFF_ code are the following:
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

## Copyrights

OFF is an open source project, it is distributed under the [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html). Anyone is interest to use, to develop or to contribute to OFF is welcome. Take a look at the [contributing guidelines](CONTRIBUTING.md) for starting to contribute to the project.

## Documentation

Detailed documentation can be found on the [GitHub pages](http://szaghi.github.com/OFF/index.html) of the project.
