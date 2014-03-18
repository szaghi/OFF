# Running Examples

The directory _examples_ into the root of the project contains some simple examples to try OFF. These examples are designed to learn OFF quickly and easily. The examples are categorized into 2 classes, namely one-dimensional and two-dimensional. The followings instructions drive you to the examples execution.

## 1D Examples: SOD Riemann Problem

The first 1D example concerns the solution of the Riemann Problem of Sod, that is a classical 1D shock tube test. A shock tube test consists in a one dimensional Riemann Problem with non reflective boundary conditions. A Riemann Problem could be defined as following. Assuming that the conservation laws

\f$ \frac{\partial }{{\partial t}}\int\limits_{{V_{ijk}}} {\overrightarrow U d{V_{ijk}}}  = R\left( {\overrightarrow U } \right) \f$

in the space domain [0,1], are coupled with the following initial conditions:

\f$ U\left( {x,t = 0} \right) = {U_0}\left( x \right) = \left\{ \begin{array}{l} {U_L}\,x < x_0\\ {U_R}\,x > x_0 \end{array} \right. \f$ \n

where \f$U_L\f$ and \f$U_R\f$ are two constant values of the conservative variables and \f$x_0\f$ is the abscissa of the interface, the Riemann Problem can be defined as a partial derivative problem coupled with a piecewise constant initial conditions having a single discontinuity. Solving a Riemann Problem consists to find the time evolution, according to the conservative laws, of the discontinuity on the initial conditions. The structure of the solution of the Riemann Problem for the Euler conservation is composed by three waves. The middle wave (traveling along \f$v\f$) is always a contact discontinuity while the left and right waves (traveling along \f$v \pm a\f$, being \f$a=\sqrt{\frac{\gamma p}{\rho}}\f$ the speed of sound) are the non-linear acoustic waves and they can be either a shocks or rarefactions. Therefore there are four possible wave patterns.

### Compiling Codes

#### Serial
    make cleanall
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes
    make OFF COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes
##### With Tecplot Inc. proprietary Library
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
    make clean
    make codes COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes

#### OpenMP
    make cleanall
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=yes MPI=no NULj=yes NULk=yes
##### With Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes

#### MPI
    make cleanall
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=no MPI=yes NULj=yes NULk=yes
##### With Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes

#### Hybrid OpenMP/MPI
    make cleanall
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=yes MPI=yes NULj=yes NULk=yes
##### With Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes

### Running Examples

#### SOD Riemann Problem
    cd examples/1D/shock_tube/SOD/

##### Building Initial and Boundary Conditions Files
    cd ibm/
    ./run_IBM.sh -out

##### Performing the Simulation
##### Serial/OpenMP
    cd ../off/
    ./run_OFF.sh -no_mpi
##### MPI
    cd ../off/
    ./run_OFF.sh -mpi

##### Postprocessing the results
    cd ../pog/
    ./run_POG.sh -tec # for Tecplot Inc. output: POG must be compiled with TECIO=yes option
    ./run_POG.sh -vtk # for VTK output
    ./run_POG.sh -gnu # for Gnuplot output

Note: the obtained results can be compared with the OFF reference ones contained into output-ref directory.

## 2D Examples: Two Dimensional Riemann Problems

### Compiling Codes

#### Serial
##### With Tecplot Inc. proprietary Library
    make cleanall
    make codes COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
    make cleanall
    make codes COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes

#### OpenMP
    make cleanall
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes
    make clean
    make OFF COMPILER=intel OPTIMIZE=yes DEBUG=no OPENMP=yes MPI=no NULk=yes
##### With Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes

#### MPI
    make cleanall
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes
    make clean
    make OFF COMPILER=intel OPTIMIZE=yes DEBUG=no OPENMP=no MPI=yes NULk=yes
##### With Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes

#### Hybrid OpenMP/MPI
    make cleanall
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes
    make clean
    make OFF COMPILER=intel OPTIMIZE=yes DEBUG=no OPENMP=yes MPI=yes NULk=yes
##### With Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
    make clean
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes

### Running Examples
Note: substitute ?? with one of 03,04,05,06,12,17

#### kt-??
    cd examples/2D/two_dimensional_riemann_problems/kt-c??/

##### Building Initial and Boundary Conditions Files
    cd ibm/
    ./run_IBM.sh -out

##### Performing the Simulation
##### Serial/OpenMP
    cd ../off/
    ./run_OFF.sh -no_mpi
##### MPI
    cd ../off/
    ./run_OFF.sh -mpi

##### Post-processing the results
    cd ../pog/
    ./run_POG.sh -tec # for Tecplot Inc. output: POG must be compiled with TECIO=yes option
    ./run_POG.sh -vtk # for VTK output
    ./run_POG.sh -gnu # for Gnuplot output
