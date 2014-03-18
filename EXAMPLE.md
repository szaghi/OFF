# Running Examples

The directory _examples_ into the root of the project contains some simple examples to try OFF. These examples are designed to learn OFF quickly and easily. The examples are categorized into 2 classes, namely one-dimensional and two-dimensional. The followings instructions drive you to the examples execution.

## 1D Examples: SOD Riemann Problem

The first 1D example concerns the solution of the Riemann Problem of Sod, that is a classical 1D shock tube test. In order to complete the SOD example, IBM, OFF and POG must be properly compiled.

### Compiling Codes

IBM and POG can be compiled only with serial support, whereas OFF has support for serial, OpenMP, MPI and OpenMP/MPI execution. POG can be compiled with or without the Tecplot Inc. proprietary IO library for saving Tecplot binary outputs.

Let us start cleaning the project root directory:

    make cleanall

Then compile codes:

##### IBM
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes
##### POG with Tecplot Inc. proprietary Library
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes TECIO=yes
##### POG without Tecplot Inc. proprietary Library
    make clean
    make codes COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes
##### OFF serial
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes
##### OFF OpenMP
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=yes MPI=no NULj=yes NULk=yes
##### OFF MPI
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=no MPI=yes NULj=yes NULk=yes
##### OFF hybrid OpenMP/MPI
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=yes MPI=yes NULj=yes NULk=yes

### Running examples

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
