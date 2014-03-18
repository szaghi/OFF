# Running Examples

The directory _examples_ into the root of the project contains some simple examples to try OFF. These examples are designed to learn OFF quickly and easily. The examples are categorized into 2 classes, namely one-dimensional and two-dimensional. The followings instructions drive you to the examples execution.

## 1D Examples: SOD Riemann Problem

The first 1D example concerns the solution of the Riemann Problem of Sod, that is a classical 1D shock tube test.

A shock tube test consists in a one dimensional Riemann Problem with non reflective boundary conditions. A Riemann Problem could be defined as following. Assuming that the conservation laws

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7B%5Cpartial%20%7D%7B%7B%5Cpartial%20t%7D%7D%5Cint%5Climits_%7B%7BV_%7Bijk%7D%7D%7D%20%7B%5Coverrightarrow%20U%20d%7BV_%7Bijk%7D%7D%7D%20%3D%20R%5Cleft%28%20%7B%5Coverrightarrow%20U%20%7D%20%5Cright%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

in the space domain [0,1], are coupled with the following initial conditions:

![equation](http://www.sciweavers.org/tex2img.php?eq=U%5Cleft%28%20%7Bx%2Ct%20%3D%200%7D%20%5Cright%29%20%3D%20%7BU_0%7D%5Cleft%28%20x%20%5Cright%29%20%3D%20%5Cleft%5C%7B%20%5Cbegin%7Barray%7D%7Bl%7D%20%7BU_L%7D%5C%2Cx%20%3C%20x_0%5C%5C%20%7BU_R%7D%5C%2Cx%20%3E%20x_0%20%5Cend%7Barray%7D%20%5Cright.&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

where ![equation](http://www.sciweavers.org/tex2img.php?eq=U_L&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0) and ![equation](http://www.sciweavers.org/tex2img.php?eq=U_R&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0) are two constant values of the conservative variables and ![equation](http://www.sciweavers.org/tex2img.php?eq=x_0&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0) is the abscissa of the interface, the Riemann Problem can be defined as a partial derivative problem coupled with a piecewise constant initial conditions having a single discontinuity. Solving a Riemann Problem consists to find the time evolution, according to the conservative laws, of the discontinuity on the initial conditions. The structure of the solution of the Riemann Problem for the Euler conservation is composed by three waves. The middle wave (traveling along _v_) is always a contact discontinuity while the left and right waves (traveling along _v +- a_, being _a_ the speed of sound) are the non-linear acoustic waves and they can be either a shocks or rarefactions. Therefore there are four possible wave patterns.

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
