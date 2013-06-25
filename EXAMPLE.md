# Running Examples

In order to run the provided examples to the following steps

## 1D Examples

### Compiling Codes

#### Serial
##### With Tecplot Inc. proprietary Library
make cleanall
make codes COMPILER=intel DEBUG=no OPENMP=no MPI=no NULj=yes NULk=yes TECIO=yes
##### Without Tecplot Inc. proprietary Library
make cleanall
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

##### Peforming the Simulation
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

## 2D Examples
