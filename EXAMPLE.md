# Examples

The directory _examples_ into the root of the project contains some simple examples to try OFF. These examples are designed to learn OFF quickly and easily. The examples are categorized into 2 classes, namely one-dimensional and two-dimensional. The followings instructions drive you to the examples execution.

## 1D Examples

### SOD Riemann Problem

The first 1D example concerns the solution of the Riemann Problem of Sod, that is a classical 1D shock tube test. In order to complete the SOD example, IBM, OFF and POG must be properly compiled.

#### Compiling codes

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

#### Running example

To run the Sod's shock tube test navigate to its working directory:

    cd examples/1D/shock_tube/SOD/

There are 3 sub-directories:
+ ibm
+ off
+ pog

Into each of the above directories there is a bash script, __run_[IBM,OFF,POG].sh__ that help you to complete the example. Moreover, there are the necessary option files to run each code. The bash scripts are designed to be executed on a workstation where you can directly execute programs without dealing with a resource manager system that schedules the jobs, such as [Torque](http://www.adaptivecomputing.com/products/open-source/torque/). For the execution of this example on a system using a jobs scheduler see the subsection [How-to deal with Torque-like scheduler](#jobscheduler).

The 3 subdirectories _ibm_, _off_ and _pog_ indicate the 3 steps that must be completed in sequence.

Firstly, navigate to _ibm_ directory and building the initial and boundary conditions and the mesh file (generating the input file for OFF):

      cd ibm

Then you have 2 options:
+ clean _ibm_ directory: `run_IBM.sh -clean`
+ run IBM for generating OFF input files: `run_IBM.sh -out`

The above commands will produce the correct input for OFF. The produced files are placed into _output_ directory inside _ibm_ one.

Now, navigate to _off_ directory

      cd ../off

Again, you have 3 options:
+ clean _off_ directory: `run_OFF.sh -clean`
+ run OFF in serial or OpenMP mode: `run_OFF.sh -no_mpi`
+ run OFF in MPI mode: `run_OFF.sh -mpi`

The above commands will perform the simulation. The produced files are placed into _output_ directory inside the _off_ one.

Finally, postprocess the result obtained. Navigate to _pog_ directory:

      cd ../pog

You have 4 options:
+ clean _pog_ directory: `run_POG.sh -clean`
+ produce Gnuplot postprocessed files: `run_POG.sh -gnu`
+ produce Tecplot postprocessed files: `run_POG.sh -tec`
+ produce VTK postprocessed files: `run_POG.sh -vtk`

The above commands will post-process the simulation results. The produced files are placed into _output_ directory inside _pog_ one. It is worth noting that the obtained results can be compared with the OFF reference ones contained into _output-ref_ directory. Moreover, into the file __exact_sod-100-t02-block-split.dat__ the exact solution of Sod's problem are reported for validation purpose.

## 2D Examples

### Two Dimensional Riemann Problems

The first 2D examples concern with 2D Riemann problems. The initial conditions are the same as the ones reported by [Kurganov and Tadmor](http://dx.doi.org/10.1002/num.10025). In order to complete the SOD example, IBM, OFF and POG must be properly compiled.

#### Compiling codes

IBM and POG can be compiled only with serial support, whereas OFF has support for serial, OpenMP, MPI and OpenMP/MPI execution. POG can be compiled with or without the Tecplot Inc. proprietary IO library for saving Tecplot binary outputs.

Let us start cleaning the project root directory:

    make cleanall

Then compile codes:

##### IBM
    make IBM COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes
##### POG with Tecplot Inc. proprietary Library
    make POG COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes TECIO=yes
##### POG without Tecplot Inc. proprietary Library
    make clean
    make codes COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes
##### OFF serial
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=no MPI=no NULk=yes
##### OFF OpenMP
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=yes MPI=no NULk=yes
##### OFF MPI
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=no MPI=yes NULk=yes
##### OFF hybrid OpenMP/MPI
    make clean
    make OFF COMPILER=intel DEBUG=no OPENMP=yes MPI=yes NULk=yes

#### Running example

To run a 2D Riemann Problem example navigate to its working directory:

    cd examples/2D/two_dimensional_riemann_problems/kt-c??/

Note: substitute ?? with one of 03,04,05,06,12,17

For each kt-?? there are 3 sub-directories:
+ ibm
+ off
+ pog

Into each of the above directories there is a bash script, __run_[IBM,OFF,POG].sh__ that help you to complete the example. Moreover, there are the necessary option files to run each code. The bash scripts are designed to be executed on a workstation where you can directly execute programs without dealing with a resource manager system that schedules the jobs, such as [Torque](http://www.adaptivecomputing.com/products/open-source/torque/). For the execution of this example on a system using a jobs scheduler see the subsection [How-to deal with Torque-like scheduler](#jobscheduler).

The 3 subdirectories _ibm_, _off_ and _pog_ indicate the 3 steps that must be completed in sequence.

Firstly, navigate to _ibm_ directory and building the initial and boundary conditions and the mesh file (generating the input file for OFF):

      cd ibm

Then you have 2 options:
+ clean _ibm_ directory: `run_IBM.sh -clean`
+ run IBM for generating OFF input files: `run_IBM.sh -out`

The above commands will produce the correct input for OFF. The produced files are placed into _output_ directory inside _ibm_ one.

Now, navigate to _off_ directory

      cd ../off

Again, you have 3 options:
+ clean _off_ directory: `run_OFF.sh -clean`
+ run OFF in serial or OpenMP mode: `run_OFF.sh -no_mpi`
+ run OFF in MPI mode: `run_OFF.sh -mpi`

The above commands will perform the simulation. The produced files are placed into _output_ directory inside the _off_ one.

Finally, postprocess the result obtained. Navigate to _pog_ directory:

      cd ../pog

You have 4 options:
+ clean _pog_ directory: `run_POG.sh -clean`
+ produce Gnuplot postprocessed files: `run_POG.sh -gnu`
+ produce Tecplot postprocessed files: `run_POG.sh -tec`
+ produce VTK postprocessed files: `run_POG.sh -vtk`

The above commands will post-process the simulation results. The produced files are placed into _output_ directory inside _pog_ one. It is worth noting that the obtained results can be compared with the OFF reference ones contained into _output-ref_ directory.

It is worth noting that all the 2D Riemann problems are set to finish at 20 time steps in order to quickly complete the examples. However, the solutions obtained have no physical meaning because the final time of integration is too small. If you want to complete one of the 2D Riemann problems reaching the same physical time of the solution of [Kurganov and Tadmor](http://dx.doi.org/10.1002/num.10025) you must edit one input file into _ibm_ directory of the kt-?? example chosen:

      cd ibm
      cd input

Into this directory there is the option file `solver.dat`: the fourth line should contain the maximum number of iterations performed, _Nmax_, just set it to 0 in order to finish the simulation at the physical time set below, _Tmax_, instead of using _Nmax_.

### <a name="jobscheduler"></a> How-to run the examples within a resource manager scheduling jobs like Torque or PBS
