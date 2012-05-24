/**
@page SOD 1D Example: SOD shock tube

@section RiemannP Riemann Problem definition

A shock tube test consists in a one dimensional Riemann Problem with non reflective boundary conditions. A Riemann Problem could be defined as following. Assuming that the conservation laws \f$ \frac{\partial }{{\partial t}}\int\limits_{{V_{ijk}}} {\overrightarrow U d{V_{ijk}}}  = R\left( {\overrightarrow U } \right) \f$ in the space domain [0,1], are coupled with the following initial conditions:

\f$ U\left( {x,t = 0} \right) = {U_0}\left( x \right) = \left\{ \begin{array}{l} {U_L}\,x < x_0\\ {U_R}\,x > x_0 \end{array} \right. \f$ \n
where \f$U_L\f$ and \f$U_R\f$ are two constant values of the conservative variables and \f$x_0\f$ is the abscissa of the interface, the Riemann Problem can be defined as a partial derivative problem coupled with a piecewise constant initial conditions having a single discontinuity. Solving a Riemann Problem consists to find the time evolution, according to the conservative laws, of the discontinuity on the initial conditions. The structure of the solution of the Riemann Problem for the Euler conservation is composed by three waves. The middle wave (traveling along \f$v\f$) is always a contact discontinuity while the left and right waves (traveling along \f$v \pm a\f$, being \f$a=\sqrt{\frac{\gamma p}{\rho}}\f$ the speed of sound) are the non-linear acoustic waves and they can be either a shocks or rarefactions. Therefore there are four possible wave patterns.

@section SodST Sod's shock tube

The Sod's shock tube is defined (in terms of primitive variables, \ref Data_Type_Primitive::Type_Primitive "see definition of primitive derived type") as:

 - Left state: \f$P_L  = \left[ {\begin{array}{*{20}{c}} \rho \\ {v}\\ {p} \end{array}} \right]=\left[ {\begin{array}{*{20}{c}} 1 \\ 0\\ 1 \end{array}} \right]\f$
 - Right state: \f$P_R  = \left[ {\begin{array}{*{20}{c}} \rho \\ {v}\\ {p} \end{array}} \right]=\left[ {\begin{array}{*{20}{c}} 0.125 \\ 0\\ 0.1 \end{array}} \right]\f$

After the rupture of the initial discontinuity the solution consists in a shock moving to the right, a rarefaction fan moving to left and a contact discontinuity, moving to the right, separating the left and right states.

@section Compiling the codes

Before running this example the IBM, @off and POG codes must be properly compiled.

For more details see \ref Compiling "Compiling Instructions". In the following subsections few guidelines are provided.
@subsection CompIBM Compiling IBM

IBM, necessary for building up the initial and boundary conditions and the mesh, could be compiled without any particular options:
@code
make IBM DEBUG=no OPTIMIZE=yes
make clean
@endcode

@subsection CompOFF Compiling OFF
@off must be compiled with <em>NULj=yes NULk=yes</em> in order to perform a correct 1D simulation. It is possible to compile serial or parallel code:

For serial code:
@code
make OFF DEBUG=no OPTIMIZE=yes NULj=yes NULk=yes
make clean
@endcode
For multi-processes/multi-threads code:
@code
make OFF DEBUG=no OPTIMIZE=yes NULj=yes NULk=yes OPENMP=yes MPI=yes
make clean
@endcode

@subsection CompPOG Compiling POG
POG must be compiled with <em>NULj=yes NULk=yes</em> in order to perform a correct 1D post-process. It is possible to compile serial or parallel (only multi-threads) code:

For serial code:
@code
make POG DEBUG=no OPTIMIZE=yes NULj=yes NULk=yes
make clean
@endcode
For multi-threads code:
@code
make OFF DEBUG=no OPTIMIZE=yes NULj=yes NULk=yes OPENMP=yes
make clean
@endcode
@note If binary Tecplot output is wanted also the option <em>TECIO=yes</em> must be passed.

@section RunSod Running the Sod's shock tube example

To run the Sod's shock tube test enter the following directory:
@code
cd examples/1D/shock_tube/SOD/
@endcode
There are 3 sub-directories:
 - ibm
 - off
 - pog

Into each of the above directories there is a bash script, <em> run_[IBM,OFF,POG].sh </em> that help you to complete the simulation.  Moreover there are the necessary option files to run each code.

To complete the simulation you must follow 3 steps in the prescribed order:
 - Building initial and boundary conditions and mesh. Into \em ibm directory there are 6 ascii input files with extension <em> .dat</em>, one symbolic link to IBM code previously compiled, and the script run_IBM.sh. The latter script has one command line argument.  If it is invoked without arguments it prints a help message. The command line argument is necessary to distinguish serial simulation from multi-processes one:
   @code
   cd ibm
   @endcode
   - For serial simulation:
     @code
     run_IBM.sh -no_mpi
     @endcode
   - For multi-processes (MPI) simulation:
     @code
     run_IBM.sh -mpi
     @endcode
   The above commands will produce the correct input for @off. The produced files are placed into \em output directory inside ibm one. Exit from ibm directory: @code cd ../@endcode
 - Running @off. off directory contains 2 symbolic link, one to <em> ../ibm/output </em> and one to OFF code previously compiled.  There are also an ascii file \em off_options.dat containing the main options for @off and a bash script run_OFF.sh. The latter script is similar to previous run_IBM.sh. If it is invoked without arguments it prints a help message. The command line argument is necessary to distinguish serial simulation from multi-processes one:
   @code
   cd off
   @endcode
   - For serial simulation:
     @code
     run_OFF.sh -no_mpi
     @endcode
   - For multi-processes (MPI) simulation:
     @code
     run_OFF.sh -mpi
     @endcode
   The above commands will perform the simulation. The produced files are placed into \em output directory inside off one. Exit from off directory: @code cd ../@endcode
 - Post-processing the simulation, running POG. The pog directory contains 2 symbolic link, one to <em> ../off/output </em> and one to POG code previously compiled.  There is also a bash script run_POG.sh. If it is invoked without arguments it prints a help message. The command line argument is necessary to distinguish VTK post-processing to Tecplot one:
   @code
   cd pog
   @endcode
   - For Tecplot post-processing:
     @code
     run_POG.sh -tec
     @endcode
   - For VTK post-processing:
     @code
     run_POG.sh -vtk
     @endcode
   The above commands will post-process the simulation results. The produced files are placed into \em output directory inside pog one. Exit from pog directory: @code cd ../@endcode


*/
