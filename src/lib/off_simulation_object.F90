#include "preprocessor_macros.h"
!< OFF simulation object definition and implementation.

module off_simulation_object
!< OFF simulation object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_error_object, only : error_object
use off_file_grid_object, only : file_grid_object
use off_file_solution_object, only : file_solution_object
use off_free_conditions_object, only : free_conditions_object
use off_mesh_object, only : mesh_object
use off_non_dimensional_numbers_object, only : non_dimensional_numbers_object
use off_os_object, only : os_object
use off_solver_object, only : solver_object
use off_time_object, only : time_object
use finer, only : file_ini
use flap, only : command_line_interface
use flow, only : eos_compressible
use foodie, only : integrand_object, integrator_multistage_object
use penf, only : I4P, I8P, MaxR8P, R8P, str

implicit none
private
public :: simulation_object

type, extends(integrand_object) :: simulation_object
   !< Simulation object class.
   !<
   !< [[simulation_object]] is a container for all simulation data for each processor/image.
   type(error_object)                   :: error                     !< Errors handler.
   type(command_line_interface)         :: cli                       !< Command line interface.
   type(os_object)                      :: os                        !< Running Operating System.
   type(file_ini)                       :: file_parameters           !< Simulation parameters file handler.
   type(file_grid_object)               :: file_grid_input           !< Grid file handler, input.
   type(file_grid_object)               :: file_grid_output          !< Grid file handler, output.
   type(file_solution_object)           :: file_initial_conditions   !< Intial conditions file handler.
   type(file_solution_object)           :: file_solution             !< Solution file handler.
   type(non_dimensional_numbers_object) :: adimensionals             !< Non dimensional numbers.
   type(free_conditions_object)         :: free_conditions           !< Free stream conditions.
   type(eos_compressible)               :: eos                       !< Equation of state.
   type(solver_object)                  :: solver                    !< solver models parameters.
   type(time_object)                    :: time                      !< Timing conditions.
   type(mesh_object)                    :: mesh                      !< Mesh data.
   logical                              :: is_cli_parsed=.false.     !< Sentinel of CLI parsing.
   logical                              :: is_output_verbose=.false. !< Verbose output.
   logical                              :: go_on_fail=.false.        !< Allow/disallow parameters loading failure.
   integer(I8P)                         :: shell_refresh_frequency=1 !< Output refresh frequency (on time steps).
   contains
      ! public methods
      procedure, pass(self) :: description                  !< Return a pretty-formatted description of the simulation.
      procedure, pass(self) :: destroy                      !< Destroy simulation data.
      procedure, pass(self) :: initialize                   !< Initialize simulation.
      procedure, pass(self) :: integrate                    !< Integrate the equations.
      procedure, pass(self) :: load_file_parameters         !< Load file parameters.
      procedure, pass(self) :: parse_command_line_interface !< Parse command line interface.
      procedure, pass(self) :: save_file_grid               !< Save file grid.
      procedure, pass(self) :: save_file_parameters         !< Save file parameters.
      procedure, pass(self) :: save_file_solution           !< Save file solution.

      ! integrand_object deferred methods
      procedure, pass(self) :: integrand_dimension !< Return integrand dimension.
      procedure, pass(self) :: t => residuals      !< Time derivative, residuals.
      ! operators
      procedure, pass(lhs) :: local_error !<`||integrand - integrand||` operator.
      ! +
      procedure, pass(lhs) :: integrand_add_integrand !< `+` operator.
      procedure, pass(lhs) :: integrand_add_real      !< `+ real` operator.
      procedure, pass(rhs) :: real_add_integrand      !< `real +` operator.
      ! *
      procedure, pass(lhs) :: integrand_multiply_integrand   !< `*` operator.
      procedure, pass(lhs) :: integrand_multiply_real        !< `* real` operator.
      procedure, pass(rhs) :: real_multiply_integrand        !< `real *` operator.
      procedure, pass(lhs) :: integrand_multiply_real_scalar !< `* real_scalar` operator.
      procedure, pass(rhs) :: real_scalar_multiply_integrand !< `real_scalar *` operator.
      ! -
      procedure, pass(lhs) :: integrand_sub_integrand !< `-` operator.
      procedure, pass(lhs) :: integrand_sub_real      !< `- real` operator.
      procedure, pass(rhs) :: real_sub_integrand      !< `real -` operator.
      ! =
      procedure, pass(lhs) :: assign_integrand !< `=` operator.
      procedure, pass(lhs) :: assign_real      !< `= real` operator.
      ! override fast operators
      procedure, pass(self) :: t_fast                              !< Time derivative, residuals, fast mode.
      procedure, pass(opr)  :: integrand_add_integrand_fast        !< `+` fast operator.
      procedure, pass(opr)  :: integrand_multiply_integrand_fast   !< `*` fast operator.
      procedure, pass(opr)  :: integrand_multiply_real_scalar_fast !< `* real_scalar` fast operator.
      procedure, pass(opr)  :: integrand_subtract_integrand_fast   !< `-` fast operator.

      ! private methods
      procedure, pass(self), private :: cli_initialize !< Initialize Command Line Interface.
      procedure, pass(self), private :: compute_dt     !< Compute the current time step by means of CFL condition.
endtype simulation_object

contains
   ! public methods
   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the simulation.
   class(simulation_object), intent(in)           :: self             !< Simulation parameters.
   character(*),             intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                  :: desc             !< Description.
   character(len=:), allocatable                  :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter                    :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//prefix_//'Grid file (input)'//NL
   desc = desc//prefix_//self%file_grid_input%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Grid file (output)'//NL
   desc = desc//prefix_//self%file_grid_output%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Initial conditions file'//NL
   desc = desc//prefix_//self%file_initial_conditions%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Solution file'//NL
   desc = desc//prefix_//self%file_solution%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Non dimensional numbers:'//NL
   desc = desc//prefix_//self%adimensionals%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Free stream conditions:'//NL
   desc = desc//prefix_//self%free_conditions%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Equation of state:'//NL
   desc = desc//prefix_//self%eos%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Solver models:'//NL
   desc = desc//prefix_//self%solver%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Time:'//NL
   desc = desc//prefix_//self%time%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Mesh:'//NL
   desc = desc//prefix_//self%mesh%description(prefix=prefix_//'  ')
   endfunction description

   _ELEMENTAL_ subroutine destroy(self)
   !< Destroy simulation data.
   class(simulation_object), intent(inout) :: self !< Simulation data.

   call self%error%destroy
   call self%cli%free
   call self%os%destroy
   call self%file_parameters%free
   call self%file_grid_input%destroy
   call self%file_grid_output%destroy
   call self%file_initial_conditions%destroy
   call self%file_solution%destroy
   call self%adimensionals%destroy
   call self%free_conditions%destroy
   call self%eos%destroy
   call self%solver%destroy
   call self%time%destroy
   call self%mesh%destroy
   self%is_cli_parsed = .false.
   self%is_output_verbose = .false.
   self%go_on_fail = .false.
   self%shell_refresh_frequency = 1
   endsubroutine destroy

   subroutine initialize(self, file_parameters)
   !< Initialize simulation.
   class(simulation_object), intent(inout)        :: self              !< Simulation data.
   character(*),             intent(in), optional :: file_parameters   !< File name of simulation parameters file.
   integer(I4P)                                   :: interfaces_number !< Number of different interfaces (level set).

   call self%destroy

   call self%error%initialize
   call self%cli_initialize
   call self%os%initialize()

   call self%file_parameters%initialize(filename=file_parameters)
   call self%file_grid_input%initialize()
   call self%file_grid_output%initialize()
   call self%file_initial_conditions%initialize()
   call self%file_solution%initialize()

   call self%adimensionals%initialize()
   call self%free_conditions%initialize()
   call self%eos%initialize()
   call self%solver%initialize()
   call self%time%initialize()

   call self%load_file_parameters(interfaces_number=interfaces_number)

   call self%mesh%initialize(eos=self%eos,                         &
                             file_grid=self%file_grid_input,       &
                             file_ic=self%file_initial_conditions, &
                             interfaces_number=interfaces_number)
   call self%mesh%immerge_immersed_boundaries(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   endsubroutine initialize

   subroutine integrate(self)
   !< Integrate the equations.
   class(simulation_object), intent(inout) :: self !< Simulation data.

   if (self%is_output_verbose) print '(A)', self%description()
   temporal_loop: do
      call self%compute_dt
      call self%time%update_dt
      associate(integrator=>self%solver%integrator)
         select type(integrator)
         class is(integrator_multistage_object)
            call integrator%integrate_fast(U=self, Dt=self%time%dt, t=self%time%t)
         endselect
      endassociate
      call self%time%print_eta(frequency=self%shell_refresh_frequency)
      call self%time%update_n
      call self%time%update_t
      call self%mesh%save_file_solution(file_solution=self%file_solution, n=self%time%n)
      if (self%time%is_the_end()) exit temporal_loop
   enddo temporal_loop
   endsubroutine integrate

   subroutine load_file_parameters(self, interfaces_number, file_name)
   !< Load file parameters.
   class(simulation_object), intent(inout)         :: self              !< Simulation object.
   integer(I4P),             intent(out), optional :: interfaces_number !< Number of different interfaces (level set).
   character(len=*),         intent(in),  optional :: file_name         !< File name.

   if (present(file_name)) call self%file_parameters%initialize(filename=file_name)
   if (self%is_output_verbose) print '(A)', 'load file "'//trim(adjustl(self%file_parameters%filename))//'"'
   call self%file_parameters%load(error=self%error%status)
   call self%error%check(message='failed to load "'//trim(adjustl(self%file_parameters%filename))// &
                         '" from procedure "simulation_object%load_file_parameters"', is_severe=.not.self%go_on_fail)
   if (self%error%status /=0) then
      write(stderr, '(A)') 'Using default simulation parameters values'
      return
   endif
   ! simulation miscellanea parameters
   call self%file_parameters%get(section_name='miscellanea', option_name='is_output_verbose', val=self%is_output_verbose, &
                                 error=self%error%status)
   call self%error%check(message='failed to load [miscellanea].(is_output_verbose)', is_severe=self%go_on_fail)
   call self%file_parameters%get(section_name='miscellanea', option_name='shell_refresh_frequency', &
                                 val=self%shell_refresh_frequency, error=self%error%status)
   call self%error%check(message='failed to load [miscellanea].(shell_refresh_frequency)', is_severe=self%go_on_fail)
   call self%file_parameters%get(section_name='miscellanea', option_name='go_on_fail', &
                                 val=self%go_on_fail, error=self%error%status)
   call self%error%check(message='failed to load [miscellanea].(go_on_fail)', is_severe=self%go_on_fail)

   ! files
   call self%file_grid_input%load_parameters_from_file(fini=self%file_parameters, options_prefix='grid_input', &
                                                       go_on_fail=self%go_on_fail)
   call self%file_grid_output%load_parameters_from_file(fini=self%file_parameters, options_prefix='grid_output', &
                                                        go_on_fail=self%go_on_fail)
   call self%file_initial_conditions%load_parameters_from_file(fini=self%file_parameters, options_prefix='initial_conditions', &
                                                               go_on_fail=self%go_on_fail)
   call self%file_solution%load_parameters_from_file(fini=self%file_parameters, options_prefix='solution', &
                                                     go_on_fail=self%go_on_fail)

   ! auxiliary
   call self%adimensionals%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)

   call self%free_conditions%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)

   call self%eos%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)

   call self%solver%load_from_file(fini=self%file_parameters, integrand_0=self, go_on_fail=self%go_on_fail)

   call self%time%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)

   ! IB
   if (present(interfaces_number)) then
      call self%file_parameters%get(section_name='immersed_boundary_bodies', option_name='files_number', &
                                    val=interfaces_number, error=self%error%status)
      if (self%error%status/=0) interfaces_number=0
   endif
   endsubroutine load_file_parameters

   subroutine parse_command_line_interface(self)
   !< Parse command line interface.
   class(simulation_object), intent(inout) :: self   !< Simulation data.
   integer(I4P)                            :: error  !< Error trapping flag.
   character(999)                          :: buffer !< Character buffer.

   call self%cli_initialize
   associate(cli=>self%cli)
      call cli%parse(error=error) ; if (error/=0) stop

      if (cli%is_passed(switch='--parameters_file')) then
         call cli%get(switch='--parameters_file', val=buffer, error=error) ; if (error/=0) stop
         call self%file_parameters%initialize(filename=buffer)
      endif

      if (cli%is_passed(switch='--go-on-fail')) then
         call cli%get(switch='--go-on-fail', val=self%go_on_fail, error=error) ; if (error/=0) stop
      endif

      if (cli%is_passed(switch='--verbose')) then
         call cli%get(switch='--verbose', val=self%is_output_verbose, error=error) ; if (error/=0) stop
      endif
   endassociate

   self%is_cli_parsed = .true.
   endsubroutine parse_command_line_interface

   subroutine save_file_grid(self, file_grid, is_parametric, file_name, ascii, metrics, off, tecplot, vtk, force)
   !< Save file grid.
   class(simulation_object), intent(inout)        :: self          !< Simulation object.
   type(file_grid_object),   intent(in), optional :: file_grid     !< File grid handler.
   logical,                  intent(in), optional :: is_parametric !< Sentinel to load grid parametric grid file.
   character(*),             intent(in), optional :: file_name     !< File name.
   logical,                  intent(in), optional :: ascii         !< Ascii/binary output.
   logical,                  intent(in), optional :: metrics       !< Save also metrics data.
   logical,                  intent(in), optional :: off           !< Save in OFF format sentinel.
   logical,                  intent(in), optional :: tecplot       !< Tecplot output format sentinel.
   logical,                  intent(in), optional :: vtk           !< VTK output format sentinel.
   logical,                  intent(in), optional :: force          !< Sentinel to force saving.

   call self%mesh%save_file_grid(file_grid=file_grid, file_name=file_name, is_parametric=is_parametric, &
                                 metrics=metrics, ascii=ascii, off=off, tecplot=tecplot, vtk=vtk, force=force)
   endsubroutine save_file_grid

   subroutine save_file_solution(self, file_solution, file_name, is_parametric, metrics, ascii, off, tecplot, vtk, n, force)
   !< Save file solution.
   class(simulation_object),   intent(inout)        :: self          !< Simulation object.
   type(file_solution_object), intent(in), optional :: file_solution !< File solution handler.
   character(*),               intent(in), optional :: file_name     !< File name.
   logical,                    intent(in), optional :: is_parametric !< Sentinel to load grid parametric grid file.
   logical,                    intent(in), optional :: metrics       !< Save metrics sentinel.
   logical,                    intent(in), optional :: ascii         !< Ascii/binary output.
   logical,                    intent(in), optional :: off           !< Save in OFF format sentinel.
   logical,                    intent(in), optional :: tecplot       !< Tecplot output format sentinel.
   logical,                    intent(in), optional :: vtk           !< VTK output format sentinel.
   integer(I8P),               intent(in), optional :: n             !< Time step.
   logical,                    intent(in), optional :: force          !< Sentinel to force saving.

   call self%mesh%save_file_solution(file_solution=file_solution, file_name=file_name, is_parametric=is_parametric, &
                                     metrics=metrics, ascii=ascii, off=off, tecplot=tecplot, vtk=vtk, n=n, force=force)
   endsubroutine save_file_solution

   subroutine save_file_parameters(self, file_name)
   !< Save file parameters.
   class(simulation_object), intent(inout) :: self      !< Simulation object.
   character(*),             intent(in)    :: file_name !< File name.
   type(file_ini)                          :: fini      !< Simulation parameters ini file handler.

   call self%adimensionals%save_into_file(fini=fini)
   call self%free_conditions%save_into_file(fini=fini)
   call self%eos%save_into_file(fini=fini)
   call self%solver%save_into_file(fini=fini)
   call self%time%save_into_file(fini=fini)
   call fini%save(filename=trim(adjustl(file_name)))
   endsubroutine save_file_parameters

   ! integrand_object deferred methods
   pure function integrand_dimension(self)
   !< Return integrand dimension.
   !<
   !< @TODO implement this.
   class(simulation_object), intent(in) :: self                !< Integrand.
   integer(I4P)                         :: integrand_dimension !< Integrand dimension.

   integrand_dimension = self%mesh%grid_dimensions%blocks_number
   endfunction integrand_dimension

   pure function residuals(self, t) result(dState_dt)
   !< Time derivative of block integrand, residuals function.
   class(simulation_object), intent(in)           :: self         !< Block.
   real(R8P),                intent(in), optional :: t            !< Time.
   real(R8P), allocatable                         :: dState_dt(:) !< Block time derivative.

   ! error residuals function to be implemented
   endfunction residuals

   pure function local_error(lhs, rhs) result(error)
   !< Estimate local truncation error between 2 integrand approximations.
   class(simulation_object), intent(in) :: lhs   !< Left hand side.
   class(integrand_object),  intent(in) :: rhs   !< Right hand side.
   real(R8P)                            :: error !< Error estimation.

   ! error local error function to be implemented
   endfunction local_error

   ! +
   pure function integrand_add_integrand(lhs, rhs) result(opr)
   !< `+` operator.
   class(simulation_object), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),  intent(in) :: rhs    !< Right hand side.
   real(R8P), allocatable               :: opr(:) !< Operator result.

   ! error add operator to be implemented
   endfunction integrand_add_integrand

   pure function integrand_add_real(lhs, rhs) result(opr)
   !< `+ real` operator.
   class(simulation_object), intent(in) :: lhs     !< Left hand side.
   real(R8P),                intent(in) :: rhs(1:) !< Right hand side.
   real(R8P), allocatable               :: opr(:)  !< Operator result.

   ! error add operator to be implemented
   endfunction integrand_add_real

   pure function real_add_integrand(lhs, rhs) result(opr)
   !< `real +` operator.
   real(R8P),                intent(in) :: lhs(1:) !< Left hand side.
   class(simulation_object), intent(in) :: rhs     !< Right hand side.
   real(R8P), allocatable               :: opr(:)  !< Operator result.

   ! error add operator to be implemented
   endfunction real_add_integrand

   ! *
   pure function integrand_multiply_integrand(lhs, rhs) result(opr)
   !< `*` operator.
   class(simulation_object), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),  intent(in) :: rhs    !< Right hand side.
   real(R8P), allocatable               :: opr(:) !< Operator result.

   ! error multiply operator to be implemented
   endfunction integrand_multiply_integrand

   pure function integrand_multiply_real(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(simulation_object),   intent(in) :: lhs     !< Left hand side.
   real(R8P),                  intent(in) :: rhs(1:) !< Right hand side.
   real(R8P), allocatable                 :: opr(:)  !< Operator result.

   ! error multiply operator to be implemented
   endfunction integrand_multiply_real

   pure function real_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   real(R8P),                intent(in) :: lhs(1:) !< Left hand side.
   class(simulation_object), intent(in) :: rhs     !< Right hand side.
   real(R8P), allocatable               :: opr(:)  !< Operator result.

   ! error multiply operator to be implemented
   endfunction real_multiply_integrand

   pure function integrand_multiply_real_scalar(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(simulation_object), intent(in) :: lhs    !< Left hand side.
   real(R8P),                intent(in) :: rhs    !< Right hand side.
   real(R8P), allocatable               :: opr(:) !< Operator result.

   ! error multiply operator to be implemented
   endfunction integrand_multiply_real_scalar

   pure function real_scalar_multiply_integrand(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   real(R8P),                intent(in) :: lhs    !< Left hand side.
   class(simulation_object), intent(in) :: rhs    !< Right hand side.
   real(R8P), allocatable               :: opr(:) !< Operator result.

   ! error multiply operator to be implemented
   endfunction real_scalar_multiply_integrand

   ! -
   pure function integrand_sub_integrand(lhs, rhs) result(opr)
   !< `-` operator.
   class(simulation_object), intent(in) :: lhs    !< Left hand side.
   class(integrand_object),  intent(in) :: rhs    !< Right hand side.
   real(R8P), allocatable               :: opr(:) !< Operator result.

   ! error subtract operator to be implemented
   endfunction integrand_sub_integrand

   pure function integrand_sub_real(lhs, rhs) result(opr)
   !< `- real` operator.
   class(simulation_object), intent(in) :: lhs     !< Left hand side.
   real(R8P),                intent(in) :: rhs(1:) !< Right hand side.
   real(R8P), allocatable               :: opr(:)  !< Operator result.

   ! error subtract operator to be implemented
   endfunction integrand_sub_real

   pure function real_sub_integrand(lhs, rhs) result(opr)
   !< `real -` operator.
   real(R8P),                intent(in) :: lhs(1:) !< Left hand side.
   class(simulation_object), intent(in) :: rhs     !< Left hand side.
   real(R8P), allocatable               :: opr(:)  !< Operator result.

   ! error subtract operator to be implemented
   endfunction real_sub_integrand

   ! =
   _PURE_ subroutine assign_integrand(lhs, rhs)
   !< `=` operator.
   class(simulation_object), intent(inout) :: lhs !< Left hand side.
   class(integrand_object),  intent(in)    :: rhs !< Right hand side.

   select type(rhs)
   class is(simulation_object)
      lhs%error                   = rhs%error
      lhs%cli                     = rhs%cli
      lhs%os                      = rhs%os
      lhs%file_parameters         = rhs%file_parameters
      lhs%file_grid_input         = rhs%file_grid_input
      lhs%file_grid_output        = rhs%file_grid_output
      lhs%file_solution           = rhs%file_solution
      lhs%adimensionals           = rhs%adimensionals
      lhs%free_conditions         = rhs%free_conditions
      lhs%eos                     = rhs%eos
      lhs%solver                  = rhs%solver
      lhs%time                    = rhs%time
      lhs%mesh                    = rhs%mesh
      lhs%is_cli_parsed           = rhs%is_cli_parsed
      lhs%is_output_verbose       = rhs%is_output_verbose
      lhs%go_on_fail              = rhs%go_on_fail
      lhs%shell_refresh_frequency = rhs%shell_refresh_frequency
   endselect
   endsubroutine assign_integrand

   pure subroutine assign_real(lhs, rhs)
   !< `= real` operator.
   class(simulation_object), intent(inout) :: lhs     !< Left hand side.
   real(R8P),                intent(in)    :: rhs(1:) !< Right hand side.

   ! error assign operator to be implemented
   endsubroutine assign_real

   ! fast operators
   subroutine t_fast(self, t)
   !< Time derivative function of integrand class, i.e. the residuals function. Fast mode acting directly on self.
   class(simulation_object), intent(inout)        :: self !< Integrand.
   real(R8P),                intent(in), optional :: t    !< Time.

   call self%mesh%impose_boundary_conditions
   call self%mesh%compute_residuals(solver=self%solver, gcu=self%solver%gcu)
   endsubroutine t_fast

   ! +
   pure subroutine integrand_add_integrand_fast(opr, lhs, rhs)
   !< `+` fast operator.
   class(simulation_object), intent(inout) :: opr !< Operator result.
   class(integrand_object),  intent(in)    :: lhs !< Left hand side.
   class(integrand_object),  intent(in)    :: rhs !< Right hand side.

   select type(lhs)
   type is(simulation_object)
      select type(rhs)
      type is(simulation_object)
         call opr%mesh%conservative_add_conservatives_fast(lhs=lhs%mesh, rhs=rhs%mesh)
      endselect
   endselect
   endsubroutine integrand_add_integrand_fast

   ! *
   pure subroutine integrand_multiply_integrand_fast(opr, lhs, rhs)
   !< `*` fast operator.
   class(simulation_object), intent(inout) :: opr !< Operator result.
   class(integrand_object),  intent(in)    :: lhs !< Left hand side.
   class(integrand_object),  intent(in)    :: rhs !< Right hand side.

   select type(lhs)
   type is(simulation_object)
      select type(rhs)
      type is(simulation_object)
         call opr%mesh%conservative_multiply_conservatives_fast(lhs=lhs%mesh, rhs=rhs%mesh)
      endselect
   endselect
   endsubroutine integrand_multiply_integrand_fast

   pure subroutine integrand_multiply_real_scalar_fast(opr, lhs, rhs)
   !< `* real_scalar` fast operator.
   class(simulation_object), intent(inout) :: opr !< Operator result.
   class(integrand_object),  intent(in)    :: lhs !< Left hand side.
   real(R8P),                intent(in)    :: rhs !< Right hand side.

   select type(lhs)
   type is(simulation_object)
      call opr%mesh%conservative_multiply_real_scalar_fast(lhs=lhs%mesh, rhs=rhs)
   endselect
   endsubroutine integrand_multiply_real_scalar_fast

   ! -
   pure subroutine integrand_subtract_integrand_fast(opr, lhs, rhs)
   !< `-` fast operator.
   class(simulation_object), intent(inout) :: opr !< Operator result.
   class(integrand_object),  intent(in)    :: lhs !< Left hand side.
   class(integrand_object),  intent(in)    :: rhs !< Right hand side.

   select type(lhs)
   type is(simulation_object)
      select type(rhs)
      type is(simulation_object)
         call opr%mesh%conservative_subtract_conservatives_fast(lhs=lhs%mesh, rhs=rhs%mesh)
      endselect
   endselect
   endsubroutine integrand_subtract_integrand_fast

   ! private methods
   subroutine cli_initialize(self)
   !< Initialize Command Line Interface.
   class(simulation_object), intent(inout) :: self  !< Simulation data.
   integer(I4P)                            :: error !< Error trapping flag.

   associate(cli=>self%cli)
      call cli%init(progname='off',                                              &
                    version='v0.0.1',                                            &
                    authors='G. Rossi, S. Zaghi',                                &
                    help='Usage: ',                                              &
                    description='Open source Finite volume Fluid dynamics code', &
                    examples=["off ..."],                                        &
                    epilog=new_line('a')//"all done")

      call cli%add(switch='--parameters_file',                &
                   switch_ab='-par',                          &
                   help='name of simulation parameters file', &
                   required=.false.,                          &
                   act='store',                               &
                   def='simulation_parameters.ini',           &
                   error=error)
      if (error/=0) stop

      call cli%add(switch='--grid_file',     &
                   switch_ab='-grd',         &
                   help='name of grid file', &
                   required=.false.,         &
                   act='store',              &
                   def='grid.grd',           &
                   error=error)
      if (error/=0) stop

      call cli%add(switch='--go-on-fail',                           &
                   switch_ab='-gof',                                &
                   help='go-on if parameters load fails somewhere', &
                   required=.false.,                                &
                   def='.false.',                                   &
                   act='store')
      if (error/=0) stop

      call cli%add(switch='--verbose',           &
                   help='enable verbose output', &
                   required=.false.,             &
                   act='store_true',             &
                   def='.false.',                &
                   error=error)
      if (error/=0) stop
   endassociate
   endsubroutine cli_initialize

   _PURE_ subroutine compute_dt(self)
   !< Compute the current time step by means of CFL condition.
   class(simulation_object), intent(inout) :: self !< Simulation data.
   integer(I4P)                            :: b    !< Counter.

   self%time%dt = MaxR8P
   do b=1, self%mesh%grid_dimensions%blocks_number
      call self%mesh%blocks(b)%compute_dt(CFL=self%time%CFL)
      self%time%dt = min(self%time%dt, self%mesh%blocks(b)%dt_min())
   enddo
   call self%time%update_dt
   endsubroutine compute_dt
endmodule off_simulation_object
