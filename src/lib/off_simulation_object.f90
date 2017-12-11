!< OFF simulation object definition and implementation.

module off_simulation_object
!< OFF simulation object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_error_object, only : error_object
use off_free_conditions_object, only : free_conditions_object
use off_mesh_object, only : mesh_object
use off_non_dimensional_numbers_object, only : non_dimensional_numbers_object
use off_os_object, only : os_object
use off_solver_object, only : solver_object
use off_time_object, only : time_object
use finer, only : file_ini
use flap, only : command_line_interface
use flow, only : eos_compressible
use foodie, only : integrand_object
use penf, only : I4P, MaxR8P, R8P, str

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
   type(non_dimensional_numbers_object) :: adimensionals             !< Non dimensional numbers.
   type(free_conditions_object)         :: free_conditions           !< Free stream conditions.
   type(eos_compressible)               :: eos                       !< Equation of state.
   type(solver_object)                  :: solver                    !< solver models parameters.
   type(time_object)                    :: time                      !< Timing conditions.
   type(mesh_object)                    :: mesh                      !< Mesh data.
   logical                              :: is_cli_parsed=.false.     !< Sentinel of CLI parsing.
   logical                              :: is_output_verbose=.false. !< Verbose output.
   logical                              :: go_on_fail=.false.        !< Allow/disallow parameters loading failure.
   contains
      ! public methods
      procedure, pass(self) :: description                  !< Return a pretty-formatted description of the simulation.
      procedure, pass(self) :: destroy                      !< Destroy simulation data.
      procedure, pass(self) :: initialize                   !< Initialize simulation.
      procedure, pass(self) :: integrate                    !< Integrate the equations.
      procedure, pass(self) :: load_file_parameters         !< Load file parameters.
      procedure, pass(self) :: parse_command_line_interface !< Parse command line interface.
      procedure, pass(self) :: save_file_parameters         !< Save file parameters.

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
      procedure, pass(self), private :: cli_initialize             !< Initialize Command Line Interface.
      procedure, pass(self), private :: compute_dt                 !< Compute the current time step by means of CFL condition.
      procedure, pass(self), private :: impose_boundary_conditions !< Impose boundary conditions.
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

   elemental subroutine destroy(self)
   !< Destroy simulation data.
   class(simulation_object), intent(inout) :: self !< Simulation data.

   call self%error%destroy
   call self%cli%free
   call self%os%destroy
   call self%file_parameters%free
   call self%adimensionals%destroy
   call self%free_conditions%destroy
   call self%eos%destroy
   call self%solver%destroy
   call self%time%destroy
   call self%mesh%destroy
   self%is_cli_parsed = .false.
   self%is_output_verbose = .false.
   self%go_on_fail = .false.
   endsubroutine destroy

   subroutine load_file_parameters(self, file_name)
   !< Load file parameters.
   class(simulation_object), intent(inout)        :: self      !< Simulation object.
   character(len=*),         intent(in), optional :: file_name !< File name.

   if (present(file_name)) call self%file_parameters%initialize(filename=file_name)
   if (self%is_output_verbose) print '(A)', 'load file "'//trim(adjustl(self%file_parameters%filename))//'"'
   call self%file_parameters%load(error=self%error%status)
   call self%error%check(message='failed to load "'//trim(adjustl(self%file_parameters%filename))//'"', &
                         is_severe=.not.self%go_on_fail)
   if (self%error%status /=0 ) then
      write(stderr, '(A)') 'Using default simulation parameters values'
      return
   endif
   call self%mesh%file_grid%load_file_name_from_file(fini=self%file_parameters, &
                                                     section_name='files', option_name='grid', go_on_fail=self%go_on_fail)
   call self%adimensionals%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   call self%free_conditions%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   call self%eos%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   call self%solver%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   call self%time%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   endsubroutine load_file_parameters

   subroutine initialize(self, os, file_parameters, adimensionals, free_conditions, eos, mesh, solver, time)
   !< Initialize simulation.
   class(simulation_object),             intent(inout)        :: self            !< Simulation data.
   type(os_object),                      intent(in), optional :: os              !< Running Operating System.
   character(*),                         intent(in), optional :: file_parameters !< File name of simulation parameters file.
   type(non_dimensional_numbers_object), intent(in), optional :: adimensionals   !< Non dimensional numbers values.
   type(free_conditions_object),         intent(in), optional :: free_conditions !< Free conditions values.
   type(eos_compressible),               intent(in), optional :: eos             !< Equation of state.
   type(mesh_object),                    intent(in), optional :: mesh            !< Mesh data.
   type(solver_object),                  intent(in), optional :: solver          !< Solver data.
   type(time_object),                    intent(in), optional :: time            !< Time data.

   call self%destroy
   call self%error%initialize
   call self%cli_initialize
   call self%os%initialize(os=os)
   call self%file_parameters%initialize(filename=file_parameters)
   call self%adimensionals%initialize(adimensionals=adimensionals)
   call self%free_conditions%initialize(free_conditions=free_conditions)
   call self%eos%initialize(eos=eos)
   call self%mesh%initialize(mesh=mesh)
   call self%solver%initialize(solver=solver)
   call self%time%initialize(time=time)
   endsubroutine initialize

   subroutine integrate(self)
   !< Integrate the equations.
   !<
   !< @TODO Implement this.
   class(simulation_object), intent(inout) :: self !< Simulation data.

   error stop 'error: simulation_object%integrate to be implemented'
   endsubroutine integrate

   subroutine parse_command_line_interface(self)
   !< Parse command line interface.
   class(simulation_object), intent(inout) :: self            !< Simulation data.
   integer(I4P)                            :: error           !< Error trapping flag.
   character(999)                          :: file_parameters !< Name of simulation parameters file.

   call self%cli_initialize
   call self%cli%parse(error=error) ; if (error/=0) stop

   call self%cli%get(switch='--parameters_file', val=file_parameters,        error=error) ; if (error/=0) stop
   call self%cli%get(switch='--go-on-fail',      val=self%go_on_fail,        error=error) ; if (error/=0) stop
   call self%cli%get(switch='--verbose',         val=self%is_output_verbose, error=error) ; if (error/=0) stop

   call self%file_parameters%initialize(filename=file_parameters)

   self%is_cli_parsed = .true.
   endsubroutine parse_command_line_interface

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
   pure subroutine assign_integrand(lhs, rhs)
   !< `=` operator.
   class(simulation_object), intent(inout) :: lhs !< Left hand side.
   class(integrand_object),  intent(in)    :: rhs !< Right hand side.

   select type(rhs)
   class is(simulation_object)
   endselect
   endsubroutine assign_integrand

   pure subroutine assign_real(lhs, rhs)
   !< `= real` operator.
   class(simulation_object), intent(inout) :: lhs     !< Left hand side.
   real(R8P),                intent(in)    :: rhs(1:) !< Right hand side.

   ! error assign operator to be implemented
   endsubroutine assign_real

   ! fast operators
   ! time derivative
   subroutine t_fast(self, t)
   !< Time derivative function of integrand class, i.e. the residuals function. Fast mode acting directly on self.
   class(simulation_object), intent(inout)        :: self !< Integrand.
   real(R8P),                intent(in), optional :: t    !< Time.

   ! error time derivate fast to be implemented
   endsubroutine t_fast

   ! +
   pure subroutine integrand_add_integrand_fast(opr, lhs, rhs)
   !< `+` fast operator.
   class(simulation_object), intent(inout) :: opr !< Operator result.
   class(integrand_object),  intent(in)    :: lhs !< Left hand side.
   class(integrand_object),  intent(in)    :: rhs !< Right hand side.

   ! error add operator fast to be implemented
   endsubroutine integrand_add_integrand_fast

   ! *
   pure subroutine integrand_multiply_integrand_fast(opr, lhs, rhs)
   !< `*` fast operator.
   class(simulation_object), intent(inout) :: opr !< Operator result.
   class(integrand_object),  intent(in)    :: lhs !< Left hand side.
   class(integrand_object),  intent(in)    :: rhs !< Right hand side.

   ! error multiply operator fast to be implemented
   endsubroutine integrand_multiply_integrand_fast

   pure subroutine integrand_multiply_real_scalar_fast(opr, lhs, rhs)
   !< `* real_scalar` fast operator.
   class(simulation_object), intent(inout) :: opr !< Operator result.
   class(integrand_object),  intent(in)    :: lhs !< Left hand side.
   real(R8P),                intent(in)    :: rhs !< Right hand side.

   ! error multiply operator fast to be implemented
   endsubroutine integrand_multiply_real_scalar_fast

   ! -
   pure subroutine integrand_subtract_integrand_fast(opr, lhs, rhs)
   !< `-` fast operator.
   class(simulation_object), intent(inout) :: opr !< Operator result.
   class(integrand_object),  intent(in)    :: lhs !< Left hand side.
   class(integrand_object),  intent(in)    :: rhs !< Right hand side.

   ! error subtract operator fast to be implemented
   endsubroutine integrand_subtract_integrand_fast

   ! private methods
   subroutine cli_initialize(self)
   !< Initialize Command Line Interface.
   class(simulation_object), intent(inout) :: self  !< Simulation data.
   integer(I4P)                            :: error !< Error trapping flag.

   associate(cli=>self%cli)
      call cli%init(progname='off',                                             &
                    version='v0.0.1',                                           &
                    authors='G. Rossi, S. Zaghi',                               &
                    help='Usage: ',                                             &
                    description='CNR-INSEAN CFD code',                          &
                    examples=["off --parameters sim_parameters.ini --verbose"], &
                    epilog=new_line('a')//"all done")

      call cli%add(switch='--parameters_file',                &
                   switch_ab='-par',                          &
                   help='name of simulation parameters file', &
                   required=.false.,                          &
                   act='store',                               &
                   def='simulation_parameters.ini',           &
                   error=error)
      if (error/=0) stop

      call cli%add(switch='--parametric_grid',          &
                   switch_ab='-pgrid',                  &
                   help='name of parametric grid file', &
                   required=.false.,                    &
                   act='store',                         &
                   def='parametric_grid.ini',           &
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

   pure subroutine compute_dt(self)
   !< Compute the current time step by means of CFL condition.
   class(simulation_object), intent(inout) :: self !< Simulation data.
   integer(I4P)                            :: b    !< Counter.

   self%time%dt = MaxR8P
   do b=1, self%mesh%grid_dimensions%blocks_number
      self%time%dt = min(self%time%dt, minval(self%mesh%blocks(b)%cell%dt))
   enddo
   call self%time%update_dt
   endsubroutine compute_dt

   subroutine impose_boundary_conditions(self)
   class(simulation_object), intent(inout) :: self !< Simulation data.
   endsubroutine impose_boundary_conditions
endmodule off_simulation_object
