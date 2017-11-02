!< OFF simulation object definition and implementation.

module off_simulation_object
!< OFF simulation object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_error_object, only : error_object
use off_free_conditions_object, only : free_conditions_object
use off_grid_object, only : grid_object
use off_non_dimensional_numbers_object, only : non_dimensional_numbers_object
use off_os_object, only : os_object
use off_solver_object, only : solver_object
use off_time_object, only : time_object
use finer, only : file_ini
use flap, only : command_line_interface
use penf, only : I4P, str

implicit none
private
public :: simulation_object

type :: simulation_object
   !< Simulation object class.
   !<
   !< [[simulation_object]] is a container for all simulation data for each processor/image.
   type(error_object)                   :: error                     !< Errors handler.
   type(command_line_interface)         :: cli                       !< Command line interface.
   type(os_object)                      :: os                        !< Running Operating System.
   type(file_ini)                       :: file_parameters           !< Simulation parameters file handler.
   type(non_dimensional_numbers_object) :: adimensionals             !< Non dimensional numbers.
   type(free_conditions_object)         :: free_conditions           !< Free stream conditions.
   type(solver_object)                  :: solver                    !< solver models parameters.
   type(time_object)                    :: time                      !< Timing conditions.
   type(grid_object)                    :: grid                      !< Grid data.
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
      ! private methods
      procedure, pass(self), private :: cli_initialize !< Initialize Command Line Interface.
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
   desc = desc//prefix_//'Solver models:'//NL
   desc = desc//prefix_//self%solver%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Time:'//NL
   desc = desc//prefix_//self%time%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Grid:'//NL
   desc = desc//prefix_//self%grid%description(prefix=prefix_//'  ')
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
   call self%solver%destroy
   call self%time%destroy
   call self%grid%destroy
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
   call self%grid%file_grid%load_file_name_from_file(fini=self%file_parameters, &
                                                     section_name='files', option_name='grid', go_on_fail=self%go_on_fail)
   call self%adimensionals%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   call self%free_conditions%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   call self%solver%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   call self%time%load_from_file(fini=self%file_parameters, go_on_fail=self%go_on_fail)
   endsubroutine load_file_parameters

   subroutine initialize(self, os, file_parameters, adimensionals, free_conditions, grid, solver, time)
   !< Initialize simulation.
   class(simulation_object),             intent(inout)        :: self            !< simulation data.
   type(os_object),                      intent(in), optional :: os              !< Running Operating System.
   character(*),                         intent(in), optional :: file_parameters !< File name of simulation parameters file.
   type(non_dimensional_numbers_object), intent(in), optional :: adimensionals   !< Non dimensional numbers values.
   type(free_conditions_object),         intent(in), optional :: free_conditions !< Free conditions values.
   type(grid_object),                    intent(in), optional :: grid            !< Grid data.
   type(solver_object),                  intent(in), optional :: solver          !< Solver data.
   type(time_object),                    intent(in), optional :: time            !< Time data.

   call self%destroy
   call self%error%initialize
   call self%cli_initialize
   call self%os%initialize(os=os)
   call self%file_parameters%initialize(filename=file_parameters)
   call self%adimensionals%initialize(adimensionals=adimensionals)
   call self%free_conditions%initialize(free_conditions=free_conditions)
   call self%grid%initialize(grid=grid)
   call self%solver%initialize(solver=solver)
   call self%time%initialize(time=time)
   endsubroutine initialize

   subroutine integrate(self)
   !< Integrate the equations.
   !<
   !< @TODO Implement this.
   class(simulation_object), intent(inout) :: self !< simulation data.

   error stop 'error: simulation_object%integrate to be implemented'
   endsubroutine integrate

   subroutine parse_command_line_interface(self)
   !< Parse command line interface.
   class(simulation_object), intent(inout) :: self            !< simulation data.
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
   call self%solver%save_into_file(fini=fini)
   call self%time%save_into_file(fini=fini)
   call fini%save(filename=trim(adjustl(file_name)))
   endsubroutine save_file_parameters

   ! private methods
   subroutine cli_initialize(self)
   !< Initialize Command Line Interface.
   class(simulation_object), intent(inout) :: self  !< simulation data.
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
endmodule off_simulation_object
