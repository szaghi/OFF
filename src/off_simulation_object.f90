!< OFF simulation object definition and implementation.

module off_simulation_object
!< OFF simulation object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_object, only : block_object
use off_error_object, only : error_object
use off_free_conditions_object, only : free_conditions_object
use off_file_grid_object, only : file_grid_object
use off_grid_dimensions_object, only : grid_dimensions_object
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
   type(file_grid_object)               :: file_grid                 !< Grid file handler.
   type(non_dimensional_numbers_object) :: adimensionals             !< Non dimensional numbers.
   type(free_conditions_object)         :: free_conditions           !< Free stream conditions.
   type(os_object)                      :: os                        !< Running Operating System.
   type(solver_object)                  :: solver                    !< solver models parameters.
   type(time_object)                    :: time                      !< Timing conditions.
   type(grid_dimensions_object)         :: grid_dimensions           !< Grid dimensions.
   type(block_object), allocatable      :: blocks(:)                 !< Blocks list.
   logical                              :: is_cli_parsed=.false.     !< Sentinel of CLI parsing.
   character(999)                       :: file_parameters=''        !< Name of simulation parameters file.
   logical                              :: is_output_verbose=.false. !< Verbose output.
   logical                              :: go_on_fail=.false.        !< Allow/disallow parameters loading failure.
   contains
      ! public methods
      procedure, pass(self) :: allocate_blocks              !< Allocate blocks accordingly to grid dimensions.
      procedure, pass(self) :: description                  !< Return a pretty-formatted description of the simulation.
      procedure, pass(self) :: destroy                      !< Destroy simulation data.
      procedure, pass(self) :: initialize                   !< Initialize simulation.
      procedure, pass(self) :: integrate                    !< Integrate the equations.
      procedure, pass(self) :: load_file_grid               !< Load grid file.
      procedure, pass(self) :: load_file_parameters         !< Load file parameters.
      procedure, pass(self) :: load_input_files             !< Load input files.
      procedure, pass(self) :: parse_command_line_interface !< Parse command line interface.
      procedure, pass(self) :: save_file_grid               !< Save grid file.
      procedure, pass(self) :: save_file_parameters         !< Save file parameters.
      ! private methods
      procedure, pass(self), private :: set_command_line_interface !< Set command line interface.
endtype simulation_object

contains
   ! public methods
   subroutine allocate_blocks(self)
   !< Allocate blocks accordingly to grid dimensions.
   class(simulation_object), intent(inout) :: self !< Simulation data.
   integer(I4P)                            :: b    !< Counter.

   if (self%grid_dimensions%blocks_number > 0) then
      if (allocated(self%blocks)) then
         call self%blocks%destroy
         deallocate(self%blocks)
      endif
      allocate(self%blocks(1:self%grid_dimensions%blocks_number))
      do b=1, self%grid_dimensions%blocks_number
         call self%blocks(b)%initialize(signature=self%grid_dimensions%block_signature(b))
      enddo
   endif
   endsubroutine allocate_blocks

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the simulation.
   class(simulation_object), intent(in)           :: self             !< Simulation parameters.
   character(*),             intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                  :: desc             !< Description.
   character(len=:), allocatable                  :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter                    :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   ! desc = desc//prefix_//'Files:'//NL
   ! desc = desc//prefix//self%files_collection%description(prefix='  ')//NL
   desc = desc//prefix_//'Non dimensional numbers:'//NL
   desc = desc//prefix_//self%adimensionals%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Free stream conditions:'//NL
   desc = desc//prefix_//self%free_conditions%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Solver models:'//NL
   desc = desc//prefix_//self%solver%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Time:'//NL
   desc = desc//prefix_//self%time%description(prefix=prefix_//'  ')//NL
   desc = desc//prefix_//'Grid dimensions:'//NL
   desc = desc//prefix_//self%grid_dimensions%description(prefix=prefix_//'  ')
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy simulation data.
   class(simulation_object), intent(inout) :: self !< Simulation data.

   call self%error%destroy
   call self%cli%free
   call self%adimensionals%destroy
   call self%file_grid%destroy
   call self%free_conditions%destroy
   call self%os%destroy
   call self%solver%destroy
   call self%time%destroy
   call self%grid_dimensions%destroy
   if (allocated(self%blocks)) then
      call self%blocks%destroy
      deallocate(self%blocks)
   endif
   self%is_cli_parsed = .false.
   self%file_parameters = ''
   self%is_output_verbose = .false.
   self%go_on_fail = .false.
   endsubroutine destroy

   subroutine initialize(self, blocks, parse_cli, load_files)
   !< Initialize simulation.
   class(simulation_object), intent(inout)        :: self        !< simulation data.
   type(block_object),       intent(in), optional :: blocks(1:)  !< Blocks structure.
   logical,                  intent(in), optional :: parse_cli   !< Enable CLI parsing.
   logical,                  intent(in), optional :: load_files  !< Enable files loading.
   logical                                        :: parse_cli_  !< Enable CLI parsing, local variable.
   logical                                        :: load_files_ !< Enable files loading, local variable.

   parse_cli_  = .false. ; if (present(parse_cli )) parse_cli_  = parse_cli
   load_files_ = .false. ; if (present(load_files)) load_files_ = load_files

   call self%destroy

   call self%error%initialize
   call self%adimensionals%initialize
   call self%file_grid%initialize
   call self%free_conditions%initialize
   call self%os%initialize
   call self%solver%initialize
   call self%time%initialize
   call self%grid_dimensions%initialize

   if (parse_cli_) call self%parse_command_line_interface

   if (load_files_) call self%load_input_files

   if (present(blocks)) then
      call self%grid_dimensions%initialize(block_signature=blocks%signature)
      allocate(self%blocks(1:size(blocks, dim=1)), source=blocks)
   endif

   if (self%is_output_verbose) print '(A)', self%description()
   endsubroutine initialize

   subroutine integrate(self)
   !< Integrate the equations.
   !<
   !< @TODO Implement this.
   class(simulation_object), intent(inout) :: self !< simulation data.

   error stop 'error: simulation_object%integrate to be implemented'
   endsubroutine integrate

   subroutine load_file_grid(self, file_basename)
   !< Load grid file.
   class(simulation_object), intent(inout) :: self          !< simulation data.
   character(*),             intent(in)    :: file_basename !< File basename.

   call self%file_grid%load_grid_dimensions_from_file(file_name=trim(adjustl(file_basename))//'.grd', &
                                                      grid_dimensions=self%grid_dimensions)
   call self%allocate_blocks
   call self%file_grid%load_nodes_from_file(file_name=trim(adjustl(file_basename))//'.grd', &
                                            grid_dimensions=self%grid_dimensions, blocks=self%blocks)
   endsubroutine load_file_grid

   subroutine load_file_parameters(self, file_name, go_on_fail)
   !< Load file parameters.
   class(simulation_object), intent(inout)        :: self        !< Simulation object.
   character(*),             intent(in)           :: file_name   !< File name.
   logical,                  intent(in), optional :: go_on_fail  !< Go on if load fails.
   logical                                        :: go_on_fail_ !< Go on if load fails, local variable.
   type(file_ini)                                 :: fini        !< Simulation parameters ini file handler.

   go_on_fail_ = .true. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call fini%load(filename=file_name, error=self%error%status)
   call self%error%check(message='failed to load "'//trim(adjustl(file_name))//'"', is_severe=.not.go_on_fail_)
   if (self%error%status /=0 ) then
      write(stderr, '(A)') 'Using default simulation parameters values'
      return
   endif

   call self%file_grid%load_file_name_from_file(fini=fini, section_name='files', option_name='grid', go_on_fail=go_on_fail)
   ! call self%files_collection%load_from_file(fini=fini, go_on_fail=go_on_fail_)
   call self%adimensionals%load_from_file(fini=fini, go_on_fail=go_on_fail_)
   call self%free_conditions%load_from_file(fini=fini, go_on_fail=go_on_fail_)
   call self%solver%load_from_file(fini=fini, go_on_fail=go_on_fail_)
   call self%time%load_from_file(fini=fini, go_on_fail=go_on_fail_)
   endsubroutine load_file_parameters

   subroutine load_input_files(self)
   !< Load from file.
   class(simulation_object), intent(inout) :: self !< simulation data.

   if (self%is_cli_parsed) then
      if (self%is_output_verbose) print '(A)', 'load file "'//trim(adjustl(self%file_parameters))//'"'
      call self%load_file_parameters(file_name=self%file_parameters, go_on_fail=self%go_on_fail)
   else
      error stop 'error: before loading input files, Command Line Interface must be parsed!'
   endif
   endsubroutine load_input_files

   subroutine parse_command_line_interface(self)
   !< Parse command line interface.
   class(simulation_object), intent(inout) :: self  !< simulation data.
   integer(I4P)                            :: error !< Error trapping flag.

   call self%set_command_line_interface
   call self%cli%parse(error=error) ; if (error/=0) stop

   call self%cli%get(switch='--parameters_file', val=self%file_parameters,   error=error) ; if (error/=0) stop
   call self%cli%get(switch='--go-on-fail',      val=self%go_on_fail,        error=error) ; if (error/=0) stop
   call self%cli%get(switch='--verbose',         val=self%is_output_verbose, error=error) ; if (error/=0) stop
   self%is_cli_parsed = .true.
   endsubroutine parse_command_line_interface

   subroutine save_file_grid(self, file_basename, ascii, metrics, off, tecplot, vtk)
   !< Save grid file.
   class(simulation_object), intent(inout)        :: self          !< simulation data.
   character(*),             intent(in)           :: file_basename !< File basename.
   logical,                  intent(in), optional :: ascii         !< Ascii/binary output.
   logical,                  intent(in), optional :: metrics       !< Save also metrics data.
   logical,                  intent(in), optional :: off           !< Save in OFF format sentinel.
   logical,                  intent(in), optional :: tecplot       !< Tecplot output format sentinel.
   logical,                  intent(in), optional :: vtk           !< VTK output format sentinel.
   integer(I4P)                                   :: b             !< Counter.
   logical                                        :: off_          !< OFF format sentinel, local variable.
   logical                                        :: vtk_          !< VTK format sentinel, local variable.
   character(len=:), allocatable                  :: file_name     !< File name buffer.

   off_ = .true.  ; if (present(off)) off_ = off
   vtk_    = .false. ; if (present(vtk   )) vtk_    = vtk

   if (off_) then
      call self%file_grid%initialize(file_name=trim(adjustl(file_basename))//'.grd')
      call self%file_grid%save_grid_dimensions_into_file(grid_dimensions=self%grid_dimensions)
      call self%file_grid%save_nodes_into_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
   endif

   if (vtk_) then
      do b=1, self%grid_dimensions%blocks_number
         file_name = trim(adjustl(file_basename))//'-block'//                 &
                     '-id_'//trim(str(n=self%blocks(b)%signature%id, no_sign=.true.))// &
                     '-lv_'//trim(str(n=self%blocks(b)%signature%level, no_sign=.true.))//'.vts'
         call self%blocks(b)%save_file_grid(file_name=file_name, ascii=ascii, metrics=metrics, vtk=vtk)
      enddo
   endif
   endsubroutine save_file_grid

   subroutine save_file_parameters(self, file_name)
   !< Save file parameters.
   class(simulation_object), intent(inout) :: self      !< Simulation object.
   character(*),             intent(in)    :: file_name !< File name.
   type(file_ini)                          :: fini      !< Simulation parameters ini file handler.

   ! call self%files_collection%save_into_file(fini=fini)
   call self%adimensionals%save_into_file(fini=fini)
   call self%free_conditions%save_into_file(fini=fini)
   call self%solver%save_into_file(fini=fini)
   call self%time%save_into_file(fini=fini)
   call fini%save(filename=trim(adjustl(file_name)))
   endsubroutine save_file_parameters

   ! private methods
   subroutine set_command_line_interface(self)
   !< Set command line interface.
   class(simulation_object), intent(inout) :: self  !< simulation data.
   integer(I4P)                            :: error !< Error trapping flag.

   call self%cli%init(progname='off',                                             &
                      version='v0.0.1',                                           &
                      authors='G. Rossi, S. Zaghi',                               &
                      help='Usage: ',                                             &
                      description='CNR-INSEAN CFD code',                          &
                      examples=["off --parameters sim_parameters.ini --verbose"], &
                      epilog=new_line('a')//"all done")

   call self%cli%add(switch='--parameters_file',                &
                     switch_ab='-par',                          &
                     help='name of simulation parameters file', &
                     required=.false.,                          &
                     act='store',                               &
                     def='simulation_parameters.ini',           &
                     error=error)
   if (error/=0) stop

   call self%cli%add(switch='--go-on-fail',                           &
                     switch_ab='-gof',                                &
                     help='go-on if parameters load fails somewhere', &
                     required=.false.,                                &
                     def='.false.',                                   &
                     act='store')
   if (error/=0) stop

   call self%cli%add(switch='--verbose',           &
                     help='enable verbose output', &
                     required=.false.,             &
                     act='store_true',             &
                     def='.false.',                &
                     error=error)
   if (error/=0) stop
   endsubroutine set_command_line_interface
endmodule off_simulation_object
