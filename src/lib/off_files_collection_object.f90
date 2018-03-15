!< OFF files collection object definition and implementation.

module off_files_collection_object
!< OFF files collection object definition and implementation.

use off_file_object, only : file_object
use finer, only : file_ini
use penf, only : I4P

implicit none
private
public :: files_collection_object

character(len=5), parameter :: INI_SECTION_NAME='files' !< INI (config) file section name containing the collection file names.

type :: files_collection_object
  !< Files collection object class.
  type(file_object) :: boundary_conditions !< Boundary conditions file.
  type(file_object) :: grid                !< Grid file.
  type(file_object) :: initial_conditions  !< Initial conditions file.
  type(file_object) :: logging             !< Logging file.
  type(file_object) :: solution            !< Solution file.
  contains
    ! public methods
    procedure, pass(self) :: description    !< Return a pretty-formatted description of the files collection.
    procedure, pass(self) :: destroy        !< Destroy files collection.
    procedure, pass(self) :: initialize     !< Initialize files collection.
    procedure, pass(self) :: load_from_file !< Load from file.
    procedure, pass(self) :: save_into_file !< Save into file.
endtype files_collection_object

contains
  ! public methods
  pure function description(self, prefix) result(desc)
  !< Return a pretty-formatted description of the files collection.
  class(files_collection_object), intent(in)           :: self             !< Files collection.
  character(*),                   intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                        :: desc             !< Description.
  character(len=:), allocatable                        :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                          :: NL=new_line('a') !< New line character.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  desc = ''
  desc = desc//prefix//'Boundary conditions file:'//NL
  desc = desc//prefix//self%boundary_conditions%description(prefix='  ')//NL
  desc = desc//prefix//'Grid file:'//NL
  desc = desc//prefix//self%grid%description(prefix='  ')//NL
  desc = desc//prefix//'Initial conditions file:'//NL
  desc = desc//prefix//self%initial_conditions%description(prefix='  ')//NL
  desc = desc//prefix//'Logging file:'//NL
  desc = desc//prefix//self%logging%description(prefix='  ')//NL
  desc = desc//prefix//'Solution file:'//NL
  desc = desc//prefix//self%solution%description(prefix='  ')
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy files collection.
  class(files_collection_object), intent(inout) :: self !< Files collection.

  call self%boundary_conditions%destroy
  call self%grid%destroy
  call self%initial_conditions%destroy
  call self%logging%destroy
  call self%solution%destroy
  endsubroutine destroy

  elemental subroutine initialize(self)
  !< Initialize Files collection.
  !<
  !< @TODO Complete this.
  class(files_collection_object), intent(inout) :: self !< Files collection.

  call self%destroy
  call self%boundary_conditions%initialize
  call self%grid%initialize
  call self%initial_conditions%initialize
  call self%logging%initialize
  call self%solution%initialize
  endsubroutine initialize

  subroutine load_from_file(self, fini, go_on_fail)
  !< Load from file.
  class(files_collection_object), intent(inout)        :: self       !< Files collection.
  type(file_ini),                 intent(in)           :: fini       !< Simulation parameters ini file handler.
  logical,                        intent(in), optional :: go_on_fail !< Go on if load fails..

  ! call self%boundary_conditions%load_file_name_from_file(fini=fini,                         &
  !                                                        section_name=INI_SECTION_NAME,     &
  !                                                        option_name='boundary_conditions', &
  !                                                        go_on_fail=go_on_fail)
  ! call self%grid%load_file_name_from_file(fini=fini,                     &
  !                                         section_name=INI_SECTION_NAME, &
  !                                         option_name='grid',            &
  !                                         go_on_fail=go_on_fail)
  ! call self%initial_conditions%load_file_name_from_file(fini=fini,                        &
  !                                                       section_name=INI_SECTION_NAME,    &
  !                                                       option_name='initial_conditions', &
  !                                                       go_on_fail=go_on_fail)
  ! call self%logging%load_file_name_from_file(fini=fini,                     &
  !                                            section_name=INI_SECTION_NAME, &
  !                                            option_name='log',             &
  !                                            go_on_fail=go_on_fail)
  ! call self%solution%load_file_name_from_file(fini=fini,                     &
  !                                             section_name=INI_SECTION_NAME, &
  !                                             option_name='solution',        &
  !                                             go_on_fail=go_on_fail)
  endsubroutine load_from_file

  subroutine save_into_file(self, fini)
  !< Save from file.
  class(files_collection_object), intent(inout) :: self !< Files collection.
  type(file_ini),                 intent(inout) :: fini !< Simulation parameters ini file handler.

  call self%boundary_conditions%save_file_name_into_file(fini=fini,                     &
                                                         section_name=INI_SECTION_NAME, &
                                                         option_name='boundary_conditions')

  call self%grid%save_file_name_into_file(fini=fini,                     &
                                          section_name=INI_SECTION_NAME, &
                                          option_name='grid')

  call self%initial_conditions%save_file_name_into_file(fini=fini,                     &
                                                        section_name=INI_SECTION_NAME, &
                                                        option_name='initial_conditions')

  call self%logging%save_file_name_into_file(fini=fini,                     &
                                             section_name=INI_SECTION_NAME, &
                                             option_name='log')

  call self%solution%save_file_name_into_file(fini=fini,                     &
                                              section_name=INI_SECTION_NAME, &
                                              option_name='solution')
  endsubroutine save_into_file
endmodule off_files_collection_object
