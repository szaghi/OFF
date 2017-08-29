!< OFF OS object definition and implementation.

module off_os_object
!< OFF OS object definition and implementation.

use off_error_object, only : error_object
use penf, only : I4P
use stringifor, only : string

implicit none
private
public :: os_object

integer(I4P), parameter :: NO_ERROR            = 0 !< No errors occurred.
integer(I4P), parameter :: ERROR_CP_FAILED     = 1 !< Failed to copy files/directories.
integer(I4P), parameter :: ERROR_MKDIR_FAILED  = 2 !< Failed to create directory.
integer(I4P), parameter :: ERROR_RM_FAILED     = 3 !< Failed to remove files/directories.
integer(I4P), parameter :: ERROR_FALLBACK_INIT = 4 !< Failed to set custom OS, used fallback system (unix).

type :: os_object
   !< OS object class.
   !<
   !< This class is designed as an helper for performing system calls, e.g. make directory, copy files, etc...
   type(error_object)            :: error           !< Error handler.
   character(len=:), allocatable :: path_separator  !< Path seperator, e.g. "/" for unix-like systems.
   character(len=:), allocatable :: cp_dir_command  !< Copy directory command.
   character(len=:), allocatable :: cp_file_command !< Copy file command.
   character(len=:), allocatable :: mkdir_command   !< Make directory command.
   character(len=:), allocatable :: rm_dir_command  !< Remove directory command.
   character(len=:), allocatable :: rm_file_command !< Remove file command.
   contains
      ! public methods
      procedure, pass(self) :: cp         !< Copy files/directories.
      procedure, pass(self) :: destroy    !< Destroy OS... not your :-)
      procedure, pass(self) :: initialize !< Initialze OS.
      procedure, pass(self) :: mkdir      !< Make directory.
      procedure, pass(self) :: rm         !< Remove files/directories.
      ! operators
      generic :: assignment(=) => os_assign_os !< Overload `=`.
      ! private methods
      procedure, pass(self) :: initialize_unix    !< Initialze OS as unix-like system.
      procedure, pass(self) :: initialize_windows !< Initialze OS as windows-like system.
      procedure, pass(lhs)  :: os_assign_os       !< Operator `=`.
endtype os_object

contains
  ! public methods
  subroutine cp(self, file_name, dir_name)
  !< Copy files/directories.
  !<
  !< @note leading and trailing white spaces are trimmed out.
  class(os_object), intent(inout)        :: self      !< OS.
  character(*),     intent(in), optional :: file_name !< File name.
  character(*),     intent(in), optional :: dir_name  !< Dir name.

  if (present(file_name)) then
     call execute_command_line(command=self%cp_file_command//' '//trim(adjustl(file_name)), exitstat=self%error%status)
  elseif (present(dir_name)) then
     call execute_command_line(command=self%cp_dir_command//' '//trim(adjustl(dir_name)), exitstat=self%error%status)
  else
     self%error%status = ERROR_CP_FAILED
  endif
  endsubroutine cp

  elemental subroutine destroy(self)
  !< Destroy OS... not your :-)
  class(os_object), intent(inout) :: self  !< OS.
  type(os_object)                 :: fresh !< Fresh instance of OS.

  self = fresh
  if (allocated(self%path_separator))  deallocate(self%path_separator)
  if (allocated(self%cp_dir_command))  deallocate(self%cp_dir_command)
  if (allocated(self%cp_file_command)) deallocate(self%cp_file_command)
  if (allocated(self%mkdir_command))   deallocate(self%mkdir_command)
  if (allocated(self%rm_dir_command))  deallocate(self%rm_dir_command)
  if (allocated(self%rm_file_command)) deallocate(self%rm_file_command)
  endsubroutine destroy

  elemental subroutine initialize(self, system)
  !< Initialize OS.
  class(os_object), intent(inout)        :: self    !< OS.
  character(*),     intent(in), optional :: system  !< System name, valid [unix, windows].
  type(string)                           :: system_ !< System name.

  self%error%status = 0
  system_ = 'UNIX'
  if (present(system)) then
     system_ = trim(adjustl(system))
     system_ = system_%upper()
  endif
  select case(system_%chars())
  case('UNIX')
     call self%initialize_unix
  case('WINDOWS')
     call self%initialize_windows
  case default
     call self%initialize_unix
     self%error%status = ERROR_FALLBACK_INIT
  endselect
  endsubroutine initialize

  subroutine mkdir(self, dir_name)
  !< Make directoriy.
  !<
  !< @note leading and trailing white spaces are trimmed out.
  class(os_object), intent(inout) :: self     !< OS.
  character(*),     intent(in)    :: dir_name !< Dir name.

  call execute_command_line(command=self%mkdir_command//' '//trim(adjustl(dir_name)), exitstat=self%error%status)
  endsubroutine mkdir

  subroutine rm(self, file_name, dir_name)
  !< Remove files/directories.
  !<
  !< @note leading and trailing white spaces are trimmed out.
  class(os_object), intent(inout)        :: self      !< OS.
  character(*),     intent(in), optional :: file_name !< File name.
  character(*),     intent(in), optional :: dir_name  !< Dir name.

  if (present(file_name)) then
     call execute_command_line(command=self%rm_file_command//' '//trim(adjustl(file_name)), exitstat=self%error%status)
  elseif (present(dir_name)) then
     call execute_command_line(command=self%rm_dir_command//' '//trim(adjustl(dir_name)), exitstat=self%error%status)
  else
     self%error%status = ERROR_RM_FAILED
  endif
  endsubroutine rm

  ! private methods
  elemental subroutine initialize_unix(self)
  !< Initialize OS as unix-like system.
  class(os_object), intent(inout) :: self !< OS.

  call self%destroy
  self%path_separator  = char(47)
  self%cp_dir_command  = 'cp -r'
  self%cp_file_command = 'cp'
  self%mkdir_command   = 'mkdir -p'
  self%rm_dir_command  = 'rm -fr'
  self%rm_file_command = 'rm -f'
  endsubroutine initialize_unix

  elemental subroutine initialize_windows(self)
  !< Initialize OS as windows-like system.
  class(os_object), intent(inout) :: self !< OS.

  call self%destroy
  self%path_separator  = char(92)
  self%cp_dir_command  = 'copy'
  self%cp_file_command = 'copy'
  self%mkdir_command   = 'mkdir'
  self%rm_dir_command  = 'del'
  self%rm_file_command = 'del'
  endsubroutine initialize_windows

   pure subroutine os_assign_os(lhs, rhs)
   !< Operator `=`.
   class(os_object), intent(inout) :: lhs !< Left hand side.
   type(os_object),  intent(in)    :: rhs !< Right hand side.

   lhs%error = rhs%error
  if (allocated(rhs%path_separator))  lhs%path_separator  = rhs%path_separator
  if (allocated(rhs%cp_dir_command))  lhs%cp_dir_command  = rhs%cp_dir_command
  if (allocated(rhs%cp_file_command)) lhs%cp_file_command = rhs%cp_file_command
  if (allocated(rhs%mkdir_command))   lhs%mkdir_command   = rhs%mkdir_command
  if (allocated(rhs%rm_dir_command))  lhs%rm_dir_command  = rhs%rm_dir_command
  if (allocated(rhs%rm_file_command)) lhs%rm_file_command = rhs%rm_file_command
   endsubroutine os_assign_os
endmodule off_os_object
