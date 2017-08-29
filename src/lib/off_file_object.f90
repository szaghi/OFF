!< OFF file object definition and implementation.

module off_file_object
!< OFF file object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_error_object, only : error_object
use finer, only : file_ini
use penf, only : I4P

implicit none
private
public :: file_object
public :: ERROR_ALREADY_CONNECTED
public :: ERROR_NOT_CONNECTED
public :: ERROR_NOT_INITIALIZED

character(len=5), parameter :: UNSET_FILE_NAME='unset'     !< Default, unset file name.
integer(I4P),     parameter :: ERROR_ALREADY_CONNECTED = 1 !< Already connected error code.
integer(I4P),     parameter :: ERROR_NOT_CONNECTED     = 2 !< Not connected error code.
integer(I4P),     parameter :: ERROR_NOT_INITIALIZED   = 3 !< Not initialized error code.

type :: file_object
   !< File object class.
   type(error_object)            :: error                  !< Errors handler.
   character(len=:), allocatable :: file_name              !< File name.
   integer(I4P)                  :: file_unit=0            !< File unit.
   logical                       :: is_initialized=.false. !< Sentinel to check if file is initialized.
   logical                       :: is_connected=.false.   !< Sentinel to check if file is connected.
   contains
      ! public methods
      procedure, pass(self) :: close                    !< Close file.
      procedure, pass(self) :: description              !< Return a pretty-formatted description of the file.
      procedure, pass(self) :: destroy                  !< Destroy file.
      procedure, pass(self) :: initialize               !< Initialize file.
      procedure, pass(self) :: load_file_name_from_file !< Load file name from file.
      procedure, pass(self) :: open                     !< Open file.
      procedure, pass(self) :: save_file_name_into_file !< Save file name into file.
      ! operators
      generic :: assignment(=) => file_assign_file !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: file_assign_file !< Operator `=`.
endtype file_object

contains
   ! public methods
   subroutine close(self)
   !< Close file.
   class(file_object), intent(inout) :: self !< File object.

   if (self%is_initialized) then
      if (self%is_connected) then
         close(unit=self%file_unit)
         self%file_unit = 0
         self%is_connected = .false.
      else
         write(stderr, '(A)') 'error: file "'//self%file_name//'" is not connected, thus its unit cannot be closed'
         self%error%status = ERROR_NOT_CONNECTED
      endif
   else
      write(stderr, '(A)') 'error: file is not initialized, thus its unit cannot be closed'
      self%error%status = ERROR_NOT_INITIALIZED
   endif
   endsubroutine close

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the file.
   class(file_object), intent(in)           :: self             !< Files collection.
   character(*),       intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable            :: desc             !< Description.
   character(len=:), allocatable            :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter              :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   if (allocated(self%file_name)) desc = desc//prefix_//'File name: '//self%file_name
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy file.
   class(file_object), intent(inout) :: self  !< File object.
   type(file_object)                 :: fresh !< Fresh instance of file object.

   self = fresh
   if (allocated(self%file_name)) deallocate(self%file_name)
   endsubroutine destroy

   elemental subroutine initialize(self, file_name)
   !< Initialize File.
   !<
   !< @note Leading and trailing white spaces are removed from file name.
   class(file_object), intent(inout)        :: self      !< File object.
   character(len=*),   intent(in), optional :: file_name !< File name.

   call self%destroy
   call self%error%initialize
   if (present(file_name)) then
      self%file_name = trim(adjustl(file_name))
   else
      self%file_name = UNSET_FILE_NAME
   endif
   self%is_initialized = .true.
   endsubroutine initialize

   subroutine load_file_name_from_file(self, fini, section_name, option_name, go_on_fail)
   !< Load file name from file.
   class(file_object), intent(inout)        :: self         !< File object.
   type(file_ini),     intent(in)           :: fini         !< Simulation parameters ini file handler.
   character(*),       intent(in)           :: section_name !< Option name into the ini file.
   character(*),       intent(in)           :: option_name  !< Option name into the ini file.
   logical,            intent(in), optional :: go_on_fail   !< Go on if load fails.
   character(999)                           :: buffer       !< Buffer string.

   call fini%get(section_name=section_name, option_name=option_name, val=buffer, error=self%error%status)
   if (present(go_on_fail)) then
      if (.not.go_on_fail) &
         call self%error%check(message='failed to load ['//section_name//'].('//option_name//')', is_severe=.not.go_on_fail)
   endif
   if (self%error%status <= 0) self%file_name = trim(adjustl(buffer))
   endsubroutine load_file_name_from_file

   subroutine open(self, file_name, format, action, access)
   !< Open file.
   class(file_object), intent(inout)        :: self       !< File object.
   character(len=*),   intent(in), optional :: file_name  !< File name.
   character(len=*),   intent(in), optional :: format     !< File format.
   character(len=*),   intent(in), optional :: action     !< File action.
   character(len=*),   intent(in), optional :: access     !< File access.
   character(len=:), allocatable            :: file_name_ !< File name, local variable.
   character(len=:), allocatable            :: format_    !< File format, local variable.
   character(len=:), allocatable            :: action_    !< File action, local variable.
   character(len=:), allocatable            :: access_    !< File access, local variable.

   format_ = 'unformatted' ; if (present(format)) format_ = format
   action_ = 'readwrite'   ; if (present(action)) action_ = action
   access_ = 'stream'      ; if (present(access)) access_ = access

   if (self%is_initialized) then
      file_name_ = self%file_name ; if (present(file_name)) file_name_ = trim(adjustl(file_name))
      if (.not.self%is_connected) then
         open(newunit=self%file_unit, file=file_name_, form=format_, action=action_, access=access_)
         self%is_connected = .true.
      else
         write(stderr, '(A)') 'error: file "'//self%file_name//'" is already connected, thus its unit cannot be re-open'
         self%error%status = ERROR_ALREADY_CONNECTED
      endif
   else
      write(stderr, '(A)') 'error: file is not initialized, thus it cannot be open'
      self%error%status = ERROR_NOT_INITIALIZED
   endif
   endsubroutine open

   subroutine save_file_name_into_file(self, fini, section_name, option_name)
   !< Save file name into file.
   class(file_object), intent(inout)        :: self         !< File object.
   type(file_ini),     intent(inout)        :: fini         !< Simulation parameters ini file handler.
   character(*),       intent(in)           :: section_name !< Option name into the ini file.
   character(*),       intent(in)           :: option_name  !< Option name into the ini file.

   call fini%add(section_name=section_name, option_name=option_name, val=self%file_name, error=self%error%status)
   endsubroutine save_file_name_into_file

   ! private methods
   pure subroutine file_assign_file(lhs, rhs)
   !< Operator `=`.
   class(file_object), intent(inout) :: lhs !< Left hand side.
   type(file_object),  intent(in)    :: rhs !< Right hand side.

                                 lhs%error          = rhs%error
   if (allocated(rhs%file_name)) lhs%file_name      = rhs%file_name
                                 lhs%file_unit      = rhs%file_unit
                                 lhs%is_initialized = rhs%is_initialized
                                 lhs%is_connected   = rhs%is_connected
   endsubroutine file_assign_file
endmodule off_file_object
