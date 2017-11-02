!< OFF free conditions object definition and implementation.

module off_free_conditions_object
!< OFF free conditions object definition and implementation.

use off_error_object, only : error_object
use finer, only : file_ini
use penf, only : I4P, R8P, str
use vecfor, only : ex, ey, ez, vector

implicit none
private
public :: free_conditions_object

character(len=15), parameter :: INI_SECTION_NAME='free_conditions' !< INI (config) file section name containing the
                                                                   !< free stream conditions.

type :: free_conditions_object
   !< Free conditions object class.
   !<
   !< Define the conditions of the **free stream**.
   type(error_object) :: error    !< Errors handler.
   type(vector)       :: velocity !< Velocity.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return a pretty-formatted description of the free conditions.
      procedure, pass(self) :: destroy        !< Destroy free conditions.
      procedure, pass(self) :: initialize     !< Initialize free conditions.
      procedure, pass(self) :: load_from_file !< Load from file.
      procedure, pass(self) :: save_into_file !< Save into file.
      ! operators
      generic :: assignment(=) => free_assign_free !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: free_assign_free !< Operator `=`.
endtype free_conditions_object

contains
   ! public methods
   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the free conditions.
   class(free_conditions_object), intent(in)           :: self             !< Free conditions.
   character(*),                  intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                       :: desc             !< Description.
   character(len=:), allocatable                       :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter                         :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//prefix_//'velocity : '//trim(str([self%velocity%x, self%velocity%y, self%velocity%z]))
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy free conditions.
   class(free_conditions_object), intent(inout) :: self  !< Free conditions.
   type(free_conditions_object)                 :: fresh !< Fresh free conditions.

   self = fresh
   endsubroutine destroy

   elemental subroutine initialize(self, free_conditions)
   !< Initialize free conditions.
   class(free_conditions_object), intent(inout)        :: self            !< Free conditions.
   type(free_conditions_object),  intent(in), optional :: free_conditions !< Free conditions values.

   call self%destroy
   if (present(free_conditions)) self = free_conditions
   endsubroutine initialize

   subroutine load_from_file(self, fini, go_on_fail)
   !< Load from file.
   class(free_conditions_object), intent(inout)        :: self          !< Free conditions.
   type(file_ini),                intent(in)           :: fini          !< Simulation parameters ini file handler.
   logical,                       intent(in), optional :: go_on_fail    !< Go on if load fails.
   real(R8P)                                           :: velocity(1:3) !< Velocity loading buffer.

   call fini%get(section_name=INI_SECTION_NAME, option_name='velocity', val=velocity, error=self%error%status)
   if (present(go_on_fail)) then
      if (.not.go_on_fail) &
         call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(velocity)', is_severe=.not.go_on_fail)
   endif
   if (self%error%status <= 0) self%velocity = ex * velocity(1) + ey * velocity(2) + ez * velocity(3)
   endsubroutine load_from_file

   subroutine save_into_file(self, fini)
   !< Save into file.
   class(free_conditions_object), intent(inout) :: self !< Free conditions.
   type(file_ini),                intent(inout) :: fini !< Simulation parameters ini file handler.

   call fini%add(section_name=INI_SECTION_NAME, option_name='velocity', val=[self%velocity%x, self%velocity%y, self%velocity%z], &
                 error=self%error%status)
   endsubroutine save_into_file

   ! private methods
   pure subroutine free_assign_free(lhs, rhs)
   !< Operator `=`.
   class(free_conditions_object), intent(inout) :: lhs !< Left hand side.
   type(free_conditions_object),  intent(in)    :: rhs !< Right hand side.

   lhs%error = rhs%error
   lhs%velocity = rhs%velocity
   endsubroutine free_assign_free
endmodule off_free_conditions_object
