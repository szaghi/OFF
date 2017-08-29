!< OFF error object definition and implementation.

module off_error_object
!< OFF error object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I4P, str

implicit none
private
public :: error_object

type :: error_object
   !< Error object class.
   !<
   !< Handler of errors/excetpions.
   integer(I4P)                  :: status=0_I4P !< Error status.
   character(len=:), allocatable :: message      !< Error message.
   contains
      ! public methods
      procedure, pass(self) :: check      !< Check error status.
      procedure, pass(self) :: destroy    !< Destroy error.
      procedure, pass(self) :: initialize !< Initialize error.
      ! operators
      generic :: assignment(=) => err_assign_err !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: err_assign_err !< Operator `=`.
endtype error_object

contains
   ! public methods
   subroutine check(self, message, is_severe)
   !< Check error status.
   !<
   !< If errors occurred, a warn is printed.
   class(error_object), intent(inout)        :: self       !< Error object.
   character(*),        intent(in), optional :: message    !< Error message.
   logical,             intent(in), optional :: is_severe  !< Enable severe error.
   logical                                   :: is_severe_ !< Enable severe error, local variable.

   if (self%status /= 0) then
      is_severe_ = .false. ; if (present(is_severe)) is_severe_ = is_severe
      if (present(message)) then
         self%message = 'error ['//trim(str(self%status))//']: '//message
      else
         self%message = 'error ['//trim(str(self%status))//']'
      endif
      write(stderr, '(A)') self%message
      if (is_severe_) stop
   endif
   endsubroutine check

   elemental subroutine destroy(self)
   !< Destroy error.
   class(error_object), intent(inout) :: self  !< Error object.
   type(error_object)                 :: fresh !< Fresh instance of error object.

   self = fresh
   if (allocated(self%message)) deallocate(self%message)
   endsubroutine destroy

   elemental subroutine initialize(self)
   !< Initialize error.
   class(error_object), intent(inout) :: self !< Error object.

   call self%destroy
   endsubroutine initialize

   ! private methods
   pure subroutine err_assign_err(lhs, rhs)
   !< Operator `=`.
   class(error_object), intent(inout) :: lhs !< Left hand side.
   type(error_object),  intent(in)    :: rhs !< Right hand side.

   lhs%status = rhs%status
   if (allocated(rhs%message)) lhs%message = rhs%message
   endsubroutine err_assign_err
endmodule off_error_object
