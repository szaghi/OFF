!< OFF boundary conditions object definition and implementation.

module off_bc_object
!< OFF boundary conditions object definition and implementation.

use penf, only : I1P
use stringifor, only : string

implicit none
private
public :: bc_object
public :: BC_FREE
public :: BC_WALL
public :: BC_PERIODIC
public :: BC_EXTRAPOLATED

integer(I1P), parameter :: BC_FREE         = 0_I1P !< Cell free boundary conditions.
integer(I1P), parameter :: BC_WALL         = 1_I1P !< Solid (inviscid) wall boundary conditions.
integer(I1P), parameter :: BC_PERIODIC     = 2_I1P !< Periodic (circular) boundary conditions, special case of adjacency.
integer(I1P), parameter :: BC_EXTRAPOLATED = 3_I1P !< Extrapolated boundary conditions.

type :: bc_object
   !< Boundary conditions object class.
   integer(I1P) :: id = BC_FREE !< Boundary conditions id.
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy bc.
      procedure, pass(self) :: initialize !< Initialize bc.
      ! operators
      generic :: assignment(=) => bc_assign_bc, bc_assign_I1P !< Overload `=`.
      ! private methods
      procedure, pass(lhs)  :: bc_assign_bc  !< Operator `=`.
      procedure, pass(lhs)  :: bc_assign_I1P !< Operator `bc = integer(I1P)`.
      procedure, pass(self) :: set_from_code !< Set boundary conditions from character code.
endtype bc_object

contains
   ! public methods
   elemental subroutine destroy(self)
   !< Destroy bc.
   class(bc_object), intent(inout) :: self  !< BC object.
   type(bc_object)                 :: fresh !< Fresh instance of bc object.

   self = fresh
   endsubroutine destroy

   pure subroutine initialize(self, id, code, bc)
   !< Initialize bc.
   class(bc_object), intent(inout)        :: self !< BC object.
   integer(I1P),     intent(in), optional :: id   !< BC id.
   character(*),     intent(in), optional :: code !< Character code.
   type(bc_object),  intent(in), optional :: bc   !< BC object.

   call self%destroy
   if (present(id))   self%id = id
   if (present(code)) call self%set_from_code(code=code)
   if (present(bc))   self%id = bc%id
   endsubroutine initialize

   ! private methods
   ! `=` operator
   pure subroutine bc_assign_bc(lhs, rhs)
   !< Operator `=`.
   class(bc_object), intent(inout) :: lhs !< Left hand side.
   type(bc_object),  intent(in)    :: rhs !< Right hand side.

   lhs%id = rhs%id
   endsubroutine bc_assign_bc

   pure subroutine bc_assign_I1P(lhs, rhs)
   !< Operator `bc = integer(I1P)`.
   class(bc_object), intent(inout) :: lhs !< Left hand side.
   integer(I1P),     intent(in)    :: rhs !< Right hand side.

   lhs%id = rhs
   endsubroutine bc_assign_I1P

   ! auxiliary methods
   pure subroutine set_from_code(self, code)
   !< Set boundary conditions from character code.
   class(bc_object), intent(inout) :: self     !< BC object.
   character(*),     intent(in)    :: code     !< Character code.
   type(string)                    :: str_code !< String code.

   str_code = trim(adjustl(code))
   str_code = str_code%upper()
   select case(str_code%chars())
   case('FREE')
      self%id = BC_FREE
   case('WALL')
      self%id = BC_WALL
   case('PERIODIC')
      self%id = BC_PERIODIC
   case('EXTRAPOLATED')
      self%id = BC_EXTRAPOLATED
   endselect
   endsubroutine set_from_code
endmodule off_bc_object

