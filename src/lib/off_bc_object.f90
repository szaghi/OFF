!< OFF boundary conditions object definition and implementation.

module off_bc_object
!< OFF boundary conditions object definition and implementation.

use flow, only : conservative_compressible
use penf, only : I1P, I4P
use stringifor, only : string

implicit none
private
public :: bc_object
public :: BC_FREE
public :: BC_WALL
public :: BC_PERIODIC
public :: BC_EXTRAPOLATED
public :: BC_ADJACENT
public :: BC_INLET_SUPERSONIC

integer(I1P), parameter :: BC_FREE             = 0_I1P !< Cell free boundary conditions.
integer(I1P), parameter :: BC_WALL             = 1_I1P !< Solid (inviscid) wall boundary conditions.
integer(I1P), parameter :: BC_PERIODIC         = 2_I1P !< Periodic (circular) boundary conditions, special case of adjacency.
integer(I1P), parameter :: BC_EXTRAPOLATED     = 3_I1P !< Extrapolated boundary conditions.
integer(I1P), parameter :: BC_ADJACENT         = 4_I1P !< Adjacent boundary conditions.
integer(I1P), parameter :: BC_INLET_SUPERSONIC = 5_I1P !< Supersonic inlet boundary conditions.

type :: bc_object
   !< Boundary conditions object class.
   integer(I1P)                                 :: id = BC_FREE !< Boundary conditions id.
   integer(I4P),                    allocatable :: adj(:)       !< Indexes of adjacent cell (b,i,j,k), if any.
   type(conservative_compressible), allocatable :: U            !< Conservatives variable, somehow imposed.
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy bc.
      procedure, pass(self) :: initialize !< Initialize bc.
      procedure, pass(self) :: is         !< Check is bc is of a given type.
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

   pure subroutine initialize(self, id, code, adj, U, bc)
   !< Initialize bc.
   class(bc_object),                intent(inout)        :: self   !< BC object.
   integer(I1P),                    intent(in), optional :: id     !< BC id.
   character(*),                    intent(in), optional :: code   !< Character code.
   integer(I4P),                    intent(in), optional :: adj(:) !< Indexes of adjactent cell.
   type(conservative_compressible), intent(in), optional :: U      !< Conservative variables, somehow imposed.
   type(bc_object),                 intent(in), optional :: bc     !< BC object.

   call self%destroy
   if (present(id))   self%id = id
   if (present(code)) call self%set_from_code(code=code, adj=adj)
   if (present(adj)) self%adj = adj
   if (present(U)) self%U = U
   if (present(bc))   self%id = bc%id
   endsubroutine initialize

   pure function is(self, id) result(is_of)
   !< Check if bc is of a given type.
   class(bc_object), intent(in) :: self  !< BC object.
   integer(I1P),     intent(in) :: id    !< BC id.
   logical                      :: is_of !< Result of the inquire.

   is_of = self%id == id
   endfunction is

   ! private methods
   ! `=` operator
   pure subroutine bc_assign_bc(lhs, rhs)
   !< Operator `=`.
   class(bc_object), intent(inout) :: lhs !< Left hand side.
   type(bc_object),  intent(in)    :: rhs !< Right hand side.

   lhs%id = rhs%id
   if (allocated(rhs%adj)) then
      lhs%adj = rhs%adj
   else
      if (allocated(lhs%adj)) deallocate(lhs%adj)
   endif
   if (allocated(rhs%U)) then
      lhs%U = rhs%U
   else
      if (allocated(lhs%U)) deallocate(lhs%U)
   endif
   endsubroutine bc_assign_bc

   pure subroutine bc_assign_I1P(lhs, rhs)
   !< Operator `bc = integer(I1P)`.
   class(bc_object), intent(inout) :: lhs !< Left hand side.
   integer(I1P),     intent(in)    :: rhs !< Right hand side.

   lhs%id = rhs
   endsubroutine bc_assign_I1P

   ! auxiliary methods
   pure subroutine set_from_code(self, code, adj)
   !< Set boundary conditions from character code.
   class(bc_object), intent(inout)        :: self     !< BC object.
   character(*),     intent(in)           :: code     !< Character code.
   integer(I4P),     intent(in), optional :: adj(:)   !< Indexes of adjactent cell.
   type(string)                           :: str_code !< String code.

   str_code = trim(adjustl(code))
   str_code = str_code%upper()
   select case(str_code%chars())
   case('FRE', 'FREE')
      self%id = BC_FREE
   case('WAL', 'WALL')
      self%id = BC_WALL
   case('PER', 'PERIODIC')
      self%id = BC_PERIODIC
   case('EXT', 'EXTRAPOLATED')
      self%id = BC_EXTRAPOLATED
   case('ADJ', 'ADJACENT')
      self%id = BC_ADJACENT
      if (.not.allocated(self%adj)) allocate(self%adj(4))
   case('INS', 'INLET_SUPERSONIC')
      self%id = BC_INLET_SUPERSONIC
      if (.not.allocated(self%U)) allocate(self%U)
   endselect
   if (present(adj)) self%adj = adj
   endsubroutine set_from_code
endmodule off_bc_object
