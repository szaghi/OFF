!< OFF cell object definition and implementation.

module off_cell_object
!< OFF cell object definition and implementation.

use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: cell_object

type :: cell_object
   !< Cell object class.
   type(vector) :: center        !< Cell center.
   real(R8P)    :: volume=0._R8P !< Cell volume.
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy cell.
      procedure, pass(self) :: initialize !< Initialize cell.
      ! operators
      generic :: assignment(=) => cell_assign_cell !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: cell_assign_cell !< Operator `=`.
endtype cell_object

contains
   ! public methods
   elemental subroutine destroy(self)
   !< Destroy cell.
   class(cell_object), intent(inout) :: self !< Cell object.
   type(cell_object)                 :: fresh !< Fresh instance of cell object.

   self = fresh
   endsubroutine destroy

   elemental subroutine initialize(self)
   !< Initialize cell.
   class(cell_object), intent(inout) :: self !< Cell object.

   call self%destroy
   endsubroutine initialize

   ! private methods
   pure subroutine cell_assign_cell(lhs, rhs)
   !< Operator `=`.
   class(cell_object), intent(inout) :: lhs !< Left hand side.
   type(cell_object),  intent(in)    :: rhs !< Right hand side.

   lhs%center = rhs%center
   lhs%volume = rhs%volume
   endsubroutine cell_assign_cell
endmodule off_cell_object
