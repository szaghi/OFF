!< OFF cell object definition and implementation.

module off_cell_object
!< OFF cell object definition and implementation.

use off_level_set_object, only : level_set_object
use penf, only : I4P, R8P
use vecfor, only : vector

implicit none
private
public :: cell_object

type :: cell_object
   !< Cell object class.
   type(vector)           :: center        !< Cell center.
   real(R8P)              :: volume=0._R8P !< Cell volume.
   type(level_set_object) :: level_set     !< Level set cell data.
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

   pure subroutine initialize(self, interfaces_number, distances)
   !< Initialize cell.
   class(cell_object), intent(inout)        :: self              !< Cell object.
   integer(I4P),       intent(in), optional :: interfaces_number !< Number of different interfaces.
   real(R8P),          intent(in), optional :: distances(:)      !< Distance from all interfaces.

   call self%destroy
   call self%level_set%initialize(interfaces_number=interfaces_number, distances=distances)
   endsubroutine initialize

   ! private methods
   pure subroutine cell_assign_cell(lhs, rhs)
   !< Operator `=`.
   class(cell_object), intent(inout) :: lhs !< Left hand side.
   type(cell_object),  intent(in)    :: rhs !< Right hand side.

   lhs%center = rhs%center
   lhs%volume = rhs%volume
   lhs%level_set = rhs%level_set
   endsubroutine cell_assign_cell
endmodule off_cell_object
