!< OFF level set object definition and implementation.

module off_level_set_object
!< OFF level set object definition and implementation.

use penf, only : I4P, R8P

implicit none
private
public :: level_set_object

type :: level_set_object
   !< Level set object class.
   !<
   !< Level set technique is used to track interfaces by means of computed distances with sign.
   real(R8P)              :: distance=0._R8P !< Distance from the closest interface.
   real(R8P), allocatable :: distances(:)    !< Distance from all interfaces.
   contains
      ! public methods
      procedure, pass(self) :: destroy         !< Destroy cell.
      procedure, pass(self) :: initialize      !< Initialize cell.
      procedure, pass(self) :: update_distance !< Update distance.
      ! operators
      generic :: assignment(=) => level_set_assign_level_set !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: level_set_assign_level_set !< Operator `=`.
endtype level_set_object

contains
   ! public methods
   elemental subroutine destroy(self)
   !< Destroy object.
   class(level_set_object), intent(inout) :: self  !< Level set object.
   type(level_set_object)                 :: fresh !< Fresh instance of level set object.

   self = fresh
   endsubroutine destroy

   pure subroutine initialize(self, interfaces_number, distances)
   !< Initialize object.
   class(level_set_object), intent(inout)        :: self              !< Level set object.
   integer(I4P),            intent(in), optional :: interfaces_number !< Number of different interfaces.
   real(R8P),               intent(in), optional :: distances(:)      !< Distance from all interfaces.

   call self%destroy
   if (present(interfaces_number)) then
      allocate(self%distances(1:interfaces_number))
      self%distances = 0._R8P
   endif
   if (present(distances)) self%distances = distances
   call self%update_distance
   endsubroutine initialize

   elemental subroutine update_distance(self)
   !< Update distance.
   class(level_set_object), intent(inout) :: self !< Level set object.

   if (allocated(self%distances)) then
      self%distance = minval(self%distances)
   endif
   endsubroutine update_distance

   ! private methods
   pure subroutine level_set_assign_level_set(lhs, rhs)
   !< Operator `=`.
   class(level_set_object), intent(inout) :: lhs !< Left hand side.
   type(level_set_object),  intent(in)    :: rhs !< Right hand side.

   lhs%distance = rhs%distance
   if (allocated(rhs%distances)) then
      lhs%distances = rhs%distances
   else
      if (allocated(lhs%distances)) deallocate(lhs%distances)
   endif
   endsubroutine level_set_assign_level_set
endmodule off_level_set_object
