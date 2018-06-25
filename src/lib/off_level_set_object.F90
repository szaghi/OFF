#include "preprocessor_macros.h"
!< OFF level set object definition and implementation.

module off_level_set_object
!< OFF level set object definition and implementation.

use penf, only : I4P, R8P
use vecfor, only : vector

implicit none
private
public :: level_set_object
public :: level_set_normal

type :: level_set_object
   !< Level set object class.
   !<
   !< Level set technique is used to track interfaces by means of computed distances with sign.
   real(R8P)              :: distance=1e10_R8P !< Distance from the closest interface.
   real(R8P), allocatable :: distances(:)      !< Distance from all interfaces.
   type(vector)           :: normal            !< Normal of level set function.
   contains
      ! public methods
      procedure, pass(self) :: compute_normal  !< Compute normal.
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
   pure subroutine compute_normal(self, cell, dx, dy, dz)
   !< Compute normal.
   !<
   !< @note Level set values must be given as a cross-cells stencil:
   !< + cell(1) => i min
   !< + cell(2) => i max
   !< + cell(3) => j min
   !< + cell(4) => j max
   !< + cell(5) => k min
   !< + cell(6) => k max
   class(level_set_object), intent(inout) :: self       !< Level set object.
   type(level_set_object),  intent(in)    :: cell(1:6)  !< Level set in the 6 cross-cells.
   real(R8P),               intent(in)    :: dx, dy, dz !< Space steps.

   self%normal = level_set_normal(cell=cell, dx=dx, dy=dy, dz=dz)
   endsubroutine compute_normal

   elemental subroutine destroy(self)
   !< Destroy object.
   class(level_set_object), intent(inout) :: self  !< Level set object.
   type(level_set_object)                 :: fresh !< Fresh instance of level set object.

   self = fresh
   endsubroutine destroy

   _PURE_ subroutine initialize(self, interfaces_number, distances)
   !< Initialize object.
   class(level_set_object), intent(inout)        :: self              !< Level set object.
   integer(I4P),            intent(in), optional :: interfaces_number !< Number of different interfaces.
   real(R8P),               intent(in), optional :: distances(:)      !< Distance from all interfaces.

   call self%destroy
   if (present(interfaces_number)) then
      if (interfaces_number>0) then
         allocate(self%distances(1:interfaces_number))
         self%distances = 1e10_R8P
      endif
   endif
   if (present(distances)) self%distances = distances
   call self%update_distance
   endsubroutine initialize

   elemental subroutine update_distance(self)
   !< Update distance.
   class(level_set_object), intent(inout) :: self !< Level set object.
   integer(I4P)                           :: d    !< Counter.

   if (allocated(self%distances)) then
      if (any(self%distances<0._R8P)) then
         self%distance = -1e10_R8P
         do d=1, size(self%distances, dim=1)
            if (self%distances(d)<0._R8P) then
               self%distance = max(self%distance, self%distances(d))
            endif
         enddo
      else
         self%distance = minval(self%distances)
      endif
   endif
   endsubroutine update_distance

   ! non TBP
   pure function level_set_normal(cell, dx, dy, dz) result(normal)
   !< Compute normal of level set function.
   !<
   !< @note Level set values must be given as a cross-cells stencil:
   !< + cell(1) => i min
   !< + cell(2) => i max
   !< + cell(3) => j min
   !< + cell(4) => j max
   !< + cell(5) => k min
   !< + cell(6) => k max
   type(level_set_object), intent(in) :: cell(1:6)  !< Level set in the 6 cross-cells.
   real(R8P),              intent(in) :: dx, dy, dz !< Space steps.
   type(vector)                       :: normal     !< Level set function normal.

   normal%x = (cell(2)%distance - cell(1)%distance) / dx
   normal%y = (cell(4)%distance - cell(3)%distance) / dy
   normal%z = (cell(6)%distance - cell(5)%distance) / dz
   call normal%normalize
   endfunction level_set_normal

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
   lhs%normal = rhs%normal
   endsubroutine level_set_assign_level_set
endmodule off_level_set_object
