#include "preprocessor_macros.h"
!< OFF cell object definition and implementation.

module off_cell_object
!< OFF cell object definition and implementation.

use off_bc_object, only : bc_object
use off_level_set_object, only : level_set_object
use flow, only : conservative_compressible, primitive_compressible
use penf, only : I4P, R8P
use vecfor, only : vector

implicit none
private
public :: cell_object

type :: cell_object
   !< Cell object class.
   type(bc_object), allocatable    :: bc            !< Cell boundary conditions (if not allocated it is a free cell).
   type(vector)                    :: center        !< Cell center.
   real(R8P)                       :: volume=0._R8P !< Cell volume.
   real(R8P)                       :: Dt=0._R8P     !< Local time step.
   type(level_set_object)          :: level_set     !< Level set cell data.
   type(primitive_compressible)    :: P             !< Primitive variables.
   type(conservative_compressible) :: U             !< Conservative variables, Integrand (state) variables.
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
   class(cell_object), intent(inout) :: self  !< Cell object.
   type(cell_object)                 :: fresh !< Fresh instance of cell object.

   self = fresh
   endsubroutine destroy

   _PURE_ subroutine initialize(self, bc, Dt, interfaces_number, distances, P, U)
   !< Initialize cell.
   class(cell_object),              intent(inout)        :: self              !< Cell object.
   type(bc_object),                 intent(in), optional :: bc                !< Boundary conditions.
   real(R8P),                       intent(in), optional :: Dt                !< Local time step.
   integer(I4P),                    intent(in), optional :: interfaces_number !< Number of different interfaces.
   real(R8P),                       intent(in), optional :: distances(:)      !< Distance from all interfaces.
   type(primitive_compressible),    intent(in), optional :: P                 !< Primitive variables.
   type(conservative_compressible), intent(in), optional :: U                 !< Conservative variables.

   call self%destroy
   if (present(bc)) then
      if (.not.allocated(self%bc)) allocate(self%bc)
      self%bc = bc
   endif
   if (present(Dt)) self%Dt = Dt
   call self%level_set%initialize(interfaces_number=interfaces_number, distances=distances)
   if (present(P)) self%P = P
   if (present(U)) self%U = U
   endsubroutine initialize

   ! private methods
   pure subroutine cell_assign_cell(lhs, rhs)
   !< Operator `=`.
   class(cell_object), intent(inout) :: lhs !< Left hand side.
   type(cell_object),  intent(in)    :: rhs !< Right hand side.

   if (allocated(rhs%bc)) then
      if (.not.allocated(lhs%bc)) allocate(lhs%bc)
      lhs%bc = rhs%bc
   else
      if (allocated(lhs%bc)) deallocate(lhs%bc)
   endif
   lhs%center = rhs%center
   lhs%volume = rhs%volume
   lhs%Dt = rhs%Dt
   lhs%level_set = rhs%level_set
   lhs%P = rhs%P
   lhs%U = rhs%U
   endsubroutine cell_assign_cell
endmodule off_cell_object
