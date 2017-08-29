!< OFF face object definition and implementation.

module off_face_object
!< OFF face object definition and implementation.

use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: face_object

type :: face_object
   !< Face object class.
   real(R8P)    :: area=0._R8P !< Area.
   type(vector) :: normal      !< Normal versor.
   contains
      ! public methods
      procedure, pass(self) :: compute_metrics !< Compute face metrics.
      procedure, pass(self) :: destroy         !< Destroy face.
      procedure, pass(self) :: initialize      !< Initialize face.
      ! operators
      generic :: assignment(=) => face_assign_face !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: face_assign_face !< Operator `=`.
endtype face_object

contains
   ! public methods
   elemental subroutine compute_metrics(self, pt1, pt2, pt3, pt4, signd)
   !< Compute face metrics.
   class(face_object), intent(inout) :: self  !< Face.
   type(vector),       intent(in)    :: pt1   !< Point 1 of face.
   type(vector),       intent(in)    :: pt2   !< Point 2 of face.
   type(vector),       intent(in)    :: pt3   !< Point 3 of face.
   type(vector),       intent(in)    :: pt4   !< Point 4 of face.
   real(R8P),          intent(in)    :: signd !< Sign of direction along normal coordinate.

   self%normal = self%normal%face_normal4(pt1 = pt1, pt2 = pt2, pt3 = pt3, pt4 = pt4) * signd
   self%area = self%normal%normL2()
   call self%normal%normalize
   endsubroutine compute_metrics

   elemental subroutine destroy(self)
   !< Destroy face.
   class(face_object), intent(inout) :: self  !< Face object.
   type(face_object)                 :: fresh !< Fresh instance of face object.

   self = fresh
   endsubroutine destroy

   pure subroutine initialize(self)
   !< Initialize face.
   class(face_object), intent(inout) :: self !< Face object.

   call self%destroy
   endsubroutine initialize

   ! private methods
   pure subroutine face_assign_face(lhs, rhs)
   !< Operator `=`.
   class(face_object), intent(inout) :: lhs !< Left hand side.
   type(face_object),  intent(in)    :: rhs !< Right hand side.

   lhs%area   = rhs%area
   lhs%normal = rhs%normal
   endsubroutine face_assign_face
endmodule off_face_object
