!< OFF node object definition and implementation.

module off_node_object
!< OFF node object definition and implementation.

use vecfor, only : vector

implicit none
private
public :: node_object

type :: node_object
   !< Node object class.
   type(vector) :: vertex !< Vertex coordinates.
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy node.
      procedure, pass(self) :: initialize !< Initialize node.
      ! operators
      generic :: assignment(=) => node_assign_node !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: node_assign_node !< Operator `=`.
endtype node_object

contains
   ! public methods
   elemental subroutine destroy(self)
   !< Destroy node.
   class(node_object), intent(inout) :: self  !< Node object.
   type(node_object)                 :: fresh !< Fresh instance of node object.

   self = fresh
   endsubroutine destroy

   pure subroutine initialize(self)
   !< Initialize node.
   class(node_object), intent(inout) :: self !< Node object.

   call self%destroy
   endsubroutine initialize

   ! private methods
   pure subroutine node_assign_node(lhs, rhs)
   !< Operator `=`.
   class(node_object), intent(inout) :: lhs !< Left hand side.
   type(node_object),  intent(in)    :: rhs !< Right hand side.

   lhs%vertex = rhs%vertex
   endsubroutine node_assign_node
endmodule off_node_object
