!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_ForestDerivedType Data_Type_Forest
!> Module definition of Type_Forest
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_ForestPublicProcedure Data_Type_Forest
!> Module definition of Type_Forest
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_ForestPrivateProcedure Data_Type_Forest
!> Module definition of Type_Forest
!> @}

!> @brief Module Data_Type_Forest contains the definition of Type_Forest type and useful procedures for its handling.
!> Type_Forest a set of trees for handling complex hierarchical data structure. All the features of Type_Tree are extended to a set
!> trees that can be interconnected.
!> As an example a forest of 2 trees can be represented (in 2D) as:
!> @code
!>
!>                     Face 2/1
!> Y                     |
!>/|\                   \|/
!> |
!> | +-------------------+-------------------+
!> | |         |         | 16 | 17 |         |
!> | |         |         |    |    |         |
!> | |   4     |   5     |----4----|   5     |
!> | |         |         | 14 | 15 |         |
!> | |         |         |    |    |         |
!> | |------ tree 1 -----|------ tree 2 -----|
!> | |  8 |  9 |         |         |         |
!> | |    |    |         |         |         |
!> | |----2----|   3     |    2    |   3     |
!> | |  6 |  7 |         |         |         |
!> | |    |    |         |         |         |
!> | +-------------------+-------------------+
!> |
!> |
!> o----------------------------------------------->X
!> @endcode
!> Note that in the above example the nodes of each tree have been indicated with the local (to the tree) ID. To define the
!> connection between different trees of a forest it is necessary to define the face/line along witch they are connected and their
!> relative orientation. In the above example the tree 1 and 2 are connected along the face 2 for tree 1 and face 1 for tree 1 and
!> their relative orientation is " x y z" meaning that the tree2 1 and 2 have x, y and z axis with the same orientation.
!> On the contrary the following example has inverted orientations:
!> @code
!>
!>                     Face 2/4
!> Y1                    |
!>/|\                   \|/
!> |                                           /|\ X2
!> | +-------------------+-------------------+  |
!> | |         |         | 21 | 19 |         |  |
!> | |         |         |    |    |         |  |
!> | |   4     |   5     |----5----|   3     |  |
!> | |         |         | 20 | 18 |         |  |
!> | |         |         |    |    |         |  |
!> | |------ tree 1 -----|------ tree 2 -----|  |
!> | |  8 |  9 |         |         |         |  |
!> | |    |    |         |         |         |  |
!> | |----2----|   3     |    4    |   2     |  |
!> | |  6 |  7 |         |         |         |  |
!> | |    |    |         |         |         |  |
!> | +-------------------+-------------------+  |
!> |                                            |
!> |                  Y2<-----------------------o
!> |
!> o----------------------------------------------->X1
!> @endcode
!> In the above example tree 1 has oX1Y1 reference system, while tree 2 has oX2Y2 system. The two trees are connected along the face
!> 2 for tree 1 and face 4 for tree 2. Their relative orientation is "-y x z" meaning that x axis of tree 1 is parallel to axis y of
!> tree 2 with inverted direction, y axis of tree 1 is parallel to x axis of tree 2 with the same direction, and z axis are the
!> same. Once the relative orientation of connected trees is defined it is possible to compute the siblings nodes of node connected
!> to the nodes of a different tree of the forest.
module Data_Type_Forest
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision   ! Integers and reals precision definition.
USE Data_Type_Tree ! Definition of Type_Tree.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type defining a Forest of trees. It is based on Type_Tree.
!> @ingroup Data_Type_ForestDerivedType
type, public:: Type_Forest
  integer(I4P)::                 Nt = 0_I4P !< Number of trees of the forest.
  type(Type_Tree), allocatable:: tree(:)    !< The Forest trees array [1:Nt].
  contains
    procedure:: free     => free_forest     ! Procedure for freeing (destroying) the forest.
    final::     finalize => finalize_forest ! Procedure for freeing dynamic memory when finalizing.
    procedure:: loopID   => loopID_forest   ! Procedure for performing a while loop returning the ID of trees.
    ! private procedures
    procedure:: assign_forest ! Procedure for assignment.
endtype Type_Forest
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_ForestPrivateProcedure
  !> @{
  !> @brief Procedure for freeing (destroying) the forest.
  elemental subroutine free_forest(forest)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Forest), intent(INOUT):: forest !< Forest.
  integer(I4P)::                      t      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(forest%tree)) then
    do t=lbound(forest%tree,dim=1),ubound(forest%tree,dim=1)
      call forest%tree(t)%free
    enddo
    deallocate(forest%tree)
    forest%Nt = 0_I4P
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_forest

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize_forest(forest)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Forest), intent(INOUT):: forest !< Forest.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call forest%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize_forest

  !> @brief Procedure for performing a while loop returning the ID value of all nodes of all trees of the foreset (for ID looping).
  function loopID_forest(forest,ID) result(again)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Forest), intent(IN)::  forest !< Forest.
  integer(I4P),       intent(OUT):: ID     !< Tree ID value.
  logical::                         again  !< Flag continuing the loop.
  integer(I4P), save::              t = 0  !< ID counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  again = .false.
  if (allocated(forest%tree)) then
    if (t==0) then
      t  = lbound(forest%tree,dim=1)
      ID =  t
      again = .true.
    elseif (t<ubound(forest%tree,dim=1)) then
      t  = t + 1
      ID = t
      again = .true.
    else
      t  = 0
      ID = 0
      again = .false.
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction loopID_forest

  ! Assignment (=)
  !> @brief Procedure for assignment between two forest variables.
  elemental subroutine assign_forest(forest1,forest2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Forest), intent(INOUT):: forest1 !< LHS.
  type(Type_Forest),  intent(IN)::    forest2 !< RHS.
  integer(I4P)::                      t       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call forest1%free
  if (allocated(forest2%tree)) then
    forest1%Nt = size(forest2%tree,dim=1)
    allocate(forest1%tree(lbound(forest2%tree,dim=1):ubound(forest2%tree,dim=1)))
    do t=lbound(forest2%tree,dim=1),ubound(forest2%tree,dim=1)
      forest1%tree(t) = forest2%tree(t)
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_forest
  !> @}
endmodule Data_Type_Forest
