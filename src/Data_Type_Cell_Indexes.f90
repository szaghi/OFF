!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Cell_IndexesDerivedType Data_Type_Cell_Indexes
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_Cell_IndexesInterface Data_Type_Cell_Indexes
!> Module definition of Type_Cell_Indexes
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Cell_IndexesPrivateProcedure Data_Type_Cell_Indexes
!> Module definition of Type_Cell_Indexes
!> @}

!> @brief This module contains the definition of Type_Cell_Indexes and its procedures.
!> Type_Cell_Indexes is a derived type containing adjacent boundary condition informations.
module Data_Type_Cell_Indexes
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing adjacent boundary condition.
!> @ingroup Data_Type_Cell_IndexesDerivedType
type, public:: Type_Cell_Indexes
  !integer(I4P):: l = 0_I4P !< l index of adjacent block.
  integer(I4P):: b = 0_I4P !< b index of adjacent block.
  integer(I4P):: i = 0_I4P !< i index of adjacent cell in the (b,l) block.
  integer(I4P):: j = 0_I4P !< j index of adjacent cell in the (b,l) block.
  integer(I4P):: k = 0_I4P !< k index of adjacent cell in the (b,l) block.
  contains
    ! operators overloading
    generic:: assignment(=) => assign_bc_adjacent
    ! private procedures
    procedure, pass(bc1), private:: assign_bc_adjacent
endtype Type_Cell_Indexes
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Cell_IndexesPrivateProcedure
  !> @{
  ! Assignment (=)
  !> @brief Procedure for assignment between two adjacent boundary conditions variables.
  elemental subroutine assign_bc_adjacent(bc1,bc2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell_Indexes), intent(INOUT):: bc1
  type(Type_Cell_Indexes),  intent(IN)::    bc2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !bc1%l = bc2%l
  bc1%b = bc2%b
  bc1%i = bc2%i
  bc1%j = bc2%j
  bc1%k = bc2%k
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_bc_adjacent
  !> @}
endmodule Data_Type_Cell_Indexes
