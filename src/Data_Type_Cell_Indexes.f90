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
!> Type_Cell_Indexes is a derived type containing cell indexes.
module Data_Type_Cell_Indexes
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing cell indexes
!> @ingroup Data_Type_Cell_IndexesDerivedType
type, public:: Type_Cell_Indexes
  integer(I8P):: ID = 0_I8P !< Block ID index.
  integer(I4P):: i  = 0_I4P !< i index in the block arrays.
  integer(I4P):: j  = 0_I4P !< j index in the block arrays.
  integer(I4P):: k  = 0_I4P !< k index in the block arrays.
  contains
    ! operators overloading
    generic:: assignment(=) => assign_cell_indexes
    ! private procedures
    procedure, pass(self1), private:: assign_cell_indexes
endtype Type_Cell_Indexes
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Cell_IndexesPrivateProcedure
  !> @{
  ! Assignment (=)
  !> @brief Procedure for assignment between two cell indexes variables.
  elemental subroutine assign_cell_indexes(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell_Indexes), intent(INOUT):: self1 !< LHS.
  type(Type_Cell_Indexes),  intent(IN)::    self2 !< RHS.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self1%ID = self2%ID
  self1%i  = self2%i
  self1%j  = self2%j
  self1%k  = self2%k
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_cell_indexes
  !> @}
endmodule Data_Type_Cell_Indexes
