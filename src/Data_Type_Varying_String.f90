!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Varying_StringDerivedType Data_Type_Varying_String
!> Module definition of Type_Varying_String
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Varying_StringPrivateProcedure Data_Type_Varying_String
!> Module definition of Type_Varying_String
!> @}

!> @brief Module Data_Type_Varying_String contains the definition of Type_Varying_String and useful procedures for its handling.
module Data_Type_Varying_String
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type defining Type_Varying_String, a useful type for implementing ISO_VARYING_STRING by means of fortran 2003
!> features.
!> @ingroup Data_Type_Varying_StringDerivedType
type, public:: Type_Varying_String
  character(len=:), allocatable:: vs !< Deferred (vaying) length string.
  contains
    procedure:: free     ! Procedure for freeing dynamic memory.
    final::     finalize ! Procedure for freeing dynamic memory when finalizing.
endtype Type_Varying_String
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_XML_TagPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free(vstring)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Varying_String), intent(INOUT):: vstring !< Varying (lenght) string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(vstring%vs)) deallocate(vstring%vs)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(vstring)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Varying_String), intent(INOUT):: vstring !< XML tag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call vstring%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
  !> @}
endmodule Data_Type_Varying_String
