!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_FaceDerivedType Data_Type_Face
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_FacePublicProcedure Data_Type_Face
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_FacePrivateProcedure Data_Type_Face
!> @}

!> This module contains the definition of Type_Face and its procedures.
!> @todo \b DocWriteRead: Complete the documentation of write and read functions
module Data_Type_Face
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision     ! Integers and reals precision definition.
USE Data_Type_BC     !< Definition of Type_BC.
USE Data_Type_Vector !< Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing (inter)face-level data.
!> (Inter)face-level type contains data (mesh and boundary conditions) of interfaces numerical grid.
!> @ingroup Data_Type_FaceDerivedType
type, public:: Type_Face
  type(Type_Vector):: N          !< Face normal versor.
  real(R8P)::         S = 0._R8P !< Face area.
  type(Type_BC)::     BC         !< Boundary conditions.
  contains
    procedure, non_overridable:: init => init_face ! Procedure for initilizing allocatable variables.
    procedure, non_overridable:: free => free_face ! Procedure for freeing the memory of allocatable variables.
endtype Type_Face
!> @brief Pointer of Type_Face for creating array of pointers of Type_SFace.
!> @ingroup Data_Type_FaceDerivedType
type, public:: Type_Face_Ptr
  type(Type_Face),  pointer:: p => null()
endtype Type_Face_Ptr
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_FacePrivateProcedure
  !> @{
  !> Subroutine for freeing dynamic data of Type_Face variables.
  elemental subroutine free_face(face)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Face), intent(INOUT):: face !< face data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call face%BC%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_face

  !> Subroutine for initializing dynamic data of Type_Face variables.
  elemental subroutine init_face(face,bc0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Face), intent(INOUT)::        face !< Cell data.
  type(Type_BC),    intent(IN), optional:: bc0  !< Boundary conditions inizialization data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(bc0)) then
    call face%BC%init(bc0=bc0)
  else
    call face%BC%init
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_face
  !> @}
endmodule Data_Type_Face
