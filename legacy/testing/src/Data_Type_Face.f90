!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_FaceDerivedType Data_Type_Face
!> Module definition of Type_Face
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_FacePrivateProcedure Data_Type_Face
!> Module definition of Type_Face
!> @}

!> This module contains the definition of Type_Face and its procedures.
!> @todo \b DocWriteRead: Complete the documentation of write and read functions
module Data_Type_Face
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                        ! Integers and reals precision definition.
USE Data_Type_BC,     only: Type_BC     ! Definition of Type_BC.
USE Data_Type_Vector, only: Type_Vector ! Definition of Type_Vector.
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
    procedure:: free  => free_face  ! Procedure for freeing dynamic memory.
    procedure:: alloc => alloc_face ! Procedure for allocating dynamic memory.
    procedure:: load  => load_face  ! Procedure for loading Type_Face data.
    procedure:: save  => save_face  ! Procedure for saving Type_Face data.
    final::     finalize            ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
    generic:: assignment(=) => assign_face
    ! private procedures
    procedure, pass(face1), private:: assign_face
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
  !> @brief Procedure for freeing dynamic data of Type_Face variables.
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

  !> @brief Subroutine for freeing dynamic memory when finalizing.
  elemental subroutine finalize(face)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Face), intent(INOUT):: face !< Cell data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call face%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic data of Type_Face variables.
  elemental subroutine alloc_face(face)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Face), intent(INOUT):: face !< Cell data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call face%BC%alloc
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_face

  !> @brief Procedure for loading Type_Face data.
  subroutine load_face(face,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Face),       intent(INOUT):: face    !< Face data.
  integer(I8P), optional, intent(IN)::    pos     !< Position specifier.
  integer(I4P), optional, intent(OUT)::   iostat  !< IO error.
  character(*), optional, intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),           intent(IN)::    unit    !< Logic unit.
  integer(I4P)::                          iostatd !< IO error.
  character(500)::                        iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    call face%N%load( unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
    read(             unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)face%S
    call face%BC%load(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
  else
    call face%N%load( unit=unit,        iostat=iostatd,iomsg=iomsgd)
    read(             unit=unit,        iostat=iostatd,iomsg=iomsgd)face%S
    call face%BC%load(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_face

  !> @brief Procedure for saving Type_Face data.
  subroutine save_face(face,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Face),       intent(IN)::  face    !< Face data.
  integer(I8P), optional, intent(IN)::  pos     !< Position specifier.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    call face%N%save( unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
    write(            unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)face%S
    call face%BC%save(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
  else
    call face%N%save( unit=unit,        iostat=iostatd,iomsg=iomsgd)
    write(            unit=unit,        iostat=iostatd,iomsg=iomsgd)face%S
    call face%BC%save(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_face

  ! Assignment (=)
  !> @brief Procedure for assignment between two faces variables.
  elemental subroutine assign_face(face1,face2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Face), intent(INOUT):: face1
  type(Type_Face),  intent(IN)::    face2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  face1%N  = face2%N
  face1%S  = face2%S
  face1%BC = face2%BC
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_face
  !> @}
endmodule Data_Type_Face
