!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_ProbePublicProcedure Data_Type_Probe
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_ProbePrivateProcedure Data_Type_Probe
!> @}

!> This module contains the definition of Type_Probe and its procedures.
module Data_Type_Probe
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: write_probe,read_probe
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing probe informations.
!> @ingroup DerivedType
type, public:: Type_Probe
  integer(I_P):: b = 1_I_P !< Block (global map) index.
  integer(I_P):: i = 0_I_P !< I direction index.
  integer(I_P):: j = 0_I_P !< J direction index.
  integer(I_P):: k = 0_I_P !< K direction index.
  contains
    procedure, non_overridable:: set ! Procedure for setting members of Type_Probe.
endtype Type_Probe
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_ProbePublicProcedure
  !> @{
  !> @brief Function for writing Type_Probe data.
  !> The probe data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function write_probe(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Probe), intent(IN), optional:: scalar         !< Scalar probe data.
  type(Type_Probe), intent(IN), optional:: array1D(:)     !< One dimensional array probe data.
  type(Type_Probe), intent(IN), optional:: array2D(:,:)   !< Two dimensional array probe data.
  type(Type_Probe), intent(IN), optional:: array3D(:,:,:) !< Three dimensional array probe data.
  character(*),     intent(IN), optional:: format         !< Format specifier.
  integer(I4P),     intent(IN)::           unit           !< Logic unit.
  integer(I_P)::                           err            !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        write(unit,*,iostat=err)scalar
      elseif (present(array1D)) then
        write(unit,*,iostat=err)array1D
      elseif (present(array2D)) then
        write(unit,*,iostat=err)array2D
      elseif (present(array3D)) then
        write(unit,*,iostat=err)array3D
      endif
    case default
      if (present(scalar)) then
        write(unit,adjustl(trim(format)),iostat=err)scalar
      elseif (present(array1D)) then
        write(unit,adjustl(trim(format)),iostat=err)array1D
      elseif (present(array2D)) then
        write(unit,adjustl(trim(format)),iostat=err)array2D
      elseif (present(array3D)) then
        write(unit,adjustl(trim(format)),iostat=err)array3D
      endif
    endselect
  else
    if (present(scalar)) then
      write(unit,iostat=err)scalar
    elseif (present(array1D)) then
      write(unit,iostat=err)array1D
    elseif (present(array2D)) then
      write(unit,iostat=err)array2D
    elseif (present(array3D)) then
      write(unit,iostat=err)array3D
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction write_probe

  !> @brief Function for reading Type_Probe data.
  !> The probe data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function read_probe(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Probe), intent(INOUT), optional:: scalar         !< Scalar probe data.
  type(Type_Probe), intent(INOUT), optional:: array1D(:)     !< One dimensional array probe data.
  type(Type_Probe), intent(INOUT), optional:: array2D(:,:)   !< Two dimensional array probe data.
  type(Type_Probe), intent(INOUT), optional:: array3D(:,:,:) !< Three dimensional array probe data.
  character(*),     intent(IN),    optional:: format         !< Format specifier.
  integer(I4P),     intent(IN)::              unit           !< Logic unit.
  integer(I_P)::                              err            !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        read(unit,*,iostat=err)scalar
      elseif (present(array1D)) then
        read(unit,*,iostat=err)array1D
      elseif (present(array2D)) then
        read(unit,*,iostat=err)array2D
      elseif (present(array3D)) then
        read(unit,*,iostat=err)array3D
      endif
    case default
      if (present(scalar)) then
        read(unit,adjustl(trim(format)),iostat=err)scalar
      elseif (present(array1D)) then
        read(unit,adjustl(trim(format)),iostat=err)array1D
      elseif (present(array2D)) then
        read(unit,adjustl(trim(format)),iostat=err)array2D
      elseif (present(array3D)) then
        read(unit,adjustl(trim(format)),iostat=err)array3D
      endif
    endselect
  else
    if (present(scalar)) then
      read(unit,iostat=err)scalar
    elseif (present(array1D)) then
      read(unit,iostat=err)array1D
    elseif (present(array2D)) then
      read(unit,iostat=err)array2D
    elseif (present(array3D)) then
      read(unit,iostat=err)array3D
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction read_probe
  !> @}

  !> @ingroup Data_Type_ProbePrivateProcedure
  !> @{
  !> Subroutine for setting members of Type_Probe.
  elemental subroutine set(probe,b,i,j,k)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Probe), intent(INOUT)::        probe !< Probe data.
  integer(I_P),      intent(IN), optional:: b     !< Block (global map) index.
  integer(I_P),      intent(IN), optional:: i     !< I direction index.
  integer(I_P),      intent(IN), optional:: j     !< J direction index.
  integer(I_P),      intent(IN), optional:: k     !< K direction index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) probe%b = b
  if (present(i)) probe%i = i
  if (present(j)) probe%j = j
  if (present(k)) probe%k = k
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set
  !> @}
endmodule Data_Type_Probe
