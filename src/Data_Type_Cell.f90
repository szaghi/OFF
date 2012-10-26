!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_Cell Data_Type_Cell
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_CellPublicProcedure Data_Type_Cell
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_CellPrivateProcedure Data_Type_Cell
!> @}

!> This module contains the definition of Type_Cell and its procedures.
!> @todo \b DocWriteRead: Complete the documentation of write and read functions
module Data_Type_Cell
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                             ! Integers and reals precision definition.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: cell_unset,cell_unset_str
public:: cell_active,cell_active_str
public:: cell_passive,cell_passive_str
public:: cell_probe,cell_probe_str
public:: cell_isoss,cell_isoss_str
public:: Nid
public:: cell_list,cell_list_str
public:: write_cell,read_cell
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_Cell
!> @{
character(3), parameter:: cell_unset_str   = 'UNS'  !< String corresponding to unset cell identifier.
integer(I1P), parameter:: cell_unset       = -1_I1P !< Unset cell identifier.
character(3), parameter:: cell_active_str  = 'ACT'  !< String corresponding to active cell identifier.
integer(I1P), parameter:: cell_active      =  0_I1P !< Active cell identifier.
character(3), parameter:: cell_passive_str = 'PAS'  !< String corresponding to passive identifier.
integer(I1P), parameter:: cell_passive     = 1_I1P  !< Passive cell identifier.
character(3), parameter:: cell_probe_str   = 'PRB'  !< String corresponding to probe cell identifier.
integer(I1P), parameter:: cell_probe       = 2_I1P  !< Probe cell identifier.
character(3), parameter:: cell_isoss_str   = 'ISS'  !< String corresponding to isotropic cell identifier.
integer(I1P), parameter:: cell_isoss       = 3_I1P  !< Isotropic cell identifier.
integer(I_P), parameter:: Nid = 5                   !< Number of possible cell identifiers.
character(3), parameter:: cell_list_str(1:Nid) = &
                                               (/ cell_unset_str,   &
                                                  cell_active_str,  &
                                                  cell_passive_str, &
                                                  cell_probe_str,   &
                                                  cell_isoss_str    &
                                                /)  !< Cell identifiers string list.
integer(I1P), parameter:: cell_list(1:Nid) = &
                                           (/ cell_unset,   &
                                              cell_active,  &
                                              cell_passive, &
                                              cell_probe,   &
                                              cell_isoss    &
                                            /)      !< Cell identifiers list.
!> @}
!> Derived type containing cell informations.
!> @ingroup DerivedType
type, public:: Type_Cell
  integer(I1P):: id = cell_active !< Id identifier.
  integer(I1P):: oc = 0_I1P       !< Octree refinement: oc>0 cell is refined oc-times.
  contains
    procedure, non_overridable:: str2id ! Procedure for setting integer id from string id.
    procedure, non_overridable:: id2str ! Procedure for converting integer id to string id.
endtype Type_Cell
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_CellPublicProcedure
  !> @{
  !> @brief Function for writing Type_Cell data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function write_cell(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Cell), intent(IN), optional:: scalar         !< Scalar vector data.
  type(Type_Cell), intent(IN), optional:: array1D(:)     !< One dimensional array vector data.
  type(Type_Cell), intent(IN), optional:: array2D(:,:)   !< Two dimensional array vector data.
  type(Type_Cell), intent(IN), optional:: array3D(:,:,:) !< Three dimensional array vector data.
  character(*),    intent(IN), optional:: format         !< Format specifier.
  integer(I4P),    intent(IN)::           unit           !< Logic unit.
  integer(I_P)::                          err            !< Error trapping flag: 0 no errors, >0 error occurs.
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
  endfunction write_cell

  !> @brief Function for reading Type_Cell data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function read_cell(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Cell), intent(INOUT), optional:: scalar         !< Scalar vector data.
  type(Type_Cell), intent(INOUT), optional:: array1D(:)     !< One dimensional array vector data.
  type(Type_Cell), intent(INOUT), optional:: array2D(:,:)   !< Two dimensional array vector data.
  type(Type_Cell), intent(INOUT), optional:: array3D(:,:,:) !< Three dimensional array vector data.
  character(*),    intent(IN),    optional:: format         !< Format specifier.
  integer(I4P),    intent(IN)::              unit           !< Logic unit.
  integer(I_P)::                             err            !< Error trapping flag: 0 no errors, >0 error occurs.
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
  endfunction read_cell
  !> @}

  !> @ingroup Data_Type_CellPrivateProcedure
  !> @{
  !> @brief Subroutine for setting integer id from string id.
  elemental subroutine str2id(cell,id_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell), intent(INOUT):: cell   !< Cell data.
  character(3),     intent(IN)::    id_str !< String of id.
  integer(I_P)::                    i      !< Cell id counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=1,Nid
    if (adjustl(trim(id_str))==adjustl(trim(cell_list_str(i)))) then
      cell%id = cell_list(i)
      exit
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine str2id

  !> @brief Function for converting integer id to string id.
  !> @return \b bc_str character(3) variable.
  elemental function id2str(cell) result(id_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Cell), intent(IN):: cell   !< Cell data.
  character(3)::                 id_str !< String of id.
  integer(I_P)::                 i      !< Cell id counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=1,Nid
    if (cell%id==cell_list(i)) then
      id_str = adjustl(trim(cell_list_str(i)))
      exit
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction id2str
  !> @}
endmodule Data_Type_Cell
