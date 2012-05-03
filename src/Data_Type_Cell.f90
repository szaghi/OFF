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
public:: init,set,get
public:: write,read
public:: get_cell_id_str,get_cell_id
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! To handle the conditions of cell there are integers and characters parameters. Each integer parameter has a
! corresponding character parameter.
! Id parameters:
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
!> Derived type containing cell informations.
type, public:: Type_Cell
  sequence
  integer(I1P):: id = cell_active !< Id identifier.
  integer(I1P):: oc = 0_I1P       !< Octree refinement: oc>0 cell is refined oc-times.
endtype Type_Cell
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Write overloading of Type_Cell variable.
!> This is a generic interface to 8 functions: there are 2 functions (one binary and another ascii) for writing scalar variables,
!> 1D/2D or 3D arrays. The functions return an error integer code. The calling signatures are:
!> @code ...
!> integer(I4P):: err,unit
!> character(1):: format="*"
!> type(Type_BC):: cell_scal,cell_1D(10),cell_2D(10,2),cell_3D(10,2,3)
!> ...
!> ! formatted writing of cell_scal, cell_1D, cell_2D and cell_3D
!> err = write(unit,format,cell_scal)
!> err = write(unit,format,cell_1D)
!> err = write(unit,format,cell_2D)
!> err = write(unit,format,cell_3D)
!> ! binary writing of cell_scal, cell_1D, cell_2D and cell_3D
!> err = write(unit,cell_scal)
!> err = write(unit,cell_1D)
!> err = write(unit,cell_2D)
!> err = write(unit,cell_3D)
!> ... @endcode
interface write
  module procedure Write_Bin_Scalar, Write_Ascii_Scalar
  module procedure Write_Bin_Array1D,Write_Ascii_Array1D
  module procedure Write_Bin_Array2D,Write_Ascii_Array2D
  module procedure Write_Bin_Array3D,Write_Ascii_Array3D
endinterface
!> @brief Read overloading of Type_Cell variable.
!> This is a generic interface to 8 functions: there are 2 functions (one binary and another ascii) for reading scalar variables,
!> 1D/2D or 3D arrays. The functions return an error integer code. The calling signatures are:
!> @code ...
!> integer(I4P):: err,unit
!> character(1):: format="*"
!> type(Type_BC):: cell_scal,cell_1D(10),cell_2D(10,2),cell_3D(10,2,3)
!> ...
!> ! formatted reading of cell_scal, cell_1D, cell_2D and cell_3D
!> err = read(unit,format,cell_scal)
!> err = read(unit,format,cell_1D)
!> err = read(unit,format,cell_2D)
!> err = read(unit,format,cell_3D)
!> ! binary reading of cell_scal, cell_1D, cell_2D and cell_3D
!> err = read(unit,cell_scal)
!> err = read(unit,cell_1D)
!> err = read(unit,cell_2D)
!> err = read(unit,cell_3D)
!> ... @endcode
interface read
  module procedure Read_Bin_Scalar, Read_Ascii_Scalar
  module procedure Read_Bin_Array1D,Read_Ascii_Array1D
  module procedure Read_Bin_Array2D,Read_Ascii_Array2D
  module procedure Read_Bin_Array3D,Read_Ascii_Array3D
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !>Function for initializing Type_Cell variable.
  !> @return \b cell Type_Cell variable.
  elemental function init(id,oc) result(cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN), optional:: id   !< Id identifier.
  integer(I1P), intent(IN), optional:: oc   !< Octree refinement: oc>0 cell is refined oc-times.
  type(Type_Cell)::                    cell !< Cell informations data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(id)) cell%id = id
  if (present(oc)) cell%oc = oc
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction init

  !>Subroutine for assignment Type_Cell variable.
  !> @return \b cell Type_Cell variable.
  elemental subroutine set(id,oc,cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),    intent(IN), optional:: id   !< Id identifier.
  integer(I1P),    intent(IN), optional:: oc   !< Octree refinement: oc>0 cell is refined oc-times.
  type(Type_Cell), intent(INOUT)::        cell !< Cell informations data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(id)) cell%id = id
  if (present(oc)) cell%oc = oc
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> Subroutine for extracting Type_Cell variable components.
  elemental subroutine get(id,oc,cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),    intent(OUT), optional:: id   !< Id identifier.
  integer(I1P),    intent(OUT), optional:: oc   !< Octree refinement: oc>0 cell is refined oc-times.
  type(Type_Cell), intent(IN)::            cell !< Cell informations data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(id)) id = cell%id
  if (present(oc)) oc = cell%oc
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get

  ! write
  function Write_Bin_Scalar(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, scalar) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN):: unit  ! logic unit
  type(Type_Cell), intent(IN):: cell
  integer(I_P)::                err   ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err) cell
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Scalar

  function Write_Ascii_Scalar(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, scalar) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN):: unit   ! logic unit
  character(*),    intent(IN):: format ! format specifier
  type(Type_Cell), intent(IN):: cell
  integer(I_P)::                err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err) cell
  case default
    write(unit,adjustl(trim(format)),iostat=err) cell
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Scalar

  function Write_Bin_Array1D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, array 1D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN):: unit     ! logic unit
  type(Type_Cell), intent(IN):: cell(:)
  integer(I_P)::                err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err) cell
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array1D

  function Write_Ascii_Array1D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, array 1D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN):: unit    ! logic unit
  character(*),    intent(IN):: format  ! format specifier
  type(Type_Cell), intent(IN):: cell(:)
  integer(I_P)::                err     ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err) cell
  case default
    write(unit,adjustl(trim(format)),iostat=err) cell
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array1D

  function Write_Bin_Array2D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, array 2D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN):: unit     ! logic unit
  type(Type_Cell), intent(IN):: cell(:,:)
  integer(I_P)::                err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err) cell
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array2D

  function Write_Ascii_Array2D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, array 2D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN):: unit    ! logic unit
  character(*),    intent(IN):: format  ! format specifier
  type(Type_Cell), intent(IN):: cell(:,:)
  integer(I_P)::                err     ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err) cell
  case default
    write(unit,adjustl(trim(format)),iostat=err) cell
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array2D

  function Write_Bin_Array3D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, array 3D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN):: unit     ! logic unit
  type(Type_Cell), intent(IN):: cell(:,:,:)
  integer(I_P)::                err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err) cell
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array3D

  function Write_Ascii_Array3D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, array 3D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN):: unit    ! logic unit
  character(*),    intent(IN):: format  ! format specifier
  type(Type_Cell), intent(IN):: cell(:,:,:)
  integer(I_P)::                err     ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err) cell
  case default
    write(unit,adjustl(trim(format)),iostat=err) cell
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array3D

  ! read
  function Read_Bin_Scalar(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, scalar) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN)::    unit  ! logic unit
  type(Type_Cell), intent(INOUT):: cell
  integer(I_P)::                   err   ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err) cell
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Scalar

  function Read_Ascii_Scalar(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, scalar) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN)::    unit   ! logic unit
  character(*),    intent(IN)::    format ! format specifier
  type(Type_Cell), intent(INOUT):: cell
  integer(I_P)::                   err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err) cell
  case default
    read(unit,adjustl(trim(format)),iostat=err) cell
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Scalar

  function Read_Bin_Array1D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, array 1D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN)::    unit  ! logic unit
  type(Type_Cell), intent(INOUT):: cell(:)
  integer(I_P)::                   err   ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err) cell
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array1D

  function Read_Ascii_Array1D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, array 1D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN)::    unit     ! logic unit
  character(*),    intent(IN)::    format   ! format specifier
  type(Type_Cell), intent(INOUT):: cell(:)
  integer(I_P)::                   err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err) cell
  case default
    read(unit,adjustl(trim(format)),iostat=err) cell
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array1D

  function Read_Bin_Array2D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, array 2D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN)::    unit  ! logic unit
  type(Type_Cell), intent(INOUT):: cell(:,:)
  integer(I_P)::                   err   ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err) cell
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array2D

  function Read_Ascii_Array2D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, array 2D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN)::    unit     ! logic unit
  character(*),    intent(IN)::    format   ! format specifier
  type(Type_Cell), intent(INOUT):: cell(:,:)
  integer(I_P)::                   err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err) cell
  case default
    read(unit,adjustl(trim(format)),iostat=err) cell
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array2D

  function Read_Bin_Array3D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, array 3D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN)::    unit  ! logic unit
  type(Type_Cell), intent(INOUT):: cell(:,:,:)
  integer(I_P)::                   err   ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err) cell
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array3D

  function Read_Ascii_Array3D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, array 3D) Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),    intent(IN)::    unit     ! logic unit
  character(*),    intent(IN)::    format   ! format specifier
  type(Type_Cell), intent(INOUT):: cell(:,:,:)
  integer(I_P)::                   err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err) cell
  case default
    read(unit,adjustl(trim(format)),iostat=err) cell
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array3D

  !> Function for getting integer id from the corresponding string id.
  !> @return \b id integer(I1P) variable.
  function get_cell_id(id_str) result(id)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: id_str !< String id.
  integer(I1P)::             id     !< Integer id.
  integer(I_P)::             i      !< Identifier counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=1,Nid
    if (adjustl(trim(id_str))==adjustl(trim(cell_list_str(i)))) then
      id = cell_list(i)
      exit
    elseif (i==Nid) then
      write(stderr,'(A)')' Attention!'
      write(stderr,'(A)')' The cell identifier:'
      write(stderr,'(A)')' '//adjustl(trim(id_str))
      write(stderr,'(A)')' is unknown!'
      stop
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction get_cell_id

  !> Function for getting string id from the corresponding integer id.
  !> @return \b id_str character(3) variable.
  function get_cell_id_str(id) result(id_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: id     !< Integer id.
  character(3)::             id_str !< String id.
  integer(I_P)::             i      !< Identifier counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=1,Nid
    if (id==cell_list(i)) then
      id_str = adjustl(trim(cell_list_str(i)))
      exit
    elseif (i==Nid) then
      write(stderr,'(A)')' Attention!'
      write(stderr,'(A)')' The cell identifier:'
      write(stderr,FI_P) id
      write(stderr,'(A)')' is unknown!'
      stop
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction get_cell_id_str
endmodule Data_Type_Cell
