module Data_Type_Cell
!-----------------------------------------------------------------------------------------------------------------------------------
!!This module contains the definition of Type_Cell and its functions and subroutines.
!-----------------------------------------------------------------------------------------------------------------------------------

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
public:: Type_Cell
public:: init,set,get
public:: write,read
public:: get_cell_id_str,get_cell_id
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! To handle the conditions of cell there are integers and characters parameters. Each integer parameter has a
! corresponding character parameter.
! Id parameters:
character(3), parameter:: cell_unset_str   = 'UNS'  ! String corresponding to unset cell identifier.
integer(I1P), parameter:: cell_unset       = -1_I1P ! Unset cell identifier.
character(3), parameter:: cell_active_str  = 'ACT'  ! String corresponding to active cell identifier.
integer(I1P), parameter:: cell_active      =  0_I1P ! Active cell identifier.
character(3), parameter:: cell_passive_str = 'PAS'  ! String corresponding to passive identifier.
integer(I1P), parameter:: cell_passive     = 1_I1P  ! Passive cell identifier.
character(3), parameter:: cell_probe_str   = 'PRB'  ! String corresponding to probe cell identifier.
integer(I1P), parameter:: cell_probe       = 2_I1P  ! Probe cell identifier.
character(3), parameter:: cell_isoss_str   = 'ISS'  ! String corresponding to isentropic cell identifier.
integer(I1P), parameter:: cell_isoss       = 3_I1P  ! Isentropic cell identifier.
integer(I_P), parameter:: Nid = 5                   ! Number of possible cell identifiers.
character(3), parameter:: cell_list_str(1:Nid) = (/ cell_unset_str,   &
                                                    cell_active_str,  &
                                                    cell_passive_str, &
                                                    cell_probe_str,   &
                                                    cell_isoss_str /) ! Cell identifiers string list.
integer(I1P), parameter:: cell_list(1:Nid) = (/ cell_unset,   &
                                                cell_active,  &
                                                cell_passive, &
                                                cell_probe,   &
                                                cell_isoss /) ! Cell identifiers list.
! Type_Cell definition:
type:: Type_Cell
  sequence
  integer(I1P):: id = cell_active ! Id identifier.
  integer(I1P):: oc = 0_I1P       ! Octree refinemnt: oc>0 cell is refined.
endtype Type_Cell
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Write overloading.
interface write
  module procedure Write_Bin_Scalar,        & ! binary scalar
                   Write_Ascii_Scalar,      & ! ascii scalar
                   Write_Bin_Vectorial1D,   & ! binary vectorial 1D
                   Write_Ascii_Vectorial1D, & ! ascii vectorial 1D
                   Write_Bin_Vectorial2D,   & ! binary vectorial 2D
                   Write_Ascii_Vectorial2D, & ! ascii vectorial 2D
                   Write_Bin_Vectorial3D,   & ! binary vectorial 3D
                   Write_Ascii_Vectorial3D    ! ascii vectorial 3D
endinterface
!!Read overloading.
interface read
  module procedure Read_Bin_Scalar,        & ! binary scalar
                   Read_Ascii_Scalar,      & ! ascii scalar
                   Read_Bin_Vectorial1D,   & ! binary vectorial 1D
                   Read_Ascii_Vectorial1D, & ! ascii vectorial 1D
                   Read_Bin_Vectorial2D,   & ! binary vectorial 2D
                   Read_Ascii_Vectorial2D, & ! ascii vectorial 2D
                   Read_Bin_Vectorial3D,   & ! binary vectorial 3D
                   Read_Ascii_Vectorial3D    ! ascii vectorial 3D
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function init(id,oc) result(cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for initializing Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN), optional:: id
  integer(I1P), intent(IN), optional:: oc
  type(Type_Cell)::                    cell
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(id)) cell%id = id
  if (present(oc)) cell%oc = oc
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction init

  elemental subroutine set(id,oc,cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment Type_Cell.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),    intent(IN), optional:: id
  integer(I1P),    intent(IN), optional:: oc
  type(Type_Cell), intent(INOUT)::        cell
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(id)) cell%id = id
  if (present(oc)) cell%oc = oc
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  elemental subroutine get(id,oc,cell)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for extraction Type_Cell attributes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),    intent(OUT), optional:: id
  integer(I1P),    intent(OUT), optional:: oc
  type(Type_Cell), intent(IN)::            cell
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

  function Write_Bin_Vectorial1D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, vectorial 1D) Type_Cell.
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
  endfunction Write_Bin_Vectorial1D

  function Write_Ascii_Vectorial1D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, vectorial 1D) Type_Cell.
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
  endfunction Write_Ascii_Vectorial1D

  function Write_Bin_Vectorial2D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, vectorial 2D) Type_Cell.
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
  endfunction Write_Bin_Vectorial2D

  function Write_Ascii_Vectorial2D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, vectorial 2D) Type_Cell.
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
  endfunction Write_Ascii_Vectorial2D

  function Write_Bin_Vectorial3D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, vectorial 3D) Type_Cell.
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
  endfunction Write_Bin_Vectorial3D

  function Write_Ascii_Vectorial3D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, vectorial 3D) Type_Cell.
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
  endfunction Write_Ascii_Vectorial3D

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

  function Read_Bin_Vectorial1D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, vectorial 1D) Type_Cell.
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
  endfunction Read_Bin_Vectorial1D

  function Read_Ascii_Vectorial1D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, vectorial 1D) Type_Cell.
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
  endfunction Read_Ascii_Vectorial1D

  function Read_Bin_Vectorial2D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, vectorial 2D) Type_Cell.
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
  endfunction Read_Bin_Vectorial2D

  function Read_Ascii_Vectorial2D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, vectorial 2D) Type_Cell.
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
  endfunction Read_Ascii_Vectorial2D

  function Read_Bin_Vectorial3D(unit,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, vectorial 3D) Type_Cell.
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
  endfunction Read_Bin_Vectorial3D

  function Read_Ascii_Vectorial3D(unit,format,cell) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, vectorial 3D) Type_Cell.
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
  endfunction Read_Ascii_Vectorial3D

  ! get_cell_id
  function get_cell_id(id_str) result(id)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for getting integer id from the corresponding character id.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: id_str ! String id.
  integer(I1P)::             id     ! Integer id.
  integer(I_P)::             i      ! Identifier counter.
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

  function get_cell_id_str(id) result(id_str)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for getting character id from the corresponding integer id.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: id     ! Integer id.
  character(3)::             id_str ! String id.
  integer(I_P)::             i      ! Identifier counter.
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
