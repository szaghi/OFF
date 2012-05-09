!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_ProbePrivateProcedure Data_Type_Probe
!> @}

!>This module contains the definition of Type_Probe and its procedures.
!> @todo \b DocComplete: Complete the documentation of internal procedures
module Data_Type_Probe
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: set
public:: write,read
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing probe informations.
!> @ingroup DerivedType
type, public:: Type_Probe
  sequence
  integer(I_P):: b = 1_I_P !< Block (global map) index.
  integer(I_P):: i = 0_I_P !< I direction index.
  integer(I_P):: j = 0_I_P !< J direction index.
  integer(I_P):: k = 0_I_P !< K direction index.
endtype Type_Probe
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Write overloading of Type_Probe variable.
!> This is a generic interface to 2 functions (one binary and another ascii) for writing scalar variables.
!> The functions return an error integer code. The calling signatures are:
!> @code ...
!> integer(I4P):: err,unit
!> character(1):: format="*"
!> type(Type_Probe):: probe
!> ...
!> ! formatted writing
!> err = write(unit,format,probe)
!> ! binary writing
!> err = write(unit,probe)
!> ... @endcode
!> @ingroup Interface
interface write
  module procedure Write_Bin,Write_Ascii
endinterface
!> @brief Read overloading of Type_Probe variable.
!> This is a generic interface to 2 functions (one binary and another ascii) for reading scalar variables.
!> The functions return an error integer code. The calling signatures are:
!> @code ...
!> integer(I4P):: err,unit
!> character(1):: format="*"
!> type(Type_Probe):: probe
!> ...
!> ! formatted reading
!> err = read(unit,format,probe)
!> ! binary reading
!> err = read(unit,probe)
!> ... @endcode
!> @ingroup Interface
interface read
  module procedure Read_Bin,Read_Ascii
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !>Function for setting components of Type_Probe variable.
  !> @return \b probe type(Type_Probe) variable.
  elemental subroutine set(b,i,j,k,probe)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),     intent(IN), optional:: b     !< Block (global map) index.
  integer(I_P),     intent(IN), optional:: i     !< I direction index.
  integer(I_P),     intent(IN), optional:: j     !< J direction index.
  integer(I_P),     intent(IN), optional:: k     !< K direction index.
  type(Type_Probe), intent(INOUT)::        probe !< Probe data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(b)) probe%b = b
  if (present(i)) probe%i = i
  if (present(j)) probe%j = j
  if (present(k)) probe%k = k
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> @ingroup Data_Type_ProbePrivateProcedure
  !> @{
  ! write
  function Write_Bin(unit,probe) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (binary, scalar) Type_Probe.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),     intent(IN):: unit  ! logic unit
  type(Type_Probe), intent(IN):: probe
  integer(I_P)::                 err   ! Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)probe
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin

  function Write_Ascii(unit,format,probe) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing (ascii, scalar) Type_Probe.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),     intent(IN):: unit   ! logic unit
  character(*),     intent(IN):: format ! format specifier
  type(Type_probe), intent(IN):: probe
  integer(I_P)::                 err    ! Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)probe
  case default
    write(unit,adjustl(trim(format)),iostat=err)probe
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii

  ! read
  function Read_Bin(unit,probe) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (binary, scalar) Type_Probe.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),     intent(IN)::    unit  ! logic unit
  type(Type_Probe), intent(INOUT):: probe
  integer(I_P)::                    err   ! Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)probe
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin

  function Read_Ascii(unit,format,probe) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading (ascii, scalar) Type_Probe.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),     intent(IN)::    unit   ! logic unit
  character(*),     intent(IN)::    format ! format specifier
  type(Type_Probe), intent(INOUT):: probe
  integer(I_P)::                    err    ! Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)probe
  case default
    read(unit,adjustl(trim(format)),iostat=err)probe
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii
  !> @}
endmodule Data_Type_Probe
