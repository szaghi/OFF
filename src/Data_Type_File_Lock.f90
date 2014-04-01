!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_LockDerivedType Data_Type_File_Lock
!> Module definition of Type_File_Lock
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_LockPrivateProcedure Data_Type_File_Lock
!> Module definition of Type_File_Lock
!> @}

!> @brief Module Data_Type_File_Lock contains the definition of Type_File_Lock.
module Data_Type_File_Lock
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                          ! Integers and reals precision definition.
USE Data_Type_File_Base, only: Type_File_Base, err_locked ! Definition of Type_File_Base.
USE Data_Type_Time,      only: Get_Date_String            ! Prcedure for getting current data.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Lock.
!> @ingroup Data_Type_File_LockDerivedType
type, public, extends(Type_File_Base):: Type_File_Lock
  contains
    procedure:: lock   ! Procedure for locking, namely creating a lockfile.
    procedure:: unlock ! Procedure for unlocking, namely deleting a previously created lockfile.
endtype Type_File_Lock
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_LockPrivateProcedure
  !> @{
  !> @brief Procedure for locking, namely creating a lockfile.
  subroutine lock(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Lock), intent(INOUT):: file_d  !< File data.
  logical::                              is_file !< Inquiring flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%set(name='lockfile',path_in='',path_out='')

  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    inquire(file=file_d%name,exist=is_file,iostat=iostat)
    if (is_file) then
      call file_d%raise_error(errtype=err_locked)
    else
      call file_d%open(ascii=.true.,action='WRITE') ; if (iostat/=0) return
      write(unit,'(A)')' Simulation started on'
      write(unit,'(A)')' '//Get_Date_String()
    endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine lock

  !> @brief Procedure for unlocking, namely deleting a previously created lockfile.
  subroutine unlock(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Lock), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%close(delete=.true.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine unlock
  !> @}
endmodule Data_Type_File_Lock
