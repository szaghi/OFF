!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_OS Data_Type_OS
!> @}

!> This module contains the definition of Type_OS and its procedures.
!> This derived type has useful parameters for performing system calls.
!> @bug <b>MS Windows directory separator's documentation</b>: \n The documentation of variable "win_sep" containing MS Windows
!>                                                              directory separator must be skipped because doxygen produces error
!>                                                              with "\" character.
module Data_Type_OS
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                    ! Integers and reals precision definition.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: uix_id, c_uix_id, uix_sep, uix_remove, uix_copy, uix_mkdir
public:: win_id, c_win_id, win_sep, win_remove, win_copy, win_mkdir
public:: init,set
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_OS
!> @{
character(3), parameter:: c_uix_id   = "UIX"   !< Unix/Linux string identifier.
integer(I1P), parameter:: uix_id     = 1_I1P   !< Unix/Linux identifier.
character(1), parameter:: uix_sep    = "/"     !< Unix/Linux directories separator.
character(2), parameter:: uix_remove = "rm"    !< Unix/Linux remove command.
character(2), parameter:: uix_copy   = "cp"    !< Unix/Linux copy command.
character(5), parameter:: uix_mkdir  = "mkdir" !< Unix/Linux make dir command.
character(3), parameter:: c_win_id   = "WIN"   !< MS Windows string identifier.
integer(I1P), parameter:: win_id     = 2_I1P   !< MS Windows identifier.
#ifndef DOXYGEN_SKIP
character(1), parameter:: win_sep    = "\"     !< MS Windows directories separator.
#endif
character(3), parameter:: win_remove = "del"   !< MS Windows remove command.
character(4), parameter:: win_copy   = "copy"  !< MS Windows copy command.
character(5), parameter:: win_mkdir  = "mkdir" !< MS Windows make dir command.
!> @}
!> Derived type contains useful parameters for performing portable system calls.
!> @ingroup DerivedType
type, public:: Type_OS
 sequence
 integer(I1P):: id     = uix_id     !< OS id.
 character(1):: sep    = uix_sep    !< OS directories separator.
 character(3):: remove = uix_remove !< OS remove command.
 character(4):: copy   = uix_copy   !< OS copy command.
 character(5):: mkdir  = uix_mkdir  !< OS make dir command.
endtype Type_OS
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !>Function for initializing Type_OS.
  function init(id,c_id,sep,remove,copy,mkdir) result(OS)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN), optional:: id     !< OS id parameter (integer).
  character(*), intent(IN), optional:: c_id   !< OS id parameter (string).
  character(*), intent(IN), optional:: sep    !< Directories separator.
  character(*), intent(IN), optional:: remove !< Remove file command.
  character(*), intent(IN), optional:: copy   !< Copy file command.
  character(*), intent(IN), optional:: mkdir  !< Make directory command.
  type(Type_OS)::                      OS     !< Output OS.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(id))     call Set(id=id,OS=OS)
  if (present(c_id))   call Set(c_id=c_id,OS=OS)
  if (present(sep))    OS%sep    = trim(sep)
  if (present(remove)) OS%remove = trim(remove)
  if (present(copy))   OS%copy   = trim(copy)
  if (present(mkdir))  OS%mkdir  = trim(mkdir)
  select case(OS%id)
  case(uix_id)
    write(stdout,'(A)')' The OS has been set as *nix like'
  case(win_id)
    write(stdout,'(A)')' The OS has been set as windows like'
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction init

  !>Subroutine for setting Type_OS (id).
  subroutine set(id,c_id,sep,remove,copy,mkdir,OS)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),  intent(IN), optional:: id     !< OS id parameter (integer).
  character(*),  intent(IN), optional:: c_id   !< OS id parameter (string).
  character(*),  intent(IN), optional:: sep    !< Directories separator.
  character(*),  intent(IN), optional:: remove !< Remove file command.
  character(*),  intent(IN), optional:: copy   !< Copy file command.
  character(*),  intent(IN), optional:: mkdir  !< Make directory command.
  type(Type_OS), intent(INOUT)::        OS     !< Output OS.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(id)) then
    OS%id     = id
    select case(id)
    case(uix_id)
      OS%sep    = uix_sep
      OS%remove = uix_remove
      OS%copy   = uix_copy
      OS%mkdir  = uix_mkdir
    case(win_id)
      OS%sep    = win_sep
      OS%remove = win_remove
      OS%copy   = win_copy
      OS%mkdir  = win_mkdir
    case default
      write(stderr,'(A)')' OS id not recognized!'
      write(stderr,'(A)')' Valid integer id are:'
      write(stderr,'(A)')' "'//trim(str(.true.,uix_id))//'" for *nix OS'
      write(stderr,'(A)')' "'//trim(str(.true.,win_id))//'" for Windows OS'
      stop
    endselect
  endif
  if (present(c_id)) then
    select case(c_id)
    case(c_uix_id)
      OS%id     = uix_id
      OS%sep    = uix_sep
      OS%remove = uix_remove
      OS%copy   = uix_copy
      OS%mkdir  = uix_mkdir
    case(c_win_id)
      OS%id     = win_id
      OS%sep    = win_sep
      OS%remove = win_remove
      OS%copy   = win_copy
      OS%mkdir  = win_mkdir
    case default
      write(stderr,'(A)')' OS id not recognized!'
      write(stderr,'(A)')' Valid charcter id are:'
      write(stderr,'(A)')' "'//c_uix_id//'" for *nix OS'
      write(stderr,'(A)')' "'//c_win_id//'" for Windows OS'
      stop
    endselect
  endif
  if (present(sep))    OS%sep    = trim(sep)
  if (present(remove)) OS%remove = trim(remove)
  if (present(copy))   OS%copy   = trim(copy)
  if (present(mkdir))  OS%mkdir  = trim(mkdir)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set
endmodule Data_Type_OS
