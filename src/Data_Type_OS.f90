!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_OSDerivedType Data_Type_OS
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_OSGlobalVarPar Data_Type_OS
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_OSPrivateProcedure Data_Type_OS
!> @}

!> This module contains the definition of Type_OS and its procedures.
!> This derived type has useful parameters for performing system calls.
!> @bug <b>MS Windows directory separator's documentation</b>: \n The documentation of variable "win_sep" containing MS Windows
!>                                                              directory separator must be skipped because doxygen produces error
!>                                                              with "\" character.
module Data_Type_OS
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                             ! Integers and reals precision definition.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stderr => ERROR_UNIT ! Standard output/error logical units.
#ifdef MPI2
USE MPI                                                      ! MPI runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: uix_id, c_uix_id, uix_sep, uix_remove, uix_copy, uix_mkdir
public:: win_id, c_win_id, win_sep, win_remove, win_copy, win_mkdir
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_OSGlobalVarPar
!> @{
character(3), parameter:: c_uix_id   = "UIX"    !< Unix/Linux string identifier.
integer(I1P), parameter:: uix_id     = 1_I1P    !< Unix/Linux identifier.
character(1), parameter:: uix_sep    = "/"      !< Unix/Linux directories separator.
character(2), parameter:: uix_remove = "rm"     !< Unix/Linux remove command.
character(2), parameter:: uix_copy   = "cp"     !< Unix/Linux copy command.
character(5), parameter:: uix_mkdir  = "mkdir"  !< Unix/Linux make dir command.
character(3), parameter:: c_win_id   = "WIN"    !< MS Windows string identifier.
integer(I1P), parameter:: win_id     = 2_I1P    !< MS Windows identifier.
#ifndef DOXYGEN_SKIP
character(1), parameter:: win_sep    = char(92) !< MS Windows directories separator.
#endif
character(3), parameter:: win_remove = "del"    !< MS Windows remove command.
character(4), parameter:: win_copy   = "copy"   !< MS Windows copy command.
character(5), parameter:: win_mkdir  = "mkdir"  !< MS Windows make dir command.
!> @}
!> @brief Derived type contains useful parameters for performing portable system calls.
!> @note By default the Unix/Linux initalization is used.
!> @ingroup Data_Type_OSDerivedType
type, public:: Type_OS
 integer(I1P):: id     = uix_id     !< OS id.
 character(1):: sep    = uix_sep    !< OS directories separator.
 character(3):: remove = uix_remove !< OS remove command.
 character(4):: copy   = uix_copy   !< OS copy command.
 character(5):: mkdir  = uix_mkdir  !< OS make dir command.
 contains
   procedure:: init ! Procedure for initializing Type_OS.
endtype Type_OS
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Subroutine for initializing Type_OS.
  !> @ingroup Data_Type_OSPrivateProcedure
  subroutine init(OS,myrank,id,c_id)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS), intent(INOUT)::        OS     !< Output OS.
  integer(I_P),   intent(IN), optional:: myrank !< Actual rank process.
  integer(I1P),   intent(IN), optional:: id     !< OS id parameter (integer).
  character(*),   intent(IN), optional:: c_id   !< OS id parameter (string).
  character(DI_P)::                      rks    !< String containing myrank.
#ifdef MPI2
  integer(I_P)::                         err    !< Error trapping flag: 0 no errors, >0 error occurs.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,0_I_P)) ; if (present(myrank)) rks = 'rank'//trim(str(.true.,myrank))
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
      write(stderr,'(A)')trim(rks)//' OS id not recognized!'
      write(stderr,'(A)')trim(rks)//' Valid integer id are:'
      write(stderr,'(A)')trim(rks)//' "'//trim(str(.true.,uix_id))//'" for *nix OS'
      write(stderr,'(A)')trim(rks)//' "'//trim(str(.true.,win_id))//'" for Windows OS'
#ifdef MPI2
      call MPI_FINALIZE(err)
#endif
      stop
    endselect
  elseif (present(c_id)) then
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
      write(stderr,'(A)')trim(rks)//' OS id not recognized!'
      write(stderr,'(A)')trim(rks)//' Valid charcter id are:'
      write(stderr,'(A)')trim(rks)//' "'//c_uix_id//'" for *nix OS'
      write(stderr,'(A)')trim(rks)//' "'//c_win_id//'" for Windows OS'
#ifdef MPI2
      call MPI_FINALIZE(err)
#endif
      stop
    endselect
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init
endmodule Data_Type_OS
