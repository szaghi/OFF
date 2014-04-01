!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_OSDerivedType Data_Type_OS
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_OSGlobalVarPar Data_Type_OS
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_OSPublicProcedure Data_Type_OS
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_OSPrivateProcedure Data_Type_OS
!> @}

!> This module contains the definition of Type_OS and its procedures.
!> This derived type has useful parameters for performing system calls.
module Data_Type_OS
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                    ! Integers and reals precision definition.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_OSGlobalVarPar
!> @{
character(3), parameter:: c_uix_id   = "UIX"      !< Unix/Linux string identifier.
integer(I1P), parameter:: uix_id     = 1_I1P      !< Unix/Linux identifier.
character(1), parameter:: uix_sep    = char(47)   !< Unix/Linux directories separator.
character(2), parameter:: uix_remove = "rm"       !< Unix/Linux remove command.
character(2), parameter:: uix_copy   = "cp"       !< Unix/Linux copy command.
character(8), parameter:: uix_mkdir  = "mkdir -p" !< Unix/Linux make dir command.
character(3), parameter:: c_win_id   = "WIN"      !< MS Windows string identifier.
integer(I1P), parameter:: win_id     = 2_I1P      !< MS Windows identifier.
character(1), parameter:: win_sep    = char(92)   !< MS Windows directories separator.
character(3), parameter:: win_remove = "del"      !< MS Windows remove command.
character(4), parameter:: win_copy   = "copy"     !< MS Windows copy command.
character(5), parameter:: win_mkdir  = "mkdir"    !< MS Windows make dir command.
!> @}
!> @brief Derived type contains list of OS commands.
!> @note By default the OS is set by _SYSTEMxxx preprocessing flag.
!> @ingroup Data_Type_OSDerivedType
type, public:: Type_OS_Command
#if defined _OSYSTEMuix
  ! Unix/Linux OS
  character(3):: rmv = uix_remove !< OS remove command.
  character(4):: cpy = uix_copy   !< OS copy command.
  character(8):: mkd = uix_mkdir  !< OS make dir command.
#elif _OSYSTEMwin
  ! MS Windows OS
  character(3):: rmv = win_remove !< OS remove command.
  character(4):: cpy = win_copy   !< OS copy command.
  character(5):: mkd = win_mkdir  !< OS make dir command.
#else
  ! OS unknown: diabling OS features
  character(3):: rmv = '' !< OS remove command.
  character(4):: cpy = '' !< OS copy command.
  character(5):: mkd = '' !< OS make dir command.
#endif
endtype Type_OS_Command
!> @brief Derived type contains useful parameters for performing portable system calls.
!> @note By default the OS is set by _SYSTEMxxx preprocessing flag.
!> @ingroup Data_Type_OSDerivedType
type, public:: Type_OS
#if defined _OSYSTEMuix
  ! Unix/Linux OS
  character(3):: c_id   = c_uix_id !< OS string id.
  integer(I1P):: id     = uix_id   !< OS id.
  character(1):: sep    = uix_sep  !< OS directories separator.
#elif _OSYSTEMwin
  ! MS Windows OS
  character(3):: c_id   = c_win_id !< OS string id.
  integer(I1P):: id     = win_id   !< OS id.
  character(1):: sep    = win_sep  !< OS directories separator.
#else
  ! OS unknown: diabling OS features, but separator is set as Uix one
  character(3):: c_id   = ''      !< OS string id.
  integer(I1P):: id     = 0_I1P   !< OS id.
  character(1):: sep    = uix_sep !< OS directories separator.
#endif
  type(Type_OS_Command):: cmd !< OS commands list.
  contains
    procedure::         print => print_os    ! Procedure for printing current OS informations.
    procedure::         init                 ! Procedure for initializing (at runtime) Type_OS.
    procedure::         basedir              ! Procedure for stripping last component from file name providing base directory name.
    procedure::         basename             ! Procedure for stripping directory name from file name.
    procedure::         string_separator_fix ! Procedure for fixing the character seperator.
    procedure, nopass:: execute_cmd          ! Procedure for executing shell commands.
    procedure::         copy_file            ! Procedure for copying file.
    procedure::         remove_file          ! Procedure for removing file.
    procedure::         make_dir             ! Procedure for making directory.
endtype Type_OS
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_OSPublicProcedure
  !> @{
  !> @}

  !> @ingroup Data_Type_OSPrivateProcedure
  !> @{
  !> @brief Subroutine for printing current OS informations.
  subroutine print_os(OS,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS),         intent(INOUT):: OS      !< OS.
  character(*), optional, intent(IN)::    pref    !< Prefixing string.
  integer(I4P), optional, intent(OUT)::   iostat  !< IO error.
  character(*), optional, intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),           intent(IN)::    unit    !< Logic unit.
  character(len=:), allocatable::         prefd   !< Prefixing string.
  integer(I4P)::                          iostatd !< IO error.
  character(500)::                        iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  select case(OS%id)
  case(uix_id,win_id)
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' The OS is '//trim(OS%c_id)
  case default
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' The OS is unknown'
  endselect
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   ID: '//trim(str(.true.,OS%id))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   SEP: '//trim(OS%sep)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Commands'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'     remove file   : '//trim(OS%cmd%rmv)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'     copy   file   : '//trim(OS%cmd%cpy)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'     make directory: '//trim(OS%cmd%mkd)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_os

  !> @brief Subroutine for initializing (at runtime) Type_OS.
  subroutine init(OS,myrank,id,c_id)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS), intent(INOUT)::        OS     !< OS.
  integer(I4P),   intent(IN), optional:: myrank !< Actual rank process necessary for concurrent multi-processes calls.
  integer(I1P),   intent(IN), optional:: id     !< OS id parameter.
  character(*),   intent(IN), optional:: c_id   !< OS id parameter (string).
  character(DI4P)::                      rks    !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(id)) then
    select case(id)
    case(uix_id)
      OS%c_id    = c_uix_id
      OS%id      = uix_id
      OS%sep     = uix_sep
      OS%cmd%rmv = uix_remove
      OS%cmd%cpy = uix_copy
      OS%cmd%mkd = uix_mkdir
    case(win_id)
      OS%c_id    = c_win_id
      OS%id      = win_id
      OS%sep     = win_sep
      OS%cmd%rmv = win_remove
      OS%cmd%cpy = win_copy
      OS%cmd%mkd = win_mkdir
    case default
      rks = 'rank'//trim(str(.true.,0_I4P)) ; if (present(myrank)) rks = 'rank'//trim(str(.true.,myrank))
      write(stderr,'(A)')trim(rks)//' Attention: OS id not recognized!'
      write(stderr,'(A)')trim(rks)//' Valid charcter id are:'
      write(stderr,'(A)')trim(rks)//' "'//c_uix_id//'" for *nix OS'
      write(stderr,'(A)')trim(rks)//' "'//c_win_id//'" for Windows OS'
      write(stderr,'(A)')trim(rks)//' The OS is set to unknow. All OS-related calls are disabled!'
    endselect
  elseif (present(c_id)) then
    select case(c_id)
    case(c_uix_id)
      OS%id      = uix_id
      OS%sep     = uix_sep
      OS%cmd%rmv = uix_remove
      OS%cmd%cpy = uix_copy
      OS%cmd%mkd = uix_mkdir
    case(c_win_id)
      OS%id      = win_id
      OS%sep     = win_sep
      OS%cmd%rmv = win_remove
      OS%cmd%cpy = win_copy
      OS%cmd%mkd = win_mkdir
    case default
      rks = 'rank'//trim(str(.true.,0_I4P)) ; if (present(myrank)) rks = 'rank'//trim(str(.true.,myrank))
      write(stderr,'(A)')trim(rks)//' Attention: OS id not recognized!'
      write(stderr,'(A)')trim(rks)//' Valid charcter id are:'
      write(stderr,'(A)')trim(rks)//' "'//c_uix_id//'" for *nix OS'
      write(stderr,'(A)')trim(rks)//' "'//c_win_id//'" for Windows OS'
      write(stderr,'(A)')trim(rks)//' The OS is set to unknow. All OS-related calls are disabled!'
    endselect
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  !> @brief The basedir function strip last component from file name providing the base directory name.
  !> @note The leading and trealing spaces are removed from the file name.
  elemental function basedir(OS,filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS),   intent(IN):: OS       !< OS.
  character(len=*), intent(IN):: filename !< File name.
  character(len=len(filename)):: basedir  !< Name of base directory.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (OS%id/=0) then
    basedir = trim(adjustl(filename(1:index(filename,OS%sep,back=.true.))))
  else
    basedir = filename
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction basedir

  !> @brief The basename function strip directory name from file name.
  !> @note The leading and trealing spaces are removed from the file name.
  elemental function basename(OS,filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS),   intent(IN):: OS       !< OS.
  character(len=*), intent(IN):: filename !< File name.
  character(len=len(filename)):: basename !< Base name of file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (OS%id/=0) then
    basename = trim(adjustl(filename(index(filename,OS%sep)+1:)))
  else
    basename = filename
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction basename

  !> @brief Function that fixes the separator character eventually present into a string converting to the correct OS's separators.
  function string_separator_fix(OS,myrank,string) result(string_fixed)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS),   intent(IN)::           OS           !< OS.
  integer(I4P),     intent(IN), optional:: myrank       !< Actual rank process necessary for concurrent multi-processes calls.
  character(len=*), intent(IN)::           string       !< String to be converted.
  character(len=len(string))::             string_fixed !< Converted string.
  integer(I_P)::                           n1           !< Characters counter.
  character(DI4P)::                        rks          !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string_fixed = string
  select case(OS%id)
  case(uix_id)
    do n1=1,len_trim(string)
      if (string(n1:n1)==win_sep) then
        string_fixed(n1:n1)=uix_sep ! conversion
      endif
    enddo
  case(win_id)
    do n1=1,len_trim(string)
      if (string(n1:n1)==uix_sep) then
        string_fixed(n1:n1)=win_sep ! conversion
      endif
    enddo
  case default
    rks = 'rank'//trim(str(.true.,0_I4P)) ; if (present(myrank)) rks = 'rank'//trim(str(.true.,myrank))
    write(stderr,'(A)')trim(rks)//' Attention: OS unknown. No seperators fix is performed!'
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction string_separator_fix

  !> @brief Function for executing shell commands. This is a tentative to make a portable procedure until the Fortran 2008 standard
  !> is completely implemented.
  function execute_cmd(command) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: command !< Name of file that must be copied.
  integer(I4P)::             err     !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef FSTD2008
  call execute_command_line(command=command,exitstat=err) ! Fortran 2008 standard procedure.
#else
  call system(command) ; err = 0_I4P ! Non standard (but widely supported) system call; it is assumed a successful execution.
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction execute_cmd

  !> @brief Function for copying a file.
  !> @note The leading and trealing spaces are removed from the files names.
  function copy_file(OS,myrank,Nproc,source_file,target_file) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS), intent(IN)::           OS          !< OS.
  integer(I4P),   intent(IN), optional:: myrank      !< Actual rank process necessary for concurrent multi-processes calls.
  integer(I4P),   intent(IN), optional:: Nproc       !< Number of MPI processes used.
  character(*),   intent(IN)::           source_file !< Name of file that must be copied.
  character(*),   intent(IN)::           target_file !< Destination path.
  integer(I4P)::                         err         !< Error trapping flag: 0 no errors, >0 error occurs.
  logical::                              exist       !< Inquiring flag.
  character(DI4P)::                      rks         !< String containing myrank.
  integer(I4P)::                         rank,Np     !< Dummy temporary variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (OS%id/=0) then
    inquire(file=adjustl(trim(source_file)),exist=exist,iostat=err)                       ! Verifing file presence.
    if (exist) then                                                                       ! Source file found.
      inquire(file=adjustl(trim(target_file)),exist=exist,iostat=err)                     ! Verifing target file.
      if (exist) then                                                                     ! Target file present.
        err=3_I4P                                                                         ! Updating err.
      else                                                                                ! Target file not found.
        err = OS%execute_cmd(trim(OS%cmd%cpy)//" "//&
                             adjustl(trim(source_file))//" "//adjustl(trim(target_file))) ! Coping file.
        inquire(file=adjustl(trim(target_file)),exist=exist,iostat=err)                   ! Verifing copied file.
      endif
    endif
  else
    err=1_I4P
    rank = 0 ; if (present(myrank)) rank = myrank ; Np = 1 ; if (present(Nproc)) Np = Nproc ; rks = 'rank'//trim(strz(Np,rank))
    write(stderr,'(A)')trim(rks)//' Attention: OS unknown. The copy_file statement cannot be executed!'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction copy_file

  !> @brief Function for removing a file.
  !> @note The leading and trealing spaces are removed from the file name.
  function remove_file(OS,myrank,Nproc,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS), intent(IN)::           OS       !< OS.
  integer(I4P),   intent(IN), optional:: myrank   !< Actual rank process necessary for concurrent multi-processes calls.
  integer(I4P),   intent(IN), optional:: Nproc    !< Number of MPI processes used.
  character(*),   intent(IN)::           filename !< Name of file that must be copied.
  integer(I4P)::                         err      !< Error trapping flag: 0 no errors, >0 error occurs.
  logical::                              exist    !< Inquiring flag.
  character(DI4P)::                      rks      !< String containing myrank.
  integer(I4P)::                         rank,Np     !< Dummy temporary variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (OS%id/=0) then
    inquire(file=adjustl(trim(filename)),exist=exist,iostat=err)           ! Verifing the presence of file.
    if (exist) then                                                        ! File found.
      err = OS%execute_cmd(trim(OS%cmd%rmv)//" "//adjustl(trim(filename))) ! Removing file.
      inquire(file=adjustl(trim(filename)),exist=exist,iostat=err)         ! Verifing the presence of file.
      if (exist) then                                                      ! File not removed.
        err=2_I4P                                                          ! Updating err.
      else                                                                 ! File removed.
        err=0_I4P                                                          ! Updating err.
      endif
    endif
  else
    err=1_I4P
    rank = 0 ; if (present(myrank)) rank = myrank ; Np = 1 ; if (present(Nproc)) Np = Nproc ; rks = 'rank'//trim(strz(Np,rank))
    write(stderr,'(A)')trim(rks)//' Attention: OS unknown. The remove_file statement cannot be executed!'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction remove_file

  !> @brief Function for creating a directory if it doesn't exist.
  !> @note The leading and trealing spaces are removed from the directory name.
  function make_dir(OS,myrank,Nproc,directory) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_OS), intent(IN)::           OS        !< OS.
  integer(I4P),   intent(IN), optional:: myrank    !< Actual rank process necessary for concurrent multi-processes calls.
  integer(I4P),   intent(IN), optional:: Nproc     !< Number of MPI processes used.
  character(*),   intent(IN)::           directory !< Name of the directory that must be created.
  integer(I4P)::                         err       !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DI4P)::                      rks       !< String containing myrank.
  integer(I4P)::                         rank,Np   !< Dummy temporary variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (OS%id/=0) then
    err = OS%execute_cmd(OS%cmd%mkd//" "//adjustl(trim(directory)))
  else
    err = 1_I4P
    rank = 0 ; if (present(myrank)) rank = myrank ; Np = 1 ; if (present(Nproc)) Np = Nproc ; rks = 'rank'//trim(strz(Np,rank))
    write(stderr,'(A)')trim(rks)//' Attention: OS unknown. The make_dir statement cannot be executed!'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction make_dir
  !> @}
endmodule Data_Type_OS
