!> @ingroup Library
!> @{
!> @defgroup Lib_IO_MiscLibrary Lib_IO_Misc
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Lib_IO_MiscGlobalVarPar Lib_IO_Misc
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_IO_MiscPublicProcedure Lib_IO_Misc
!> @}

!> This module contains miscellanea procedures for input/output and strings operations.
!> This is a library module.
!> @ingroup Lib_IO_MiscLibrary
module Lib_IO_Misc
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                    ! Integers and reals precision definition.
USE Data_Type_OS                                                                    ! Definition of Type_OS.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
#ifdef MPI2
USE MPI                                                                             ! MPI runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: OS
public:: stdout,stderr
public:: Get_Unit
public:: remove_file
public:: copy_file
public:: lc_file
public:: File_Not_Found
public:: basedir,basename
public:: make_dir
public:: Upper_Case
public:: Lower_Case
public:: string_OS_sep
public:: tokenize
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Lib_IO_MiscGlobalVarPar
type(Type_OS):: OS !< OS definition (see \ref Data_Type_OS::Type_OS "definition").
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_IO_MiscPublicProcedure
  !> @{
  !> @brief The Get_Unit function returns a free logic unit for opening a file. The unit value is returned by the function, and also
  !> by the optional argument "Free_Unit". This allows the function to be used directly in an open statement like:
  !> open(unit=Get_Unit(myunit),...) ; read(myunit)...
  !> If no units are available, -1 is returned.
  integer function Get_Unit(Free_Unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer, intent(OUT), optional:: Free_Unit !< Free logic unit.
  integer::                        n1        !< Counter.
  integer::                        ios       !< Inquiring flag.
  logical::                        lopen     !< Inquiring flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Get_Unit = -1
  n1=1
  do
    if ((n1/=stdout).AND.(n1/=stderr)) then
      inquire (unit=n1,opened=lopen,iostat=ios)
      if (ios==0) then
        if (.NOT.lopen) then
          Get_Unit = n1 ; if (present(Free_Unit)) Free_Unit = Get_Unit
          return
        endif
      endif
    endif
    n1=n1+1
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Get_Unit

  !> @brief Function for creating a directory if it doesn't exist.
  !> @return \b err integer(I_P) variable for error trapping.
  function make_dir(myrank,directory) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN), optional:: myrank    !< Actual rank process necessary for concurrent multi-processes calls.
  character(*),  intent(IN)::           directory !< Name of the directory that must be created.
  integer(I4P)::                        err       !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                        UnitFree  !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Creating a test file in the directory.
  if (present(myrank)) then
    open(unit=Get_Unit(UnitFree),file=trim(directory)//'ExiSt.p'//trim(strz(3,myrank)),iostat=err)
  else
    open(unit=Get_Unit(UnitFree),file=trim(directory)//'ExiSt.p'//trim(strz(3,0     )),iostat=err)
  endif
  close(UnitFree)
  if (err/=0_I_P) then                                                     ! Directory not found.
    call system(OS%mkdir//" "//directory(1:len(trim(directory))-1))        ! Directory creation.
    err = 0_I4P                                                            ! Updating err.
  else                                                                     ! Directory present; removing test file.
    if (present(myrank)) then
      err = remove_file(trim(directory)//'ExiSt.p'//trim(strz(3,myrank)))  ! Removing file test.
    else
      err = remove_file(trim(directory)//'ExiSt.p'//trim(strz(3,0     )))  ! Removing file test.
    endif
    err = 0_I4P                                                            ! Updating err.
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction make_dir

  !> @brief Function for removing a file.
  !> @return \b err integer(I_P) variable for error trapping.
  function remove_file(filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),  intent(IN):: filename !< Name of file that must be copied.
  integer(I4P)::              err      !< Error trapping flag: 0 no errors, >0 error occurs.
  logical(4)::                exist    !< Inquiring flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=exist)            ! Verifing the presence of file.
  if (.NOT.exist) then                                         ! File not found.
    err=1_I_P                                                  ! Updating err.
  else                                                         ! File found.
    call system(trim(OS%remove)//" "//adjustl(trim(filename))) ! Removing file.
    inquire(file=adjustl(trim(filename)),exist=exist)          ! Verifing the presence of file.
    if (exist) then                                            ! File not removed.
      err=2_I4P                                                ! Updating err.
    else                                                       ! File removed.
      err=0_I4P                                                ! Updating err.
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction remove_file

  !> @brief Function for coping a file.
  !> @return \b err integer(I_P) variable for error trapping.
  function copy_file(source_file,target_file) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),  intent(IN):: source_file !< Name of file that must be copied.
  character(*),  intent(IN):: target_file !< Destination path.
  integer(I4P)::              err         !< Error trapping flag: 0 no errors, >0 error occurs.
  logical(4)::                exist       !< Inquiring flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(source_file)),exist=exist)                                             ! Verifing the presence of file.
  if (.NOT.exist) then                                                                             ! Source file not found.
    err=2_I4P                                                                                      ! Updating err.
  else                                                                                             ! Source file found.
    inquire(FILE=adjustl(trim(target_file)),EXIST=exist)                                           ! Verifing target file.
    if (exist) then                                                                                ! Target file already present.
      err=3_I4P                                                                                    ! Updating err.
    else                                                                                           ! Target file not found.
      call system(trim(OS%copy)//" "//adjustl(trim(source_file))//" "//adjustl(trim(target_file))) ! Coping file.
      inquire(FILE=adjustl(trim(target_file)),EXIST=exist)                                         ! Verifing copied file.
      if (.NOT.exist) then                                                                         ! File not copied.
        err=1_I4P                                                                                  ! Updating err.
      else                                                                                         ! File copied.
        err=0_I4P                                                                                  ! Updating err.
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction copy_file

  !> @brief Function for calculating the number of lines (records) of a sequential file.
  !>@return \b n integer(I4P) variable
  function lc_file(filename) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: filename ! File name.
  integer(I4P)::             n        ! Number of lines (records).
  logical(4)::               is_file  ! Inquiring flag.
  character(11)::            fileform ! File format: FORMATTED or UNFORMATTED.
  integer(I4P)::             unitfile ! Logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,form=fileform) ! Verifing the presence of the file.
  if (.NOT.is_file) then                                            ! File not found.
    n = -1_I4P                                                      ! Returning -1.
  else                                                              ! File found.
    n = 0_I4P                                                       ! Initializing number of records.
    open(unit   = Get_Unit(unitfile),      &
         file   = adjustl(trim(filename)), &
         action = 'READ',                  &
         form   = adjustl(trim(Upper_Case(fileform))))              ! Opening file.
    select case(adjustl(trim(Upper_Case(fileform))))
    case('FORMATTED')
      do
        read(unitfile,*,end=10)                                     ! Reading record.
        n = n + 1_I4P                                               ! Updating number of records.
      enddo
    case('UNFORMATTED')
      do
        read(unitfile,end=10)                                       ! Reading record.
        n = n + 1_I4P                                               ! Updating number of records.
      enddo
    endselect
    10 continue
    close(unitfile)                                                 ! Closing file.
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lc_file

  !> @brief The basedir function strip last component from file name providing the base directory name.
  !>@return \b basedir character(*) variable
  function basedir(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: filename !< File name.
  character(len=len(filename)):: basedir  !< Name of base directory.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  basedir = adjustl(trim(filename(1:index(filename,OS%sep,back=.true.))))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction basedir

  !> @brief The basename function strip directory from file name.
  !>@return \b basename character(*) variable
  function basename(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: filename !< File name.
  character(len=len(filename)):: basename !< Base name of file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  basename = adjustl(trim(filename(index(filename,OS%sep)+1:)))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction basename

  !> @brief The Upper_Case function converts the lower case characters of a string to upper case one.
  !>@return \b Upper_Case character(*) variable
  function Upper_Case(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: string     !< String to be converted.
  character(len=len(string))::   Upper_Case !< Converted string.
  integer::                      n1         !< Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Upper_Case = string
  do n1=1,len(string)
    select case(ichar(string(n1:n1)))
    case(97:122)
      Upper_Case(n1:n1)=char(ichar(string(n1:n1))-32) ! upper case conversion
    endselect
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Upper_Case

  !> @brief The Lower_Case function converts the upper case characters of a string to lower case one.
  !>@return \b Lower_Case character(*) variable
  function Lower_Case(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: string     !< String to be converted.
  character(len=len(string))::   Lower_Case !< Converted string.
  integer::                      n1         !< Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Lower_Case = string
  do n1=1,len_trim(string)
    select case(ichar(string(n1:n1)))
    case(65:90)
      Lower_Case(n1:n1)=char(ichar(string(n1:n1))+32) ! lower case conversion
    endselect
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Lower_Case

  !> @brief The string_OS_sep function converts directories separators into a string to the correct OS's separators.
  !>@return \b string_OS_sep character(*) variable
  function string_OS_sep(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: string        ! String to be converted.
  character(len=len(string))::   string_OS_sep ! Converted string.
  integer(I_P)::                 n1            ! Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string_OS_sep = string
  select case(OS%id)
  case(uix_id)
    do n1=1,len_trim(string)
      if (string(n1:n1)==win_sep) then
        string_OS_sep(n1:n1)=uix_sep ! conversion
      endif
    enddo
  case(win_id)
    do n1=1,len_trim(string)
      if (string(n1:n1)==uix_sep) then
        string_OS_sep(n1:n1)=win_sep ! conversion
      endif
    enddo
  endselect
  !---------------------------------------------------------------------------------------------------------------------------------
  return
  endfunction string_OS_sep

  !> @brief Subroutine for tokenizing a string in order to parse it.
  !> @note The dummy array containing tokens must allocatable and its character elements must have the same lenght of the input
  !> string. If the length of the delimiter is higher than the input string one then the output tokens array is allocated with
  !> only one element set to char(0).
  subroutine tokenize(strin,delimiter,Nt,toks)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*),          intent(IN)::               strin     !< String to be tokenized.
  character(len=*),          intent(IN)::               delimiter !< Delimiter of tokens.
  integer(I4P),              intent(OUT), optional::    Nt        !< Number of tokens.
  character(len=len(strin)), intent(OUT), allocatable:: toks(:)   !< Tokens.
  character(len=len(strin))::                           strsub    !< Temporary string.
  integer(I4P)::                                        dlen      !< Delimiter length.
  integer(I4P)::                                        c,n,t     !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initialization
  strsub = strin
  dlen = len(delimiter)
  if (dlen>len(strin)) then
    if (allocated(toks)) deallocate(toks) ; allocate(toks(1:1)) ; toks(1) = char(0) ; if (present(Nt)) Nt = 1
    return
  endif
  ! computing the number of tokens
  n = 1
  do c=1,len(strsub)-dlen ! loop over string characters
    if (strsub(c:c+dlen-1)==delimiter) n = n + 1
  enddo
  if (allocated(toks)) deallocate(toks) ; allocate(toks(1:n))
  ! tokenization
  do t=1,n ! loop over tokens
    c = index(strsub,delimiter)
    if (c>0) then
      toks(t) = strsub(1:c-1)
      strsub = strsub(c+dlen:)
    else
      toks(t) = strsub
    endif
  enddo
  if (present(Nt)) Nt = n
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine tokenize

  !> @brief Subroutine for printing to stderr a "file not found error".
  subroutine File_Not_Found(myrank,filename,cpn)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN), optional:: myrank   !< Actual rank process.
  character(*), intent(IN)::           filename !< Name of file where option variables are saved.
  character(*), intent(IN)::           cpn      !< Calling procedure name.
  character(DI_P)::                    rks      !< String containing myrank.
  integer(I_P)::                       err      !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,0_I_P)) ; if (present(myrank)) rks = 'rank'//trim(str(.true.,myrank))
  write(stderr,'(A)',iostat=err)trim(rks)//' File '//adjustl(trim(filename))//' Not Found!'
  write(stderr,'(A)',iostat=err)trim(rks)//' Calling procedure "'//adjustl(trim(cpn))//'"'
#ifdef MPI2
    call MPI_FINALIZE(err)
#endif
  stop
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine File_Not_Found
  !> @}
endmodule Lib_IO_Misc
