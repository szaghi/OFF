!> @ingroup GlobalVarPar
!> @{
!> @defgroup Lib_IO_Misc Lib_IO_Misc
!> @}

!> This module contains miscellanea procedures for input/output and string operations.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Library
module Lib_IO_Misc
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
USE Data_Type_OS ! Definition of Type_OS.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: OS
public:: stdout,stderr
public:: Get_Unit
public:: make_dir
public:: remove_file
public:: copy_file
public:: lc_file
public:: basedir,basename
public:: Upper_Case
public:: Lower_Case
public:: string_OS_sep
public:: File_Not_Found
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Lib_IO_Misc
type(Type_OS):: OS !< OS definition (see \ref Data_Type_OS::Type_OS "definition").
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief The GetUnit function is used for getting a free logic unit.
  !>@return Free_Unit
  function Get_Unit() result(Free_Unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer:: Free_Unit !< Free logic unit.
  integer:: n1        !< Counter.
  integer:: ios       !< Inquiring flag.
  logical:: lopen     !< Inquiring flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Free_Unit = -1
  n1=1
  do
    if ((n1/=stdout).AND.(n1/=stderr)) then
      inquire (unit=n1,opened=lopen,iostat=ios)
      if (ios==0) then
        if (.NOT.lopen) then
          Free_Unit = n1
          return
        endif
      endif
    endif
    n1=n1+1
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Get_Unit

  function make_dir(myrank,directory) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!The make_dir function creates a directory if it doesn't exist.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),  intent(IN), optional:: myrank    ! Actual rank process necessary for concurrent multi-processes calls.
  character(*),  intent(IN)::           directory ! Name of the directory that must be created.
  integer(I4P)::                        err       ! Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                        UnitFree  ! Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  UnitFree = Get_Unit()                                                                  ! Getting free logic unit.
  if (present(myrank)) then
    open(unit=UnitFree,file=trim(directory)//'ExiSt.p'//trim(strz(3,myrank)),iostat=err) ! Creating a test file in the directory.
  else
    open(unit=UnitFree,file=trim(directory)//'ExiSt.p'//trim(strz(3,0     )),iostat=err) ! Creating a test file in the directory.
  endif
  close(UnitFree)                                                                        ! Close the test file.
  if (err/=0_I_P) then                                                                   ! Directory not found.
    call system(OS%mkdir//" "//directory(1:len(trim(directory))-1))                      ! Directory creation.
    err = 0_I4P                                                                          ! Updating err.
  else                                                                                   ! Directory present; removing test file.
    if (present(myrank)) then
      err = remove_file(trim(directory)//'ExiSt.p'//trim(strz(3,myrank)))                ! Removing file test.
    else
      err = remove_file(trim(directory)//'ExiSt.p'//trim(strz(3,0     )))                ! Removing file test.
    endif
    err = 0_I4P                                                                          ! Updating err.
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction make_dir

  function remove_file(filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!The remove_file function removes file in a portable way.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),  intent(IN):: filename ! Name of file that must be copied.
  integer(I4P)::              err      ! Error trapping flag: 0 no errors, >0 error occurs.
  logical(4)::                exist    ! Inquiring flag.
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

  function copy_file(source_file,target_file) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!The copy_file function copies file to a destination path in a portable way.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),  intent(IN):: source_file ! Name of file that must be copied.
  character(*),  intent(IN):: target_file ! Destination path.
  integer(I4P)::              err         ! Error trapping flag: 0 no errors, >0 error occurs.
  logical(4)::                exist       ! Inquiring flag.
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

  function lc_file(filename,fileform) result(n)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! The lc_file function calculates the number of lines (records) of a sequential file.
  !---------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN):: filename ! File name.
  character(*), intent(IN):: fileform ! File format: FORMATTED or UNFORMATTED.
  integer(I4P)::             n        ! Number of lines (records).
  logical(4)::               is_file  ! Inquiring flag.
  integer(I4P)::             unitfile ! File logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file)  ! Verifing the presence of the file.
  if (.NOT.is_file) then                               ! File not found.
    n = -1_I4P                                         ! Returning -1.
  else                                                 ! File found.
    n = 0_I4P                                          ! Initializing number of records.
    unitfile = Get_Unit()                              ! Getting free logic unit.
    open(unit   = unitfile,                &
         file   = adjustl(trim(filename)), &
         action = 'READ',                  &
         form   = adjustl(trim(Upper_Case(fileform)))) ! Opening file.
    select case(adjustl(trim(Upper_Case(fileform))))
    case('FORMATTED')
      do
        read(unitfile,*,end=10)                        ! Reading record.
        n = n + 1_I4P                                  ! Updating number of records.
      enddo
    case('UNFORMATTED')
      do
        read(unitfile,end=10)                          ! Reading record.
        n = n + 1_I4P                                  ! Updating number of records.
      enddo
    endselect
    10 continue
    close(unitfile)                                    ! Closing file.
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lc_file

  function basedir(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! The basedir function strip last component from file name providing the base directory name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: filename ! File name.
  character(len=len(filename)):: basedir  ! Name of base directory.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  basedir = adjustl(trim(filename(1:index(filename,OS%sep,back=.true.))))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction basedir

  function basename(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! The basename function strip directory from file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: filename ! File name.
  character(len=len(filename)):: basename ! Base name of file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  basename = adjustl(trim(filename(index(filename,OS%sep)+1:)))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction basename

  !> @brief The Upper_Case function converts the lower case characters of a string to upper case one.
  !>@return Upper_Case
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
      Upper_Case(n1:n1)=char(ichar(string(n1:n1))-32) ! Upper case conversion
    endselect
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Upper_Case

  function Lower_Case(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!The Lower_Case function converts the upper case characters of a string to lower case one. Use this function in
  !!order to achieve case-insensitive.
  !---------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=*), intent(IN):: string     ! String to be converted.
  character(len=len(string))::   Lower_Case ! Converted string.
  integer::                      n1         ! Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !!The following is the code snippet of Upper_Case function.
  !!
  !(\doc)codesnippet
  Lower_Case = string
  do n1=1,len_trim(string)
    select case(ichar(string(n1:n1)))
    case(65:90)
      Lower_Case(n1:n1)=char(ichar(string(n1:n1))+32) ! conversion
    endselect
  enddo
  return
  !(doc/)codesnippet
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Lower_Case

  function string_OS_sep(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!The string_OS_sep function converts directories separators into a string to the correct OS's separators.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: string        ! String to be converted.
  character(len=len(string))::   string_OS_sep ! Converted string.
  integer(I_P)::                 n1            ! Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !!The following is the code snippet of string_OS_sep function.
  !!
  !(\doc)codesnippet
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
  !(doc/)codesnippet
  !---------------------------------------------------------------------------------------------------------------------------------
  return
  endfunction string_OS_sep

  subroutine File_Not_Found(myrank,filename,cpn)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for printing to stderr a "file not found error".
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  myrank   ! Actual rank process.
  character(*), intent(IN)::  filename ! Name of file where option variables are saved.
  character(*), intent(IN)::  cpn      ! Calling procedure name.
  integer(I_P)::              err      ! Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stderr,'(A,'//FI_P//')',iostat=err)' My RANK is: ',myrank
  write(stderr,'(A)',           iostat=err)' File '//adjustl(trim(filename))//' Not Found!'
  write(stderr,'(A)',           iostat=err)' Calling procedure "'//adjustl(trim(cpn))//'"'
#ifdef MPI2
    call MPI_FINALIZE(err)
#endif
  stop
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine File_Not_Found
endmodule Lib_IO_Misc
