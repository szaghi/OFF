!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_BaseDerivedType Data_Type_File_Base
!> Module definition of Type_File_Base
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_BasePrivateProcedure Data_Type_File_Base
!> Module definition of Type_File_Base
!> @}

!> @brief Module Data_Type_File_Base contains the definition of Type_File_Base, that is an abstract base type upon which the file
!> structures adopted by OFF project are constructed.
!> @note Type_File_Base is an abstract type, thus it cannot be directly used. It must be extendend in order to implement a real
!> file structure.
module Data_Type_File_Base
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                   ! Integers and reals precision definition.
USE Data_Type_OS, only: Type_OS    ! Definition of Type_OS.
USE Lib_IO_Misc                    ! Procedures for IO and strings operations.
USE Lib_Strings,  only: Upper_Case ! Library for strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Public parameters.
character(1), public, parameter:: cr                   = achar(10) !< Carriage-return character for stream-record finalizing.
integer(I4P), public, parameter:: err_unkw             = 10100     !< Unknown error.
integer(I4P), public, parameter:: err_not_LE           = 10101     !< Not Little Endian architecture error ID.
integer(I4P), public, parameter:: err_not_BE           = 10102     !< Not Big Endian architecture error ID.
integer(I4P), public, parameter:: err_not_connected    = 10103     !< File not connected error ID.
integer(I4P), public, parameter:: err_connected        = 10104     !< File connected error ID.
integer(I4P), public, parameter:: err_unkw_open_action = 10105     !< Unknown open action.
integer(I4P), public, parameter:: err_locked           = 10106     !< Present lockfile.
integer(I4P), public, parameter:: err_Nproc_unmatch    = 10107     !< Number of MPI processes inconsistent.
integer(I4P), public, parameter:: err_not_tecio        = 10108     !< Tecplot IO library unavailable.
integer(I4P), public, parameter:: err_gnu_binary       = 10109     !< Impossible to save Gnuplot binary file.
integer(I4P), public, parameter:: err_unkw_OS          = 10110     !< Unknown OS: impossible backup file.
!> @brief Derived type containing the definition of Type_File_Base.
!> @ingroup Data_Type_File_BaseDerivedType
type, public, abstract:: Type_File_Base
  ! public members
  character(len=:),    public, allocatable::  path_in               !< Input path.
  character(len=:),    public, allocatable::  path_out              !< Output path.
  character(len=:),    public, allocatable::  name                  !< Name of file.
  integer(I4P),        public::               unit      = 0_I4P     !< Logic unit of file.
  integer(I4P),        public::               iostat    = 0_I4P     !< IO error.
  character(500),      public::               iomsg                 !< IO error message.
  character(len=:),    public,  allocatable:: errpref               !< Prefixing string for stderr outputs.
  character(len=:),    public,  allocatable:: outpref               !< Prefixing string for stdout outputs.
  procedure(), nopass, public, pointer::      final_call => null()  !< Customized procedure to be finally called.
  ! private members
  logical,           private::              connected = .false.   !< Flag for inquiring the file connection (open) to unit file.
  logical,           private::              hd_passed = .false.   !< Flag for inquiring the (over)passing of the file header.
  integer(I8P),      private::              fout      = 1_I8P     !< Output frequency.
  character(len=:),  private, allocatable:: final_act             !< Final action to be executed as fallback.
  logical,           private::              initialized = .false. !< Flag for checking if the file data have been initialized.
  contains
    procedure:: raise_error                        ! Procedure for raising an error.
    procedure:: free_base  => free_file_base       ! Procedure for freeing dynamic memory.
    procedure:: init_base  => init_file_base       ! Procedure for initializing file structure.
    procedure:: set        => set_file_base        ! Procedure for setting file structure.
    procedure:: open       => open_file_base       ! Procedure for opening file.
    procedure:: close      => close_file_base      ! Procedure for closing file.
    procedure:: backup     => backup_file_base     ! Procedure for creating a backup copy of file.
    procedure:: get_stream => get_stream_file_base ! Procedure for getting contents of file as single stream.
    procedure:: fallback   => fallback_file_base   ! Procedure for a sicure fallback when severe error occurs.
    procedure:: is_to_save                         ! Procedure for testing if is time to save file accordingly to fout.
    procedure:: is_severe_error                    ! Procedure for testing if a severe error is occurred.
endtype Type_File_Base
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_BasePrivateProcedure
  !> @{
  !> @brief Procedure raising an error.
  subroutine raise_error(file_d,errtype)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base), intent(INOUT):: file_d   !< File data.
  integer(I4P),          intent(IN)::    errtype  !< Type of error arised.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(file_d%errpref)) file_d%errpref = ''
  select case(errtype)
  case(err_not_LE)
    file_d%iomsg  = file_d%errpref//' Attention: trying to load Little Endian file on Big Endian architecture!'
    file_d%iostat = err_not_LE
  case(err_not_BE)
    file_d%iomsg  = file_d%errpref//' Attention: trying to load Big Endian file on Little Endian architecture!'
    file_d%iostat = err_not_BE
  case(err_not_connected)
    file_d%iomsg  = file_d%errpref//' Attention: file "'//trim(adjustl(file_d%name))//'" is not connected!'
    file_d%iostat = err_not_connected
  case(err_connected)
    file_d%iomsg  = file_d%errpref//' Attention: file "'//trim(adjustl(file_d%name))//'" is already connected!'
    file_d%iostat = err_connected
  case(err_unkw_open_action)
    file_d%iomsg  = file_d%errpref//' Attention: unknown open action, valid actions are "READ" and "WRITE"!'
    file_d%iostat = err_unkw_open_action
  case(err_locked)
    file_d%iomsg  = file_d%errpref//' Attention: lockfile has been found into the working directory!'
    file_d%iostat = err_locked
  case(err_Nproc_unmatch)
    file_d%iomsg  = file_d%errpref//' Attention: Nproc defined into procmap file differs from the one used!'
    file_d%iostat = err_Nproc_unmatch
  case(err_not_tecio)
    file_d%iomsg  = file_d%errpref//' Attention: impossible to save binary tecplot file without Tecplot tecio library!'
    file_d%iostat = err_not_tecio
  case(err_gnu_binary)
    file_d%iomsg  = file_d%errpref//' Attention: impossible to save Gnuplot binary file!'
    file_d%iostat = err_gnu_binary
  case(err_unkw_OS)
    file_d%iomsg  = file_d%errpref//' Attention: Unknown OS type. Some operations are impossible (e.g. backup file)!'
    file_d%iostat = err_unkw_OS
  case default
    file_d%iomsg  = file_d%errpref//' Attention: unknown error!'
    file_d%iostat = err_unkw
  endselect
#ifdef DEBUG
  write(stderr,'(A)')file_d%errpref//trim(file_d%iomsg)
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine raise_error

  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_file_base(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(file_d%path_in )) deallocate(file_d%path_in )
  if (allocated(file_d%path_out)) deallocate(file_d%path_out)
  if (allocated(file_d%name    )) deallocate(file_d%name    )
  if (allocated(file_d%errpref )) deallocate(file_d%errpref )
  if (allocated(file_d%outpref )) deallocate(file_d%outpref )
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_file_base

  !> @brief Procedure for setting Type_File_Base.
  elemental subroutine init_file_base(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base),  intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  file_d%path_in     = ''
  file_d%path_out    = ''
  file_d%errpref     = ''
  file_d%outpref     = ''
  file_d%final_act   = 'STOP'
  file_d%initialized = .true.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_file_base

  !> @brief Procedure for setting Type_File_Base.
  subroutine set_file_base(file_d,path_in,path_out,unit,connected,hd_passed,fout,errpref,outpref,final_act,iostat,iomsg,name,&
                           final_call)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base),           intent (INOUT):: file_d     !< File data.
  character(*), optional,          intent (IN)::    path_in    !< Input path.
  character(*), optional,          intent (IN)::    path_out   !< Output path.
  integer(I4P), optional,          intent (IN)::    unit       !< Logic unit of file.
  logical,      optional,          intent (IN)::    connected  !< Flag for inquiring the file connection (open) to unit file.
  logical,      optional,          intent (IN)::    hd_passed  !< Flag for inquiring the (over)passing of the file header.
  integer(I8P), optional,          intent (IN)::    fout       !< Output frequency.
  character(*), optional,          intent (IN)::    errpref    !< Prefixing string for stderr outputs.
  character(*), optional,          intent (IN)::    outpref    !< Prefixing string for stdout outputs.
  character(*), optional,          intent (IN)::    final_act  !< Final action to be executed as fallback.
  integer(I4P), optional,          intent (IN)::    iostat     !< IO error.
  character(*), optional,          intent (IN)::    iomsg      !< IO error message.
  character(*), optional,          intent (IN)::    name       !< Name of file.
  procedure(),  optional, pointer, intent (IN)::    final_call !< Customized procedure to be finally called.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.file_d%initialized) call file_d%init_base
  if (present(name      )) file_d%name       =  name
  if (present(path_in   )) file_d%path_in    =  path_in
  if (present(path_out  )) file_d%path_out   =  path_out
  if (present(unit      )) file_d%unit       =  unit
  if (present(connected )) file_d%connected  =  connected
  if (present(hd_passed )) file_d%hd_passed  =  hd_passed
  if (present(fout      )) file_d%fout       =  fout
  if (present(errpref   )) file_d%errpref    =  errpref
  if (present(outpref   )) file_d%outpref    =  outpref
  if (present(final_act )) file_d%final_act  =  final_act
  if (present(iostat    )) file_d%iostat     =  iostat
  if (present(iomsg     )) file_d%iomsg      =  iomsg
  if (present(final_call)) file_d%final_call => final_call
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_file_base

  !> @brief Procedure for opening file.
  !> @note Accordingly to the passed arguments different files are open:
  !>   * the trim(path_in/out)//name if neither "n" nor "flip" are passed;
  !>   * the trim(path_in/out)//name//'-N_\#n' if "n" is passed;
  !>   * the trim(path_in/out)//name//'-F_\#flip' if "flip" is passed.
  subroutine open_file_base(file_d,n,flip,ascii,append,replace,action)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base),  intent(INOUT):: file_d  !< File data.
  integer(I8P), optional, intent(IN)::    n       !< Time step number.
  integer(I1P), optional, intent(IN)::    flip    !< Flip-flop number.
  logical,      optional, intent(IN)::    ascii   !< Flag for inquiring the file format.
  logical,      optional, intent(IN)::    append  !< Flag for inquiring the file open position.
  logical,      optional, intent(IN)::    replace !< Flag for inquiring the file open replace status.
  character(*),           intent(IN)::    action  !< Action to be perfomed on file: "READ" or "WRITE".
  character(len=:), allocatable::         fname   !< Temporary variable for creating file name.
  character(len=:), allocatable::         acc     !< File access.
  character(len=:), allocatable::         frm     !< File format.
  character(len=:), allocatable::         pos     !< File position.
  character(len=:), allocatable::         stat    !< File status.
  character(len=:), allocatable::         act     !< File action.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  acc  = 'STREAM'      ; if (present(ascii  )) acc  = 'SEQUENTIAL'
  frm  = 'UNFORMATTED' ; if (present(ascii  )) frm  = 'FORMATTED'
  pos  = 'REWIND'      ; if (present(append )) pos  = 'APPEND'
  stat = 'UNKNOWN'     ; if (present(replace)) stat = 'REPLACE'
  act = trim(adjustl(Upper_Case(action)))
  if (file_d%connected) then
    call file_d%raise_error(errtype=err_connected)
    return
  endif
  ! building up complete file name according to the action
  select case(act)
  case('READ')
    fname = trim(adjustl(file_d%path_in))//trim(adjustl(file_d%name))
  case('WRITE')
    fname = trim(adjustl(file_d%path_out))//trim(adjustl(file_d%name))
  case default
    call file_d%raise_error(errtype=err_unkw_open_action)
  endselect
  ! building up complete file name according to the passed arguments
  if (present(n)) then
    fname = trim(fname)//'-N_'//trim(strz(10,n))
  elseif (present(flip)) then
    fname = trim(fname)//'-F_'//trim(strz(1,flip))
  endif
  ! in the case of READ action inquiring the existance of the file
  if (act=='READ') then
    call inquire_file(cpn='open_file_base',pref=file_d%errpref,iostat=file_d%iostat,iomsg=file_d%iomsg,file=trim(fname))
    if (file_d%iostat/=0) return
  endif
  ! opening the file
  open(unit=Get_Unit(file_d%unit),file=trim(fname),&
       access=acc,position=pos,status=stat,form=frm,action=act,iostat=file_d%iostat,iomsg=file_d%iomsg)
  file_d%connected = .true.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine open_file_base

  !> @brief Procedure for closing file.
  subroutine close_file_base(file_d,delete)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base), intent(INOUT):: file_d !< File data.
  logical, optional,     intent(IN)::    delete !< Status flag for deleting file when closing.
  character(len=:), allocatable::        status !< Status of closing.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  status = 'KEEP' ; if (present(delete)) status = 'DELETE'
  if (file_d%connected) then
    close(unit=file_d%unit,status=status,iostat=file_d%iostat,iomsg=file_d%iomsg)
    call file_d%set(unit=0_I4P,connected=.false.,hd_passed=.false.)
  else
    call file_d%raise_error(errtype=err_not_connected)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine close_file_base

  !> @brief Procedure for creating a backup copy of file. The orginal file is copied into output path.
  subroutine backup_file_base(file_d,OS)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base), intent(INOUT):: file_d !< File data.
  type(Type_OS),         intent(IN)::    OS     !< Running architecture.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (OS%id/=0) then ! OS statements are available
    file_d%iostat = OS%copy_file(source_file=trim(adjustl(file_d%path_in ))//trim(adjustl(file_d%name)),&
                                 target_file=trim(adjustl(file_d%path_out))//trim(adjustl(file_d%name)))
  else
    call file_d%raise_error(errtype = err_unkw_OS)
    call file_d%fallback
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine backup_file_base

  !> @brief Procedure for getting contents of file as single stream.
  subroutine  get_stream_file_base(file_d,delimiter_start,delimiter_end,stream)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base),         intent(INOUT):: file_d          !< File data.
  character(*), optional,        intent(IN)::    delimiter_start !< Delimiter from which start the stream.
  character(*), optional,        intent(IN)::    delimiter_end   !< Delimiter to which end the stream.
  character(len=:), allocatable, intent(OUT)::   stream          !< Output string containing the file data as a single stream.
  character(len=:), allocatable::                fname           !< Temporary variable for creating file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  fname = trim(adjustl(file_d%path_in))//trim(adjustl(file_d%name))
  if (present(delimiter_start).and.present(delimiter_end)) then
    call read_file_as_stream(pref=file_d%errpref,iostat=file_d%iostat,iomsg=file_d%iomsg,filename=trim(fname),stream=stream, &
                             delimiter_start=delimiter_start,delimiter_end=delimiter_end)
  else
    call read_file_as_stream(pref=file_d%errpref,iostat=file_d%iostat,iomsg=file_d%iomsg,filename=trim(fname),stream=stream)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get_stream_file_base

  !> @brief Procedure for a sicure fallback when severe error occurs.
  subroutine fallback_file_base(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base),          intent(IN):: file_d     !< File data.
#ifdef _MPI
  integer(I4P)::                               err        !< Error traping flag.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (file_d%is_severe_error()) then
    write(stderr,'(A)')file_d%errpref//file_d%iomsg
    if (associated(file_d%final_call)) then
      call file_d%final_call
    else
      select case(file_d%final_act)
      case('STOP')
#ifdef _MPI
        call MPI_FINALIZE(err)
#endif
        stop
      case('CONTINUE')
        ! do nothing
      case default
      endselect
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine fallback_file_base

  !> @brief Procedure for testing if is time to save file accordingly to fout.
  elemental function is_to_save(file_d,n) result(yes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base), intent(IN):: file_d !< File data.
  integer(I8P),          intent(IN):: n      !< Current time step.
  logical::                           yes    !< Is to save or not.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  yes=(mod(n,file_d%fout)==0)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_to_save

  !> @brief Procedure for testing if a severe error is occurred.
  elemental function is_severe_error(file_d) result(yes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Base), intent(IN):: file_d !< File data.
  logical::                           yes    !< Is to save or not.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  yes=((file_d%iostat/=0).and.(file_d%iostat/=iostat_end).and.(file_d%iostat/=iostat_eor))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_severe_error
  !> @}
endmodule Data_Type_File_Base
