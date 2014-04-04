!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Command_Line_ArgumentDerivedType Data_Type_Command_Line_Argument
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_Command_Line_ArgumentInterface Data_Type_Command_Line_Argument
!> Module definition of Type_Command_Line_Argument
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Command_Line_ArgumentPrivateProcedure Data_Type_Command_Line_Argument
!> Module definition of Type_Command_Line_Argument
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_Command_Line_ArgumentPublicProcedure Data_Type_Command_Line_Argument
!> Module definition of Type_Command_Line_Argument
!> @}

!> @brief This module contains the definition of Type_Command_Line_Argument and its procedures.
!> Type_Command_Line_Argument (CLA) is a derived type containing the useful data for handling command line arguments in order to
!> easy implement flexible a Command Line Interface (CLI).
!> @note Presently there is no support for positional CLAs, but only for named ones.
!> @note Presently there is no support for multiple valued CLAs, but only for single valued ones (or without any value, i.e. logical
!> CLA).
!> @todo Add support for positional CLAs.
!> @todo Add support for multiple valued (list of values) CLAs.
module Data_Type_Command_Line_Argument
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
USE Lib_IO_Misc  ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: cli_parser
public:: cli_print
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
integer(I4P),  parameter:: max_val_len        = 1000          !< Maximum number of characters of CLA value.
character(5),  parameter:: action_store       = 'STORE'       !< CLA that stores a value associated to its switch.
character(10), parameter:: action_store_true  = 'STORE_TRUE'  !< CLA that stores .true. without the necessity of a value.
character(11), parameter:: action_store_false = 'STORE_FALSE' !< CLA that stores .false. without the necessity of a value.
!> Derived type containing the useful data for handling command line arguments in order to easy implement flexible a Command Line
!> Interface (CLI).
!> @note If not otherwise declared the action on CLA value is set to "store" a value that must be passed after the switch name.
!> @ingroup Data_Type_Command_Line_ArgumentDerivedType
type, public:: Type_Command_Line_Argument
  character(len=:), allocatable:: switch           !< Switch name.
  character(len=:), allocatable:: switch_ab        !< Abbreviated switch name.
  character(len=:), allocatable:: help             !< Help message describing the CLA.
  logical::                       required=.false. !< Flag for set required argument.
  logical::                       passed  =.false. !< Flag for checking if CLA has been passed to CLI.
  character(len=:), allocatable:: act              !< CLA value action.
  character(len=:), allocatable:: def              !< Default value.
  character(len=:), allocatable:: val              !< CLA value.
  contains
    procedure:: free     ! Procedure for freeing dynamic memory.
    procedure:: init     ! Procedure for initializing CLA.
    procedure:: get      ! Procedure for getting CLA value after it has been corrected parsed from CLI parser.
    final::     finalize ! Procedure for freeing dynamic memory when finalizing.
endtype Type_Command_Line_Argument
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Command_Line_ArgumentPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free(cla)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Command_Line_Argument), intent(INOUT):: cla !< CLA data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cla%switch   )) deallocate(cla%switch   )
  if (allocated(cla%switch_ab)) deallocate(cla%switch_ab)
  if (allocated(cla%help     )) deallocate(cla%help     )
  if (allocated(cla%act      )) deallocate(cla%act      )
  if (allocated(cla%def      )) deallocate(cla%def      )
  if (allocated(cla%val      )) deallocate(cla%val      )
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(cla)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Command_Line_Argument), intent(INOUT):: cla !< CLA data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call cla%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for initializing CLA.
  !> @note If not otherwise declared the action on CLA value is set to "store" a value that must be passed after the switch name.
  subroutine init(cla,pref,switch_ab,help,required,act,def,switch)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Command_Line_Argument), intent(INOUT):: cla       !< CLA data.
  character(*), optional,            intent(IN)::    pref      !< Prefixing string.
  character(*), optional,            intent(IN)::    switch_ab !< Abbreviated switch name.
  character(*), optional,            intent(IN)::    help      !< Help message describing the CLA.
  logical,      optional,            intent(IN)::    required  !< Flag for set required argument.
  character(*), optional,            intent(IN)::    act       !< CLA value action.
  character(*), optional,            intent(IN)::    def       !< Default value.
  character(*),                      intent(IN)::    switch    !< Switch name.
  character(len=:), allocatable::                    prefd     !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
                         cla%switch   =switch
                         cla%switch_ab=switch                 ;if (present(switch_ab)) cla%switch_ab=switch_ab
                         cla%help     ='Undocumented argument';if (present(help     )) cla%help     =help
  if (present(required)) cla%required =required
                         cla%act      =action_store           ;if (present(act      )) cla%act      =trim(adjustl(Upper_Case(act)))
                         cla%def      =''                     ;if (present(def      )) cla%def      =def
  if ((.not.cla%required).and.(.not.present(def))) then
    prefd = '' ; if (present(pref)) prefd = pref
    write(stderr,'(A)')prefd//' Error: the CLA "'//cla%switch//'" is not set as "required" but no default value has been passed!'
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  !> @brief Procedure for initializing CLA.
  !> @note For logical type CLA the value is directly read without any robust error trapping.
  subroutine get(cla,pref,val)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Command_Line_Argument), intent(INOUT):: cla      !< CLA data.
  character(*), optional,            intent(IN)::    pref     !< Prefixing string.
  class(*),                          intent(INOUT):: val      !< CLA value.
  integer::                                          val_kind !< Kind value for numeric value.
  character(len=:), allocatable::                    prefd    !< Prefixing string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  select type(val)
  type is(real)
    val_kind = kind(val)
    if (cla%passed) then
      ! cla%val has a value, cla being passed to CLI
      if (cla%act==action_store) then
        if     (val_kind==R16P) then
          val = cton(str=trim(adjustl(cla%val)),knd=1._R16P)
        elseif (val_kind==R8P) then
          val = cton(str=trim(adjustl(cla%val)),knd=1._R8P)
        elseif (val_kind==R4P) then
          val = cton(str=trim(adjustl(cla%val)),knd=1._R4P)
        endif
      endif
    else
      ! cla%val has not a value
      if (cla%required) then
        ! cla%def has not a value, cla being required
        write(stderr,'(A)')prefd//' Error: CLA "'//trim(adjustl(cla%switch))//'" is required by CLI but it has not been passed!'
      else
        ! cla%def has a value
        if (cla%act==action_store) then
          if     (val_kind==R16P) then
            val = cton(str=trim(adjustl(cla%def)),knd=1._R16P)
          elseif (val_kind==R8P) then
            val = cton(str=trim(adjustl(cla%def)),knd=1._R8P)
          elseif (val_kind==R4P) then
            val = cton(str=trim(adjustl(cla%def)),knd=1._R4P)
          endif
        endif
      endif
    endif
  type is(integer)
    val_kind = kind(val)
    if (cla%passed) then
      ! cla%val has a value, cla being passed to CLI
      if (cla%act==action_store) then
        if     (val_kind==I8P) then
          val = cton(str=trim(adjustl(cla%val)),knd=1_I8P)
        elseif (val_kind==I4P) then
          val = cton(str=trim(adjustl(cla%val)),knd=1_I4P)
        elseif (val_kind==I2P) then
          val = cton(str=trim(adjustl(cla%val)),knd=1_I2P)
        elseif (val_kind==I1P) then
          val = cton(str=trim(adjustl(cla%val)),knd=1_I1P)
        endif
      endif
    else
      ! cla%val has not a value
      if (cla%required) then
        ! cla%def has not a value, cla being required
        write(stderr,'(A)')prefd//' Error: CLA "'//trim(adjustl(cla%switch))//'" is required by CLI but it has not been passed!'
      else
        ! cla%def has a value
        if (cla%act==action_store) then
          if     (val_kind==I8P) then
            val = cton(str=trim(adjustl(cla%def)),knd=1_I8P)
          elseif (val_kind==I4P) then
            val = cton(str=trim(adjustl(cla%def)),knd=1_I4P)
          elseif (val_kind==I2P) then
            val = cton(str=trim(adjustl(cla%def)),knd=1_I2P)
          elseif (val_kind==I1P) then
            val = cton(str=trim(adjustl(cla%def)),knd=1_I1P)
          endif
        endif
      endif
    endif
  type is(logical)
    if (cla%passed) then
      ! cla%val has a value, cla being passed to CLI
      if (cla%act==action_store) then
        read(cla%val,*)val
      elseif (cla%act==action_store_true) then
        val = .true.
      elseif (cla%act==action_store_false) then
        val = .false.
      endif
    else
      ! cla%val has not a value
      if (cla%required) then
        ! cla%def has not a value, cla being required
        write(stderr,'(A)')prefd//' Error: CLA "'//trim(adjustl(cla%switch))//'" is required by CLI but it has not been passed!'
      else
        ! cla%def has a value
        read(cla%def,*)val
      endif
    endif
  type is(character(*))
    if (cla%passed) then
      ! cla%val has a value, cla being passed to CLI
      if (cla%act==action_store) then
        val = cla%val
      endif
    else
      ! cla%val has not a value
      if (cla%required) then
        ! cla%def has not a value, cla being required
        write(stderr,'(A)')prefd//' Error: CLA "'//trim(adjustl(cla%switch))//'" is required by CLI but it has not been passed!'
      else
        ! cla%def has a value
        val = cla%def
      endif
    endif
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get
  !> @}

  !> @ingroup Data_Type_Command_Line_ArgumentPublicProcedure
  !> @{
  !> @brief Procedure for parsing Command Line Arguments by means of a previously initialized CLA list.
  !> @note The leading and trailing white spaces are removed from CLA values.
  subroutine cli_parser(pref,help,examples,progname,cla_list,error)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), optional,           intent(IN)::    pref           !< Prefixing string.
  character(*), optional,           intent(IN)::    help           !< Help message describing the Command Line Interface.
  character(*), optional,           intent(IN)::    examples(1:)   !< Examples of correct usage.
  character(*),                     intent(IN)::    progname       !< Program name.
  type(Type_Command_Line_Argument), intent(INOUT):: cla_list(1:)   !< CLA list.
  integer(I4P),                     intent(OUT)::   error          !< Error trapping flag.
  integer(I4P)::                                    Nca            !< Number of command line arguments passed.
  integer(I4P)::                                    Nca_required   !< Number of command line arguments that CLI requires.
  integer(I4P)::                                    Nca_optional   !< Number of command line arguments that are optional for CLI.
  character(len=:), allocatable::                   ca_switch      !< Switch name.
  character(max_val_len)::                          ca_val         !< Switch value.
  integer(I4P)::                                    max_switch_len !< Number of command line arguments passed.
  character(len=:), allocatable::                   cli_help       !< Dummy variable for CLI help.
  logical::                                         found          !< Flag for checking if ca_switch has been found in cla_list.
  character(len=:), allocatable::                   prefd          !< Prefixing string.
  integer(I4P)::                                    c,cc           !< Counter for command line arguments.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  error = 0
  prefd = '' ; if (present(pref)) prefd = pref
  ! setting the general CLI help message
  if (present(help)) then
    cli_help = help
  else
    cli_help = ' The Command Line Interface (CLI) has the following options'
  endif
  ! counting the required and optional CLA and computing the maximum length of switch name
  Nca_required   = 0_I4P
  Nca_optional   = 0_I4P
  max_switch_len = 0_I4P
  do c=1,size(cla_list,dim=1)
    if (cla_list(c)%required) then
      Nca_required = Nca_required + 1_I4P
    else
      Nca_optional = Nca_optional + 1_I4P
    endif
    max_switch_len = max(max_switch_len,len_trim(adjustl(cla_list(c)%switch)))
  enddo
  ca_switch = repeat(' ',max_switch_len)
  ! counting the passed CLA
  Nca = command_argument_count()
  if (Nca<Nca_required) then
    write(stderr,'(A)')prefd//' Error: the Command Line Interface requires at least '//trim(str(.true.,Nca_required))//&
                              ' arguments to be passed whereas only '//trim(str(.true.,Nca))//' have been!'
    call print_usage
    error = 1
    return
  else
    ! parsing switch
    c = 0
    do while (c<Nca)
      c = c + 1
      call get_command_argument(c,ca_switch)
      found = .false.
      do cc=1,size(cla_list,dim=1)
        if (trim(adjustl(cla_list(cc)%switch   ))==trim(adjustl(ca_switch)).or.&
            trim(adjustl(cla_list(cc)%switch_ab))==trim(adjustl(ca_switch))) then
          if (cla_list(cc)%act==action_store) then
            c = c + 1
            call get_command_argument(c,ca_val)
            cla_list(cc)%val = trim(adjustl(ca_val))
          endif
          cla_list(cc)%passed = .true.
          found = .true.
        endif
      enddo
      if (.not.found) then
        write(stderr,'(A)')prefd//' Error: switch "'//trim(adjustl(ca_switch))//'" is unknown!'
        call print_usage
        error = 2
        return
      endif
    enddo
  endif
  ! checking if all required CLAs have been passed
  do c=1,size(cla_list,dim=1)
    if (cla_list(c)%required) then
      if (.not.cla_list(c)%passed) then
        write(stderr,'(A)')prefd//' Error: CLA "'//trim(adjustl(cla_list(c)%switch))//&
                                  '" is required by CLI but it has not been passed!'
        call print_usage
        error = 3
        return
      endif
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Procedure for printing the correct use Command Line Interface accordingly to the cla_list passed.
    subroutine print_usage
    !-------------------------------------------------------------------------------------------------------------------------------
    character(len=:), allocatable:: cla_list_sign !< Complete signature of CLA list.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    cla_list_sign = '   '//progname//' '
    do c=1,size(cla_list,dim=1)
      if (cla_list(c)%act==action_store) then
        if (cla_list(c)%required) then
          cla_list_sign = trim(cla_list_sign)//' '//trim(adjustl(cla_list(c)%switch))//' value'
        else
          cla_list_sign = trim(cla_list_sign)//' ['//trim(adjustl(cla_list(c)%switch))//' value]'
        endif
      else
        if (cla_list(c)%required) then
          cla_list_sign = trim(cla_list_sign)//' '//trim(adjustl(cla_list(c)%switch))
        else
          cla_list_sign = trim(cla_list_sign)//' ['//trim(adjustl(cla_list(c)%switch))//']'
        endif
      endif
    enddo
    write(stdout,'(A)')prefd//cli_help
    write(stdout,'(A)')prefd//cla_list_sign
    write(stdout,'(A)')prefd//' Each Command Line Arguments (CLA) has the following meaning:'
    do c=1,size(cla_list,dim=1)
      if (cla_list(c)%act==action_store) then
        if (trim(adjustl(cla_list(c)%switch))/=trim(adjustl(cla_list(c)%switch_ab))) then
          write(stdout,'(A)')prefd//'   ['//trim(adjustl(cla_list(c)%switch))//' value] or ['//&
                                            trim(adjustl(cla_list(c)%switch_ab))//' value]'
        else
          write(stdout,'(A)')prefd//'   ['//trim(adjustl(cla_list(c)%switch))//' value]'
        endif
      else
        if (trim(adjustl(cla_list(c)%switch))/=trim(adjustl(cla_list(c)%switch_ab))) then
          write(stdout,'(A)')prefd//'   ['//trim(adjustl(cla_list(c)%switch))//'] or ['//trim(adjustl(cla_list(c)%switch_ab))//']'
        else
          write(stdout,'(A)')prefd//'   ['//trim(adjustl(cla_list(c)%switch))//']'
        endif
      endif
      write(stdout,'(A)')prefd//'     '//trim(adjustl(cla_list(c)%help))
      if (cla_list(c)%required) then
        write(stdout,'(A)')prefd//'     It is a non optional CLA thus must be passed to CLI'
      else
        write(stdout,'(A)')prefd//'     It is a optional CLA which default value is "'//trim(adjustl(cla_list(c)%def))//'"'
      endif
    enddo
    if (present(examples)) then
      write(stdout,'(A)')prefd//' Usage examples:'
      do c=1,size(examples,dim=1)
        write(stdout,'(A)')prefd//'   -) '//trim(examples(c))
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine print_usage
  endsubroutine cli_parser

  !> @brief Procedure for printing Command Line Arguments being previously parsed.
  subroutine cli_print(pref,cla_list)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), optional,           intent(IN)::  pref           !< Prefixing string.
  type(Type_Command_Line_Argument), intent(IN)::  cla_list(1:)   !< CLA list.
  character(len=:), allocatable::                 prefd          !< Prefixing string.
  integer(I4P)::                                  c              !< Counter for command line arguments.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(stdout,'(A)')prefd//' CLAs parsed:'
  do c=1,size(cla_list,dim=1)
    if (cla_list(c)%passed) then
      if (cla_list(c)%act==action_store) then
        write(stdout,'(A)')prefd//'   '//cla_list(c)%switch//' '//trim(adjustl(cla_list(c)%val))
      else
        write(stdout,'(A)')prefd//'   '//cla_list(c)%switch//' '//trim(adjustl(cla_list(c)%act))
      endif
    else
      write(stdout,'(A)')prefd//'   '//cla_list(c)%switch//' '//trim(adjustl(cla_list(c)%def))
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine cli_print
  !> @}
endmodule Data_Type_Command_Line_Argument
