!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_OFF_OptionsDerivedType Data_Type_File_OFF_Options
!> Module definition of Type_File_OFF_Options
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_File_OFF_OptionsInterface Data_Type_File_OFF_Options
!> Module definition of Type_File_OFF_Options
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_File_OFF_OptionsPublicProcedure Data_Type_File_OFF_Options
!> Module definition of Type_File_OFF_Options
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_OFF_OptionsPrivateProcedure Data_Type_File_OFF_Options
!> Module definition of Type_File_OFF_Options
!> @}

!> @brief Module Data_Type_File_OFF_Options contains the definition of Type_File_OFF_Options, that is the main OFF options file.
module Data_Type_File_OFF_Options
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                              ! Integers and reals precision definition.
USE Data_Type_File_Base, only: Type_File_Base ! Definition of Type_File_Base.
USE Data_Type_OS,        only: Type_OS        ! Definition of Type_OS.
USE Data_Type_XML_Tag,   only: Type_XML_Tag   ! Definition of Type_XML_Tag.
USE Lib_IO_Misc,         only: iostat_eor     ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_OFF_Options.
!> @ingroup Data_Type_File_OFF_OptionsDerivedType
type, public, extends(Type_File_Base):: Type_File_OFF_Options
  character(len=:), allocatable:: fn_solv         !< File name of solver options file.
  character(len=:), allocatable:: fn_mesh         !< File name of mesh file.
  character(len=:), allocatable:: fn_bc           !< File name of boundary conditions file.
  character(len=:), allocatable:: fn_init         !< File name of initial conditions file.
  character(len=:), allocatable:: fn_proc         !< File name of processes/blocks map file.
  character(len=:), allocatable:: fn_sol          !< File name of solution file.
  character(len=:), allocatable:: fn_prof         !< File name of profiling file.
  integer(I8P)::                  shl_out = 1_I8P !< Shell output refresh frequency.
  integer(I8P)::                  sol_out = 1_I8P !< Solution output saving frequency.
  integer(I8P)::                  rst_out = 1_I8P !< Solution restart saving frequency.
  integer(I8P)::                  prb_out = 1_I8P !< Probes output saving frequency.
  contains
    procedure:: free  => free_off_options  ! Procedure for freeing dynamic memory.
    procedure:: load  => load_off_options  ! Procedure for loading OFF options file.
    procedure:: save  => save_off_options  ! Procedure for saving OFF options file.
    procedure:: print => print_off_options ! Procedure for printing OFF options with a pretty format.
    final::     finalize_off_options       ! Procedure for freeing dynamic memory when finalizing.
endtype Type_File_OFF_Options
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_OFF_OptionsPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_off_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_OFF_Options), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free_base
  if (allocated(file_d%fn_solv)) deallocate(file_d%fn_solv)
  if (allocated(file_d%fn_mesh)) deallocate(file_d%fn_mesh)
  if (allocated(file_d%fn_bc  )) deallocate(file_d%fn_bc  )
  if (allocated(file_d%fn_init)) deallocate(file_d%fn_init)
  if (allocated(file_d%fn_proc)) deallocate(file_d%fn_proc)
  if (allocated(file_d%fn_sol )) deallocate(file_d%fn_sol )
  if (allocated(file_d%fn_prof)) deallocate(file_d%fn_prof)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_off_options

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize_off_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_File_OFF_Options), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize_off_options

  !> @brief Procedure for loading OFF options file.
  !> @note This options file is a very simple XML file that drives the running of OFF program.
  !> @note This options file is formatted by means of XML syntax. A validated file should be formatted as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <OFF_Options>
  !>     <inputs>
  !>       <path> ... </path>
  !>       <solver_file_name> ... </solver_file_name>
  !>       <mesh_file_name>  ... </mesh_file_name>
  !>       <bc_file_name>    ... </bc_file_name>
  !>       <init_file_name>  ... </init_file_name>
  !>       <procmap_file_name>  ... </procmap_file_name>
  !>     </inputs>
  !>     <outputs>
  !>       <path> ... </path>
  !>       <solution_file_name>  ... </solution_file_name>
  !>       <profiling_file_name>  ... </profiling_file_name>
  !>       <shell_out_frequency> ... </shell_out_frequency>
  !>       <solution_out_frequency> ... </solution_out_frequency>
  !>       <solution_restart_frequency> ... </solution_restart_frequency>
  !>       <probe_out_frequency> ... </probe_out_frequency>
  !>     </outputs>
  !>   </OFF_Options>
  !> @endcode
  !> The main tag is 'OFF_Options' that contains the 'inputs' and 'outputs' tags. The tags can appear in any order. The 'inputs'
  !> tag is composed by 6 nested tags, namely 'path', 'solver_file_name', 'mesh_file_name', 'bc_file_name', 'init_file_name', and
  !> 'procmap_file_name', whereas 'outputs' tag has 7 nested tags, namely 'path', 'solution_file_name', 'profiling_file_name',
  !> 'shell_out_frequency', 'solution_out_frequency', 'solution_restart_frequency'. and 'probe_out_frequency'.
  !> The OFF options file contains the following data:
  !>   - file_d%path_in;
  !>   - file_d%path_out;
  !>   - file_d%fn_solv;
  !>   - file_d%fn_mesh;
  !>   - file_d%fn_bc;
  !>   - file_d%fn_init;
  !>   - file_d%fn_proc;
  !>   - file_d%fn_sol;
  !>   - file_d%fn_prof;
  !>   - file_d%shl_out;
  !>   - file_d%sol_out;
  !>   - file_d%rst_out;
  !>   - file_d%prb_out.
  !> It is worth nothing that the syntax of tag names and attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_off_options(file_d,OS)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_OFF_Options), intent(INOUT):: file_d !< File data.
  type(Type_OS),                intent(IN)::    OS     !< Running architecture.
  character(len=:), allocatable::               line   !< Dummy string for parsing option file.
  type(Type_XML_Tag)::                          tag    !< Tag parsing data.
  integer(I4P)::                                t,c    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg,&
            shl_out=>file_d%shl_out,sol_out=>file_d%sol_out,rst_out=>file_d%rst_out,prb_out=>file_d%prb_out)
    call file_d%open(ascii=.true.,action='READ') ; if (iostat/=0) return
    allocate(character(len=1000):: tag%string%vs,line)
    Read_Loop: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)line ; if (iostat/=0.and.iostat/=iostat_eor) return
      if     (index(string=line,substring='<inputs')>0) then
        do t=1,6
          read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)tag%string%vs ; if (iostat/=0.and.iostat/=iostat_eor) return
          if     (index(string=tag%string%vs,substring='<path')>0) then
            call tag%set(tag_name='path') ; call tag%get_value
            tag%tag_val%vs=OS%string_separator_fix(string=trim(tag%tag_val%vs))
            c = len_trim(tag%tag_val%vs) ; if (tag%tag_val%vs(c:c)/=OS%sep) tag%tag_val%vs = trim(tag%tag_val%vs)//OS%sep
            call file_d%set(path_in=tag%tag_val%vs)
          elseif (index(string=tag%string%vs,substring='<solver_file_name')>0) then
            call tag%set(tag_name='solver_file_name') ; call tag%get_value
            file_d%fn_solv = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<mesh_file_name')>0) then
            call tag%set(tag_name='mesh_file_name') ; call tag%get_value
            file_d%fn_mesh = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<bc_file_name')>0) then
            call tag%set(tag_name='bc_file_name') ; call tag%get_value
            file_d%fn_bc = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<init_file_name')>0) then
            call tag%set(tag_name='init_file_name') ; call tag%get_value
            file_d%fn_init = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<procmap_file_name')>0) then
            call tag%set(tag_name='procmap_file_name') ; call tag%get_value
            file_d%fn_proc = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          endif
        enddo
      elseif (index(string=line,substring='<outputs')>0) then
        do t=1,7
          read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)tag%string%vs ; if (iostat/=0.and.iostat/=iostat_eor) return
          if     (index(string=tag%string%vs,substring='<path')>0) then
            call tag%set(tag_name='path') ; call tag%get_value
            tag%tag_val%vs=OS%string_separator_fix(string=trim(tag%tag_val%vs))
            c = len_trim(tag%tag_val%vs) ; if (tag%tag_val%vs(c:c)/=OS%sep) tag%tag_val%vs = trim(tag%tag_val%vs)//OS%sep
            call file_d%set(path_out=tag%tag_val%vs)
          elseif (index(string=tag%string%vs,substring='<solution_file_name')>0) then
            call tag%set(tag_name='solution_file_name') ; call tag%get_value
            file_d%fn_sol = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<profiling_file_name')>0) then
            call tag%set(tag_name='profiling_file_name') ; call tag%get_value
            file_d%fn_prof = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<shell_out_frequency')>0) then
            call tag%set(tag_name='shell_out_frequency') ; call tag%get_value
            shl_out = cton(str=tag%tag_val%vs,knd=1_I8P)
          elseif (index(string=tag%string%vs,substring='<solution_out_frequency')>0) then
            call tag%set(tag_name='solution_out_frequency') ; call tag%get_value
            sol_out = cton(str=tag%tag_val%vs,knd=1_I8P)
          elseif (index(string=tag%string%vs,substring='<solution_restart_frequency')>0) then
            call tag%set(tag_name='solution_restart_frequency') ; call tag%get_value
            rst_out = cton(str=tag%tag_val%vs,knd=1_I8P)
          elseif (index(string=tag%string%vs,substring='<probe_out_frequency')>0) then
            call tag%set(tag_name='probe_out_frequency') ; call tag%get_value
            prb_out = cton(str=tag%tag_val%vs,knd=1_I8P)
          endif
        enddo
      endif
    enddo Read_Loop
    10 call file_d%close
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_off_options

  !> @brief Procedure for saving OFF options file.
  !> @note The definition of the OFF options file syntax is described in the 'load_off_options' documentation.
  subroutine save_off_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_OFF_Options), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    call file_d%open(ascii=.true.,action='WRITE') ; if (iostat/=0) return
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<?xml version="1.0"?>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<OFF_Options>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <inputs>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <path>'//             trim(file_d%path_in)//'</path>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <solver_file_name>'// trim(file_d%fn_solv)//'</solver_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <mesh_file_name>'//   trim(file_d%fn_mesh)//'</mesh_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <bc_file_name>'//     trim(file_d%fn_bc  )//'</bc_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <init_file_name>'//   trim(file_d%fn_init)//'</init_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <procmap_file_name>'//trim(file_d%fn_proc)//'</procmap_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  </inputs>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <outputs>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <path>'//                      trim(file_d%path_out)//'</path>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <solution_file_name>'//        trim(file_d%fn_sol  )//&
                                                           '</solution_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <profiling_file_name>'//       trim(file_d%fn_prof )//&
                                                           '</profiling_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <shell_out_frequency>'//       trim(str(n=file_d%shl_out ))//&
                                                           '</shell_out_frequency>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <solution_out_frequency>'//    trim(str(n=file_d%sol_out ))//&
                                                           '</solution_out_frequency>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <solution_restart_frequency>'//trim(str(n=file_d%rst_out ))//&
                                                           '</solution_restart_frequency>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <probe_out_frequency>'//       trim(str(n=file_d%prb_out ))//&
                                                           '</probe_out_frequency>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  </outputs>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'</OFF_Options>'
    call file_d%close
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_off_options

  !> @brief Procedure for printing OFF options with a pretty format.
  subroutine print_off_options(file_d,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_OFF_Options), intent(INOUT):: file_d  !< File data.
  character(*), optional,       intent(IN)::    pref    !< Prefixing string.
  integer(I4P), optional,       intent(OUT)::   iostat  !< IO error.
  character(*), optional,       intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),                 intent(IN)::    unit    !< Logic unit.
  character(len=:), allocatable::               prefd   !< Prefixing string.
  integer(I4P)::                                iostatd !< IO error.
  character(500)::                              iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Inputs'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   path: '//file_d%path_in
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   solver file name: '//file_d%fn_solv
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   mesh file name: '//file_d%fn_mesh
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   bc file name: '//file_d%fn_bc
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   init file name: '//file_d%fn_init
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   procmap file name: '//file_d%fn_proc
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Outputs'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   path: '//file_d%path_out
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   solution file name: '//file_d%fn_mesh
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   profiling file name: '//file_d%fn_prof
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   shell output frequency: '//    trim(str(n=file_d%shl_out))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   solution output frequency: '// trim(str(n=file_d%sol_out))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   solution restart frequency: '//trim(str(n=file_d%rst_out))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   probe output frequency: '//    trim(str(n=file_d%prb_out))
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_off_options
  !> @}
endmodule Data_Type_File_OFF_Options
