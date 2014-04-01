!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_IBM_OptionsDerivedType Data_Type_File_IBM_Options
!> Module definition of main fluid dynamic file structures, Type_File_IBM_Options
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_IBM_OptionsPrivateProcedure Data_Type_File_IBM_Options
!> Module definition of main fluid dynamic file structures, Type_File_IBM_Options
!> @}

!> @brief Module Data_Type_File_IBM_Options contains the definition of Type_File_IBM_Options, that is the main IBM options file.
module Data_Type_File_IBM_Options
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                        ! Integers and reals precision definition.
USE Data_Type_File_Base,      only: Type_File_Base      ! Definition of Type_File_Base.
USE Data_Type_OS,             only: Type_OS             ! Definition of Type_OS.
USE Data_Type_Varying_String, only: Type_Varying_String ! Definition of Type_XML_Tag.
USE Data_Type_XML_Tag,        only: Type_XML_Tag        ! Definition of Type_XML_Tag.
USE Lib_IO_Misc                                         ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_IBM_Options.
!> @ingroup Data_Type_File_IBM_OptionsDerivedType
type, public, extends(Type_File_Base):: Type_File_IBM_Options
  character(len=:), allocatable::          fn_spec        !< File name of species file.
  character(len=:), allocatable::          fn_mesh        !< File name of mesh file.
  character(len=:), allocatable::          fn_bc          !< File name of boundary conditions file.
  character(len=:), allocatable::          fn_init        !< File name of initial conditions file.
  character(len=:), allocatable::          desc_type      !< Blocks description type: 'BLOCKS' for simple (Cartesian descriptions),
                                                          !< 'ICEMCFD' for complex grids.
  type(Type_Varying_String), allocatable:: fn_blk_desc(:) !< File names of blocks description.
  integer(I4P)::                           Nl             !< Number of grid levels.
  contains
    procedure:: free  => free_ibm_options  ! Procedure for freeing dynamic memory.
    procedure:: alloc => alloc_ibm_options ! Procedure for allocating dynamic memory.
    procedure:: load  => load_ibm_options  ! Procedure for loading IBM options file.
    procedure:: save  => save_ibm_options  ! Procedure for saving IBM options file.
    procedure:: print => print_ibm_options ! Procedure for printing IBM options with a pretty format.
    final::     finalize_ibm_options       ! Procedure for freeing dynamic memory when finalizing.
endtype Type_File_IBM_Options
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_IBM_OptionsPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_ibm_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_IBM_Options), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free_base
  if (allocated(file_d%fn_spec    )) deallocate(file_d%fn_spec  )
  if (allocated(file_d%fn_mesh    )) deallocate(file_d%fn_mesh  )
  if (allocated(file_d%fn_bc      )) deallocate(file_d%fn_bc    )
  if (allocated(file_d%fn_init    )) deallocate(file_d%fn_init  )
  if (allocated(file_d%desc_type  )) deallocate(file_d%desc_type)
  if (allocated(file_d%fn_blk_desc)) then
    call file_d%fn_blk_desc%free ; deallocate(file_d%fn_blk_desc)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_ibm_options

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize_ibm_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_File_IBM_Options), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize_ibm_options

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc_ibm_options(file_d,Nf)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_IBM_Options), intent(INOUT):: file_d !< File data.
  integer(I4P),                 intent(IN)::    Nf     !< Number of blocks description files.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(file_d%fn_blk_desc)) deallocate(file_d%fn_blk_desc) ; allocate(file_d%fn_blk_desc(1:Nf))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_ibm_options

  !> @brief Procedure for loading IBM options.
  !> @note This options file is a very simple XML file that drives the running of IBM program.
  !> @note This options file is formatted by means of XML syntax. A validated file should be formatted as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <IBM_Options>
  !>     <inputs>
  !>       <path> ... </path>
  !>       <species_file_name> ... </species_file_name>
  !>       <blocks_description type="blocks(or icemcfd)"> ... </blocks_description>
  !>     </inputs>
  !>     <outputs>
  !>       <path> ... </path>
  !>       <grid_levels> ... </grid_levels>
  !>       <mesh_file_name> ... </mesh_file_name>
  !>       <bc_file_name> ... </bc_file_name>
  !>       <init_file_name> ... </init_file_name>
  !>     </outputs>
  !>   </IBM_Options>
  !> @endcode
  !> The main tag is 'IBM_Options' that contains the the 'inputs' and 'outputs' tags. The tags can appear in any order. The 'inputs'
  !> tag is composed by 3 nested tags, namely 'path', 'species_file_name' and 'blocks_description', whereas 'outputs' tag has 5
  !> nested tags, namely 'path', 'grid_levels', 'mesh_file_name', 'bc_file_name' and 'init_file_name'.
  !> The IBM options file contains the following data:
  !>   - file_d%path_in;
  !>   - file_d%path_out;
  !>   - file_d%fn_spec;
  !>   - file_d%fn_mesh;
  !>   - file_d%fn_bc;
  !>   - file_d%fn_init;
  !>   - file_d%desc_type;
  !>   - file_d%fn_blk_desc;
  !>   - file_d%Nl.
  !> It is worth nothing that the syntax of tag names and attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_ibm_options(file_d,OS)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_IBM_Options), intent(INOUT):: file_d   !< File data.
  type(Type_OS),                intent(IN)::    OS       !< Running architecture.
  character(len=:), allocatable::               line     !< Dummy string for parsing option file.
  type(Type_XML_Tag)::                          tag      !< Tag parsing data.
  integer(I4P)::                                t,f,c,Nf !< Counters.
  !--------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg,Nl=>file_d%Nl)
    call file_d%open(ascii=.true.,action='READ') ; if (iostat/=0) return
    allocate(character(len=1000):: tag%string%vs,line)
    Read_Loop: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)line ; if (iostat/=0.and.iostat/=iostat_eor) return
      if     (index(string=line,substring='<inputs')>0) then
        do t=1,3
          read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)tag%string%vs ; if (iostat/=0.and.iostat/=iostat_eor) return
          if     (index(string=tag%string%vs,substring='<path')>0) then
            call tag%set(tag_name='path') ; call tag%get_value
            tag%tag_val%vs=OS%string_separator_fix(string=trim(tag%tag_val%vs))
            c = len_trim(tag%tag_val%vs) ; if (tag%tag_val%vs(c:c)/=OS%sep) tag%tag_val%vs = trim(tag%tag_val%vs)//OS%sep
            call file_d%set(path_in=tag%tag_val%vs)
          elseif (index(string=tag%string%vs,substring='<species_file_name')>0) then
            call tag%set(tag_name='species_file_name') ; call tag%get_value
            file_d%fn_spec = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<blocks_description')>0) then
            call tag%set(tag_name='blocks_description',att_name=['type']) ; call tag%get_value ; call tag%get_attributes
            file_d%desc_type = Upper_Case(trim(tag%att_val(1)%vs))
            tag%tag_val%vs = trim(unique(string=tag%tag_val%vs,substring=' '))
            Nf = count(string=trim(tag%tag_val%vs),substring=' ')+1
            call file_d%alloc(Nf=Nf)
            do f=1,Nf
              c = index(string=tag%tag_val%vs,substring=' ')
              if (c>0) then
                file_d%fn_blk_desc(f)%vs = OS%string_separator_fix(string=trim(tag%tag_val%vs(1:c)))
              else
                file_d%fn_blk_desc(f)%vs = OS%string_separator_fix(string=trim(tag%tag_val%vs))
              endif
              tag%tag_val%vs = trim(tag%tag_val%vs(c+1:))
            enddo
          endif
        enddo
      elseif (index(string=line,substring='<outputs')>0) then
        do t=1,5
          read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)tag%string%vs ; if (iostat/=0.and.iostat/=iostat_eor) return
          if     (index(string=tag%string%vs,substring='<path')>0) then
            call tag%set(tag_name='path') ; call tag%get_value
            tag%tag_val%vs=OS%string_separator_fix(string=trim(tag%tag_val%vs))
            c = len_trim(tag%tag_val%vs) ; if (tag%tag_val%vs(c:c)/=OS%sep) tag%tag_val%vs = trim(tag%tag_val%vs)//OS%sep
            call file_d%set(path_out=tag%tag_val%vs)
          elseif (index(string=tag%string%vs,substring='<grid_levels')>0) then
            call tag%set(tag_name='grid_levels') ; call tag%get_value
            Nl = cton(str=tag%tag_val%vs,knd=1_I4P)
          elseif (index(string=tag%string%vs,substring='<mesh_file_name')>0) then
            call tag%set(tag_name='mesh_file_name') ; call tag%get_value
            file_d%fn_mesh = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<bc_file_name')>0) then
            call tag%set(tag_name='bc_file_name') ; call tag%get_value
            file_d%fn_bc = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          elseif (index(string=tag%string%vs,substring='<init_file_name')>0) then
            call tag%set(tag_name='init_file_name') ; call tag%get_value
            file_d%fn_init = OS%string_separator_fix(string=trim(tag%tag_val%vs))
          endif
        enddo
      endif
    enddo Read_Loop
    10 call file_d%close
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_ibm_options

  !> @brief Procedure for saving IBM options file.
  !> @note The definition of the IBM options file syntax is described in the 'load_ibm_options' documentation.
  subroutine save_ibm_options(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_IBM_Options), intent(INOUT):: file_d !< File data.
  character(len=:), allocatable::               string !< Dummy string.
  integer(I4P)::                                f      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg,Nl=>file_d%Nl)
    call file_d%open(ascii=.true.,action='WRITE') ; if (iostat/=0) return
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<?xml version="1.0"?>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'<IBM_Options>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <inputs>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <path>'//trim(file_d%path_in)//'</path>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <species_file_name>'//trim(file_d%fn_spec)//'</species_file_name>'
    string = '    <blocks_description type="'//trim(Upper_Case(file_d%desc_type))//'">'
    do f=1,size(file_d%fn_blk_desc,dim=1)
      string = trim(string)//' '//trim(file_d%fn_blk_desc(f)%vs)
    enddo
    string = trim(string)//'</blocks_description>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)trim(string)
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  </inputs>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  <outputs>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <path>'//trim(file_d%path_out)//'</path>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <grid_levels>'//trim(str(n=Nl))//'</grid_levels>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <mesh_file_name>'//trim(file_d%fn_mesh)//'</mesh_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <bc_file_name>'//trim(file_d%fn_bc)//'</bc_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'    <init_file_name>'//trim(file_d%fn_init)//'</init_file_name>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'  </outputs>'
    write(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg)'</IBM_Options>'
    call file_d%close
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_ibm_options

  !> @brief Procedure for printing IBM options with a pretty format.
  subroutine print_ibm_options(file_d,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_IBM_Options), intent(INOUT):: file_d  !< File data.
  character(*), optional,       intent(IN)::    pref    !< Prefixing string.
  integer(I4P), optional,       intent(OUT)::   iostat  !< IO error.
  character(*), optional,       intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),                 intent(IN)::    unit    !< Logic unit.
  character(len=:), allocatable::               prefd   !< Prefixing string.
  integer(I4P)::                                iostatd !< IO error.
  character(500)::                              iomsgd  !< Temporary variable for IO error message.
  character(len=:), allocatable::               string  !< Dummy string.
  integer(I4P)::                                f       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Inputs'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   path: '//file_d%path_in
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   species file name: '//file_d%fn_spec
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   blocks description type: '//Upper_Case(file_d%desc_type)
  string = ''
  do f=1,size(file_d%fn_blk_desc,dim=1)
    string = trim(string)//' '//trim(file_d%fn_blk_desc(f)%vs)
  enddo
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   blocks description files: '//string
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Outputs'
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   path: '//file_d%path_out
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   mesh file name: '//file_d%fn_mesh
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   bc file name: '//file_d%fn_bc
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   init file name: '//file_d%fn_init
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_ibm_options
  !> @}
endmodule Data_Type_File_IBM_Options
