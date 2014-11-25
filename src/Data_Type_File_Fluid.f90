!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_FluidDerivedType Data_Type_File_Fluid
!> Module definition of main fluid dynamic file structures, Type_File_Fluid
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_FluidPrivateProcedure Data_Type_File_Fluid
!> Module definition of main fluid dynamic file structures, Type_File_Fluid
!> @}

!> @brief Module Data_Type_File_Fluid contains the definition of Type_File_Fluid, that is the main fluid dynamic file structures
!> adopted by OFF project.
!> @todo \b encodin: Implement ascii and base64 encoding besides raw one
module Data_Type_File_Fluid
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                 ! Integers and reals precision definition.
USE Data_Type_File_Base, only: Type_File_Base,cr ! Definition of Type_File_Base.
USE Data_Type_Global,    only: Type_Global       ! Definition of Type_Global.
USE Data_Type_Parallel,  only: Type_Parallel     ! Definition of Type_Parallel.
USE Data_Type_Species,   only: Type_Species      ! Definition of Type_Species.
USE Data_Type_SBlock,    only: Type_SBlock       ! Definition of Type_SBlock.
USE Data_Type_Time_Step, only: Type_Time_Step    ! Definition of Type_Time_Step.
USE Data_Type_XML_Tag,   only: Type_XML_Tag      ! Definition of Type_XML_Tag.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Fluid.
!> @ingroup Data_Type_File_FluidDerivedType
type, public, extends(Type_File_Base):: Type_File_Fluid
  contains
    procedure:: load_header => load_header_file_fluid      ! Procedure for loading header from fluid file.
    procedure:: save_header => save_header_file_fluid      ! Procedure for saving header into fluid file.
    procedure:: load_block  => load_block_cells_file_fluid ! Procedure for loading blocks cells fluid dynamic data from fluid file.
    procedure:: save_block  => save_block_cells_file_fluid ! Procedure for saving blocks cells fluid dynamic data into fluid file.
    procedure:: load        => load_file_fluid             ! Procedure for loading fluid file.
    procedure:: save        => save_file_fluid             ! Procedure for saving fluid file.
endtype Type_File_Fluid
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_FluidPrivateProcedure
  !> @{
  !> @brief Procedure for loading header of fluid file. This procedure loads the header of fluid file that contains the fluid
  !> dynamic data of all blocks (not only the ones belonging to myrank).
  !> @note This header filecontains the following data:
  !>   - species0;
  !>   - time_step.
  !> The fluid file is an unformatted file with STREAM access. However it is formatted by means of XML syntax, thus
  !> string containing a tag must be terminated by means of a carriage-return (cr=achar(10)) character.
  !> The header of fluid file is composed as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Fluid_Header>
  !>     <Initial_Species>
  !>       <Specie s="1"> <specie .../></Specie>
  !>       <Specie s="2"> <specie .../></Specie>
  !>       ...
  !>       <Specie s="Ns"><specie .../></Specie>
  !>     </Initial_Species>
  !>     <Time>...</Time>
  !>     <Time_Step>...</Time_Step>
  !>   </Fluid_Header>
  !> @endcode
  !> The fluid file header is composed by means of 3 tags:
  !>   * 'Initial_Species' defining the initial species; it contains 'Ns' 'Specie' tags that are parsed accordingly the definition
  !>     of Type_Specie; it is no necessary to explicitely define the number of initial species 'Ns' that is automatically computed
  !>     parsing the header; the tag 'Specie' is an 'overloading' tag of 'specie' one actually defining the specie componenents; the
  !>     overloading tag 'Specie' introduces the attribute 's' defining the numeration of each specie into the global numeration;
  !>     the total number of species defined is automatically computed parsing the whole file, thus it is not necessary to
  !>     explicitely define it; the tag 'specie' is defined into 'Data_Type_Specie' module;
  !>   * 'Time' defining the actual time of the saved fluid dynamic data;
  !>   * 'Time_Step' defining the actual time-step of the saved fluid dynamic data.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_header_file_fluid(file_d,species0,time_step)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Fluid), intent(INOUT):: file_d    !< File data.
  type(Type_Species),     intent(INOUT):: species0  !< Initial species.
  type(Type_Time_Step),   intent(INOUT):: time_step !< Time step data.
  type(Type_XML_Tag)::                    tag       !< Dummy XML tag for parsing file.
  integer(I4P)::                          s         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(t => time_step%t,n => time_step%n)
    allocate(character(len=1000):: tag%string%vs)
    call file_d%open(ascii=.true.,action='READ') ; if (file_d%iostat/=0) return
    species0%Ns = 0
    Species_Count: do
      read(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)tag%string%vs ; if (file_d%iostat/=0) return
      if (index(string=tag%string%vs,substring='<Specie ')>0) then
        species0%Ns = species0%Ns + 1
      endif
    enddo Species_Count
    10 rewind(unit=file_d%unit)
    call species0%alloc
    call species0%compute_Npc
    Read_Loop: do
      read(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)tag%string%vs ; if (file_d%iostat/=0) return
      if     (index(string=tag%string%vs,substring='<Specie ')>0) then
        call tag%set(tag_name='Specie',att_name=['s']) ; call tag%get_attributes ; call tag%get_value
        s = cton(str=trim(adjustl(tag%att_val(1)%vs)),knd=1_I4P)
        call species0%heats(s)%load_str_xml(str_xml=tag%tag_val%vs)
      elseif (index(string=tag%string%vs,substring='<Time>')>0) then
        call tag%set(tag_name='Time') ; call tag%get_value
        t = cton(str=trim(adjustl(tag%tag_val%vs)),knd=1._R8P)
      elseif (index(string=tag%string%vs,substring='<Time_Step>')>0) then
        call tag%set(tag_name='Time_Step') ; call tag%get_value
        n = cton(str=trim(adjustl(tag%tag_val%vs)),knd=1_I8P)
      elseif (index(string=tag%string%vs,substring='</Fluid_Header')>0) then
        exit Read_Loop
      endif
    enddo Read_Loop
    20 call file_d%close
  endassociate
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_header_file_fluid

  !> @brief Procedure for saving header of fluid file.
  !> @note The definition of the header fluid file syntax is described in the 'load_header_file_fluid' documentation.
  subroutine save_header_file_fluid(file_d,n,flip,species0,time_step)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Fluid), intent(INOUT):: file_d    !< File data.
  integer(I8P), optional, intent(IN)::    n         !< Time step number.
  integer(I1P), optional, intent(IN)::    flip      !< Flip-flop number.
  type(Type_Species),     intent(IN)::    species0  !< Initial species.
  type(Type_Time_Step),   intent(IN)::    time_step !< Time step data.
  integer(I4P)::                          s         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if     (present(n)) then
    call file_d%open(replace=.true.,action='WRITE',n=n) ; if (file_d%iostat/=0) return
  elseif (present(flip)) then
    call file_d%open(replace=.true.,action='WRITE',flip=flip) ; if (file_d%iostat/=0) return
  else
    call file_d%open(replace=.true.,action='WRITE') ; if (file_d%iostat/=0) return
  endif
  associate(t => time_step%t,n => time_step%n)
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<?xml version="1.0"?>'//cr
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<Fluid_Header>'//cr
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'  <Initial_Species>'//cr
    do s=1,species0%Ns
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'    <Specie s="'//trim(str(n=s))//'">'//&
        species0%heats(s)%save_str_xml()//'</Specie>'//cr
    enddo
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'  </Initial_Species>'//cr
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'  <Time>'//     trim(str(n=t))//'</Time>'//cr
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'  <Time_Step>'//trim(str(n=n))//'</Time_Step>'//cr
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'</Fluid_Header>'//cr
    call file_d%close
  endassociate
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_header_file_fluid

  !> @brief Procedure for loading blocks cells fluid dyamic data from fluid file.
  !> @note After the header (see 'load_header_file_fluid' documentation) the fluid file contains the cells fluid dynamic data. These
  !> are saved inside the 'Fluid_Cells' tag. The tag syntax is the following:
  !> @code
  !>   <Fluid_Cells>
  !>     <block_fluid_cells ID="1" encoding="..."> fluid_data </block_fluid_cells>
  !>     <block_fluid_cells ID="2" encoding="..."> fluid_data </block_fluid_cells>
  !>     ...
  !>     <block_fluid_cells ID="Nb_tot" encoding="..."> fluid_data </block_fluid_cells>
  !>   </Fluid_Cells>
  !> @endcode
  !> Each 'block_fluid_cells' tag has 2 attributes, 'ID' indicating the block ID-key and 'encoding' that can
  !> be set to 3 different values: 'ascii', 'raw' or 'base64' in the case fluid_data are saved in asciii, raw-binary or Base64-ascii
  !> encoded format, respectively.
  !> Note that it is not necessary to explicitely declare the global blocks number, Nb_tot, and grid levels that are computed
  !> parsing the mesh file that must be loaded before the present file.
  !> The 'block_fluid_cells' tag contains the following data:
  !>   - block%C.
  !> It is worth noting that the 'block_fluid_cells' tags can appear in any order.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_block_cells_file_fluid(file_d,block,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Fluid), intent(INOUT):: file_d !< File data.
  type(Type_SBlock),      intent(INOUT):: block  !< Block data.
  integer(I8P),           intent(IN)::    ID     !< ID-key of block to be loaded.
  type(Type_XML_Tag)::                    tag    !< Dummy XML tag for parsing file.
  character(len=1)::                      c1     !< Dummy string for parsing file.
  character(len=:), allocatable::         c2     !< Dummy string for parsing file.
  integer(I8P)::                          ib     !< Counter.
  integer(I4P)::                          i,j,k  !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%open(action='READ') ; if (file_d%iostat/=0) return
  ! reading file until 'Fluid_Cells' is reached
  Fluid_Cells_Tag_Search: do
    read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return
    if (index(string=c1,substring='<')>0) then
      c2 = c1
      Tag_End_Serch: do
        read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return
        c2 = c2//c1
        if (index(string=c1,substring='>')>0) exit Tag_End_Serch
      enddo Tag_End_Serch
      if (index(string=c2,substring='<Fluid_Cells>')>0) then
        read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return ! reading cr char
        exit Fluid_Cells_Tag_Search
      endif
    endif
  enddo Fluid_Cells_Tag_Search
  10 continue
  call tag%set(tag_name='block_fluid_cells')
  call tag%alloc(att_name=.true.,Na=2)
  tag%att_name(1)%vs = 'ID'
  tag%att_name(2)%vs = 'encoding'
  Tag_Block_Fluid_Cells_Parsing: do
    read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)c1 ; if (file_d%iostat/=0) return
    if (c1=='<') then
      c2 = c1
      Tag_Block_Fluid_Cells_End_Serch: do
        read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)c1 ; if (file_d%iostat/=0) return
        c2 = c2//c1
        if (c1=='>') exit Tag_Block_Fluid_Cells_End_Serch
      enddo Tag_Block_Fluid_Cells_End_Serch
      if     (index(string=c2,substring='<block_fluid_cells')>0) then
        call tag%set(string=c2) ; call tag%get_attributes
        ib = cton(str=tag%att_val(1)%vs,knd=1_I8P)
        if (ib==ID) then
          associate(gc => block%dims%gc,Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
            do k=1-gc(5),Nk+gc(6)
              do j=1-gc(3),Nj+gc(4)
                do i=1-gc(1),Ni+gc(2)
                  call block%C(i,j,k)%load(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
                enddo
              enddo
            enddo
          endassociate
          exit Tag_Block_Fluid_Cells_Parsing
        else
          cycle Tag_Block_Fluid_Cells_Parsing
        endif
      elseif (index(string=c2,substring='</Fluid_Cells')>0) then
        exit Tag_Block_Fluid_Cells_Parsing
      endif
    endif
  enddo Tag_Block_Fluid_Cells_Parsing
  20 call file_d%close
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_block_cells_file_fluid

  !> @brief Procedure for saving blocks cells fluid dynamic data from mesh file.
  !> @note The definition of the cells fluid dynamic data syntax is described in the 'load_cells_file_mesh' documentation.
  subroutine save_block_cells_file_fluid(file_d,n,flip,fluid_cells_tag_open,fluid_cells_tag_close,block,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Fluid),        intent(INOUT):: file_d                !< File data.
  integer(I8P),        optional, intent(IN)::    n                     !< Time step number.
  integer(I1P),        optional, intent(IN)::    flip                  !< Flip-flop number.
  logical,             optional, intent(IN)::    fluid_cells_tag_open  !< Switch for opening the 'Fluid_Cells' tag.
  logical,             optional, intent(IN)::    fluid_cells_tag_close !< Switch for closing the 'Fluid_Cells' tag.
  type(Type_SBlock),   optional, intent(IN)::    block                 !< Block data.
  integer(I8P),        optional, intent(IN)::    ID                    !< ID-key of block to be saved.
  integer(I4P)::                                 i,j,k                 !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if     (present(n)) then
    call file_d%open(append=.true.,action='WRITE',n=n) ; if (file_d%iostat/=0) return
  elseif (present(flip)) then
    call file_d%open(append=.true.,action='WRITE',flip=flip) ; if (file_d%iostat/=0) return
  else
    call file_d%open(append=.true.,action='WRITE') ; if (file_d%iostat/=0) return
  endif
  if     (present(fluid_cells_tag_open )) then
    if (fluid_cells_tag_open) then
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<Fluid_Cells>'//cr
    endif
  elseif (present(fluid_cells_tag_close)) then
    if (fluid_cells_tag_close) then
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'</Fluid_Cells>'//cr
    endif
  elseif (present(block).and.present(ID)) then
    associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
      write(unit=unit,iostat=iostat,iomsg=iomsg)'<block_fluid_cells ID="'//trim(str(n=ID))//'" encoding="raw">'
      associate(gc => block%dims%gc,Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
        do k=1-gc(5),Nk+gc(6)
          do j=1-gc(3),Nj+gc(4)
            do i=1-gc(1),Ni+gc(2)
              call block%C(i,j,k)%save(unit=unit,iostat=iostat,iomsg=iomsg)
            enddo
          enddo
        enddo
      endassociate
      write(unit=unit,iostat=iostat,iomsg=iomsg)'</block_fluid_cells>'//cr
    endassociate
  endif
  call file_d%close
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block_cells_file_fluid

  !> @brief Procedure for loading fluid file.
  !> @note The fluid file is an unformatted file with STREAM access. However it is formatted by means of XML syntax, thus
  !> string containing a tag must be terminated by means of a carriage-return (cr=char(10)) character.
  !> The fluid file contains the following data:
  !>   - global%species0;
  !>   - global%time_step.
  !>   - global%block%C for all blocks.
  !> The fluid file is composed as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Fluid_Header>
  !>   ...
  !>   </Fluid_Header>
  !>   <Fluid_Cells>
  !>   ...
  !>   </Fluid_Cells>
  !> @endcode
  !> For 'Fluid_Header' and 'Fluid_Cells' tags see the documentation of 'load_header_file_fluid' and 'load_block_cells_file_fluid',
  !> respectively.
  !> The blocks loaded are only the ones belonging to myrank. This is accomplished in a transparent way. It is worth nothing that
  !> the syntax of tag attributes is case sensitive, whereas the syntax of their values is case insensitive.
  subroutine load_file_fluid(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Fluid), intent(INOUT):: file_d !< File data.
  type(Type_Global),      intent(INOUT):: global !< Global data.
  type(Type_SBlock), pointer::            block  !< Block pointer for scanning global%block tree.
  integer(I8P)::                          ID     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%load_header(species0=global%species0,time_step=global%time_step)
  call global%set_Ns_from_species0
  call global%set_Nrk_from_rk_ord
  call global%blocks_dims_update
  do while(global%block%loopID(ID=ID))
    block => global%block%dat(ID=ID)
    call block%alloc(members=.true.)
    call file_d%load_block(block=block,ID=ID)
    call block%update_primitive(species0=global%species0)
    call block%primitive2conservative
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_file_fluid

  !> @brief Procedure for saving mesh file.
  !> @note The definition of fluid file syntax is described in the 'load_file_fluid' documentation.
  subroutine save_file_fluid(file_d,n,flip,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Fluid), intent(INOUT):: file_d !< File data.
  integer(I8P), optional, intent(IN)::    n      !< Time step number.
  integer(I1P), optional, intent(IN)::    flip   !< Flip-flop number.
  type(Type_Global),      intent(IN)::    global !< Global data.
  type(Type_SBlock), pointer::            block  !< Block pointer for scanning global%block tree.
  integer(I8P)::                          ID     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (global%parallel%is_master()) then
    if     (present(n)) then
      call file_d%save_header(species0=global%species0,time_step=global%time_step,n=n)
      call file_d%save_block(fluid_cells_tag_open=.true.,n=n)
    elseif (present(flip)) then
      call file_d%save_header(species0=global%species0,time_step=global%time_step,flip=flip)
      call file_d%save_block(fluid_cells_tag_open=.true.,flip=flip)
    else
      call file_d%save_header(species0=global%species0,time_step=global%time_step)
      call file_d%save_block(fluid_cells_tag_open=.true.)
    endif
  endif
#ifdef _MPI
#else
  do while(global%block%loopID(ID=ID))
    block => global%block%dat(ID=ID)
    if     (present(n)) then
      call file_d%save_block(block=block,ID=ID,n=n)
    elseif (present(flip)) then
      call file_d%save_block(block=block,ID=ID,flip=flip)
    else
      call file_d%save_block(block=block,ID=ID)
    endif
  enddo
#endif
  if (global%parallel%is_master()) then
    if     (present(n)) then
      call file_d%save_block(fluid_cells_tag_close=.true.,n=n)
    elseif (present(flip)) then
      call file_d%save_block(fluid_cells_tag_close=.true.,flip=flip)
    else
      call file_d%save_block(fluid_cells_tag_close=.true.)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_file_fluid
  !> @}
endmodule Data_Type_File_Fluid
