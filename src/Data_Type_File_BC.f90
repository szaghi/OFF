!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_BCDerivedType Data_Type_File_BC
!> Module definition of main mesh file structures, Type_File_BC
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_BCPrivateProcedure Data_Type_File_BC
!> Module definition of main mesh file structures, Type_File_BC
!> @}

!> @brief Module Data_Type_File_BC contains the definition of Type_File_BC, that is the main mesh file structures adopted by
!> OFF project.
!> @todo \b encodin: Implement ascii and base64 encoding besides raw one
module Data_Type_File_BC
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                 ! Integers and reals precision definition.
USE Data_Type_File_Base, only: Type_File_Base,cr ! Definition of Type_File_Base.
USE Data_Type_Global,    only: Type_Global       ! Definition of Type_Global.
USE Data_Type_Parallel,  only: Type_Parallel     ! Definition of Type_Parallel.
USE Data_Type_SBlock,    only: Type_SBlock       ! Definition of Type_SBlock.
USE Data_Type_XML_Tag,   only: Type_XML_Tag      ! Definition of Type_XML_Tag.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_BC.
!> @ingroup Data_Type_File_BCDerivedType
type, public, extends(Type_File_Base):: Type_File_BC
  contains
    procedure:: load_block => load_block_bc_file_bc ! Procedure for loading block boundary conditions from BC file.
    procedure:: save_block => save_block_bc_file_bc ! Procedure for saving block boundary conditions into BC file.
    procedure:: load       => load_file_bc          ! Procedure for loading boundary conditions file.
    procedure:: save       => save_file_bc          ! Procedure for saving boundary conditions file.
endtype Type_File_BC
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_BCPrivateProcedure
  !> @{
  !> @brief Procedure for loading boundary conditions file.
  !> @note The boundary conditions file is an unformatted file with STREAM access. However it is formatted by means of XML syntax,
  !> thus string containing a tag must be terminated by means of a carriage-return (cr=char(10)) character.
  !> @note The boundary conditions file is composed as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Boundary_Conditions>
  !>     <block_bc b="1"      l="1"  encoding="..."> bc_data </block_bc>
  !>     <block_bc b="2"      l="1"  encoding="..."> bc_data </block_bc>
  !>     ...
  !>     <block_bc b="Nb_tot" l="1"  encoding="..."> bc_data </block_bc>
  !>     <block_bc b="1"      l="2"  encoding="..."> bc_data </block_bc>
  !>     <block_bc b="2"      l="2"  encoding="..."> bc_data </block_bc>
  !>     ...
  !>     <block_bc b="Nb_tot" l="2"  encoding="..."> bc_data </block_bc>
  !>     ...
  !>     <block_bc b="1"      l="Nl" encoding="..."> bc_data </block_bc>
  !>     <block_bc b="2"      l="Nl" encoding="..."> bc_data </block_bc>
  !>     ...
  !>     <block_bc b="Nb_tot" l="Nl" encoding="..."> bc_data </block_bc>
  !>   </Boundary_Conditions>
  !> @endcode
  !> @note Note that it is not necessary to explicitely declare the global blocks number, Nb_tot, and grid levels that are computed
  !> parsing the mesh file that must be loaded before the present file.
  !> Each 'block_bc' tag has 3 attributes, 'b' and 'l' indicating the block number and grid level, and 'encoding' that can
  !> be set to 3 different values: 'ascii', 'raw' or 'base64' in the case fluid_data are saved in asciii, raw-binary or Base64-ascii
  !> encoded format, respectively.
  !> @note The 'block_bc' tag contains the following data:
  !>   - block(b,l)%Fi%BC.
  !>   - block(b,l)%Fj%BC.
  !>   - block(b,l)%Fk%BC.
  !> It is worth noting that the 'block_bc' tags can appear in any order.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_block_bc_file_bc(file_d,block,b,l,parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_BC), intent(INOUT):: file_d      !< File data.
  type(Type_SBlock),   intent(INOUT):: block       !< Block data.
  integer(I4P),        intent(IN)::    b           !< Local index of block to be loaded.
  integer(I4P),        intent(IN)::    l           !< Level of block to be loaded.
  type(Type_Parallel), intent(IN)::    parallel    !< Parallel data.
  type(Type_XML_Tag)::                 tag         !< Dummy XML tag for parsing file.
  character(len=1)::                   c1          !< Dummy string for parsing file.
  character(len=:), allocatable::      c2          !< Dummy string for parsing file.
  integer(I4P)::                       bb,ll,i,j,k !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(blockmap => parallel%blockmap)
    call file_d%open(action='READ') ; if (file_d%iostat/=0) return
    ! reading file until 'Boundary_Conditions' is reached
    Boundary_Conditions_Tag_Search: do
      read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return
      if (index(string=c1,substring='<')>0) then
        c2 = c1
        Tag_End_Serch: do
          read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return
          c2 = c2//c1
          if (index(string=c1,substring='>')>0) exit Tag_End_Serch
        enddo Tag_End_Serch
        if (index(string=c2,substring='<Boundary_Conditions>')>0) then
          read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return ! reading cr char
          exit Boundary_Conditions_Tag_Search
        endif
      endif
    enddo Boundary_Conditions_Tag_Search
    10 continue
    call tag%set(tag_name='block_bc')
    call tag%alloc(att_name=.true.,Na=3)
    tag%att_name(1)%vs = 'b'
    tag%att_name(2)%vs = 'l'
    tag%att_name(3)%vs = 'encoding'
    Tag_Block_BC_Parsing: do
      read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)c1 ; if (file_d%iostat/=0) return
      if (c1=='<') then
        c2 = c1
        Tag_Block_BC_End_Serch: do
          read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)c1 ; if (file_d%iostat/=0) return
          c2 = c2//c1
          if (c1=='>') exit Tag_Block_BC_End_Serch
        enddo Tag_Block_BC_End_Serch
        if     (index(string=c2,substring='<block_bc')>0) then
          call tag%set(string=c2) ; call tag%get_attributes
          bb = cton(str=tag%att_val(1)%vs,knd=1_I4P) ; ll = cton(str=tag%att_val(2)%vs,knd=1_I4P)
          if (bb==blockmap(b).and.ll==l) then
            associate(gc => block%dims%gc,Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
              do k=1-gc(5),Nk+gc(6)
                do j=1-gc(3),Nj+gc(4)
                  do i=0-gc(1),Ni+gc(2)
                    call block%Fi(i,j,k)%BC%load(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
                  enddo
                enddo
              enddo
              do k=1-gc(5),Nk+gc(6)
                do j=0-gc(3),Nj+gc(4)
                  do i=1-gc(1),Ni+gc(2)
                    call block%Fj(i,j,k)%BC%load(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
                  enddo
                enddo
              enddo
              do k=0-gc(5),Nk+gc(6)
                do j=1-gc(3),Nj+gc(4)
                  do i=1-gc(1),Ni+gc(2)
                    call block%Fk(i,j,k)%BC%load(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
                  enddo
                enddo
              enddo
            endassociate
            exit Tag_Block_BC_Parsing
          else
            cycle Tag_Block_BC_Parsing
          endif
        elseif (index(string=c2,substring='</Boundary_Conditions')>0) then
          exit Tag_Block_BC_Parsing
        endif
      endif
    enddo Tag_Block_BC_Parsing
    20 call file_d%close
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_block_bc_file_bc

  !> @brief Procedure for saving boundary conditions file.
  !> @note The definition of the boundary conditions file syntax is described in the 'load_file_bc' documentation.
  !> @note For a valid file, this procedure must be called a first time for opening the 'Boundary_Conditions' tag with only the
  !> argument "bc_tag_open=.true." and finally, after all blocks bc data have been saved, with only the argument
  !> "bc_tag_close=.true." for closing the main tag.
  subroutine save_block_bc_file_bc(file_d,bc_tag_open,bc_tag_close,block,b,l,parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_BC),           intent(INOUT):: file_d       !< File data.
  logical,             optional, intent(IN)::    bc_tag_open  !< Switch for opening the 'Boundary_Conditions' tag (and init file).
  logical,             optional, intent(IN)::    bc_tag_close !< Switch for closing the 'Boundary_Conditions' tag.
  type(Type_SBlock),   optional, intent(IN)::    block        !< Block data.
  integer(I4P),        optional, intent(IN)::    b            !< Local index of block to be saved.
  integer(I4P),        optional, intent(IN)::    l            !< Level of block to be saved.
  type(Type_Parallel), optional, intent(IN)::    parallel     !< Parallel data.
  integer(I4P)::                                 i,j,k        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if     (present(bc_tag_open)) then
    if (bc_tag_open) then
      call file_d%open(replace=.true.,action='WRITE') ; if (file_d%iostat/=0) return
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<?xml version="1.0"?>'//cr
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<Boundary_Conditions>'//cr
      call file_d%close
    endif
  elseif (present(bc_tag_close)) then
    if (bc_tag_close) then
      call file_d%open(append=.true.,action='WRITE') ; if (file_d%iostat/=0) return
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'</Boundary_Conditions>'//cr
      call file_d%close
    endif
  elseif (present(block).and.present(b).and.present(l).and.present(parallel)) then
    associate(blockmap => parallel%blockmap)
      call file_d%open(append=.true.,action='WRITE') ; if (file_d%iostat/=0) return
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<block_bc b="'//trim(str(n=blockmap(b)))//'" l="'//&
        trim(str(n=l))//'" encoding="raw">'
      associate(gc => block%dims%gc,Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
        do k=1-gc(5),Nk+gc(6)
          do j=1-gc(3),Nj+gc(4)
            do i=0-gc(1),Ni+gc(2)
              call block%Fi(i,j,k)%BC%save(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
            enddo
          enddo
        enddo
        do k=1-gc(5),Nk+gc(6)
          do j=0-gc(3),Nj+gc(4)
            do i=1-gc(1),Ni+gc(2)
              call block%Fj(i,j,k)%BC%save(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
            enddo
          enddo
        enddo
        do k=0-gc(5),Nk+gc(6)
          do j=1-gc(3),Nj+gc(4)
            do i=1-gc(1),Ni+gc(2)
              call block%Fk(i,j,k)%BC%save(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
            enddo
          enddo
        enddo
      endassociate
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'</block_bc>'//cr
      call file_d%close
    endassociate
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block_bc_file_bc

  !> @brief Procedure for loading boundary conditions file.
  !> @note The boundary conditions file is an unformatted file with STREAM access. However it is formatted by means of XML syntax,
  !> thus string containing a tag must be terminated by means of a carriage-return (cr=char(10)) character.
  !> @note The boundary conditions file is composed as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Boundary_Conditions>
  !>   ...
  !>   </Boundary_Conditions>
  !> @endcode
  !> For 'Boundary_Conditions' tag and its value see the documentation of 'load_block_bc_file_bc'.
  !> The boundary conditions file contains the following data:
  !>   - global%block(b,l)%Fi%BC for all b in [1,Nb_tot] and l in [1,Nl].
  !>   - global%block(b,l)%Fj%BC for all b in [1,Nb_tot] and l in [1,Nl].
  !>   - global%block(b,l)%Fk%BC for all b in [1,Nb_tot] and l in [1,Nl].
  !> The blocks data are local to myrank, namely the blocks loaded are only the ones belonging to myrank. This is accomplished in
  !> a transparent way by means of the global%parallel%blockmap. It is worth nothing that the syntax of tag attributes is case
  !> sensitive, whereas the syntax of their values is case insensitive.
  subroutine load_file_bc(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_BC), intent(INOUT):: file_d !< File data.
  type(Type_Global),   intent(INOUT):: global !< Global data.
  integer(I4P)::                       b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,global%mesh_dims%Nl
    do b=1,global%mesh_dims%Nb
      call file_d%load_block(block=global%block(b,l),b=b,l=l,parallel=global%parallel)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_file_bc

  !> @brief Procedure for saving boundary conditions file.
  !> @note The definition of boundary conditions file syntax is described in the 'load_file_bc' documentation.
  !> @note It is worth noting that the header and the 'Boundary_Conditions' tag are saved only by myrank==0, while all processes
  !> (in an parallel MPI framework) save their own blocks nodes.
  subroutine save_file_bc(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_BC), intent(INOUT):: file_d !< File data.
  type(Type_Global),   intent(IN)::    global !< Global data.
  integer(I4P)::                       b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (global%parallel%myrank==0) then
    call file_d%save_block(bc_tag_open=.true.)
  endif
#ifdef _MPI
#else
  do l=1,global%mesh_dims%Nl
    do b=1,global%mesh_dims%Nb
      call file_d%save_block(block=global%block(b,l),b=b,l=l,parallel=global%parallel)
    enddo
  enddo
#endif
  if (global%parallel%myrank==0) call file_d%save_block(bc_tag_close=.true.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_file_bc
  !> @}
endmodule Data_Type_File_BC
