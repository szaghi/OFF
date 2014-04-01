!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_MeshDerivedType Data_Type_File_Mesh
!> Module definition of main mesh file structures, Type_File_Mesh
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_MeshPrivateProcedure Data_Type_File_Mesh
!> Module definition of main mesh file structures, Type_File_Mesh
!> @}

!> @brief Module Data_Type_File_Mesh contains the definition of Type_File_Mesh, that is the main mesh file structures adopted by
!> OFF project.
!> @todo \b encodin: Implement ascii and base64 encoding besides raw one
module Data_Type_File_Mesh
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                          ! Integers and reals precision definition.
USE Data_Type_File_Base,       only: Type_File_Base,cr    ! Definition of Type_File_Base.
USE Data_Type_Global,          only: Type_Global          ! Definition of Type_Global.
USE Data_Type_Mesh_Dimensions, only: Type_Mesh_Dimensions ! Definition of Type_Mesh_Dimensions.
USE Data_Type_Parallel,        only: Type_Parallel        ! Definition of Type_Parallel.
USE Data_Type_SBlock,          only: Type_SBlock          ! Definition of Type_SBlock.
USE Data_Type_XML_Tag,         only: Type_XML_Tag         ! Definition of Type_XML_Tag.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Mesh.
!> @ingroup Data_Type_File_MeshDerivedType
type, public, extends(Type_File_Base):: Type_File_Mesh
  contains
    procedure:: load_header => load_header_file_mesh      ! Procedure for loading header from mesh file.
    procedure:: save_header => save_header_file_mesh      ! Procedure for saving header into mesh file.
    procedure:: load_block  => load_block_nodes_file_mesh ! Procedure for loading block nodes coordinates from mesh file.
    procedure:: save_block  => save_block_nodes_file_mesh ! Procedure for saving bloack nodes coordinates into mesh file.
    procedure:: load        => load_file_mesh             ! Procedure for loading mesh file.
    procedure:: save        => save_file_mesh             ! Procedure for saving mesh file.
endtype Type_File_Mesh
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_MeshPrivateProcedure
  !> @{
  !> @brief Procedure for loading header of mesh file. This procedure loads the header of mesh file that contains the mesh
  !> dimensions, namely the dimensions of all blocks (not only the ones belonging to myrank).
  !> @note The header of this file contains the following data:
  !>   - mesh_dim.
  !> @note The mesh file is an unformatted file with STREAM access. However it is formatted by means of XML syntax, thus
  !> string containing a tag must be terminated by means of a carriage-return (cr=char(10)) character.
  !> @note The header of mesh file is composed as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Mesh_Header>
  !>     <block_dims b="1" l="1"><dimensions .../></block_dims>
  !>     <block_dims b="2" l="1"><dimensions .../></block_dims>
  !>     <block_dims b="3" l="1"><dimensions .../></block_dims>
  !>     ...
  !>     <block_dims b="Nb_tot" l="1"><dimensions .../></block_dims>
  !>     <block_dims b="1" l="2"><dimensions .../></block_dims>
  !>     <block_dims b="2" l="2"><dimensions .../></block_dims>
  !>     <block_dims b="3" l="2"><dimensions .../></block_dims>
  !>     ...
  !>     <block_dims b="Nb_tot" l="2"><dimensions .../></block_dims>
  !>     ...
  !>     <block_dims b="1" l="Nl"><dimensions .../></block_dims>
  !>     <block_dims b="2" l="Nl"><dimensions .../></block_dims>
  !>     <block_dims b="3" l="Nl"><dimensions .../></block_dims>
  !>     ...
  !>     <block_dims b="Nb_tot" l="Nl"><dimensions .../></block_dims>
  !>   </Mesh_Header>
  !> @endcode
  !> @note Note that it is not necessary to explicitely declare the global blocks number, Nb_tot, and number of grids levels that
  !> are computed parsing the header. The 'block_dims' tags can appear in any order.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_header_file_mesh(file_d,mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh),      intent(INOUT):: file_d    !< File data.
  type(Type_Mesh_Dimensions), intent(INOUT):: mesh_dims !< File data.
  type(Type_XML_Tag)::                        tag       !< Dummy XML tag for parsing file.
  integer(I4P)::                              b,l       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(character(len=1000):: tag%string%vs)
  call file_d%open(ascii=.true.,action='READ') ; if (file_d%iostat/=0) return
  Nl_Nb_tot_Read: do
    read(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)tag%string%vs ; if (file_d%iostat/=0) return
    if (index(string=tag%string%vs,substring='<block_dims')>0) then
      call tag%set(tag_name='block_dims',att_name=['b','l']) ; call tag%get_attributes
      b = cton(str=tag%att_val(1)%vs,knd=1_I4P) ; l = cton(str=tag%att_val(2)%vs,knd=1_I4P)
      mesh_dims%Nb_tot = max(b,mesh_dims%Nb_tot) ; mesh_dims%Nl = max(l,mesh_dims%Nl)
    elseif (index(string=tag%string%vs,substring='</Mesh_Header')>0) then
      exit Nl_Nb_tot_Read
    endif
  enddo Nl_Nb_tot_Read
  10 rewind(unit=file_d%unit)
  call mesh_dims%alloc
  Mesh_Dims_Read: do
    read(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)tag%string%vs ; if (file_d%iostat/=0) return
    if     (index(string=tag%string%vs,substring='<block_dims')>0) then
      call tag%set(tag_name='block_dims',att_name=['b','l']) ; call tag%get_attributes ; call tag%get_value
      b = cton(str=tag%att_val(1)%vs,knd=1_I4P) ; l = cton(str=tag%att_val(2)%vs,knd=1_I4P)
      call mesh_dims%block_dims(b,l)%load_str_xml(str_xml=tag%tag_val%vs)
    elseif (index(string=tag%string%vs,substring='</Mesh_Header')>0) then
      exit Mesh_Dims_Read
    endif
  enddo Mesh_Dims_Read
  20 call file_d%close
  mesh_dims%set = .true.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_header_file_mesh

  !> @brief Procedure for saving header of mesh file. This procedure saves the header of mesh file that contains the mesh
  !> dimensions, namely the dimensions of all blocks (not only the ones belonging to myrank).
  !> @note The definition of the header mesh file syntax is described in the 'load_header_file_mesh' documentation.
  subroutine save_header_file_mesh(file_d,mesh_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh),      intent(INOUT):: file_d    !< File data.
  type(Type_Mesh_Dimensions), intent(IN)::    mesh_dims !< File data.
  integer(I4P)::                              b,l       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%open(replace=.true.,action='WRITE') ; if (file_d%iostat/=0) return
  write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<?xml version="1.0"?>'//cr
  write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<Mesh_Header>'//cr
  do l=1,mesh_dims%Nl ; do b=1,mesh_dims%Nb_tot
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'  <block_dims b="'//trim(str(n=b))//&
      '" l="'//trim(str(n=l))//'">'//mesh_dims%block_dims(b,l)%save_str_xml()//'</block_dims>'//cr
  enddo ; enddo
  write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'</Mesh_Header>'//cr
  call file_d%close
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_header_file_mesh

  !> @brief Procedure for loading block nodes coordinates from mesh file.
  !> @note After the header (see 'load_header_file_mesh' documentation) the mesh file contains the nodes coordinates of the whole
  !> blocks. Each block has its own tag containing its own nodes coordinates: this procedure load the data of a single block
  !> providing its indexes pair b-l in the local (to myrank) numeration. The blocks tags are saved inside the 'Mesh_Nodes' tag.
  !> @note The 'block_nodes' tag contains the following data:
  !>   - block%node.
  !> The tag syntax is the following:
  !>   <Mesh_Nodes>
  !>     <block_nodes b="..." l="..." encoding="..."> nodes_data </block_nodes>
  !>   </Mesh_Nodes>
  !> Each 'block_nodes' tag has 3 attributes, 'b' and 'l' indicating the block number and grid level, and 'encoding' that can be
  !> set to 3 different values: 'ascii', 'raw' or 'base64' in the case nodes_data are saved in asciii, raw-binary or Base64-ascii
  !> encoded format, respectively.
  !> It is worth noting that the 'block_nodes' tags can appear in any order.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_block_nodes_file_mesh(file_d,block,b,l,parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh), intent(INOUT):: file_d      !< File data.
  type(Type_SBlock),     intent(INOUT):: block       !< Block data.
  integer(I4P),          intent(IN)::    b,l         !< Local index of block to be loaded.
  type(Type_Parallel),   intent(IN)::    parallel    !< Parallel data.
  type(Type_XML_Tag)::                   tag         !< Dummy XML tag for parsing file.
  character(len=1)::                     c1          !< Dummy string for parsing file.
  character(len=:), allocatable::        c2          !< Dummy string for parsing file.
  integer(I4P)::                         bb,ll,i,j,k !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(blockmap => parallel%blockmap)
    call file_d%open(action='READ') ; if (file_d%iostat/=0) return
    ! reading file until 'Mesh_Nodes' is reached
    Mesh_Nodes_Tag_Search: do
      read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return
      if (index(string=c1,substring='<')>0) then
        c2 = c1
        Tag_End_Serch: do
          read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return
          c2 = c2//c1
          if (index(string=c1,substring='>')>0) exit Tag_End_Serch
        enddo Tag_End_Serch
        if (index(string=c2,substring='<Mesh_Nodes>')>0) then
          read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)c1 ; if (file_d%iostat/=0) return ! reading cr char
          exit Mesh_Nodes_Tag_Search
        endif
      endif
    enddo Mesh_Nodes_Tag_Search
    10 continue
    call tag%set(tag_name='block_nodes')
    call tag%alloc(att_name=.true.,Na=3)
    tag%att_name(1)%vs = 'b'
    tag%att_name(2)%vs = 'l'
    tag%att_name(3)%vs = 'encoding'
    Tag_Block_Nodes_Parsing: do
      read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)c1 ; if (file_d%iostat/=0) return
      if (c1=='<') then
        c2 = c1
        Tag_Block_Nodes_End_Serch: do
          read(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)c1 ; if (file_d%iostat/=0) return
          c2 = c2//c1
          if (c1=='>') exit Tag_Block_Nodes_End_Serch
        enddo Tag_Block_Nodes_End_Serch
        if     (index(string=c2,substring='<block_nodes')>0) then
          call tag%set(string=c2) ; call tag%get_attributes
          bb = cton(str=tag%att_val(1)%vs,knd=1_I4P) ; ll = cton(str=tag%att_val(2)%vs,knd=1_I4P)
          if (bb==blockmap(b).and.ll==l) then
            associate(gc => block%dims%gc,Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
              do k=0-gc(5),Nk+gc(6)
                do j=0-gc(3),Nj+gc(4)
                  do i=0-gc(1),Ni+gc(2)
                    call block%node(i,j,k)%load(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
                  enddo
                enddo
              enddo
            endassociate
            exit Tag_Block_Nodes_Parsing
          else
            cycle Tag_Block_Nodes_Parsing
          endif
        elseif (index(string=c2,substring='</Mesh_Nodes')>0) then
          exit Tag_Block_Nodes_Parsing
        endif
      endif
    enddo Tag_Block_Nodes_Parsing
    20 call file_d%close
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_block_nodes_file_mesh

  !> @brief Procedure for saving blocks nodes coordinates from mesh file.
  !> @note The definition of the nodes coordinates syntax is described in the 'load_block_nodes_file_mesh' documentation.
  !> @note For a valid file, this procedure must be called a first time for opening the 'Mesh_Nodes' tag with only the argument
  !> "mesh_nodes_tag_open=.true." and finally, after all blocks nodes have been saved, with only the argument
  !> "mesh_nodes_tag_close=.true." for closing the main tag.
  subroutine save_block_nodes_file_mesh(file_d,mesh_nodes_tag_open,mesh_nodes_tag_close,block,b,l,parallel)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh),         intent(INOUT):: file_d               !< File data.
  logical,             optional, intent(IN)::    mesh_nodes_tag_open  !< Switch for opening the 'Mesh_Nodes' tag.
  logical,             optional, intent(IN)::    mesh_nodes_tag_close !< Switch for closing the 'Mesh_Nodes' tag.
  type(Type_SBlock),   optional, intent(IN)::    block                !< Block data.
  integer(I4P),        optional, intent(IN)::    b,l                  !< Local index of block to be saved.
  type(Type_Parallel), optional, intent(IN)::    parallel             !< Parallel data.
  integer(I4P)::                                 i,j,k                !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if     (present(mesh_nodes_tag_open )) then
    if (mesh_nodes_tag_open) then
      call file_d%open(append=.true.,action='WRITE') ; if (file_d%iostat/=0) return
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<Mesh_Nodes>'//cr
      call file_d%close
    endif
  elseif (present(mesh_nodes_tag_close)) then
    if (mesh_nodes_tag_close) then
      call file_d%open(append=.true.,action='WRITE') ; if (file_d%iostat/=0) return
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'</Mesh_Nodes>'//cr
      call file_d%close
    endif
  elseif (present(block).and.present(b).and.present(l).and.present(parallel)) then
    associate(blockmap => parallel%blockmap)
      call file_d%open(append=.true.,action='WRITE') ; if (file_d%iostat/=0) return
      write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<block_nodes b="'//trim(str(n=blockmap(b)))//&
        '" l="'//trim(str(n=l))//'" encoding="raw">'
        associate(gc => block%dims%gc,Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
          do k=0-gc(5),Nk+gc(6)
            do j=0-gc(3),Nj+gc(4)
              do i=0-gc(1),Ni+gc(2)
                call block%node(i,j,k)%save(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)
              enddo
            enddo
          enddo
        endassociate
        write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'</block_nodes>'//cr
      call file_d%close
    endassociate
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block_nodes_file_mesh

  !> @brief Procedure for loading mesh file.
  !> @note The mesh file is an unformatted file with STREAM access. However it is formatted by means of XML syntax, thus
  !> string containing a tag must be terminated by means of a carriage-return (cr=char(10)) character.
  !> @note The mesh file contains the following data:
  !>   - global%mesh_dim;
  !>   - global%block(b,l)%node for all b in [1,Nb_tot] and l in [1,Nl].
  !> @note The mesh file is composed as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Mesh_Header>
  !>   ...
  !>   </Mesh_Header>
  !>   <Mesh_Nodes>
  !>   ...
  !>   </Mesh_Nodes>
  !> @endcode
  !> For 'Mesh_Header' and 'Mesh_Nodes' tags see the documentation of 'load_header_file_mesh' and 'load_nodes_file_mesh',
  !> respectively.
  !> Note that the header information must be loaded by all processes into a parallel MPI framework thus there is no control
  !> relying on myrank. On the contrary, the blocks data are local to myrank, namely the blocks loaded are only the ones belonging
  !> to myrank. This is accomplished in a transparent way by means of the global%parallel%blockmap.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_file_mesh(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh), intent(INOUT):: file_d !< File data.
  type(Type_Global),     intent(INOUT):: global !< Global data.
  integer(I4P)::                         b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%load_header(mesh_dims=global%mesh_dims)
  call global%alloc
  call global%blocks_dims_update
  call global%block%alloc
  do l=1,global%mesh_dims%Nl
    do b=1,global%mesh_dims%Nb
      call file_d%load_block(block=global%block(b,l),b=b,l=l,parallel=global%parallel)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_file_mesh

  !> @brief Procedure for saving mesh file.
  !> @note The definition of mesh file syntax is described in the 'load_file_mesh' documentation.
  !> @note It is worth noting that the header and the 'Mesh_Nodes' tag are saved only by myrank==0, while all processes (in an
  !> parallel MPI framework) save their own blocks nodes.
  subroutine save_file_mesh(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh), intent(INOUT):: file_d !< File data.
  type(Type_Global),     intent(IN)::    global !< Global data.
  integer(I4P)::                         b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (global%parallel%myrank==0) then
    call file_d%save_header(mesh_dims=global%mesh_dims)
    call file_d%save_block(mesh_nodes_tag_open=.true.)
  endif
#ifdef _MPI
#else
  do l=1,global%mesh_dims%Nl
    do b=1,global%mesh_dims%Nb
      call file_d%save_block(block=global%block(b,l),b=b,l=l,parallel=global%parallel)
    enddo
  enddo
#endif
  if (global%parallel%myrank==0) call file_d%save_block(mesh_nodes_tag_close=.true.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_file_mesh
  !> @}
endmodule Data_Type_File_Mesh
