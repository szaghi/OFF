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
USE IR_Precision                                            ! Integers and reals precision definition.
USE Data_Type_Block_Dimensions, only: Type_Block_Dimensions ! Definition of Type_Block_Dimensions.
USE Data_Type_File_Base,        only: Type_File_Base,cr     ! Definition of Type_File_Base.
USE Data_Type_Global,           only: Type_Global           ! Definition of Type_Global.
USE Data_Type_Parallel,         only: Type_Parallel         ! Definition of Type_Parallel.
USE Data_Type_SBlock,           only: Type_SBlock           ! Definition of Type_SBlock.
USE Data_Type_Tree,             only: Type_Tree             ! Definition of Type_Tree.
USE Data_Type_XML_Tag,          only: Type_XML_Tag          ! Definition of Type_XML_Tag.
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
  !> The mesh file is an unformatted file with STREAM access. However it is formatted by means of XML syntax, thus
  !> string containing a tag must be terminated by means of a carriage-return (cr=achar(10)) character.
  !> The header of mesh file is composed as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Mesh_Header ancestors_number="Na" refinement_ratio="ref_ratio">
  !>     <block_dims ID="1"><dimensions .../></block_dims>
  !>     <block_dims ID="2"><dimensions .../></block_dims>
  !>     <block_dims ID="3"><dimensions .../></block_dims>
  !>     ...
  !>     <block_dims ID="Nb_tot"><dimensions .../></block_dims>
  !>   </Mesh_Header>
  !> @endcode
  !> The 'Mesh_Header' tag has 2 attributes 'ancestors_number' and 'refinement_ratio': these attributes are necessary for AMR
  !> (Adaptive Mesh Refinement) mesh and they define the number of original blocks (at root refinement level, level=1) and the
  !> AMR ratio (8 for octree and 4 for quadtree meshes), respectively. In case they are not defined the number of ancestor blocks
  !> is set equal to the number of blocks actually stored into the file, namely 'Nb_tot' and the AMR ratio is set to one, namely
  !> no AMR. Each 'block_dims' tag is identified by means of a unique ID-key, defined by the attribute 'ID', that is the block
  !> number for no AMR mesh and the Morton-number for AMR one. The number of blocks actually stored, 'Nb_tot', is the number of
  !> ancestor blocks (Na) for no AMR mesh or Nb_tot >= Na for refined mesh.
  !> Note that it is not necessary to explicitly declare the global blocks number, Nb_tot, that is computed parsing the header.
  !> The 'block_dims' tags can appear in any order.
  !> This file contains the following data:
  !>   - global%block_dims;
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_header_file_mesh(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh), intent(INOUT):: file_d              !< File data.
  type(Type_Global),     intent(INOUT):: global              !< Global data.
  type(Type_XML_Tag)::                   tag                 !< Dummy XML tag for parsing file.
  type(Type_Block_Dimensions)::          blkdims             !< Block dimensions.
  integer(I4P)::                         Nb_tot,Na,ref_ratio !< Counters.
  integer(I8P)::                         ID                  !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  allocate(character(len=1000):: tag%string%vs)
  call file_d%open(ascii=.true.,action='READ') ; if (file_d%iostat/=0) return
  Nb_tot    = 0
  Na        = 0
  ref_ratio = 1
  Mesh_Header: do
    read(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg,end=10)tag%string%vs ; if (file_d%iostat/=0) return
    if (index(string=tag%string%vs,substring='<Mesh_Header')>0) then
      call tag%set(tag_name='Mesh_Header',att_name=['ancestors_number','refinement_ratio']) ; call tag%get_attributes
      if (allocated(tag%att_val(1)%vs)) Na        = cton(str=tag%att_val(1)%vs,knd=1_I4P)
      if (allocated(tag%att_val(2)%vs)) ref_ratio = cton(str=tag%att_val(2)%vs,knd=1_I4P)
      ! if Na is not defined it is computed counting the number of blocks stored
      ! if ref_ratio is not defined it is left equal to one, namely no AMR
    elseif (index(string=tag%string%vs,substring='<block_dims')>0) then
      Nb_tot = Nb_tot + 1
    elseif (index(string=tag%string%vs,substring='</Mesh_Header')>0) then
      exit Mesh_Header
    endif
  enddo Mesh_Header
  10 rewind(unit=file_d%unit)
  if (Na==0) Na = Nb_tot
  call global%block_dims%init(Na=Na,ref_ratio=ref_ratio)
  Mesh_Dims_Read: do
    read(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg,end=20)tag%string%vs ; if (file_d%iostat/=0) return
    if     (index(string=tag%string%vs,substring='<block_dims')>0) then
      call tag%set(tag_name='block_dims',att_name=['ID']) ; call tag%get_attributes ; call tag%get_value
      ID = cton(str=tag%att_val(1)%vs,knd=1_I8P)
      call blkdims%load_str_xml(str_xml=tag%tag_val%vs)
      call global%block_dims%put(ID=ID,d=blkdims)
    elseif (index(string=tag%string%vs,substring='</Mesh_Header')>0) then
      exit Mesh_Dims_Read
    endif
  enddo Mesh_Dims_Read
  20 call file_d%close
  call global%block_dims%update
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_header_file_mesh

  !> @brief Procedure for saving header of mesh file. This procedure saves the header of mesh file that contains the mesh
  !> dimensions, namely the dimensions of all blocks (not only the ones belonging to myrank).
  !> @note The definition of the header mesh file syntax is described in the 'load_header_file_mesh' documentation.
  subroutine save_header_file_mesh(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh), intent(INOUT):: file_d  !< File data.
  type(Type_Global),     intent(IN)::    global  !< Global data.
  type(Type_Block_Dimensions), pointer:: blkdims !< Block dimensions.
  integer(I8P)::                         ID      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%open(replace=.true.,action='WRITE') ; if (file_d%iostat/=0) return
  write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<?xml version="1.0"?>'//cr
  write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<Mesh_Header '//global%block_dims%str_Na_ref_ratio()//'>'//cr
  do while(global%block_dims%loopID(ID=ID))
    blkdims => global%block_dims%dat(ID=ID)
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'  <block_dims ID="'//trim(str(n=ID))//'">'//&
      blkdims%save_str_xml()//'</block_dims>'//cr
  enddo
  write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'</Mesh_Header>'//cr
  call file_d%close
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_header_file_mesh

  !> @brief Procedure for loading block nodes coordinates from mesh file.
  !> @note After the header (see 'load_header_file_mesh' documentation) the mesh file contains the nodes coordinates of the whole
  !> blocks. Each block has its own tag containing its own nodes coordinates: this procedure load the data of a single block
  !> providing its indexes pair b-l in the local (to myrank) numeration. The blocks tags are saved inside the 'Mesh_Nodes' tag.
  !> The 'block_nodes' tag contains the following data:
  !>   - block%node.
  !> The tag syntax is the following:
  !> @code
  !>   <Mesh_Nodes>
  !>     <block_nodes ID="..." encoding="..."> nodes_data </block_nodes>
  !>   </Mesh_Nodes>
  !> @endcode
  !> Each 'block_nodes' tag has 2 attributes, 'ID' indicating the block ID-key, and 'encoding' that can be set to 3 different
  !> values: 'ascii', 'raw' or 'base64' in the case nodes_data are saved in asciii, raw-binary or Base64-ascii
  !> encoded format, respectively.
  !> It is worth noting that the 'block_nodes' tags can appear in any order.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_block_nodes_file_mesh(file_d,block,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh), intent(INOUT):: file_d !< File data.
  type(Type_SBlock),     intent(INOUT):: block  !< Block data.
  integer(I8P),          intent(IN)::    ID     !< ID-key of block to be loaded.
  type(Type_XML_Tag)::                   tag    !< Dummy XML tag for parsing file.
  character(len=1)::                     c1     !< Dummy string for parsing file.
  character(len=:), allocatable::        c2     !< Dummy string for parsing file.
  integer(I8P)::                         ib     !< Counter.
  integer(I4P)::                         i,j,k  !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
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
  tag%att_name(1)%vs = 'ID'
  tag%att_name(2)%vs = 'encoding'
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
        ib = cton(str=tag%att_val(1)%vs,knd=1_I8P)
        if (ib==ID) then
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_block_nodes_file_mesh

  !> @brief Procedure for saving blocks nodes coordinates from mesh file.
  !> @note The definition of the nodes coordinates syntax is described in the 'load_block_nodes_file_mesh' documentation.
  !> For a valid file, this procedure must be called a first time for opening the 'Mesh_Nodes' tag with only the argument
  !> "mesh_nodes_tag_open=.true." and finally, after all blocks nodes have been saved, with only the argument
  !> "mesh_nodes_tag_close=.true." for closing the main tag.
  subroutine save_block_nodes_file_mesh(file_d,mesh_nodes_tag_open,mesh_nodes_tag_close,block,ID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh),       intent(INOUT):: file_d               !< File data.
  logical,           optional, intent(IN)::    mesh_nodes_tag_open  !< Switch for opening the 'Mesh_Nodes' tag.
  logical,           optional, intent(IN)::    mesh_nodes_tag_close !< Switch for closing the 'Mesh_Nodes' tag.
  type(Type_SBlock), optional, intent(IN)::    block                !< Block data.
  integer(I8P),      optional, intent(IN)::    ID                   !< ID-key of block to be saved.
  integer(I4P)::                               i,j,k                !< Counters.
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
  elseif (present(block).and.present(ID)) then
    call file_d%open(append=.true.,action='WRITE') ; if (file_d%iostat/=0) return
    write(unit=file_d%unit,iostat=file_d%iostat,iomsg=file_d%iomsg)'<block_nodes ID="'//trim(str(n=ID))//'" encoding="raw">'
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
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block_nodes_file_mesh

  !> @brief Procedure for loading mesh file.
  !> @note The mesh file is an unformatted file with STREAM access. However it is formatted by means of XML syntax, thus
  !> string containing a tag must be terminated by means of a carriage-return (cr=char(10)) character.
  !> The mesh file contains the following data:
  !>   - global%block_dims;
  !>   - global%block%node for all blocks.
  !> The mesh file is composed as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Mesh_Header ancestors_number="Na" refinement_ratio="ref_ratio">
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
  !> to myrank. This is accomplished in a transparent way.
  !> It is worth nothing that the syntax of tag attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_file_mesh(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh), intent(INOUT):: file_d !< File data.
  type(Type_Global),     intent(INOUT):: global !< Global data.
  type(Type_SBlock), pointer::           block  !< Block pointer for scanning global%block tree.
  integer(I8P)::                         ID     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%load_header(global=global)
  call global%parallel%compute_BPmap(block_dims=global%block_dims)
  call global%blocks_init
  do while(global%block%loopID(ID=ID))
    block => global%block%dat(ID=ID)
    call file_d%load_block(block=block,ID=ID)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_file_mesh

  !> @brief Procedure for saving mesh file.
  !> @note The definition of mesh file syntax is described in the 'load_file_mesh' documentation.
  !> It is worth noting that the header and the 'Mesh_Nodes' tag are saved only by myrank==0, while all processes (in an
  !> parallel MPI framework) save their own blocks nodes.
  subroutine save_file_mesh(file_d,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Mesh), intent(INOUT):: file_d !< File data.
  type(Type_Global),     intent(IN)::    global !< Global data.
  type(Type_SBlock), pointer::           block  !< Block pointer for scanning global%block tree.
  integer(I8P)::                         ID     !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (global%parallel%is_master()) then
    call file_d%save_header(global=global)
    call file_d%save_block(mesh_nodes_tag_open=.true.)
  endif
#ifdef _MPI
#else
  do while(global%block%loopID(ID=ID))
    block => global%block%dat(ID=ID)
    call file_d%save_block(block=block,ID=ID)
  enddo
#endif
  if (global%parallel%is_master()) call file_d%save_block(mesh_nodes_tag_close=.true.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_file_mesh
  !> @}
endmodule Data_Type_File_Mesh
