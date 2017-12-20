!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_File_Blocks_CartesianDerivedType Data_Type_File_Blocks_Cartesian
!> Module definition of Type_File_Blocks_Cartesian
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_File_Blocks_CartesianPrivateProcedure Data_Type_File_Blocks_Cartesian
!> Module definition of Type_File_Blocks_Cartesian
!> @}

!> @brief Module Data_Type_File_Blocks_Cartesian contains the definition of Type_File_Blocks_Cartesian, that is
!> a simple ascii file describing a Cartesian structured blocks.
module Data_Type_File_Blocks_Cartesian
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                 ! Integers and reals precision definition.
USE Data_Type_File_Base, only: Type_File_Base    ! Definition of Type_File_Base.
USE Data_Type_Global,    only: Type_Global       ! Definition of Type_Global.
USE Data_Type_Region,    only: Type_Region       ! Definition of Type_Region.
USE Data_Type_SBlock,    only: Type_SBlock       ! Definition of Type_SBlock.
USE Data_Type_XML_Tag,   only: Type_XML_Tag      ! Definition of Type_XML_Tag.
USE Lib_IO_Misc,         only: stderr,iostat_eor ! IO error utilities.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Blocks_Cartesian.
!> @note This options file is a very simple XML file and it is useful only for Cartesian grids.
!> @note This options file is formatted by means of XML syntax. A validated file should be formatted as following:
!> @code
!>   <?xml version="1.0"?>
!>   <IBM_Blocks>
!>     <Blocks_Cartesian>
!>       <block b="1">
!>         <dimensions .../>
!>         <extents .../>
!>         <boundaryconditions ../>
!>         <primitives .../>
!>       </block>
!>       <block b="2">
!>         <dimensions .../>
!>         <extents .../>
!>         <boundaryconditions ../>
!>         <primitives .../>
!>       </block>
!>       ...
!>     </Blocks_Cartesian>
!>     <Regions>
!>       <region r="1">
!>         <shape .../>
!>         <primitives .../>
!>       </region>
!>       <region r="2">
!>         <shape .../>
!>         <primitives .../>
!>       </region>
!>       ...
!>     </Regions>
!>   </IBM_Blocks>
!> @endcode
!> The main tag is named 'IBM_Blocks'. It can contain two nested tags, 'Blocks_Cartesian' that is mandatory and 'Regions' that is
!> optional. The tag 'Blocks_Cartesian' contains the blocks description. The blocks description can appear
!> in any order, the blocks numeration is made by means of the attribute 'b' of the 'block' tag. Each block tag must be formatted
!> with exactly xml tags as following:
!> @code
!>   <block b="#block_number">
!>     <dimensions gc="#imin #imax #jmin #jmax #kmin #kmax" Ni="#cell-i" Nj="#cell-j" Nk="#cell-k"/>
!>     <extents xmin="#xmin" ymin="#ymin" zmin="#zmin" xmax="#xmax" ymax="#ymax" zmax="#zmax"/>
!>     <boundaryconditions imin="#imin [d]" imax="#imax [d]" jmin="#jmin [d]" jmax="#jmax [d]" kmin="#kmin [d]" kmax="#kmax [d]"/>
!>     <primitives r="rho#1 rho#2 ... rho#Ns" u="x-velocity" v="y-velocity" w="z-velocity" p="pressure"/>
!>   </block>
!> @endcode
!> The parsing of block description tags is performed by means of proper parser that is defined within the corresponding type,
!> namely the types 'Type_Block_Dimensions', 'Type_Block_Extents', 'Type_Block_BC' and 'Type_Primitives'. The above 4 tags can
!> appear in any order as well as their attributes. The tag 'Regions' is optional and it can be used to set particular regions
!> of the blocks to a different initial conditions. Each region tag must be formatted with exactly xml tags as following:
!> @code
!>   <region r="#region_number">
!>     <primitives r="rho#1 rho#2 ... rho#Ns" u="x-velocity" v="y-velocity" w="z-velocity" p="pressure"/>
!>     <shape tp="type_of_region" center="#x #y #z" radius="#radius" height="#height"/>
!>   </region>
!> @endcode
!> The 'primitives' tag is parsed as in the 'block' tag, while the 'shape' tag is directly parsed here. Presently, only the
!> 'cylinder' shape type is implemented. The above 2 tags can appear in any order as well as their attributes. The regions
!> description can appear in any order, the regions numeration is made by means of the attribute 'r' of the 'region' tag. The
!> total number of regions defined is automatically computed parsing the whole file
!> @note This file contains the following data:
!>   - file_d%block.
!> It is worth nothing that the syntax of tag names and attributes is case sensitive, whereas the syntax of their values is case
!> insensitive.
!> @ingroup Data_Type_File_Blocks_CartesianDerivedType
type, public, extends(Type_File_Base):: Type_File_Blocks_Cartesian
  integer(I4P)::                   Nb = 0_I4P !< Number of blocks (IDs into AMR mesh).
  type(Type_SBlock), allocatable:: block(:)   !< Blocks input data [1:mesh_dims%Nb].
  integer(I4P)::                   Nr = 0_I4P !< Number of regions.
  type(Type_Region), allocatable:: region(:)  !< Regions input data [1:Nr].
  contains
    procedure:: free  => free_blocks   ! Procedure for freeing dynamic memory.
    procedure:: alloc => alloc_blocks  ! Procedure for allocating dynamic memory.
    procedure:: load  => load_blocks   ! Procedure for loading IBM blocks description.
    procedure:: save  => save_blocks   ! Procedure for saving IBM blocks description.
    procedure:: compute_mesh_dims      ! Procedure for computing mesh dimensions.
    final::     finalize               ! Procedure for freeing dynamic memory when finalizing.
endtype Type_File_Blocks_Cartesian
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_Blocks_CartesianPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_blocks(file_d,also_base)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Blocks_Cartesian), intent(INOUT):: file_d    !< File data.
  logical, optional,                 intent(IN)::    also_base !< Flag for freeing also file_base dynamic memory.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(also_base)) then
    if (also_base) call file_d%free_base
  endif
  if (allocated(file_d%block))  then
    call file_d%block%free ; deallocate(file_d%block)
  endif
  if (allocated(file_d%region))  then
    call file_d%region%free ; deallocate(file_d%region)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_blocks

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_File_Blocks_Cartesian), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free(also_base=.true.)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc_blocks(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Blocks_Cartesian), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free
  if (file_d%Nb>0) allocate(file_d%block( 1:file_d%Nb))
  if (file_d%Nr>0) allocate(file_d%region(1:file_d%Nr))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_blocks

  !> @brief Procedure for loading Cartesian blocks description file.
  subroutine load_blocks(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Blocks_Cartesian), intent(INOUT):: file_d         !< File data.
  character(len=:), allocatable::                    stream         !< String containing the file data as a single stream.
  type(Type_XML_Tag)::                               tag            !< Main XML tag.
  integer(I4P)::                                     t,tt,b,bb,r,rr !< Counters.
  !--------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    call file_d%get_stream(stream=stream,delimiter_start='<IBM_Blocks',delimiter_end='</IBM_Blocks>')!; if (iostat/=0) return
    call tag%set(tag_name='IBM_Blocks',search_into_string=stream)
    call tag%get_nested
    if (associated(tag%nested)) then
      ! computing the number of blocks and regions
      do t=1,size(tag%nested)
        call tag%nested(t)%get_nested
        select case(tag%nested(t)%tag_name%vs)
        case('Blocks_Cartesian')
          file_d%Nb = size(tag%nested(t)%nested)
        case('Regions')
          file_d%Nr = size(tag%nested(t)%nested)
        endselect
      enddo
      call file_d%alloc
      ! loading all data
      do t=1,size(tag%nested)
        select case(tag%nested(t)%tag_name%vs)
        case('Blocks_Cartesian')
          do b=1,size(tag%nested(t)%nested)
            call tag%nested(t)%nested(b)%get_nested
            bb = cton(str=tag%nested(t)%nested(b)%att_val(1)%vs,knd=1_I4P)
            associate(block => file_d%block(bb))
              do tt=1,size(tag%nested(t)%nested(b)%nested)
                select case(tag%nested(t)%nested(b)%nested(tt)%tag_name%vs)
                case('dimensions')
                  call block%dims%load_str_xml(str_xml=tag%nested(t)%nested(b)%nested(tt)%string%vs)
                case('extents')
                  call block%exts%load_str_xml(str_xml=tag%nested(t)%nested(b)%nested(tt)%string%vs)
                case('boundaryconditions')
                  call block%BC%load_str_xml(str_xml=tag%nested(t)%nested(b)%nested(tt)%string%vs)
                case('primitives')
                  call block%IC%load_str_xml(str_xml=tag%nested(t)%nested(b)%nested(tt)%string%vs)
                endselect
              enddo
            endassociate
          enddo
        case('Regions')
          do r=1,size(tag%nested(t)%nested)
            call tag%nested(t)%nested(r)%get_nested
            rr = cton(str=tag%nested(t)%nested(r)%att_val(1)%vs,knd=1_I4P)
            associate(region => file_d%region(rr))
              do tt=1,size(tag%nested(t)%nested(r)%nested)
                select case(tag%nested(t)%nested(r)%nested(tt)%tag_name%vs)
                case('shape')
                  call region%load_shape_str_xml(str_xml=tag%nested(t)%nested(r)%nested(tt)%string%vs)
                case('primitives')
                  call region%prim%load_str_xml(str_xml=tag%nested(t)%nested(r)%nested(tt)%string%vs)
                endselect
              enddo
            endassociate
          enddo
        endselect
      enddo
    endif
  endassociate
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_blocks

  !> @brief Procedure for saving IBM blocks description file.
  subroutine save_blocks(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Blocks_Cartesian), intent(INOUT):: file_d !< File data.
  integer(I4P)::                                     b,r    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%open(ascii=.true.,action='WRITE') ; if (file_d%iostat/=0) return
  write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'<?xml version="1.0"?>'
  write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'<Blocks_Cartesian>'
  if (file_d%Nb > 0) then
    do b=1,file_d%Nb
      associate(block => file_d%block(b))
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'<block b="'//trim(str(n=b))//'">'
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)block%dims%save_str_xml()
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)block%exts%save_str_xml()
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)block%BC%save_str_xml()
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)block%IC%save_str_xml()
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'</block>'
      endassociate
    enddo
  endif
  write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'</Blocks_Cartesian>'
  if (file_d%Nr > 0) then
    do b=1,file_d%Nr
      associate(region => file_d%region(r))
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'<region r="'//trim(str(n=r))//'">'
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)region%save_shape_str_xml()
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)region%prim%save_str_xml()
        write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'</region>'
      endassociate
    enddo
  endif
  call file_d%close
  call file_d%fallback
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_blocks

  !> @brief Procedure for computing the mesh dimensions associated to the Cartesian blocks loaded.
  subroutine compute_mesh_dims(file_d,global,ref_ratio)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Blocks_Cartesian), intent(INOUT):: file_d    !< File data.
  type(Type_Global),                 intent(INOUT):: global    !< Global data.
  integer(I4P),                      intent(IN)::    ref_ratio !< Refinement ratio.
  integer(I4P)::                                     b         !< Counter.
  !--------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(block=>file_d%block,Nb_tot=>file_d%Nb)
    call global%block_dims%init(Na=Nb_tot,ref_ratio=ref_ratio)
    do b=1,Nb_tot
      call global%block_dims%put(ID=int(b,kind=I8P),d=block(b)%dims)
    enddo
    call global%block_dims%update
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_mesh_dims
  !> @}
endmodule Data_Type_File_Blocks_Cartesian
