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
USE IR_Precision                                          ! Integers and reals precision definition.
USE Data_Type_File_Base,       only: Type_File_Base       ! Definition of Type_File_Base.
USE Data_Type_Mesh_Dimensions, only: Type_Mesh_Dimensions ! Definition of Type_Mesh_Dimensions.
USE Data_Type_Region,          only: Type_Region          ! Definition of Type_Region.
USE Data_Type_SBlock,          only: Type_SBlock          ! Definition of Type_SBlock.
USE Lib_IO_Misc,               only: stderr,iostat_eor    ! IO error utilities.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_File_Blocks_Cartesian.
!> @ingroup Data_Type_File_Blocks_CartesianDerivedType
type, public, extends(Type_File_Base):: Type_File_Blocks_Cartesian
  type(Type_Mesh_Dimensions)::     mesh_dims  !< Mesh dimensions associated to the file.
  type(Type_SBlock), allocatable:: block(:)   !< Blocks input data [1:mesh_dims%Nb].
  integer(I4P)::                   Nr = 0_I4P !< Number of regions.
  type(Type_Region), allocatable:: region(:)  !< Regions input data [1:Nr].
  contains
    procedure:: free  => free_blocks   ! Procedure for freeing dynamic memory.
    procedure:: alloc => alloc_blocks  ! Procedure for allocating dynamic memory.
    procedure:: load  => load_blocks   ! Procedure for loading IBM blocks description.
    procedure:: save  => save_blocks   ! Procedure for saving IBM blocks description.
    procedure:: compute_mesh_dims      ! Procedure for saving IBM blocks into scratch files.
    final::     finalize               ! Procedure for freeing dynamic memory when finalizing.
endtype Type_File_Blocks_Cartesian
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_File_Blocks_CartesianPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  elemental subroutine free_blocks(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Blocks_Cartesian), intent(INOUT):: file_d !< File data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call file_d%free_base
  if (allocated(file_d%block))  then
    call file_d%block%free ; deallocate(file_d%block)
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
  call file_d%free
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
  if (allocated(file_d%block))  then
    call file_d%block%free ; deallocate(file_d%block)
  endif
  allocate(file_d%block(1:file_d%mesh_dims%Nb))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_blocks

  !> @brief Procedure for loading Cartesian blocks description file.
  !> @note This options file is a very simple XML file and it is useful only for Cartesian grids.
  !> @note This options file is formatted by means of XML syntax. A validated file should be formatted as following:
  !> @code
  !>   <?xml version="1.0"?>
  !>   <Main_Tag>
  !>     <Initial_Species>
  !>       <Specie s="1">...</Specie>
  !>       <Specie s="2">...</Specie>
  !>       ...
  !>       <Specie s="Ns">...</Specie>
  !>     </Initial_Species>
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
  !>   </Main_Tag>
  !> @endcode
  !> The main tag 'Main_Tag' is optional and can be set to any name, whereas the tags 'Initial_Species' and 'Blocks_Cartesian' are
  !> mandatory. The tag 'Initial_Species' contains the definition of each initial specie. Each specie if defined by a 'Specie' tag
  !> that is an 'overloading' tag of 'specie' one actually defining the specie components. The overloading tag 'Specie' introduces
  !> the attribute 's' defining the numeration of each specie into the global numeration. The total number of species defined is
  !> automatically computed parsing the whole file, thus it is not necessary to explicitly define it. The tag 'specie' is defined
  !> into 'Data_Type_Specie' module. The tag 'Blocks_Cartesian' contains the blocks description. The blocks description can appear
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
  !>   - file_d%mesh_dims%Nb;
  !>   - file_d%block.
  !> It is worth nothing that the syntax of tag names and attributes is case sensitive, whereas the syntax of their values is case
  !> insensitive.
  subroutine load_blocks(file_d)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Blocks_Cartesian), intent(INOUT):: file_d      !< File data.
  character(len=:), allocatable::                    line,tag    !< Dummy strings for parsing option file.
  integer(I4P)::                                     b,c1,c2,r,t !< Counters.
  !--------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(unit=>file_d%unit,iostat=>file_d%iostat,iomsg=>file_d%iomsg)
    call file_d%open(ascii=.true.,action='READ') ; if (file_d%iostat/=0) return
    allocate(character(len=1000):: line,tag)
    file_d%mesh_dims%Nb = 0
    Blocks_Count: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=10)line ; if (iostat/=0.and.iostat/=iostat_eor) return
      if (index(string=line,substring='<block')>0) file_d%mesh_dims%Nb = file_d%mesh_dims%Nb + 1
    enddo Blocks_Count
    10 rewind(unit=unit)
    file_d%mesh_dims%Nb_tot = file_d%mesh_dims%Nb
    call file_d%alloc
    file_d%Nr = 0
    Regions_Count: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=20)line ; if (iostat/=0.and.iostat/=iostat_eor) return
      if (index(string=line,substring='<region')>0) file_d%Nr = file_d%Nr + 1
    enddo Regions_Count
    20 rewind(unit=unit)
    if (file_d%Nr>0) then
      if (allocated(file_d%region)) deallocate(file_d%region) ; allocate(file_d%region(1:file_d%Nr))
    endif
    Read_Loop: do
      read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=30)line ; if (iostat/=0.and.iostat/=iostat_eor) return
      if (index(string=line,substring='<block')>0) then
        c1=index(string=line,substring='b="')+3
        c2=index(string=line(c1:),substring='"')-1
        b = cton(str=line(c1:c1+c2-1),knd=1_I4P)
        associate(block => file_d%block(b))
          Block_Tag_Loop: do t=1,4
            read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=30)tag ; if (iostat/=0.and.iostat/=iostat_eor) return
            if     (index(string=tag,substring='<dimensions')>0) then
              call block%dims%load_str_xml(str_xml=trim(tag))
            elseif (index(string=tag,substring='<extents')>0) then
              call block%exts%load_str_xml(str_xml=trim(tag))
            elseif (index(string=tag,substring='<boundaryconditions')>0) then
              call block%BC%load_str_xml(str_xml=trim(tag))
            elseif (index(string=tag,substring='<primitives')>0) then
              call block%IC%load_str_xml(str_xml=trim(tag))
            else
              write(stderr,'(A)')' Attention: the description of block '//trim(str(.true.,b))//' is not complete!'
              exit
            endif
          enddo Block_Tag_Loop
        endassociate
      elseif (index(string=line,substring='<region')>0) then
        c1=index(string=line,substring='r="')+3
        c2=index(string=line(c1:),substring='"')-1
        r = cton(str=line(c1:c1+c2-1),knd=1_I4P)
        associate(region => file_d%region(r))
          Region_Tag_Loop: do t=1,2
            read(unit=unit,fmt='(A)',iostat=iostat,iomsg=iomsg,end=30)tag ; if (iostat/=0.and.iostat/=iostat_eor) return
            if     (index(string=tag,substring='<shape')>0) then
              call region%load_shape_str_xml(str_xml=trim(tag))
            elseif (index(string=tag,substring='<primitives')>0) then
              call region%prim%load_str_xml(str_xml=trim(tag))
            else
              write(stderr,'(A)')' Attention: the description of region '//trim(str(.true.,r))//' is not complete!'
              exit
            endif
          enddo Region_Tag_Loop
        endassociate
      endif
    enddo Read_Loop
    30 call file_d%close
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_blocks

  !> @brief Procedure for saving IBM blocks description file.
  !> @note The definition of the IBM blocks description file syntax is described in the 'load_blocks' documentation.
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
  do b=1,file_d%mesh_dims%Nb
    associate(block => file_d%block(b))
      write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'<block b="'//trim(str(n=b))//'">'
      write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)block%dims%save_str_xml()
      write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)block%exts%save_str_xml()
      write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)block%BC%save_str_xml()
      write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)block%IC%save_str_xml()
      write(unit=file_d%unit,fmt='(A)',iostat=file_d%iostat,iomsg=file_d%iomsg)'</block>'
    endassociate
  enddo
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_blocks

  !> @brief Procedure for computing the mesh dimensions associated to the Cartesian blocks loaded.
  subroutine compute_mesh_dims(file_d,Nl)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Blocks_Cartesian), intent(INOUT):: file_d !< File data.
  integer(I4P),                      intent(IN)::    Nl     !< Number of grid levels.
  integer(I4P)::                                     b,l    !< Counters.
  !--------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(mesh_dims=>file_d%mesh_dims,block=>file_d%block)
    mesh_dims%Nl = Nl
    call mesh_dims%alloc
    do l=1,mesh_dims%Nl ; do b=1,mesh_dims%Nb
      if     (mod(block(b)%dims%Ni,(2**(l-1)))/=0) then
        write(stderr,'(A)')' Attention the number of grid levels used is not consistent with the number of cells'
        write(stderr,'(A)')' Impossible to compute grid level '//trim(str(.true.,l))//' of block '//trim(str(.true.,b))
        write(stderr,'(A)')' Inconsistent direction i, Ni '//trim(str(.true.,block(b)%dims%Ni))
        stop
      elseif (mod(block(b)%dims%Nj,(2**(l-1)))/=0) then
        write(stderr,'(A)')' Attention the number of grid levels used is not consistent with the number of cells'
        write(stderr,'(A)')' Impossible to compute grid level '//trim(str(.true.,l))//' of block '//trim(str(.true.,b))
        write(stderr,'(A)')' Inconsistent direction j, Nj '//trim(str(.true.,block(b)%dims%Nj))
        stop
      elseif (mod(block(b)%dims%Nk,(2**(l-1)))/=0) then
        write(stderr,'(A)')' Attention the number of grid levels used is not consistent with the number of cells'
        write(stderr,'(A)')' Impossible to compute grid level '//trim(str(.true.,l))//' of block '//trim(str(.true.,b))
        write(stderr,'(A)')' Inconsistent direction k, Nk '//trim(str(.true.,block(b)%dims%Nk))
        stop
      endif
      mesh_dims%block_dims(b,l)%gc =      block(b)%dims%gc
      mesh_dims%block_dims(b,l)%Ni =      block(b)%dims%Ni/(2**(l-1))
      mesh_dims%block_dims(b,l)%Nj =      block(b)%dims%Nj/(2**(l-1))
      mesh_dims%block_dims(b,l)%Nk =      block(b)%dims%Nk/(2**(l-1))
      mesh_dims%block_dims(b,l)%Ns = size(block(b)%IC%r,dim=1)
      call mesh_dims%block_dims(b,l)%compute_NpNc
    enddo ; enddo
    mesh_dims%set = .true.
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_mesh_dims
  !> @}
endmodule Data_Type_File_Blocks_Cartesian
