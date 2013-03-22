!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_AMRBlockDerivedType Data_Type_AMRBlock
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_AMRBlockPublicProcedure Data_Type_AMRBlock
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_AMRBlockPrivateProcedure Data_Type_AMRBlock
!> @}

!> @brief Module Data_Type_AMRBlock contains the definition of Type_AMRBlock type and useful procedures for its handling.
module Data_Type_AMRBlock
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision        ! Integers and reals precision definition.
USE Data_Type_Cell      ! Definition of Type_Cell.
USE Data_Type_Global    ! Definition of Type_Global.
USE Data_Type_HashID    ! Definition of Type_HashID.
USE Data_Type_HashTCell ! Definition of Type_HashTCell.
USE Data_Type_HashTFace ! Definition of Type_HashTFace.
USE Data_Type_HashTNode ! Definition of Type_HashTNode.
USE Data_Type_Vector    ! Definition of Type_Vector.
USE Lib_IO_Misc         ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: level
public:: firstID
public:: lastID
public:: siblings
public:: child
public:: directparent
public:: parentatlevel
public:: directchildren
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing AMR (hierarchical, dynamic) grid data.
!> AMR grid type contains data (mesh, boundary conditions and fluid dynamic data) of Adaptive Mesh Refinement hierarchical (implicit
!> connectivity) numerical grid. It is specular to a Structured block: the 'base' cells are structured cells with the ability to
!> self refine adaptively. The hierarchical, dynamic data structured (dynamic octree) is based on hash tables. The hash tables are
!> linearized in order to facilitate parallelization. 'Space filling curve' techniques are used to this purpose.
!> @note The space filling technique defines a unique key (ID) for each cell at each refinement level.
!> Considering a parent cell (named 0), each child cell that can be attained within the refinement hierarchy from the parent cell,
!> can be identified by a unique 'path', as shown in the following code for 1 refinement level.
!> The ordering of the children of each cell (that constitute an 'octant') is as following:
!> @code
!>                  +----+----+
!>                 /|   /|   /|
!>                / |  6 |  7 |
!>               /  +----*----+
!>              /  /|/  /|/  /|
!>             /  / |  4 |  5 |
!>            /  / /+----+----+
!>           /  / // / // /  /
!>          +----+----+/ /  /
!>         /| / /| / /| /  /
!>        / |/ //|/ //|/  /
!>       /  +----*----+  /
!>      /  /|// /|// /| /
!>     /  / |/ / |/ / |/
!>    /  / /+----+----+
!>   /  / // / // /  /
!>  +----+----+/ /  /       y(j)  z(k)
!>  | /  | /  | /  /        ^    ^
!>  |/ 2/|/ 3/|/  /         |   /
!>  +----*----+  /          |  /
!>  | /  | /  | /           | /
!>  |/ 0 |/ 1 |/            |/
!>  +----+----+             +------->x(i)
!> @endcode
!> The path is read from right to left, from lower to higher levels. The child numbers of the path are called 'digits'.
!> Given a cell of level 'l', only 'l digits' associated with the lowest levels are of interest. As a consequence the cell
!> identifier is constructed as following:
!> @code
!> ID = (block,bcell,level,path,reserved_bits)
!> 'block' are the bits identifying block;
!> 'bcell' are the bits identifying the (parent) base cell;
!> 'level' are the bits identifying the level;
!> 'path' are the bits describing the path;
!> 'reserved_bits' are the bits reserved for other purposes.
!> @endcode
!> Because the digits of the path range in [0,1,2...,7] each digit needs 4 bits. The maximum level of refinement is arbitrarily
!> chosen to be 12, thus the level information needs 4 bits. The path needs digit_bits*level_max = 4 bits*12 = 48 bits. The base
!> cells are aritrarly set to a maximum of 16777216 needing 24 bits to be reprensented (this choice allows, for example, block of
!> base cells of the dimensions of 256x256x256). The maximum number of blocks is arbitrarily chosen to be 65536 needing 16 bits.
!> The non reserved bits needed for building the cell identifier are 16+24+4+48 = 92.
!> @ingroup Data_Type_AMRBlockDerivedType
type, public:: Type_AMRBlock
  ! Structured, base grid informations.
  type(Type_Global), pointer:: global                  !< Pointer to global-level data.
  integer(I1P)::               gc(1:6)=&
                                      (/1_I1P,1_I1P, & ! gc(1) => left i, gc(2) => right i.
                                        1_I1P,1_I1P, & ! gc(3) => left j, gc(4) => right j.
                                        1_I1P,1_I1P  & ! gc(5) => left k, gc(6) => right k.
                                        /)             !< Number of ghost cells for the 6 faces of the block.
  integer(I2P)::               Ni = 0_I2P              !< Number of (base) cells in i direction.
  integer(I2P)::               Nj = 0_I2P              !< Number of (base) cells in j direction.
  integer(I2P)::               Nk = 0_I2P              !< Number of (base) cells in k direction.
  ! Hierarchical (dynamic) data structure based on linearized octree and hash tables.
  type(Type_HashTNode):: node !< (Hash table of) Nodes coordinates.
  type(Type_HashTFace):: Fi   !< (Hash table of) Faces i data.
  type(Type_HashTFace):: Fj   !< (Hash table of) Faces j data.
  type(Type_HashTFace):: Fk   !< (Hash table of) Faces k data.
  type(Type_HashTCell):: C    !< (Hash table of) Cells data.
  contains
    procedure, non_overridable:: init => init_amrblock ! Procedure for initializing the AMR grid.
    procedure, non_overridable:: load_basemesh_dims  ! Procedure for loading the base mesh data dimensions.
endtype Type_AMRBlock
  !integer(I4P)::              Nc = 0_I4P ! Number of cells on the current process
  !integer(I8P), allocatable:: IDfirst(:) ! ID of the first cell on each process [1:Nproc].
  !integer(I8P), allocatable:: IDlast(:)  ! ID of the last  cell on each process [1:Nproc].
  !integer(I8P), allocatable:: ID(:)      ! ID of each cell on current process [1:Nc].
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_AMRPrivateProcedure
  !> @{
  !> @brief Subroutine for initializing the AMR grid.
  subroutine init_amrblock(block,b)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_AMRBlock), intent(INOUT):: block !< AMR block.
  integer(I2P),         intent(IN)::    b     !< Block number (in the global ordering).
  integer(I2P)::                        i,j,k !< Counters.
  type(Type_HashID)::                   ID    !< ID value.
  type(Type_Cell)::                     cell  !< Cell data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing hash tables
  call block%node%init
  call block%Fi%init
  call block%Fj%init
  call block%Fk%init
  call block%C%init
  ! inserting data into the hash tables
  do k=1-block%gc(5),block%Nk+block%gc(6)
    do j=1-block%gc(3),block%Nj+block%gc(4)
      do i=1-block%gc(1),block%Ni+block%gc(2)
        call ID%build(b=b,i=i,j=j,k=k,l=0_I1P,p=0_I8P)
        call block%C%put(ID,cell)
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_amrblock

  !> Function for loading the base mesh data dimensions of block from the mesh file 'filename'.
  !> @return \b err integer(I_P) variable.
  function load_basemesh_dims(block,ascii,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_AMRBlock), intent(INOUT):: block    !< Block level data.
  logical, optional,    intent(IN)::    ascii    !< Flag for ascii file.
  character(*),         intent(IN)::    filename !< Name of file where mesh variables are saved.
  integer(I_P)::                        UnitFree !< Free logic unit.
  logical::                             is_file  !< Flag for inquiring the presence of mesh file.
  integer(I_P)::                        err      !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(block%global%myrank,filename,'load_basemesh_dims')
  if (present(ascii)) then
    open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err)block%gc
    read(UnitFree,*,iostat=err)block%Ni,block%Nj,block%Nk
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree ), file = trim(filename), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)block%gc
    read(UnitFree,iostat=err)block%Ni,block%Nj,block%Nk
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_basemesh_dims
  !> @}

  !> @ingroup Data_Type_AMRBlockPublicProcedure
  !> @{
  !> Function for computing the level correspoinding to the passed ID.
  elemental function level(ID) result(l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN) :: ID !< Input ID.
  integer(I4P)::              l  !< Level corresponding to the passed ID.
  integer(I8P)::              i  !< ID counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  l = 0
  i = ID
  do while (i/=0_I8P)
    i = (i-1_I8P)/8_I8P
    l = l + 1
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction level

  !> Function for computing the first ID of the passed level.
  elemental function firstID(level) result(fid)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN) :: level    !< Passed level.
  integer(I8P)::              fid      !< First ID of the passed level.
  integer(I4P)::              l        !< Level counter:
  integer(I8P)::              idoffset !< IDs offset.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  fid = 0_I8P
  idoffset = 1_I8P
  do l=0,level-1
    fid = fid + idoffset
    idoffset = idoffset*8
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction firstID

  !> Function for computing the last ID of the passed level.
  elemental function lastID(level) result(lid)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN) :: level    !< Passed level.
  integer(I8P)::              lid      !< First ID of the passed level.
  integer(I4P)::              l        !< Level counter:
  integer(I8P)::              idoffset !< IDs offset.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  lid = -1_I8P
  idoffset = 1_I8P
  do l=0,level
    lid = lid + idoffset
    idoffset = idoffset*8
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction lastID

  !> Function for computing the IDs of siblings of a given ID.
  pure function siblings(ID) result(sbs)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: ID       !< Input ID.
  integer(I8P)::             sbs(1:7) !< IDs of siblings.
  integer(I4P)::             myid     !< ID of input ID into siblings list.
  integer(I4P)::             s        !< Siblings counter.
  integer(I4P)::             i        !< ID counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  myid = int(ID - ((ID - 1_I8P)/8_I8P)*8_I8P)
  s = 0
  do i=1,8
    if (i/=myid) then
      s = s + 1
      sbs(s) = ((ID - 1_I8P)/8_I8P)*8_I8P + int(i,I8P)
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction siblings

  !> Function for computing the ID of child of a given (parent) ID.
  elemental function child(pID) result(cID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: pID !< Input (parent) ID.
  integer(I4P)::             cID !< Child ID.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cID = int(pID - ((pID - 1_I8P)/8_I8P)*8_I8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction child

  !> Function for computing the direct parent ID of a given ID.
  elemental function directparent(ID) result(p)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: ID !< Input ID.
  integer(I8P)::             p  !< Direct parent ID.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  p = (ID - 1_I8P)/8_I8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  end function directparent

  !> Function for computing the parent ID of a given ID at a given level.
  elemental function parentatlevel(ID,l) result(pID)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: ID  !< Input ID.
  integer(I4P), intent(IN):: l   !< Input level.
  integer(I8P)::             pID !< Parent ID of input ID at level 'l'.
  integer(I4P)::             i   !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pID = ID
  do i=1,level(ID)-l
    pID = (pID-1)/8
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction parentatlevel

  !> Function for computing the direct children IDs of a given ID.
  pure function directchildren(ID) result(children)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: ID            !< Input ID.
  integer(I8P)::             children(1:8) !< Children IDs.
  integer(I8P)::             offset        !< Children offset.
  integer(I8P)::             c             !< Children counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  children(:) = 0_I8P
  offset = ID*8_I8P
  do c=1,8
    children(c) = offset + c
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  end function directchildren
  !> @}
endmodule Data_Type_AMRBlock
