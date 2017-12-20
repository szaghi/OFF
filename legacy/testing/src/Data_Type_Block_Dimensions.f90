!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Block_DimensionsDerivedType Data_Type_Block_Dimensions
!> Module definition of Type_Block_Dimensions
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_Block_DimensionsInterface Data_Type_Block_Dimensions
!> Module definition of Type_Block_Dimensions
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Block_DimensionsPrivateProcedure Data_Type_Block_Dimensions
!> Module definition of Type_Block_Dimensions
!> @}

!> @brief Module Data_Type_Block_Dimensions contains the definition of Type_Block_Dimensions and useful procedures for its handling.
!> Type_Block_Dimensions contains all the data defining the main dimensions of block.
module Data_Type_Block_Dimensions
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                          ! Integers and reals precision definition.
USE Data_Type_XML_Tag, only: Type_XML_Tag ! Definition of Type_XML_Tag.
USE Lib_Strings,       only: unique       ! Library for strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the main dimensions of a block.
!> @ingroup Data_Type_Block_DimensionsDerivedType
type, public:: Type_Block_Dimensions
  integer(I1P):: gc(1:6)=&
                        [1_I1P,1_I1P, & ! gc(1) => left i, gc(2) => right i.
                         1_I1P,1_I1P, & ! gc(3) => left j, gc(4) => right j.
                         1_I1P,1_I1P  & ! gc(5) => left k, gc(6) => right k.
                         ]              !< Number of ghost cells for the 6 faces of the block.
  integer(I4P):: Ni     = 0_I4P         !< Number of cells in i direction.
  integer(I4P):: Nj     = 0_I4P         !< Number of cells in j direction.
  integer(I4P):: Nk     = 0_I4P         !< Number of cells in k direction.
  integer(I4P):: Np     = 7_I4P         !< Number of primitive variables    (Np = Ns + 6).
  integer(I4P):: Nc     = 5_I4P         !< Number of conservative variables (Nc = Ns + 4).
  integer(I4P):: Ns     = 1_I4P         !< Number of species.
  integer(I1P):: Nrk    = 1_I4P         !< Number of Runge-Kutta stages.
  contains
    procedure:: Nnodes => compute_Nnodes_block ! Procedure for computing the number of nodes.
    procedure:: Ncells => compute_Ncells_block ! Procedure for computing the number of cells.
    procedure:: NFi    => compute_NFi_block    ! Procedure for computing the number i-faces.
    procedure:: NFj    => compute_NFj_block    ! Procedure for computing the number j-faces.
    procedure:: NFk    => compute_NFk_block    ! Procedure for computing the number k-faces.
    procedure:: load   => load_block_dims      ! Procedure for loading dimensions from mesh file.
    procedure:: save   => save_block_dims      ! Procedure for saving dimensions into mesh file.
    procedure:: compute_NpNc                   ! Procedure for computing Np and Nc according to Ns.
    procedure:: load_str_xml                   ! Procedure for loading dimensions from a string in XML format.
    procedure:: save_str_xml                   ! Procedure for saving dimensions into a string in XML format.
    procedure:: print  => print_block_dims     ! Procedure for printing dimensions with a pretty format.
    ! operators overloading
    generic:: assignment(=) => assign_blk_dims
    ! private procedures
    procedure, pass(blk1), private:: assign_blk_dims
  endtype Type_Block_Dimensions
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Block_DimensionsPrivateProcedure
  !> @{
  !> @brief Procedure for computing the number of nodes of a block.
  elemental function compute_Nnodes_block(block_dims) result(Nnodes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(IN):: block_dims !< Block dimensions.
  integer(I4P)::                             Nnodes     !< Number of nodes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nnodes = ((block_dims%Ni + block_dims%gc(2)) - (0_I4P - block_dims%gc(1)) + 1_I4P)*& ! i direction.
           ((block_dims%Nj + block_dims%gc(4)) - (0_I4P - block_dims%gc(3)) + 1_I4P)*& ! j direction.
           ((block_dims%Nk + block_dims%gc(6)) - (0_I4P - block_dims%gc(5)) + 1_I4P)   ! k direction.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_Nnodes_block

  !> @brief Procedure for computing the number of cells of a block.
  elemental function compute_Ncells_block(block_dims) result(Ncells)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(IN):: block_dims !< Block dimensions.
  integer(I4P)::                             Ncells     !< Number of cells.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ncells = ((block_dims%Ni + block_dims%gc(2)) - (1_I4P - block_dims%gc(1)) + 1_I4P)*& ! i direction.
           ((block_dims%Nj + block_dims%gc(4)) - (1_I4P - block_dims%gc(3)) + 1_I4P)*& ! j direction.
           ((block_dims%Nk + block_dims%gc(6)) - (1_I4P - block_dims%gc(5)) + 1_I4P)   ! k direction.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_Ncells_block

  !> @brief Procedure for computing the number of i-faces of a block.
  elemental function compute_NFi_block(block_dims) result(NFi)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(IN):: block_dims !< Block dimensions.
  integer(I4P)::                             NFi        !< Number of i-faces.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  NFi = ((block_dims%Ni + block_dims%gc(2)) - (0_I4P - block_dims%gc(1)) + 1_I4P)*& ! i direction.
        ((block_dims%Nj + block_dims%gc(4)) - (1_I4P - block_dims%gc(3)) + 1_I4P)*& ! j direction.
        ((block_dims%Nk + block_dims%gc(6)) - (1_I4P - block_dims%gc(5)) + 1_I4P)   ! k direction.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_NFi_block

  !> @brief Procedure for computing the number of j-faces of a block.
  elemental function compute_NFj_block(block_dims) result(NFj)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(IN):: block_dims !< Block dimensions.
  integer(I4P)::                             NFj        !< Number of j-faces.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  NFj = ((block_dims%Ni + block_dims%gc(2)) - (1_I4P - block_dims%gc(1)) + 1_I4P)*& ! i direction.
        ((block_dims%Nj + block_dims%gc(4)) - (0_I4P - block_dims%gc(3)) + 1_I4P)*& ! j direction.
        ((block_dims%Nk + block_dims%gc(6)) - (1_I4P - block_dims%gc(5)) + 1_I4P)   ! k direction.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_NFj_block

  !> @brief Procedure for computing the number of k-faces of a block.
  elemental function compute_NFk_block(block_dims) result(NFk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(IN):: block_dims !< Block dimensions.
  integer(I4P)::                             NFk        !< Number of j-faces.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  NFk = ((block_dims%Ni + block_dims%gc(2)) - (1_I4P - block_dims%gc(1)) + 1_I4P)*& ! i direction.
        ((block_dims%Nj + block_dims%gc(4)) - (1_I4P - block_dims%gc(3)) + 1_I4P)*& ! j direction.
        ((block_dims%Nk + block_dims%gc(6)) - (0_I4P - block_dims%gc(5)) + 1_I4P)   ! k direction.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_NFk_block

  !> @brief Procedure for loading block dimensions from file.
  subroutine load_block_dims(block_dims,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(INOUT):: block_dims !< Block dimensions.
  integer(I8P), optional,       intent(IN)::    pos        !< Position specifier.
  integer(I4P), optional,       intent(OUT)::   iostat     !< IO error.
  character(*), optional,       intent(OUT)::   iomsg      !< IO error message.
  integer(I4P),                 intent(IN)::    unit       !< Logic unit.
  integer(I4P)::                                iostatd    !< IO error.
  character(500)::                              iomsgd     !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    read(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)block_dims%gc,block_dims%Ni,block_dims%Nj,block_dims%Nk,&
                                                       block_dims%Np,block_dims%Nc,block_dims%Ns,block_dims%Nrk
  else
    read(unit=unit,        iostat=iostatd,iomsg=iomsgd)block_dims%gc,block_dims%Ni,block_dims%Nj,block_dims%Nk,&
                                                       block_dims%Np,block_dims%Nc,block_dims%Ns,block_dims%Nrk
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_block_dims

  !> @brief Procedure for saving block dimensions into file.
  subroutine save_block_dims(block_dims,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(IN)::  block_dims !< Block dimensions.
  integer(I8P), optional,       intent(IN)::  pos        !< Position specifier.
  integer(I4P), optional,       intent(OUT):: iostat     !< IO error.
  character(*), optional,       intent(OUT):: iomsg      !< IO error message.
  integer(I4P),                 intent(IN)::  unit       !< Logic unit.
  integer(I4P)::                              iostatd    !< IO error.
  character(500)::                            iomsgd     !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    write(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)block_dims%gc,block_dims%Ni,block_dims%Nj,block_dims%Nk,&
                                                        block_dims%Np,block_dims%Nc,block_dims%Ns,block_dims%Nrk
  else
    write(unit=unit,        iostat=iostatd,iomsg=iomsgd)block_dims%gc,block_dims%Ni,block_dims%Nj,block_dims%Nk,&
                                                        block_dims%Np,block_dims%Nc,block_dims%Ns,block_dims%Nrk
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block_dims

  !> @brief Procedure for computing Np and Nc according to Ns.
  elemental subroutine compute_NpNc(block_dims)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(INOUT):: block_dims !< Block dimensions.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  block_dims%Np  = block_dims%Ns + 6
  block_dims%Nc  = block_dims%Ns + 4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_NpNc

  !> @brief Procedure for loading block dimensions from a string formatted in XML syntax.
  !> @note The XML string of block dimensions data must have the following syntax:
  !> @code
  !>   <dimensions gc="#imin #imax #jmin #jmax #kmin #kmax" Ni="#cell-i" Nj="#cell-j" Nk="#cell-k"/>
  !> @endcode
  !> The attributes gc defines the ghost cells number along i, j and k directions, the attributes Ni, Nj and Nk define the number of
  !> (inner) cells along i, j  and k directions, respectively.
  !> The order of attributes in not influent. It is worth nothing that the syntax of tag names and attributes is case
  !> sensitive, whereas the syntax of their values is case insensitive.
  subroutine load_str_xml(block_dims,str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(INOUT):: block_dims !< Block dimensions.
  character(*),                 intent(IN)::    str_xml    !< String containing block dimensions data in XML syntax.
  type(Type_XML_Tag)                            tag        !< XML tag data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call tag%set(string=str_xml,tag_name='dimensions',att_name=['gc','Ni','Nj','Nk']) ! setting xml tag data
  call tag%get_attributes                                                           ! parsing xml tag attributes
  tag%att_val(1)%vs = trim(unique(string=tag%att_val(1)%vs,substring=' '))          ! removing the eventual multiples white spaces
  ! loading dimensions data from xml tag attributes values
  read(tag%att_val(1)%vs,*)block_dims%gc
  block_dims%Ni = cton(str=tag%att_val(2)%vs,knd=1_I4P)
  block_dims%Nj = cton(str=tag%att_val(3)%vs,knd=1_I4P)
  block_dims%Nk = cton(str=tag%att_val(4)%vs,knd=1_I4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_str_xml

  !> @brief Procedure for saving block dimensions into a string formatted in XML syntax.
  !> @note The definition of the XML string syntax is described in the 'load_str_xml' documentation.
  pure function save_str_xml(block_dims) result(str_xml)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(IN):: block_dims !< Block dimensions.
  character(len=:), allocatable::            str_xml    !< String containing block dimensions data in XML syntax.
  integer(I4P)::                             f          !< Block faces counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  str_xml = ' '
  do f=1,6
    str_xml = trim(str_xml)//' '//trim(adjustl(str(n=block_dims%gc(f))))
  enddo
  str_xml = '<dimensions gc="'//trim(adjustl(str_xml))//                                                     &
            '" Ni="'//trim(str(n=block_dims%Ni))//'" Nj="'//trim(str(n=block_dims%Nj))//'" Nk="'//trim(str(n=block_dims%Nk))//'"/>'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_str_xml

  !> @brief Procedure for printing in a pretty ascii format the components of type Type_Block_Dimensions.
  subroutine print_block_dims(block_dims,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(IN)::  block_dims !< Block dimensions.
  character(*), optional,       intent(IN)::  pref       !< Prefixing string.
  integer(I4P), optional,       intent(OUT):: iostat     !< IO error.
  character(*), optional,       intent(OUT):: iomsg      !< IO error message.
  integer(I4P),                 intent(IN)::  unit       !< Logic unit.
  character(len=:), allocatable::             prefd      !< Prefixing string for outputs.
  integer(I4P)::                              iostatd    !< IO error.
  character(500)::                            iomsgd     !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  associate(gc=>block_dims%gc,Ni=>block_dims%Ni,Nj=>block_dims%Nj,Nk=>block_dims%Nk,Np=>block_dims%Np,Nc=>block_dims%Nc,&
            Ns=>block_dims%Ns,Nrk=>block_dims%Nrk)
    write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' gc(1)='//trim(str(n=gc(1)))//' gc(2)='//trim(str(n=gc(2)))
    write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' gc(3)='//trim(str(n=gc(3)))//' gc(4)='//trim(str(n=gc(4)))
    write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' gc(5)='//trim(str(n=gc(5)))//' gc(6)='//trim(str(n=gc(6)))
    write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Ni='//trim(str(n=Ni))//' Nj='//trim(str(n=Nj))//' Nk='//trim(str(n=Nk))
    write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Np='//trim(str(n=Np))//' Nc='//trim(str(n=Nc))//' Ns='//trim(str(n=Ns))
    write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Nrk='//trim(str(n=Nrk))
  endassociate
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_block_dims

  ! Assignment (=)
  !> @brief Procedure for assignment between two blocks dimensions variables.
  elemental subroutine assign_blk_dims(blk1,blk2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block_Dimensions), intent(INOUT):: blk1
  type(Type_Block_Dimensions),  intent(IN)::    blk2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  blk1%gc  = blk2%gc
  blk1%Ni  = blk2%Ni
  blk1%Nj  = blk2%Nj
  blk1%Nk  = blk2%Nk
  blk1%Np  = blk2%Np
  blk1%Nc  = blk2%Nc
  blk1%Ns  = blk2%Ns
  blk1%Nrk = blk2%Nrk
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_blk_dims
  !> @}
endmodule Data_Type_Block_Dimensions
