!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_OctantPublicProcedure Data_Type_Octant
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_OctantPrivateProcedure Data_Type_Octant
!> @}

!> @brief Module Data_Type_Octant contains the definition of Type_Octant type and useful procedures for its
!> handling.
module Data_Type_Octant
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision           !< Integers and reals precision definition.
!USE Data_Type_BC           !< Definition of Type_BC.
!USE Data_Type_Cell         !< Definition of Type_Cell.
!USE Data_Type_Conservative !< Definition of Type_Conservative.
USE Data_Type_Global       !< Definition of Type_Global.
!USE Data_Type_Primitive    !< Definition of Type_Primitive.
USE Data_Type_SBlock       !< Definition of Type_SBlock.
!USE Data_Type_Vector       !< Definition of Type_Vector.
!USE Lib_IO_Misc            !< Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the octant data.
!> Octant type contains data (mesh, boundary conditions and fluid dynamic data) of a refined cell. An octant is a hierarchy of 8
!> cells structured block with recursive refinement possibility of each cell as new finer octant. Each octant contains its octree
!> refinement level, informations of its parent cell, a 8 cells structured block data and 8 pointers to the eventually present
!> octants children.
!> @note All octant members are allocatable thus its memory requirement is null if the parent cell is not refined.
!> @note The 8 cells building up the octant are organized as a structured block in order to re-use all the procedures developed for
!> for a standard, non AMR, structured block.
!> @ingroup DerivedType
type, public:: Type_Octant
  integer(I1P),      allocatable:: level  !< Refinement level of the octant.
  integer(I4P),      allocatable:: ipr    !< I index of parent cell.
  integer(I4P),      allocatable:: jpr    !< J index of parent cell.
  integer(I4P),      allocatable:: kpr    !< K index of parent cell.
  type(Type_SBlock), allocatable:: sblock !< Structured block of 8 cells (plus ghost ones) building up the octant data.
  !> Pointers to children octants [1-gc(1):2+gc(2),1-gc(3):2+gc(4),1-gc(5):2+gc(6)].
  type(Type_Octant), pointer::     child(:,:,:) => null()
  contains
    procedure, non_overridable:: destroy => destroy_octant ! Procedure for destroying octant.
    procedure, non_overridable:: create => create_octant   ! Procedure for creating octant.
    procedure, non_overridable:: coarse => coarse_cell     ! Procedure for coarsening (destroy children) a cell of the octant.
    procedure, non_overridable:: refine => refine_cell     ! Procedure for refining (create a child) a cell of the octant.
endtype Type_Octant
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_OctantPrivateProcedure
  !> @{
  !> Subroutine for destroying data of Type_Octant variables.
  recursive subroutine destroy_octant(octant)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Octant), intent(INOUT):: octant  !< Octant data.
  integer(I_P)::                      gc(1:6) !< Temporary variable for storing block ghost cells number.
  integer(I4P)::                      i       !< I direction counter.
  integer(I4P)::                      j       !< J direction counter.
  integer(I4P)::                      k       !< K direction counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! refinement level
  if (allocated(octant%level)) deallocate(octant%level)
  ! parent cell informations
  if (allocated(octant%ipr)) deallocate(octant%ipr)
  if (allocated(octant%jpr)) deallocate(octant%jpr)
  if (allocated(octant%kpr)) deallocate(octant%kpr)
  ! structured block data
  if (allocated(octant%sblock)) then
    gc(1:6) = octant%sblock%gc(1:6)
    call octant%sblock%free ; deallocate(octant%sblock)
  endif
  ! children octants pointers
  if (associated(octant%child)) then
    do k=1-gc(5),2+gc(6)
      do j=1-gc(3),2+gc(4)
        do i=1-gc(1),2+gc(2)
          call octant%child(i,j,k)%destroy
        enddo
      enddo
    enddo
    deallocate(octant%child)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy_octant

  !> Subroutine for creating Type_Octant variables.
  subroutine create_octant(octant,level,ipr,jpr,kpr,pblock,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Octant), intent(INOUT):: octant   !< Octant data.
  integer(I1P),       intent(IN)::    level    !< Refinement level of the octant.
  integer(I4P),       intent(IN)::    ipr      !< I index of parent cell.
  integer(I4P),       intent(IN)::    jpr      !< J index of parent cell.
  integer(I4P),       intent(IN)::    kpr      !< K index of parent cell.
  type(Type_SBlock),  intent(IN)::    pblock   !< Parent block.
  type(Type_Global),  intent(IN)::    global   !< Global data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !! destroye octant if previously created
  !call octant%destroy
  ! refinement level
  if (.not.allocated(octant%level)) allocate(octant%level) ; octant%level = level
  ! parent cell informations
  if (.not.allocated(octant%ipr)) allocate(octant%ipr) ; octant%ipr = ipr
  if (.not.allocated(octant%jpr)) allocate(octant%jpr) ; octant%jpr = jpr
  if (.not.allocated(octant%kpr)) allocate(octant%kpr) ; octant%kpr = kpr
  ! structured block data
  if (.not.allocated(octant%sblock)) then
    allocate(octant%sblock)
    octant%sblock%gc = pblock%gc
    octant%sblock%Ni = 2
    octant%sblock%Nj = 2
    octant%sblock%Nk = 2
    call octant%sblock%alloc(global=global)
    !call init_octant_mesh()
    !call init_octant_bc()
    !call init_octant_fluid()
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create_octant

  !> Subroutine for coarsening (destroy children) a cell of the octant.
  recursive subroutine coarse_cell(octant,i,j,k)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Octant), intent(INOUT):: octant  !< Octant data.
  integer(I4P),       intent(IN)::    i       !< I index of the cell of the octant to be refined.
  integer(I4P),       intent(IN)::    j       !< J index of the cell of the octant to be refined.
  integer(I4P),       intent(IN)::    k       !< K index of the cell of the octant to be refined.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (associated(octant%child)) then
    call octant%child(i,j,k)%destroy
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine coarse_cell

  !> Subroutine for refining (create a child) a cell of the octant.
  recursive subroutine refine_cell(octant,i,j,k,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Octant), intent(INOUT):: octant  !< Octant data.
  integer(I4P),       intent(IN)::    i       !< I index of the cell of the octant to be refined.
  integer(I4P),       intent(IN)::    j       !< J index of the cell of the octant to be refined.
  integer(I4P),       intent(IN)::    k       !< K index of the cell of the octant to be refined.
  type(Type_Global),  intent(IN)::    global  !< Global data.
  integer(I_P)::                      gc(1:6) !< Temporary variable for storing block ghost cells number.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.associated(octant%child)) then
    gc(1:6) = octant%sblock%gc(1:6)
    allocate(octant%child(1-gc(1):2+gc(2),1-gc(3):2+gc(4),1-gc(5):2+gc(6)))
  endif
  call octant%child(i,j,k)%create(level=octant%level+1_I1P,ipr=i,jpr=j,kpr=k,pblock=octant%sblock,global=global)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine refine_cell
  !> @}
endmodule Data_Type_Octant
