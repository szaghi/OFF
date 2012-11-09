!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_AMRBlockPublicProcedure Data_Type_AMRBlock
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_AMRBlockPrivateProcedure Data_Type_AMRBlock
!> @}

!> @brief Module Data_Type_AMRBlock contains the definition of Type_AMRBlock type and useful procedures for its
!> handling.
module Data_Type_AMRBlock
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision     !< Integers and reals precision definition.
USE Data_Type_Global !< Definition of Type_Global.
USE Data_Type_Octant !< Definition of Type_Octant.
USE Data_Type_SBlock !< Definition of Type_SBlock.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the AMR block data.
!> AMRBlock type contains data (mesh, boundary conditions and fluid dynamic data) of an Adaptive Mesh Refinement block. It consists
!> of a level 0 structured block and a 3D array of cells defined as Type_Octant. An octant is a hierarchy of 8 cells structured
!> block with recursive refinement possibility of each cell as new finer octant.
!> @note Type_AMRBlock \b extends Type_SBlock. As a consequence each procedure working on Type_SBlock can work also on Type_AMRBlock
!> variables. The type-bound procedures of Type_SBlock are overrided: \b free and \b alloc in order to take into account also the
!> AMR octants data.
!> @ingroup DerivedType
type, public, extends(Type_SBlock):: Type_AMRBlock
  type(Type_Octant), allocatable:: octant(:,:,:) !< Cells octants [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
endtype Type_AMRBlock
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_AMRBlockPrivateProcedure
  !> @{
  !> Subroutine for freeing dynamic data of Type_AMRBlock variables.
  subroutine free_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block !< Block data.
  integer(I4P)::                      i     !< I direction counter.
  integer(I4P)::                      j     !< J direction counter.
  integer(I4P)::                      k     !< K direction counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! mesh data
  if (allocated(block%node)) deallocate(block%node)
  if (allocated(block%NFi )) deallocate(block%NFi )
  if (allocated(block%NFj )) deallocate(block%NFj )
  if (allocated(block%NFk )) deallocate(block%NFk )
  if (allocated(block%Si  )) deallocate(block%Si  )
  if (allocated(block%Sj  )) deallocate(block%Sj  )
  if (allocated(block%Sk  )) deallocate(block%Sk  )
  if (allocated(block%V   )) deallocate(block%V   )
  if (allocated(block%cent)) deallocate(block%cent)
  if (allocated(block%cell)) deallocate(block%cell)
  ! boundary conditions data
  if (allocated(block%BCi)) then
    call block%BCi%free ; deallocate(block%BCi)
  endif
  if (allocated(block%BCj)) then
    call block%BCj%free ; deallocate(block%BCj)
  endif
  if (allocated(block%BCk)) then
    call block%BCk%free ; deallocate(block%BCk)
  endif
  ! fluid dynamic data
  if (allocated(block%Dt)) deallocate(block%Dt)
  if (allocated(block%P)) then
    call block%P%free ; deallocate(block%P)
  endif
  if (allocated(block%U)) then
    call block%U%free ; deallocate(block%U)
  endif
  if (allocated(block%KS)) then
    call block%KS%free ; deallocate(block%KS)
  endif
  select type(block)
  type is(Type_SBlock)
    ! necessary for avoid "class default" branch when invoked with Type_Block variable
  class is(Type_AMRBlock)
    ! octants data
    if (allocated(block%octant)) then
      do k=1-block%gc(5),block%Nk+block%gc(6)
        do j=1-block%gc(3),block%Nj+block%gc(4)
          do i=1-block%gc(1),block%Ni+block%gc(2)
            call block%octant(i,j,k)%destroy
          enddo
        enddo
      enddo
    endif
  class default
    stop "Unexpected type of block: free memory failed!"
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_block

  !> Subroutine for allocating dynamic data of Type_AMRBlock variables.
  subroutine alloc_block(block,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block data.
  type(Type_Global),  intent(IN)::    global   !< Global data.
  integer(I_P)::                      gc(1:6)  !< Temporary variable  for storing block ghost cells number.
  integer(I_P)::                      Ni,Nj,Nk !< Temporary variables for storing block dimensions.
  integer(I_P)::                      Ns       !< Temporary variable  for storing number of species.
  integer(I_P)::                      rk_ord   !< Temporary variable  for storing rk_ord.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! free dynamic data of block if previously allocated
  call block%free
  ! storing block dimensions into temporary variables to simplify the code
  gc(1:6) = block%gc(1:6)
  Ni      = block%Ni
  Nj      = block%Nj
  Nk      = block%Nk
  Ns      = global%Ns
  rk_ord  = global%rk_ord
  ! mesh data
  allocate(block%node(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  allocate(block%NFi (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%NFj (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%NFk (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  allocate(block%Si  (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%Si = 0._R_P
  allocate(block%Sj  (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%Sj = 0._R_P
  allocate(block%Sk  (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; block%Sk = 0._R_P
  allocate(block%V   (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%V  = 0._R_P
  allocate(block%cell(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%cent(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  ! boundary conditions data
  allocate(block%BCi(0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%BCj(1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%BCk(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  ! fluid dynamic data
  allocate(block%Dt(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; block%Dt = 0._R_P
  allocate(block%P (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; call block%P%init(Ns=Ns)
  allocate(block%U (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; call block%U%init(Ns=Ns)
  allocate(block%KS(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6),1:rk_ord)) ; call block%KS%init(Ns=Ns)
  select type(block)
  type is (Type_SBlock)
    ! necessary for avoid "class default" branch when invoked with Type_Block variable
  class is (Type_AMRBlock)
    ! octants data
    allocate(block%octant(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  class default
    stop "Unexpected type of block: free memory failed!"
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_block
  !> @}
endmodule Data_Type_AMRBlock
