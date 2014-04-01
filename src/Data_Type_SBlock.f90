!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_SBlockDerivedType Data_Type_SBlock
!> Module definition of Type_SBlock
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_SBlockInterface Data_Type_SBlock
!> Module definition of Type_SBlock
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_SBlockPrivateProcedure Data_Type_SBlock
!> Module definition of Type_SBlock
!> @}

!> @brief Module Data_Type_SBlock contains the definition of Type_SBlock (structured block) type and useful procedures for its
!> handling.
module Data_Type_SBlock
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                               ! Integers and reals precision definition.
USE Data_Type_BC,                 only: Type_BC                ! Definition of Type_BC.
USE Data_Type_Block_BC,           only: Type_Block_BC          ! Definition of Type_Block_BC.
USE Data_Type_Block_Dimensions,   only: Type_Block_Dimensions  ! Definition of Type_Block_Dimensions.
USE Data_Type_Block_Extents,      only: Type_Block_Extents     ! Definition of Type_Block_Extents.
USE Data_Type_Cell,               only: Type_Cell              ! Definition of Type_Cell.
USE Data_Type_Conservative,       only: Type_Conservative      ! Definition of Type_Conservative.
USE Data_Type_Face,               only: Type_Face              ! Definition of Type_Face.
USE Data_Type_Mesh_Dimensions,    only: Type_Mesh_Dimensions   ! Definition of Type_Mesh_Dimensions.
USE Data_Type_Primitive,          only: Type_Primitive         ! Definition of Type_Primitive.
USE Data_Type_Region,             only: Type_Region            ! Definition of Type_Region.
USE Data_Type_Species,            only: Type_Species           ! Definition of Type_Species.
USE Data_Type_Space_Step,         only: Type_Space_Step        ! Definition of Type_Space_Step.
USE Data_Type_Time_Step,          only: Type_Time_Step         ! Definition of Type_Time_Step.
USE Data_Type_Vector,             only: Type_Vector,ex,ey,ez   ! Definition of Type_Vector.
USE Lib_Fluxes_Convective,        only: fluxes_convective      ! Procedure for convective fluxes.
USE Lib_IO_Misc,                  only: Upper_Case             ! Procedures for IO and strings operations.
USE Lib_Runge_Kutta,              only: rk_stage,rk_time_integ ! Runge-Kutta time integration library.
USE Lib_Thermodynamic_Laws_Ideal, only: a                      ! Procedure for computing speed of sound.
USE Lib_Variables_Conversions,    only: cons2prim              ! Pocedures for varibles set conversions.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing structured block-level data.
!> Structured block-level type contains data (mesh, boundary conditions and fluid dynamic data) of structured (implicit
!> connectivity) numerical grid. A structured block is an hexahedron with quadrilateral faces ussing the following internal
!> numeration for vertices and faces:
!> @code
!>
!> /|\Y
!>  |                            F(4)         _ F(6)
!>  |                            /|\          /!
!>  |                        7    |          /    8
!>  |                         *------------------*
!>  |                        /|   |        /    /|
!>  |                       / |   |       /    / |
!>  |                      /  |   |      /    /  |
!>  |                     /   |   |     /    /   |
!>  |                    /    |   |    +    /    |
!>  |                   /     |   |        /     |
!>  |                  /      |   +       /      |
!>  |                 /      5|          /       |6
!>  |                /        * --------/--------*
!>  |      F(1)<----/----+   /         /        /
!>  |              *------------------*    +-------->F(2)
!>  |             3|       /          |4      /
!>  |              |      /           |      /
!>  |              |     /        +   |     /
!>  |              |    /         |   |    /
!>  |              |   /      +   |   |   /
!>  |              |  /      /    |   |  /
!>  |              | /      /     |   | /
!>  |              |/      /      |   |/
!>  |              *------------------*
!>  |             1      /        |    2
!>  |                   /        \|/
!>  |   _ Z           |/_       F(3)
!>  |   /|         F(5)
!>  |  /
!>  | /
!>  |/                                                    X
!>  o----------------------------------------------------->
!>
!> @ingroup Data_Type_SBlockDerivedType
type, public:: Type_SBlock
  ! Block global data
  type(Type_Block_Extents)::       exts        !< Block extents, i.e. bounding box diagonal.
  type(Type_Block_BC)::            BC          !< Block (faces) boundary conditions.
  type(Type_Primitive)::           IC          !< Block initial conditions.
  type(Type_Block_Dimensions)::    dims        !< Block dimensions (gc,Ni,Nj,Nk,...).
  ! Block cells data
  type(Type_Vector), allocatable:: node(:,:,:) !< Nodes coord.  [0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  type(Type_Face),   allocatable:: Fi(:,:,:)   !< Faces i data  [0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Face),   allocatable:: Fj(:,:,:)   !< Faces j data  [1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Face),   allocatable:: Fk(:,:,:)   !< Faces k data  [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  type(Type_Cell),   allocatable:: C(:,:,:)    !< Cells data    [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  contains
    procedure:: free               => free_block               ! Procedure for freeing memory.
    procedure:: alloc              => alloc_block              ! Procedure for allocating memory.
    procedure:: interpolate_primitive                          ! Procedure for interpolating primitive variables at nodes.
    procedure:: metrics            => metrics_block            ! Procedure for computing block metrics.
    procedure:: metrics_correction => metrics_correction_block ! Procedure for correcting block metrics of bc cells.
    procedure:: node2center        => node2center_block        ! Procedure for computing cell center coo from cell nodes.
    procedure:: min_space_step     => min_space_step_block     ! Procedure for computing the minimum value of space step.
    procedure:: create_uniform_grid                            ! Procedure for creating a uniform grid provided block extents.
    procedure:: create_grid_from_finer                         ! Procedure for creating the grid from the grid of a finer block.
    procedure:: set_cells_bc                                   ! Procedure for setting cells boundary conditions from block ones.
    procedure:: set_region_ic                                  ! Procedure for setting initial condition in a region of block.
    procedure:: primitive2conservative                         ! Procedure for converting primitive to conservative variables.
    procedure:: conservative2primitive                         ! Procedure for converting conservative to primitive variables.
    procedure:: update_primitive                               ! Procedure for updating primitive variables of all cells of a block.
    procedure:: compute_time                                   ! Procedure for evaluating the local and global time step value.
    procedure:: residuals                                      ! Procedure for computing the residuals.
    procedure:: rk_stages_sum                                  ! Procedure for summing Runge-Kutta stages.
    procedure:: rk_time_integration                            ! Procedure for computing Runge-Kutta one time step integration.
    procedure:: load               => load_block               ! Procedure for loading block.
    procedure:: save               => save_block               ! Procedure for saving block.
    procedure:: print              => print_block              ! Procedure for printing block infos with pretty format.
    procedure:: print_info_mesh                                ! Procedure for printing block mesh infos with pretty format.
    procedure:: print_info_bc                                  ! Procedure for printing block bc infos with pretty format.
    procedure:: print_info_fluid                               ! Procedure for printing block fluid infos with pretty format.
    procedure:: mirror                                         ! Procedure for generating a mirrored block.
    final::     finalize                                       ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
    generic:: assignment(=) => assign_blk
    ! private procedures
    procedure, pass(blk1), private:: assign_blk
endtype Type_SBlock
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_SBlockPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  !> @note If this procedure is called without dummy arguments all dynamic memory variables are deallocated.
  elemental subroutine free_block(block,global_data,cells_data)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block       !< Block data.
  logical, optional,  intent(IN)::    global_data !< Switch for freeing only global block data.
  logical, optional,  intent(IN)::    cells_data  !< Switch for freeing only cells data array.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(global_data))  then
    if (global_data) call free_global_data(block)
  endif
  if (present(cells_data))  then
    if (cells_data) call free_cells_data(block)
  endif
  if ((.not.present(global_data)).and.(.not.present(cells_data)))  then
    call free_global_data(block)
    call free_cells_data(block)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Procedure for freeing dynamic memory of global data of block.
    elemental subroutine free_global_data(block)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    class(Type_SBlock), intent(INOUT):: block !< Block data.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    call block%BC%free
    call block%IC%free
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine free_global_data

    !> @brief Procedure for freeing dynamic memory of cells data of block.
    elemental subroutine free_cells_data(block)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    class(Type_SBlock), intent(INOUT):: block !< Block data.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (allocated(block%node)) deallocate(block%node)
    if (allocated(block%Fi)) then
      call block%Fi%free ; deallocate(block%Fi)
    endif
    if (allocated(block%Fj)) then
      call block%Fj%free ; deallocate(block%Fj)
    endif
    if (allocated(block%Fk)) then
      call block%Fk%free ; deallocate(block%Fk)
    endif
    if (allocated(block%C)) then
      call block%C%free ; deallocate(block%C)
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine free_cells_data
  endsubroutine free_block

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_SBlock), intent(INOUT):: block !< Block data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call block%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for allocating dynamic memory.
  elemental subroutine alloc_block(block,members)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block   !< Block data.
  logical, optional,  intent(IN)::    members !< Switch for allocating members memory (actually only C data).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(members)) then
    if (members) call block%C%alloc( Ns=block%dims%Ns,Nrk=block%dims%Nrk)
    !if (members) call block%IC%alloc(Ns=block%dims%Ns)
  else
    call block%free(cells_data=.true.)
    associate(gc => block%dims%gc(1:6),Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
      allocate(block%node(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      allocate(block%Fi  (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
      allocate(block%Fj  (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
      allocate(block%Fk  (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      allocate(block%C   (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    endassociate
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_block

  !> @brief Procedure for interpolating primitive variables at nodes.
  pure subroutine interpolate_primitive(block,primN)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),                intent(IN)::    block        !< Block data.
  type(Type_Primitive), allocatable, intent(INOUT):: primN(:,:,:) !< Nodes-interpolated primitive variables.
  type(Type_Primitive), allocatable::                primC(:,:,:) !< Cell (original) primitive variables.
  real(R8P)::                                        mf           !< Mean factor.
  integer(I4P)::                                     i,j,k        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Ni=>block%dims%Ni,Nj=>block%dims%Nj,Nk=>block%dims%Nk,Ns=>block%dims%Ns)
    if (allocated(primN)) deallocate(primN) ; allocate(primN(0:Ni,  0:Nj,  0:Nk  ))
                                              allocate(primC(0:Ni+1,0:Nj+1,0:Nk+1))
    call primN%alloc(Ns=Ns)
    primC = block%C(0:Ni+1,0:Nj+1,0:Nk+1)%P
#if !defined NULi && !defined NULj && !defined NULk
    ! 3D data
    mf = 0.125_R8P
#elif defined NULi
    primC(0   ,:,:) = 0._R8P
    primC(Ni+1,:,:) = 0._R8P
#if !defined NULj && !defined NULk
    ! 2D data
    mf = 0.25_R8P
#elif defined NULj
    ! 1D data
    mf = 0.5_R8P
    primC(:,0   ,:) = 0._R8P
    primC(:,Nj+1,:) = 0._R8P
#elif defined NULk
    ! 1D data
    mf = 0.5_R8P
    primC(:,:,0   ) = 0._R8P
    primC(:,:,Nk+1) = 0._R8P
#endif
#elif defined NULj
    primC(:,0   ,:) = 0._R8P
    primC(:,Nj+1,:) = 0._R8P
#if !defined NULi && !defined NULk
    ! 2D data
    mf = 0.25_R8P
#elif defined NULi
    ! 1D data
    mf = 0.5_R8P
    primC(0   ,:,:) = 0._R8P
    primC(Ni+1,:,:) = 0._R8P
#elif defined NULk
    ! 1D data
    mf = 0.5_R8P
    primC(:,:,0   ) = 0._R8P
    primC(:,:,Nk+1) = 0._R8P
#endif
#elif defined NULk
    primC(:,:,0   ) = 0._R8P
    primC(:,:,Nk+1) = 0._R8P
#if !defined NULi && !defined NULj
    ! 2D data
    mf = 0.25_R8P
#elif defined NULi
    ! 1D data
    mf = 0.5_R8P
    primC(0   ,:,:) = 0._R8P
    primC(Ni+1,:,:) = 0._R8P
#elif defined NULj
    ! 1D data
    mf = 0.5_R8P
    primC(:,0   ,:) = 0._R8P
    primC(:,Nj+1,:) = 0._R8P
#endif
#endif
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(i,j,k)         &
    !$OMP SHARED(Ni,Nj,Nk,primN,primC,mf)
    !$OMP DO
    do k=0,Nk
      do j=0,Nj
        do i=0,Ni
            primN(i,j,k) = primC(i+1,j+1,k+1) + primC(i,j+1,k+1) &
                         + primC(i+1,j  ,k+1) + primC(i,j,  k+1) &
                         + primC(i+1,j+1,k  ) + primC(i,j+1,k  ) &
                         + primC(i+1,j  ,k  ) + primC(i,j  ,k  )
            primN(i,j,k) = mf*primN(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL
    deallocate(primC)
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interpolate_primitive

  !> @brief Procedure for computing the metrics of structured block.
  elemental subroutine metrics_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block             !< Block-level data.
  type(Type_Vector)::                 NFS,s1,s2,nd,db   !< Dummy vector variables.
  real(R8P)::                         signi,signj,signk !< Dummy variables for checking the directions of normals.
  real(R8P)::                         Vx,Vy,Vz          !< Dummy variables for computing volume.
  real(R8P)::                         xp,yp,zp          !< Dummy variables for computing face coordinates.
  real(R8P)::                         xm,ym,zm          !< Dummy variables for computing face coordinates.
  integer(I4P)::                      i,j,k             !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing faces normals
  ! positioning at the middle of the block
  i = max(1,block%dims%Ni)
  j = max(1,block%dims%Nj)
  k = max(1,block%dims%Nk)
  ! checking the direction of i normals
  s1 = block%node(i,j  ,k) - block%node(i,j-1,k-1)
  s2 = block%node(i,j-1,k) - block%node(i,j,  k-1)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%node(i,  j,k)+block%node(i,  j-1,k)+block%node(i,  j,k-1)+block%node(i,  j-1,k-1))
  s2 = 0.25_R_P*(block%node(i-1,j,k)+block%node(i-1,j-1,k)+block%node(i-1,j,k-1)+block%node(i-1,j-1,k-1))
  db = s1 - s2
  signi = sign(1._R_P,(nd.dot.db))
  ! checking the direction of j normals
  s1 = block%node(i,j,k  ) - block%node(i-1,j,k-1)
  s2 = block%node(i,j,k-1) - block%node(i-1,j,k  )
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%node(i,j,  k)+block%node(i-1,j,  k)+block%node(i,j,  k-1)+block%node(i-1,j,  k-1))
  s2 = 0.25_R_P*(block%node(i,j-1,k)+block%node(i-1,j-1,k)+block%node(i,j-1,k-1)+block%node(i-1,j-1,k-1))
  db = s1 - s2
  signj = sign(1._R_P,(nd.dot.db))
  ! checking the direction of k normals
  s1 = block%node(i,  j,k) - block%node(i-1,j-1,k)
  s2 = block%node(i-1,j,k) - block%node(i,  j-1,k)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(block%node(i,j,k  )+block%node(i-1,j,k  )+block%node(i,j-1,k  )+block%node(i-1,j-1,k  ))
  s2 = 0.25_R_P*(block%node(i,j,k-1)+block%node(i-1,j,k-1)+block%node(i,j-1,k-1)+block%node(i-1,j-1,k-1))
  db = s1 - s2
  signk = sign(1._R_P,(nd.dot.db))
  !!$OMP PARALLEL DEFAULT(NONE)                        &
  !!$OMP PRIVATE(i,j,k,NFS,Vx,Vy,Vz,xp,yp,zp,xm,ym,zm) &
  !!$OMP SHARED(block,signi,signj,signk)
  !!$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=0,block%dims%Ni
        call NFS%face_normal4(pt1 = block%node(i,j-1,k-1), &
                              pt2 = block%node(i,j  ,k-1), &
                              pt3 = block%node(i,j  ,k  ), &
                              pt4 = block%node(i,j-1,k  ))
        NFS = NFS*signi
        block%Fi(i,j,k)%N = NFS%normalized()
        block%Fi(i,j,k)%S = NFS%normL2()
      enddo
    enddo
  enddo
  !!$OMP DO
  do k=1,block%dims%Nk
    do j=0,block%dims%Nj
      do i=1,block%dims%Ni
        call NFS%face_normal4(pt1 = block%node(i-1,j,k-1), &
                              pt2 = block%node(i-1,j,k  ), &
                              pt3 = block%node(i  ,j,k  ), &
                              pt4 = block%node(i  ,j,k-1))
        NFS = NFS*signj
        block%Fj(i,j,k)%N = NFS%normalized()
        block%Fj(i,j,k)%S = NFS%normL2()
      enddo
    enddo
  enddo
  !!$OMP DO
  do k=0,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
        call NFS%face_normal4(pt1 = block%node(i-1,j-1,k), &
                              pt2 = block%node(i  ,j-1,k), &
                              pt3 = block%node(i  ,j  ,k), &
                              pt4 = block%node(i-1,j  ,k))
        NFS = NFS*signk
        block%Fk(i,j,k)%N = NFS%normalized()
        block%Fk(i,j,k)%S = NFS%normL2()
      enddo
    enddo
  enddo
  ! computing finte volumes
  !!$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
        Vx = 0._R_P
        Vy = 0._R_P
        Vz = 0._R_P

        xp = 0.25_R_P*(block%node(i  ,j  ,k  )%x + block%node(i  ,j  ,k-1)%x + &
                       block%node(i  ,j-1,k  )%x + block%node(i  ,j-1,k-1)%x)
        yp = 0.25_R_P*(block%node(i  ,j  ,k  )%y + block%node(i  ,j  ,k-1)%y + &
                       block%node(i  ,j-1,k  )%y + block%node(i  ,j-1,k-1)%y)
        zp = 0.25_R_P*(block%node(i  ,j  ,k  )%z + block%node(i  ,j  ,k-1)%z + &
                       block%node(i  ,j-1,k  )%z + block%node(i  ,j-1,k-1)%z)
        xm = 0.25_R_P*(block%node(i-1,j  ,k  )%x + block%node(i-1,j  ,k-1)%x + &
                       block%node(i-1,j-1,k  )%x + block%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%node(i-1,j  ,k  )%y + block%node(i-1,j  ,k-1)%y + &
                       block%node(i-1,j-1,k  )%y + block%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%node(i-1,j  ,k  )%z + block%node(i-1,j  ,k-1)%z + &
                       block%node(i-1,j-1,k  )%z + block%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%Fi(i,j,k)%N%x*block%Fi(i,j,k)%S - xm*block%Fi(i-1,j,k)%N%x*block%Fi(i-1,j,k)%S
        Vy = Vy + yp*block%Fi(i,j,k)%N%y*block%Fi(i,j,k)%S - ym*block%Fi(i-1,j,k)%N%y*block%Fi(i-1,j,k)%S
        Vz = Vz + zp*block%Fi(i,j,k)%N%z*block%Fi(i,j,k)%S - zm*block%Fi(i-1,j,k)%N%z*block%Fi(i-1,j,k)%S

        xp = 0.25_R_P*(block%node(i  ,j  ,k  )%x + block%node(i  ,j  ,k-1)%x + &
                       block%node(i-1,j  ,k  )%x + block%node(i-1,j  ,k-1)%x)
        yp = 0.25_R_P*(block%node(i  ,j  ,k  )%y + block%node(i  ,j  ,k-1)%y + &
                       block%node(i-1,j  ,k  )%y + block%node(i-1,j  ,k-1)%y)
        zp = 0.25_R_P*(block%node(i  ,j  ,k  )%z + block%node(i  ,j  ,k-1)%z + &
                       block%node(i-1,j  ,k  )%z + block%node(i-1,j  ,k-1)%z)
        xm = 0.25_R_P*(block%node(i  ,j-1,k  )%x + block%node(i  ,j-1,k-1)%x + &
                       block%node(i-1,j-1,k  )%x + block%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%node(i  ,j-1,k  )%y + block%node(i  ,j-1,k-1)%y + &
                       block%node(i-1,j-1,k  )%y + block%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%node(i  ,j-1,k  )%z + block%node(i  ,j-1,k-1)%z + &
                       block%node(i-1,j-1,k  )%z + block%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%Fj(i,j,k)%N%x*block%Fj(i,j,k)%S - xm*block%Fj(i,j-1,k)%N%x*block%Fj(i,j-1,k)%S
        Vy = Vy + yp*block%Fj(i,j,k)%N%y*block%Fj(i,j,k)%S - ym*block%Fj(i,j-1,k)%N%y*block%Fj(i,j-1,k)%S
        Vz = Vz + zp*block%Fj(i,j,k)%N%z*block%Fj(i,j,k)%S - zm*block%Fj(i,j-1,k)%N%z*block%Fj(i,j-1,k)%S

        xp = 0.25_R_P*(block%node(i  ,j  ,k  )%x + block%node(i  ,j-1,k  )%x + &
                       block%node(i-1,j  ,k  )%x + block%node(i-1,j-1,k  )%x)
        yp = 0.25_R_P*(block%node(i  ,j  ,k  )%y + block%node(i  ,j-1,k  )%y + &
                       block%node(i-1,j  ,k  )%y + block%node(i-1,j-1,k  )%y)
        zp = 0.25_R_P*(block%node(i  ,j  ,k  )%z + block%node(i  ,j-1,k  )%z + &
                       block%node(i-1,j  ,k  )%z + block%node(i-1,j-1,k  )%z)
        xm = 0.25_R_P*(block%node(i  ,j  ,k-1)%x + block%node(i  ,j-1,k-1)%x + &
                       block%node(i-1,j  ,k-1)%x + block%node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(block%node(i  ,j  ,k-1)%y + block%node(i  ,j-1,k-1)%y + &
                       block%node(i-1,j  ,k-1)%y + block%node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(block%node(i  ,j  ,k-1)%z + block%node(i  ,j-1,k-1)%z + &
                       block%node(i-1,j  ,k-1)%z + block%node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*block%Fk(i,j,k)%N%x*block%Fk(i,j,k)%S - xm*block%Fk(i,j,k-1)%N%x*block%Fk(i,j,k-1)%S
        Vy = Vy + yp*block%Fk(i,j,k)%N%y*block%Fk(i,j,k)%S - ym*block%Fk(i,j,k-1)%N%y*block%Fk(i,j,k-1)%S
        Vz = Vz + zp*block%Fk(i,j,k)%N%z*block%Fk(i,j,k)%S - zm*block%Fk(i,j,k-1)%N%z*block%Fk(i,j,k-1)%S

        block%C(i,j,k)%V = max(Vx,Vy,Vz)
      enddo
    enddo
  enddo
  !!$OMP END PARALLEL
#ifdef NULi
  block%Fj(:,:,:)%N%x=0._R_P
  block%Fk(:,:,:)%N%x=0._R_P
  block%Fi(:,:,:)%N%y=0._R_P
  block%Fi(:,:,:)%N%z=0._R_P
#endif
#ifdef NULj
  block%Fi(:,:,:)%N%y=0._R_P
  block%Fk(:,:,:)%N%y=0._R_P
  block%Fj(:,:,:)%N%x=0._R_P
  block%Fj(:,:,:)%N%z=0._R_P
#endif
#ifdef NULk
  block%Fi(:,:,:)%N%z=0._R_P
  block%Fj(:,:,:)%N%z=0._R_P
  block%Fk(:,:,:)%N%x=0._R_P
  block%Fk(:,:,:)%N%y=0._R_P
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine metrics_block

  !> @brief Procedure for correcting the metrics of natural (and negative volume) boundary conditions cells of structured block.
  elemental subroutine metrics_correction_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block-level data.
  logical::                           correct  !< Flag for inquiring if metrics must be corrected.
  logical::                           wall     !< Flag for inquiring if bc is "wall-type": different corrections must be used.
  real(R8P)::                         tm       !< Tangential metrics parameter (-1 for wall-type bc).
  real(R8P)::                         sn       !< Normal metrics coefficient correction.
  integer(I4P)::                      i,j,k    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
  !!$OMP PARALLEL DEFAULT(NONE)            &
  !!$OMP PRIVATE(i,j,k,correct,wall,tm,sn) &
  !!$OMP SHARED(Ni,Nj,Nk,block)
  ! left i
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      correct = ((.not.block%Fi(0,j,k)%BC%is_adj()).or.(.not.block%Fi(0,j,k)%BC%is_per()).or.&
                 (block%C( 0,j,k)%V<(0.2_R_P*block%C(1,j,k)%V)))
      wall    = ((block%Fi(0,j,k)%BC%is_ext()).or.(block%Fi(0,j,k)%BC%is_ref()))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
        ! normal metrics
        sn=2._R_P*((block%Fi(1,j,k)%N*block%Fi(1,j,k)%S).dot.block%Fi(0,j,k)%N)
             block%Fi(-1,j,  k  )%N = -(block%Fi(1,j,k)%N*block%Fi(1,j,k)%S) + sn*block%Fi(0,j,k)%N
             block%Fi(-1,j,  k  )%S = block%Fi(-1,j,k)%N%normL2()
        call block%Fi(-1,j,  k  )%N%normalize
        ! tangential metrics
             block%Fj( 0,j  ,k  )%N = tm*block%Fj(1,j  ,k  )%N
             block%Fj( 0,j-1,k  )%N = tm*block%Fj(1,j-1,k  )%N
             block%Fj( 0,j  ,k-1)%N = tm*block%Fj(1,j  ,k-1)%N
             block%Fj( 0,j-1,k-1)%N = tm*block%Fj(1,j-1,k-1)%N
             block%Fj( 0,j  ,k  )%S = tm*block%Fj(1,j  ,k  )%S
             block%Fj( 0,j-1,k  )%S = tm*block%Fj(1,j-1,k  )%S
             block%Fj( 0,j  ,k-1)%S = tm*block%Fj(1,j  ,k-1)%S
             block%Fj( 0,j-1,k-1)%S = tm*block%Fj(1,j-1,k-1)%S

             block%Fk( 0,j  ,k  )%N = tm*block%Fk(1,j  ,k  )%N
             block%Fk( 0,j-1,k  )%N = tm*block%Fk(1,j-1,k  )%N
             block%Fk( 0,j  ,k-1)%N = tm*block%Fk(1,j  ,k-1)%N
             block%Fk( 0,j-1,k-1)%N = tm*block%Fk(1,j-1,k-1)%N
             block%Fk( 0,j  ,k  )%S = tm*block%Fk(1,j  ,k  )%S
             block%Fk( 0,j-1,k  )%S = tm*block%Fk(1,j-1,k  )%S
             block%Fk( 0,j  ,k-1)%S = tm*block%Fk(1,j  ,k-1)%S
             block%Fk( 0,j-1,k-1)%S = tm*block%Fk(1,j-1,k-1)%S
        ! volume
             block%C(  0,j,  k  )%V = block%C(    1,j,  k  )%V
      endif
    enddo
  enddo
  ! right i
  !!$OMP DO
  do k=1,Nk
    do j=1,Nj
      correct = ((.not.block%Fi(Ni+1,j,k)%BC%is_adj()).or.(.not.block%Fi(Ni+1,j,k)%BC%is_per()).or.&
                 (block%C( Ni+1,j,k)%V<(0.2_R_P*block%C(Ni,j,k)%V)))
      wall    = ((block%Fi(0,j,k)%BC%is_ext()).OR.(block%Fi(0,j,k)%BC%is_ref()))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
        ! normal metrics
        sn=2._R_P*((block%Fi(1,j,k)%N*block%Fi(1,j,k)%S).dot.block%Fi(0,j,k)%N)
             block%Fi(-1,j,  k  )%N = -(block%Fi(1,j,k)%N*block%Fi(1,j,k)%S) + sn*block%Fi(0,j,k)%N
             block%Fi(-1,j,  k  )%S = block%Fi(-1,j,k)%N%normL2()
        call block%Fi(-1,j,  k  )%N%normalize
        ! tangential metrics
             block%Fj( 0,j  ,k  )%N = tm*block%Fj(1,j  ,k  )%N
             block%Fj( 0,j-1,k  )%N = tm*block%Fj(1,j-1,k  )%N
             block%Fj( 0,j  ,k-1)%N = tm*block%Fj(1,j  ,k-1)%N
             block%Fj( 0,j-1,k-1)%N = tm*block%Fj(1,j-1,k-1)%N
             block%Fj( 0,j  ,k  )%S = tm*block%Fj(1,j  ,k  )%S
             block%Fj( 0,j-1,k  )%S = tm*block%Fj(1,j-1,k  )%S
             block%Fj( 0,j  ,k-1)%S = tm*block%Fj(1,j  ,k-1)%S
             block%Fj( 0,j-1,k-1)%S = tm*block%Fj(1,j-1,k-1)%S

             block%Fk( 0,j  ,k  )%N = tm*block%Fk(1,j  ,k  )%N
             block%Fk( 0,j-1,k  )%N = tm*block%Fk(1,j-1,k  )%N
             block%Fk( 0,j  ,k-1)%N = tm*block%Fk(1,j  ,k-1)%N
             block%Fk( 0,j-1,k-1)%N = tm*block%Fk(1,j-1,k-1)%N
             block%Fk( 0,j  ,k  )%S = tm*block%Fk(1,j  ,k  )%S
             block%Fk( 0,j-1,k  )%S = tm*block%Fk(1,j-1,k  )%S
             block%Fk( 0,j  ,k-1)%S = tm*block%Fk(1,j  ,k-1)%S
             block%Fk( 0,j-1,k-1)%S = tm*block%Fk(1,j-1,k-1)%S
        ! volume
             block%C(  0,j,  k  )%V = block%C(    1,j,  k  )%V
      endif
    enddo
  enddo
  ! left j
  !!$OMP DO
  do k=1,Nk
    do i=1,Ni
      correct = ((.not.block%Fj(i,0,k)%BC%is_adj()).or.(.not.block%Fj(i,0,k)%BC%is_per()).or.&
                 (block%C( i,0,k)%V<(0.2_R_P*block%C(i,1,k)%V)))
      wall    = ((block%Fj(i,0,k)%BC%is_ext()).OR.(block%Fj(i,0,k)%BC%is_ref()))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
        ! normal metrics
        sn=2._R_P*((block%Fj(i,1,k)%N*block%Fj(i,1,k)%S).dot.block%Fj(i,0,k)%N)
             block%Fj(i, -1,k  )%N = -(block%Fj(i,1,k)%N*block%Fj(i,1,k)%S) + sn*block%Fj(i,0,k)%N
             block%Fj(i, -1,k  )%S = block%Fj(i,-1,k)%N%normL2()
        call block%Fj(i, -1,k  )%N%normalize
        ! tangential metrics
             block%Fi(i  ,0,k  )%N = tm*block%Fi(i  ,1,k  )%N
             block%Fi(i-1,0,k  )%N = tm*block%Fi(i-1,1,k  )%N
             block%Fi(i  ,0,k-1)%N = tm*block%Fi(i  ,1,k-1)%N
             block%Fi(i-1,0,k-1)%N = tm*block%Fi(i-1,1,k-1)%N
             block%Fi(i  ,0,k  )%S = tm*block%Fi(i  ,1,k  )%S
             block%Fi(i-1,0,k  )%S = tm*block%Fi(i-1,1,k  )%S
             block%Fi(i  ,0,k-1)%S = tm*block%Fi(i  ,1,k-1)%S
             block%Fi(i-1,0,k-1)%S = tm*block%Fi(i-1,1,k-1)%S

             block%Fk(i  ,0,k  )%N = tm*block%Fk(i  ,1,k  )%N
             block%Fk(i-1,0,k  )%N = tm*block%Fk(i-1,1,k  )%N
             block%Fk(i  ,0,k-1)%N = tm*block%Fk(i  ,1,k-1)%N
             block%Fk(i-1,0,k-1)%N = tm*block%Fk(i-1,1,k-1)%N
             block%Fk(i  ,0,k  )%S = tm*block%Fk(i  ,1,k  )%S
             block%Fk(i-1,0,k  )%S = tm*block%Fk(i-1,1,k  )%S
             block%Fk(i  ,0,k-1)%S = tm*block%Fk(i  ,1,k-1)%S
             block%Fk(i-1,0,k-1)%S = tm*block%Fk(i-1,1,k-1)%S
        ! volume
             block%C( i,  0,k  )%V = block%C(    i,  1,k  )%V
      endif
    enddo
  enddo
  ! right j
  !!$OMP DO
  do k=1,Nk
    do i=1,Ni
      correct = ((.not.block%Fj(i,Nj+1,k)%BC%is_adj()).or.(.not.block%Fj(i,Nj+1,k)%BC%is_per()).or.&
                 (block%C( i,Nj+1,k)%V<(0.2_R_P*block%C(i,Nj,k)%V)))
      wall    = ((block%Fj(i,Nj+1,k)%BC%is_ext()).OR.(block%Fj(i,Nj+1,k)%BC%is_ref()))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
        ! normal metrics
        sn=2._R_P*((block%Fj(i,Nj-1,k)%N*block%Fj(i,Nj-1,k)%S).dot.block%Fj(i,Nj,k)%N)
             block%Fj(i,Nj+1,  k  )%N = -(block%Fj(i,Nj-1,k)%N*block%Fj(i,Nj-1,k)%S) + sn*block%Fj(i,Nj,k)%N
             block%Fj(i,Nj+1,  k  )%S = block%Fj(i,Nj+1,k)%N%normL2()
        call block%Fj(i,Nj+1,  k  )%N%normalize
        ! tangential metrics
             block%Fi(i  ,Nj+1,k  )%N = tm*block%Fi(i  ,Nj,k  )%N
             block%Fi(i-1,Nj+1,k  )%N = tm*block%Fi(i-1,Nj,k  )%N
             block%Fi(i  ,Nj+1,k-1)%N = tm*block%Fi(i  ,Nj,k-1)%N
             block%Fi(i-1,Nj+1,k-1)%N = tm*block%Fi(i-1,Nj,k-1)%N
             block%Fi(i  ,Nj+1,k  )%S = tm*block%Fi(i  ,Nj,k  )%S
             block%Fi(i-1,Nj+1,k  )%S = tm*block%Fi(i-1,Nj,k  )%S
             block%Fi(i  ,Nj+1,k-1)%S = tm*block%Fi(i  ,Nj,k-1)%S
             block%Fi(i-1,Nj+1,k-1)%S = tm*block%Fi(i-1,Nj,k-1)%S

             block%Fk(i  ,Nj+1,k  )%N = tm*block%Fk(i  ,Nj,k  )%N
             block%Fk(i-1,Nj+1,k  )%N = tm*block%Fk(i-1,Nj,k  )%N
             block%Fk(i  ,Nj+1,k-1)%N = tm*block%Fk(i  ,Nj,k-1)%N
             block%Fk(i-1,Nj+1,k-1)%N = tm*block%Fk(i-1,Nj,k-1)%N
             block%Fk(i  ,Nj+1,k  )%S = tm*block%Fk(i  ,Nj,k  )%S
             block%Fk(i-1,Nj+1,k  )%S = tm*block%Fk(i-1,Nj,k  )%S
             block%Fk(i  ,Nj+1,k-1)%S = tm*block%Fk(i  ,Nj,k-1)%S
             block%Fk(i-1,Nj+1,k-1)%S = tm*block%Fk(i-1,Nj,k-1)%S
        ! volume
             block%C( i,  Nj+1,k  )%V = block%C(    i,  Nj,k  )%V
      endif
    enddo
  enddo
  ! left k
  !!$OMP DO
  do j=1,Nj
    do i=1,Ni
      correct = ((.not.block%Fk(i,j,0)%BC%is_adj()).or.(.not.block%Fk(i,j,0)%BC%is_per()).or.&
                 (block%C( i,j,0)%V<(0.2_R_P*block%C(i,j,1)%V)))
      wall    = ((block%Fk(i,j,0)%BC%is_ext()).OR.(block%Fk(i,j,0)%BC%is_ref()))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
        ! normal metrics
        sn=2._R_P*((block%Fk(i,j,1)%N*block%Fk(i,j,1)%S).dot.block%Fk(i,j,0)%N)
             block%Fk(i,  j, -1)%N = -(block%Fk(i,j,1)%N*block%Fk(i,j,1)%S) + sn*block%Fk(i,j,0)%N
             block%Fk(i,  j, -1)%S = block%Fk(i,j,-1)%N%normL2()
        call block%Fk(i,  j, -1)%N%normalize
        ! tangential metrics
             block%Fi(i  ,j  ,0)%N = tm*block%Fi(i  ,j  ,1)%N
             block%Fi(i-1,j  ,0)%N = tm*block%Fi(i-1,j  ,1)%N
             block%Fi(i  ,j-1,0)%N = tm*block%Fi(i  ,j-1,1)%N
             block%Fi(i-1,j-1,0)%N = tm*block%Fi(i-1,j-1,1)%N
             block%Fi(i  ,j  ,0)%S = tm*block%Fi(i  ,j  ,1)%S
             block%Fi(i-1,j  ,0)%S = tm*block%Fi(i-1,j  ,1)%S
             block%Fi(i  ,j-1,0)%S = tm*block%Fi(i  ,j-1,1)%S
             block%Fi(i-1,j-1,0)%S = tm*block%Fi(i-1,j-1,1)%S

             block%Fj(i  ,j  ,0)%N = tm*block%Fj(i  ,j  ,1)%N
             block%Fj(i-1,j  ,0)%N = tm*block%Fj(i-1,j  ,1)%N
             block%Fj(i  ,j-1,0)%N = tm*block%Fj(i  ,j-1,1)%N
             block%Fj(i-1,j-1,0)%N = tm*block%Fj(i-1,j-1,1)%N
             block%Fj(i  ,j  ,0)%S = tm*block%Fj(i  ,j  ,1)%S
             block%Fj(i-1,j  ,0)%S = tm*block%Fj(i-1,j  ,1)%S
             block%Fj(i  ,j-1,0)%S = tm*block%Fj(i  ,j-1,1)%S
             block%Fj(i-1,j-1,0)%S = tm*block%Fj(i-1,j-1,1)%S
        ! volume
             block%C( i,  j,  0)%V = block%C(    i,  j,  1)%V
      endif
    enddo
  enddo
  ! right k
  !!$OMP DO
  do j=1,Nj
    do i=1,Ni
      correct = ((.not.block%Fk(i,j,Nk+1)%BC%is_adj()).or.(.not.block%Fk(i,j,Nk+1)%BC%is_per()).or.&
                 (block%C( i,j,Nk+1)%V<(0.2_R_P*block%C(i,j,Nk)%V)))
      wall    = ((block%Fk(i,j,Nk+1)%BC%is_ext()).OR.(block%Fk(i,j,Nk+1)%BC%is_ref()))
      tm = 1._R_P
      if (wall) tm = -1._R_P
      if (correct) then
        ! normal metrics
        sn=2._R_P*((block%Fk(i,j,Nk-1)%N*block%Fk(i,j,Nk-1)%S).dot.block%Fk(i,j,Nk)%N)
             block%Fk(i,  j,  Nk+1)%N = -(block%Fk(i,j,Nk-1)%N*block%Fk(i,j,Nk-1)%S) + sn*block%Fk(i,j,Nk)%N
             block%Fk(i,  j,  Nk+1)%S = block%Fk(i,j,Nk+1)%N%normL2()
        call block%Fk(i,  j,  Nk+1)%N%normalize
        ! tangential metrics
             block%Fi(i  ,j  ,Nk+1)%N = tm*block%Fi(i  ,j  ,Nk)%N
             block%Fi(i-1,j  ,Nk+1)%N = tm*block%Fi(i-1,j  ,Nk)%N
             block%Fi(i  ,j-1,Nk+1)%N = tm*block%Fi(i  ,j-1,Nk)%N
             block%Fi(i-1,j-1,Nk+1)%N = tm*block%Fi(i-1,j-1,Nk)%N
             block%Fi(i  ,j  ,Nk+1)%S = tm*block%Fi(i  ,j  ,Nk)%S
             block%Fi(i-1,j  ,Nk+1)%S = tm*block%Fi(i-1,j  ,Nk)%S
             block%Fi(i  ,j-1,Nk+1)%S = tm*block%Fi(i  ,j-1,Nk)%S
             block%Fi(i-1,j-1,Nk+1)%S = tm*block%Fi(i-1,j-1,Nk)%S

             block%Fj(i  ,j  ,Nk+1)%N = tm*block%Fj(i  ,j  ,Nk)%N
             block%Fj(i-1,j  ,Nk+1)%N = tm*block%Fj(i-1,j  ,Nk)%N
             block%Fj(i  ,j-1,Nk+1)%N = tm*block%Fj(i  ,j-1,Nk)%N
             block%Fj(i-1,j-1,Nk+1)%N = tm*block%Fj(i-1,j-1,Nk)%N
             block%Fj(i  ,j  ,Nk+1)%S = tm*block%Fj(i  ,j  ,Nk)%S
             block%Fj(i-1,j  ,Nk+1)%S = tm*block%Fj(i-1,j  ,Nk)%S
             block%Fj(i  ,j-1,Nk+1)%S = tm*block%Fj(i  ,j-1,Nk)%S
             block%Fj(i-1,j-1,Nk+1)%S = tm*block%Fj(i-1,j-1,Nk)%S
        ! volume
             block%C( i,  j,  Nk+1)%V = block%C(    i,  j,  Nk)%V
      endif
    enddo
  enddo
  !!$OMP END PARALLEL
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine metrics_correction_block

  !> @brief Procedure for computing cell center coordinates from cell nodes ones.
  elemental subroutine node2center_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block  !< Block-level data.
  integer(I4P)::                      i,j,k !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block)
  !$OMP DO
  do k=1-block%dims%gc(5),block%dims%Nk+block%dims%gc(6)
    do j=1-block%dims%gc(3),block%dims%Nj+block%dims%gc(4)
      do i=1-block%dims%gc(1),block%dims%Ni+block%dims%gc(2)
        block%C(i,j,k)%cent%x = (block%node(i,  j,  k  )%x + &
                                 block%node(i-1,j,  k  )%x + &
                                 block%node(i  ,j-1,k  )%x + &
                                 block%node(i  ,j  ,k-1)%x + &
                                 block%node(i-1,j-1,k-1)%x + &
                                 block%node(i  ,j-1,k-1)%x + &
                                 block%node(i-1,j  ,k-1)%x + &
                                 block%node(i-1,j-1,k  )%x)*0.125_R_P
        block%C(i,j,k)%cent%y = (block%node(i,  j,  k  )%y + &
                                 block%node(i-1,j,  k  )%y + &
                                 block%node(i  ,j-1,k  )%y + &
                                 block%node(i  ,j  ,k-1)%y + &
                                 block%node(i-1,j-1,k-1)%y + &
                                 block%node(i  ,j-1,k-1)%y + &
                                 block%node(i-1,j  ,k-1)%y + &
                                 block%node(i-1,j-1,k  )%y)*0.125_R_P
        block%C(i,j,k)%cent%z = (block%node(i,  j,  k  )%z + &
                                 block%node(i-1,j,  k  )%z + &
                                 block%node(i  ,j-1,k  )%z + &
                                 block%node(i  ,j  ,k-1)%z + &
                                 block%node(i-1,j-1,k-1)%z + &
                                 block%node(i  ,j-1,k-1)%z + &
                                 block%node(i-1,j  ,k-1)%z + &
                                 block%node(i-1,j-1,k  )%z)*0.125_R_P
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine node2center_block

  !> @brief Procedure for computing (an estimation of) the minimum value of space step.
  elemental function min_space_step_block(block) result(mss)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block !< Block-level data.
  real(R8P)::                      mss   !< Minimum space step.
  integer(I4P)::                   i,j,k !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mss = MaxR8P
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block)          &
  !$OMP REDUCTION(min: mss)
  !$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
          mss = min(mss,                                                               &
                    (0.25_R_P*(block%node(i  ,j  ,k  )%x + block%node(i  ,j-1,k  )%x +  &
                               block%node(i  ,j-1,k-1)%x + block%node(i  ,j  ,k-1)%x)-  &
                     0.25_R_P*(block%node(i-1,j  ,k  )%x + block%node(i-1,j-1,k  )%x +  &
                               block%node(i-1,j-1,k-1)%x + block%node(i-1,j  ,k-1)%x)), &
                    (0.25_R_P*(block%node(i  ,j  ,k  )%y + block%node(i-1,j  ,k  )%y +  &
                               block%node(i-1,j  ,k-1)%y + block%node(i  ,j  ,k-1)%y)-  &
                     0.25_R_P*(block%node(i  ,j-1,k  )%y + block%node(i-1,j-1,k  )%y +  &
                               block%node(i-1,j-1,k-1)%y + block%node(i  ,j-1,k-1)%y)), &
                    (0.25_R_P*(block%node(i  ,j  ,k  )%z + block%node(i-1,j-1,k  )%z +  &
                               block%node(i-1,j  ,k  )%z + block%node(i  ,j-1,k  )%z)-  &
                     0.25_R_P*(block%node(i  ,j  ,k-1)%z + block%node(i-1,j-1,k-1)%z +  &
                               block%node(i-1,j  ,k-1)%z + block%node(i  ,j-1,k-1)%z)))
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction min_space_step_block

  !> @brief Procedure for creating a uniform grid provided block extents.
  !> @note This procedure reallocate block variables calling block%alloc procedure, thus all previous array data are deleted.
  elemental subroutine create_uniform_grid(block,exts)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),                 intent(INOUT):: block   !< Block-level data.
  type(Type_Block_Extents), optional, intent(IN)::    exts    !< Block extents.
  type(Type_Vector)::                                 deltas  !< Uniform space steps.
  type(Type_Vector)::                                 Nijk    !< Number of cells along each direction in vector form.
  integer(I4P)::                                      i,j,k   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(exts)) block%exts = exts
  call block%alloc
  associate(gc => block%dims%gc, Ni => block%dims%Ni, Nj => block%dims%Nj, Nk => block%dims%Nk)
    Nijk = real(Ni,R8P)*ex+real(Nj,R8P)*ey+real(Nk,R8P)*ez
    deltas = (block%exts%emax-block%exts%emin)/Nijk
    do k=0-gc(5),Nk+gc(6)
      do j=0-gc(3),Nj+gc(4)
        do i=0-gc(1),Ni+gc(2)
          block%node(i,j,k) = block%exts%emin + (real(i,R8P)*deltas%x)*ex + (real(j,R8P)*deltas%y)*ey + (real(k,R8P)*deltas%z)*ez
        enddo
      enddo
    enddo
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create_uniform_grid

  !> @brief Procedure for creating the grid from the grid of a finer block.
  !> @note Presently, it is assumed that the finer block has a doubled number of cells along each direction, namely
  !> \f$ Ni_f = 2Ni\,  Nj_f = 2Nj\,  Nk_f = 2Nk\f$, \f$Ni(j,k)_f\f$ being the finer block dimensions.
  elemental subroutine create_grid_from_finer(block,block_f)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),       intent(INOUT):: block   !< Block-level data.
  type(Type_SBlock),        intent(IN)::    block_f !< Finer block-level data.
  integer(I4P)::                            i,j,k   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(gc => block%dims%gc, Ni => block%dims%Ni, Nj => block%dims%Nj, Nk => block%dims%Nk)
    ! inner nodes
    do k=0,Nk
      do j=0,Nj
        do i=0,Ni
          block%node(i,j,k) = block_f%node(i*2,j*2,k*2)
        enddo
      enddo
    enddo
    ! left i
    ! to be completed...
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create_grid_from_finer

  !> @brief Procedure for setting cells boundary conditions from block ones.
  !> @note This procedure is useful only for Cartesian grids.
  elemental subroutine set_cells_bc(block,mesh_dims,l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),         intent(INOUT):: block     !< Block-level data.
  type(Type_Mesh_Dimensions), intent(IN)::    mesh_dims !< Mesh dimensions.
  integer(I4P),               intent(IN)::    l         !< Grid level.
  integer(I4P)::                              i,j,k     !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(gc => block%dims%gc, Ni => block%dims%Ni, Nj => block%dims%Nj, Nk => block%dims%Nk)
    ! left i
    do k=1,Nk ; do j=1,Nj ; do i=1-gc(1),0
      block%Fi(i,j,k)%BC%tp = block%BC%F(1)%tp ; call block%Fi(i,j,k)%BC%alloc
      if (block%BC%F(1)%is_adj()) then
        block%Fi(i,j,k)%BC%adj%b = block%BC%F(1)%adj%b
        block%Fi(i,j,k)%BC%adj%i = mesh_dims%block_dims(block%BC%F(1)%adj%b,l)%Ni/(2**(l-1))+i
        block%Fi(i,j,k)%BC%adj%j = j
        block%Fi(i,j,k)%BC%adj%k = k
      elseif (block%BC%F(1)%is_in1()) then
        block%Fi(i,j,k)%BC%inf = block%BC%F(1)%inf
      endif
    enddo ; enddo ; enddo
     ! right i
    do k=1,Nk ; do j=1,Nj ; do i=Ni+1,Ni+gc(2)
      block%Fi(i-1,j,k)%BC%tp = block%BC%F(2)%tp ; call block%Fi(i-1,j,k)%BC%alloc
      if (block%BC%F(2)%is_adj()) then
        block%Fi(i-1,j,k)%BC%adj%b = block%BC%F(2)%adj%b
        block%Fi(i-1,j,k)%BC%adj%i = i-block%dims%Ni/(2**(l-1))
        block%Fi(i-1,j,k)%BC%adj%j = j
        block%Fi(i-1,j,k)%BC%adj%k = k
      elseif (block%BC%F(2)%is_in1()) then
        block%Fi(i-1,j,k)%BC%inf = block%BC%F(2)%inf
      endif
    enddo ; enddo ; enddo
    ! left j
    do k=1,Nk ; do i=1,Ni ; do j=1-gc(3),0
      block%Fj(i,j,k)%BC%tp = block%BC%F(3)%tp ; call block%Fj(i,j,k)%BC%alloc
      if (block%BC%F(3)%is_adj()) then
        block%Fj(i,j,k)%BC%adj%b = block%BC%F(3)%adj%b
        block%Fj(i,j,k)%BC%adj%i = i
        block%Fj(i,j,k)%BC%adj%j = mesh_dims%block_dims(block%BC%F(3)%adj%b,l)%Nj/(2**(l-1))+j
        block%Fj(i,j,k)%BC%adj%k = k
      elseif (block%BC%F(3)%is_in1()) then
        block%Fj(i,j,k)%BC%inf = block%BC%F(3)%inf
      endif
    enddo ; enddo ; enddo
    ! right j
    do k=1,Nk ; do i=1,Ni ; do j=Nj+1,Nj+gc(4)
      block%Fj(i,j-1,k)%BC%tp = block%BC%F(4)%tp ; call block%Fj(i,j-1,k)%BC%alloc
      if (block%BC%F(4)%is_adj()) then
        block%Fj(i,j-1,k)%BC%adj%b = block%BC%F(4)%adj%b
        block%Fj(i,j-1,k)%BC%adj%i = i
        block%Fj(i,j-1,k)%BC%adj%j = j-block%dims%Nj/(2**(l-1))
        block%Fj(i,j-1,k)%BC%adj%k = k
      elseif (block%BC%F(4)%is_in1()) then
        block%Fj(i,j-1,k)%BC%inf = block%BC%F(4)%inf
      endif
    enddo ; enddo ; enddo
    ! left k
    do j=1,Nj ; do i=1,Ni ; do k=1-gc(5),0
      block%Fk(i,j,k)%BC%tp = block%BC%F(5)%tp ; call block%Fk(i,j,k)%BC%alloc
      if (block%BC%F(5)%is_adj()) then
        block%Fk(i,j,k)%BC%adj%b = block%BC%F(5)%adj%b
        block%Fk(i,j,k)%BC%adj%i = i
        block%Fk(i,j,k)%BC%adj%j = j
        block%Fk(i,j,k)%BC%adj%k = mesh_dims%block_dims(block%BC%F(5)%adj%b,l)%Nk/(2**(l-1))+k
      elseif (block%BC%F(5)%is_in1()) then
        block%Fk(i,j,k)%BC%inf = block%BC%F(5)%inf
      endif
    enddo ; enddo ; enddo
    ! right k
    do j=1,Nj ; do i=1,Ni ; do k=Nk+1,Nk+gc(6)
      block%Fk(i,j,k-1)%BC%tp = block%BC%F(6)%tp ; call block%Fk(i,j,k-1)%BC%alloc
      if (block%BC%F(6)%is_adj()) then
        block%Fk(i,j,k-1)%BC%adj%b = block%BC%F(6)%adj%b
        block%Fk(i,j,k-1)%BC%adj%i = i
        block%Fk(i,j,k-1)%BC%adj%j = j
        block%Fk(i,j,k-1)%BC%adj%k = k-block%dims%Nk/(2**(l-1))
      elseif (block%BC%F(6)%is_in1()) then
        block%Fk(i,j,k-1)%BC%inf = block%BC%F(6)%inf
      endif
    enddo ; enddo ; enddo
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_cells_bc

  !> @brief Procedure for setting initial condition in a region of block.
  subroutine set_region_ic(block,region,species0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block       !< Block data.
  type(Type_Region),  intent(IN)::    region      !< Region data.
  type(Type_Species), intent(IN)::    species0    !< Initial species.
  type(Type_Vector)::                 radius      !< Radius from region center.
  !integer(I4P)::                      icc,jcc,kcc !< Indexes of the cell closest to the region center.
  integer(I4P)::                      i,j,k       !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call block%node2center
  associate(gc => block%dims%gc(1:6),Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
    !call find_closest_to_center(icc=icc,jcc=jcc,kcc=kcc)
    select case(trim(adjustl(Upper_Case(region%shtype))))
    case('CYLINDER-X')! x-axis aligned
    case('CYLINDER-Y')! y-axis aligned
    case('CYLINDER-Z')! z-axis aligned
      do k=1,Nk
        do j=1,Nj
          do i=1,Ni
            if (block%C(i,j,k)%cent%z>=region%center%z.and.block%C(i,j,k)%cent%z<=region%height) then
              radius = block%C(i,j,k)%cent-region%center
              if (radius%normL2()<region%radius) then
                block%C(i,j,k)%P=region%prim
                call block%C(i,j,k)%P%compute_d
                call block%C(i,j,k)%P%compute_g(species0=species0)
              endif
            endif
          enddo
        enddo
      enddo
    case default
    endselect
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  !contains
   !subroutine find_closest_to_center(icc,jcc,kcc)
   !!-------------------------------------------------------------------------------------------------------------------------------
   !implicit none
   !integer(I4P), intent(OUT):: icc,jcc,kcc !< Indexes of the cell closest to the region center.
   !type(Type_Vector)::         distance    !< Distance counter.
   !real(R8P)::                 dist_m      !< Minimum distance.
   !integer(I4P)::              i,j,k       !< Counters.
   !!-------------------------------------------------------------------------------------------------------------------------------
   !
   !!-------------------------------------------------------------------------------------------------------------------------------
   !icc = 0 ; jcc = 0 ; kcc = 0
   !dist_m = MaxR8P
   !associate(gc => block%dims%gc(1:6),Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
   !  do k=1,Nk
   !    do j=1,Nj
   !      do i=1,Ni
   !        distance = region%center-block%C(i,j,k)%cent
   !        if (distance%normL2()<=dist_m) then
   !          dist_m = distance%normL2()
   !          icc = i ; jcc = j ; kcc = k
   !        endif
   !      enddo
   !    enddo
   !  enddo
   !endassociate
   !return
   !!-------------------------------------------------------------------------------------------------------------------------------
   !endsubroutine find_closest_to_center
  endsubroutine set_region_ic

  !> @brief Procedure for converting primitive variables to conservative variables of all cells of a block.
  !> @note Only the inner cells of the block are converted.
  elemental subroutine primitive2conservative(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block !< Block-level data.
  integer(I4P)::                      i,j,k !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block)
  !$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
        call block%C(i,j,k)%prim2cons
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine primitive2conservative

  !> @brief Procedure for converting conservative variables to primitive variables of all cells of a block.
  !> @note Only the inner cells of the block are converted.
  elemental subroutine conservative2primitive(block,species0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block-level data.
  type(Type_Species), intent(IN)::    species0 !< Initial species.
  integer(I_P)::                      i,j,k    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block,species0)
  !$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
        call block%C(i,j,k)%cons2prim(species0=species0)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine conservative2primitive

  !> @brief Procedure for updating primitive variables of all cells of a block.
  !> @note Only the inner cells of the block are updated.
  elemental subroutine update_primitive(block,species0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block-level data.
  type(Type_Species), intent(IN)::    species0 !< Initial species.
  integer(I_P)::                      i,j,k    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block,species0)
  !$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
        call block%C(i,j,k)%P%compute_d
        call block%C(i,j,k)%P%compute_g(species0=species0)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update_primitive

  !> @brief Procedure for evaluating the local and global time step value by CFL condition.
  elemental subroutine compute_time(block,time_step,Dtmin)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),   intent(INOUT):: block                         !< Block-level data.
  type(Type_Time_Step), intent(IN)::    time_step                     !< Time-stepping data.
  real(R8P),            intent(OUT)::   Dtmin                         !< Minimum Dt.
  real(R8P)::                           vmax                          !< Maximum speed of waves.
  real(R8P)::                           ss                            !< Speed of sound.
  real(R8P)::                           vmiL,vmiR,vmjL,vmjR,vmkL,vmkR !< Dummy velocities.
  type(Type_Vector)::                   vm                            !< Dummy vectorial velocities.
  integer(I4P)::                        i,j,k                         !< Space counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the minimum Dt into the inner cells
  !$OMP PARALLEL DEFAULT(NONE)                                   &
  !$OMP PRIVATE(i,j,k,vmax,ss,vmiL,vmiR,vmjL,vmjR,vmkL,vmkR,vm)  &
  !$OMP SHARED(block,Dtmin,time_step)
  !$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
        ! computing the local speed of sound
        ss = a(p=block%C(i,j,k)%P%p,r=block%C(i,j,k)%P%d,g=block%C(i,j,k)%P%g)
        ! evaluating the maximum propagation speed of acoustic segnals multiplied for face area
        ! left i
        vm   = 0.5_R_P*(block%C(i-1,j,k)%P%v+block%C(i,j,k)%P%v)
        vmiL = (vm.dot.block%Fi(i-1,j,k)%N)*block%Fi(i-1,j,k)%S
        vmiL = abs(vmiL) + ss
        ! right i
        vm   = 0.5_R_P*(block%C(i,j,k)%P%v+block%C(i+1,j,k)%P%v)
        vmiR = (vm.dot.block%Fi(i,j,k)%N)*block%Fi(i,j,k)%S
        vmiR = abs(vmiR) + ss
        ! left j
        vm   = 0.5_R_P*(block%C(i,j-1,k)%P%v+block%C(i,j,k)%P%v)
        vmjL = (vm.dot.block%Fj(i,j-1,k)%N)*block%Fj(i,j-1,k)%S
        vmjL = abs(vmjL) + ss
        ! right j
        vm   = 0.5_R_P*(block%C(i,j,k)%P%v+block%C(i,j+1,k)%P%v)
        vmjR = (vm.dot.block%Fj(i,j,k)%N)*block%Fj(i,j,k)%S
        vmjR = abs(vmjR) + ss
        ! left k
        vm   = 0.5_R_P*(block%C(i,j,k-1)%P%v+block%C(i,j,k)%P%v)
        vmkL = (vm.dot.block%Fk(i,j,k-1)%N)*block%Fk(i,j,k-1)%S
        vmkL = abs(vmkL) + ss
        ! right k
        vm   = 0.5_R_P*(block%C(i,j,k)%P%v+block%C(i,j,k+1)%P%v)
        vmkR = (vm.dot.block%Fk(i,j,k)%N)*block%Fk(i,j,k)%S
        vmkR = abs(vmkR) + ss
        ! vmax
        vmax = max(vmiL,vmiR,vmjL,vmjR,vmkL,vmkR)
        block%C(i,j,k)%Dt = block%C(i,j,k)%V/vmax*time_step%CFL
      enddo
    enddo
  enddo
  ! computing minimum Dt
  !$OMP SINGLE
  Dtmin = minval(block%C(1:block%dims%Ni,1:block%dims%Nj,1:block%dims%Nk)%Dt)
  !$OMP END SINGLE
  ! ghost cells estrapolation: imposing the minum value of Dt
  ! left i frame
  !$OMP DO
  do k=1-block%dims%gc(5),block%dims%Nk+block%dims%gc(6)
    do j=1-block%dims%gc(3),block%dims%Nj+block%dims%gc(4)
      do i=1-block%dims%gc(1),0
        block%C(i,j,k)%Dt = Dtmin
      enddo
    enddo
  enddo
  ! right i frame
  !$OMP DO
  do k=1-block%dims%gc(5),block%dims%Nk+block%dims%gc(6)
    do j=1-block%dims%gc(3),block%dims%Nj+block%dims%gc(4)
      do i=block%dims%Ni+1,block%dims%Ni+block%dims%gc(2)
        block%C(i,j,k)%Dt = Dtmin
      enddo
    enddo
  enddo
  ! left j frame
  !$OMP DO
  do k=1-block%dims%gc(5),block%dims%Nk+block%dims%gc(6)
    do j=1-block%dims%gc(3),0
      do i=1-block%dims%gc(1),block%dims%Ni+block%dims%gc(2)
        block%C(i,j,k)%Dt = Dtmin
      enddo
    enddo
  enddo
  ! right j frame
  !$OMP DO
  do k=1-block%dims%gc(5),block%dims%Nk+block%dims%gc(6)
    do j=block%dims%Nj+1,block%dims%Nj+block%dims%gc(4)
      do i=1-block%dims%gc(1),block%dims%Ni+block%dims%gc(2)
        block%C(i,j,k)%Dt = Dtmin
      enddo
    enddo
  enddo
  ! left k frame
  !$OMP DO
  do k=1-block%dims%gc(5),0
    do j=1-block%dims%gc(3),block%dims%Nj+block%dims%gc(4)
      do i=1-block%dims%gc(1),block%dims%Ni+block%dims%gc(2)
        block%C(i,j,k)%Dt = Dtmin
      enddo
    enddo
  enddo
  ! right k frame
  !$OMP DO
  do k=block%dims%Nk+1,block%dims%Nk+block%dims%gc(6)
    do j=1-block%dims%gc(3),block%dims%Nj+block%dims%gc(4)
      do i=1-block%dims%gc(1),block%dims%Ni+block%dims%gc(2)
        block%C(i,j,k)%Dt = Dtmin
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_time

  !> @brief Procedure for computing the residuals. This the space operator. The residuals are stored in block%C%KS(s1) conservative
  !> variables.
  elemental subroutine residuals(block,s1,space_step,species0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),    intent(INOUT):: block      !< Block-level data.
  integer(I1P),          intent(IN)::    s1         !< Current Runge-kutta stage.
  type(Type_Space_Step), intent(IN)::    space_step !< Space stepping data.
  type(Type_Species),    intent(IN)::    species0   !< Initial species.
  type(Type_Conservative), allocatable:: Fic(:,:,:) !< I convective fluxes.
  type(Type_Conservative), allocatable:: Fjc(:,:,:) !< J convective fluxes.
  type(Type_Conservative), allocatable:: Fkc(:,:,:) !< K convective fluxes.
  type(Type_Conservative), allocatable:: Fid(:,:,:) !< I diffusive fluxes.
  type(Type_Conservative), allocatable:: Fjd(:,:,:) !< J diffusive fluxes.
  type(Type_Conservative), allocatable:: Fkd(:,:,:) !< K diffusive fluxes.
  integer(I1P)::                         gcu        !< Number of ghost cells used.
  integer(I4P)::                         i,j,k,s    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(gc=>block%dims%gc,Ni=>block%dims%Ni,Nj=>block%dims%Nj,Nk=>block%dims%Nk,Ns=>block%dims%Ns,&
            gco=>space_step%gco,inviscid=>space_step%inviscid)
    allocate(Fic(0:Ni,1:Nj,1:Nk)) ; call Fic%alloc(Ns=Ns)
    allocate(Fjc(1:Ni,0:Nj,1:Nk)) ; call Fjc%alloc(Ns=Ns)
    allocate(Fkc(1:Ni,1:Nj,0:Nk)) ; call Fkc%alloc(Ns=Ns)
    allocate(Fid(0:Ni,1:Nj,1:Nk)) ; call Fid%alloc(Ns=Ns)
    allocate(Fjd(1:Ni,0:Nj,1:Nk)) ; call Fjd%alloc(Ns=Ns)
    allocate(Fkd(1:Ni,1:Nj,0:Nk)) ; call Fkd%alloc(Ns=Ns)
    ! computing convective fluxes
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(i,j,k,s)       &
    !$OMP SHARED(s1,gc,gco,gcu,Ni,Nj,Nk,Ns,block,species0,inviscid,Fic,Fjc,Fkc,Fid,Fjd,Fkd)
#ifndef NULi
    ! i direction
    !$OMP SINGLE
    gcu = min(gco,gc(1),gc(2))
    !$OMP END SINGLE
    !$OMP DO
    do k=1,Nk
      do j=1,Nj
        call fluxes_convective(gc       = gcu,                        &
                               N        = Ni,                         &
                               species0 = species0,                   &
                               F        = block%Fi(0-gcu:Ni+gcu,j,k), &
                               C        = block%C (1-gcu:Ni+gcu,j,k), &
                               Fl       = Fic     (    0:Ni    ,j,k))
      enddo
    enddo
#endif
#ifndef NULj
    ! j direction
    !$OMP SINGLE
    gcu = min(gco,gc(3),gc(4))
    !$OMP END SINGLE
    !$OMP DO
    do k=1,Nk
      do i=1,Ni
        call fluxes_convective(gc       = gcu,                        &
                               N        = Nj,                         &
                               species0 = species0,                   &
                               F        = block%Fj(i,0-gcu:Nj+gcu,k), &
                               C        = block%C (i,1-gcu:Nj+gcu,k), &
                               Fl       = Fjc     (i,    0:Nj    ,k))
      enddo
    enddo
#endif
#ifndef NULk
    ! k direction
    !$OMP SINGLE
    gcu = min(gco,gc(5),gc(6))
    !$OMP END SINGLE
    !$OMP DO
    do j=1,Nj
      do i=1,Ni
        call fluxes_convective(gc       = gcu,                        &
                               N        = Nk,                         &
                               species0 = species0,                   &
                               F        = block%Fk(i,j,0-gcu:Nk+gcu), &
                               C        = block%C (i,j,1-gcu:Nk+gcu), &
                               Fl       = Fkc     (i,j,    0:Nk    ))
      enddo
    enddo
#endif
    ! computing diffusive fluxes
    if (.not.inviscid) then
#ifndef NULi
      ! i direction
      !$OMP DO
      do k=1,Nk
        do j=1,Nj
          do i=0,Ni
            !call fluxes_diffusive(block=block,i=i,j=j,k=k,dir='i',F=Fid(i,j,k))
          enddo
        enddo
      enddo
#endif
#ifndef NULj
      ! j direction
      !$OMP DO
      do k=1,Nk
        do j=0,Nj
          do i=1,Ni
            !call fluxes_diffusive(block=block,i=i,j=j,k=k,dir='j',F=Fjd(i,j,k))
          enddo
        enddo
      enddo
#endif
#ifndef NULk
      ! k direction
      !$OMP DO
      do k=0,Nk
        do j=1,Nj
          do i=1,Ni
            !call fluxes_diffusive(block=block,i=i,j=j,k=k,dir='k',F=Fkd(i,j,k))
          enddo
        enddo
      enddo
#endif
    endif
    ! computing the residuals
    !$OMP DO
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          ! overloaded operators form: not efficient!
          !block%KS(i,j,k,s1) =                                                                                     &
          !  (block%Si(i-1,j,  k  )*(Fic(i-1,j,  k  )+Fid(i-1,j,  k  )) - block%Si(i,j,k)*(Fic(i,j,k)+Fid(i,j,k)) + &
          !   block%Sj(i,  j-1,k  )*(Fjc(i,  j-1,k  )+Fjd(i,  j-1,k  )) - block%Sj(i,j,k)*(Fjc(i,j,k)+Fjd(i,j,k)) + &
          !   block%Sk(i,  j,  k-1)*(Fkc(i,  j,  k-1)+Fkd(i,  j,  k-1)) - block%Sk(i,j,k)*(Fkc(i,j,k)+Fkd(i,j,k))   &
          !  )/block%V(i,j,k)
          do s=1,Ns
            block%C(i,j,k)%KS(s1)%rs(s) =                                            &
            (block%Fi(i-1,j,  k  )%S*(Fic(i-1,j,  k  )%rs(s)+Fid(i-1,j,  k  )%rs(s))-&
             block%Fi(i,  j,  k  )%S*(Fic(i,  j,  k  )%rs(s)+Fid(i,  j,  k  )%rs(s))+&
             block%Fj(i,  j-1,k  )%S*(Fjc(i,  j-1,k  )%rs(s)+Fjd(i,  j-1,k  )%rs(s))-&
             block%Fj(i,  j,  k  )%S*(Fjc(i,  j,  k  )%rs(s)+Fjd(i,  j,  k  )%rs(s))+&
             block%Fk(i,  j,  k-1)%S*(Fkc(i,  j,  k-1)%rs(s)+Fkd(i,  j,  k-1)%rs(s))-&
             block%Fk(i,  j,  k  )%S*(Fkc(i,  j,  k  )%rs(s)+Fkd(i,  j,  k  )%rs(s)) &
            )/block%C(i,j,k)%V
          enddo
          block%C(i,j,k)%KS(s1)%rv =                                                                                         &
          (block%Fi(i-1,j,  k  )%S*(Fic(i-1,j,  k  )%rv+Fid(i-1,j,  k  )%rv)-block%Fi(i,j,k)%S*(Fic(i,j,k)%rv+Fid(i,j,k)%rv)+&
           block%Fj(i,  j-1,k  )%S*(Fjc(i,  j-1,k  )%rv+Fjd(i,  j-1,k  )%rv)-block%Fj(i,j,k)%S*(Fjc(i,j,k)%rv+Fjd(i,j,k)%rv)+&
           block%Fk(i,  j,  k-1)%S*(Fkc(i,  j,  k-1)%rv+Fkd(i,  j,  k-1)%rv)-block%Fk(i,j,k)%S*(Fkc(i,j,k)%rv+Fkd(i,j,k)%rv) &
          )/block%C(i,j,k)%V
          block%C(i,j,k)%KS(s1)%re =                                                                                     &
          (block%Fi(i-1,j,  k  )%S*(Fic(i-1,j,  k  )%re+Fid(i-1,j,  k  )%re)-block%Fi(i,j,k)%S*(Fic(i,j,k)%re+Fid(i,j,k)%re)+&
           block%Fj(i,  j-1,k  )%S*(Fjc(i,  j-1,k  )%re+Fjd(i,  j-1,k  )%re)-block%Fj(i,j,k)%S*(Fjc(i,j,k)%re+Fjd(i,j,k)%re)+&
           block%Fk(i,  j,  k-1)%S*(Fkc(i,  j,  k-1)%re+Fkd(i,  j,  k-1)%re)-block%Fk(i,j,k)%S*(Fkc(i,j,k)%re+Fkd(i,j,k)%re) &
          )/block%C(i,j,k)%V
        enddo
      enddo
    enddo
#ifdef NULi
    !$OMP DO
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          block%C(i,j,k)%KS(s1)%rv%x = 0._R8P
        enddo
      enddo
    enddo
#endif
#ifdef NULj
    !$OMP DO
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          block%C(i,j,k)%KS(s1)%rv%y = 0._R8P
        enddo
      enddo
    enddo
#endif
#ifdef NULk
    !$OMP DO
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          block%C(i,j,k)%KS(s1)%rv%z = 0._R8P
        enddo
      enddo
    enddo
#endif
    !$OMP END PARALLEL
    deallocate(Fic,Fjc,Fkc,Fid,Fjd,Fkd)
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine residuals

  !> @brief Procedure for summing Runge-Kutta stages for updating primitive variables (block%P).
  elemental subroutine rk_stages_sum(block,s1,species0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block    !< Block-level data.
  integer(I1P),       intent(IN)::    s1       !< Current Runge-Kutta stage.
  type(Type_Species), intent(IN)::    species0 !< Initial species.
  integer(I4P)::                      i,j,k    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(s1,block,species0)
  !$OMP DO
  do k=1_I4P,block%dims%Nk
    do j=1_I4P,block%dims%Nj
      do i=1_I4P,block%dims%Ni
        call rk_stages_sum_backend(s1=s1,Ns=block%dims%Ns,Nc=block%dims%Nc, &
                                   species0=species0,                       &
                                   Dt=block%C(i,j,k)%Dt,U=block%C(i,j,k)%U,KS=block%C(i,j,k)%KS(1:s1-1),P=block%C(i,j,k)%P)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure subroutine rk_stages_sum_backend(s1,Ns,Nc,species0,Dt,U,KS,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I1P),            intent(IN)::    s1               !< Current Runge-Kutta stage.
    integer(I_P),            intent(IN)::    Ns               !< Number of species.
    integer(I_P),            intent(IN)::    Nc               !< Number of conservative variables.
    type(Type_Species),      intent(IN)::    species0         !< Initial species.
    real(R_P),               intent(IN)::    Dt               !< Current time step.
    type(Type_Conservative), intent(IN)::    U                !< Conservative variables.
    type(Type_Conservative), intent(IN)::    KS(1:)           !< Runge-Kutta stages.
    type(Type_Primitive),    intent(INOUT):: P                !< Primitive variables.
    real(R_P)::                              KSa(1:Nc,1:s1-1) !< Dummy Runge-Kutta stages.
    type(Type_Conservative)::                Ud               !< Dummy conservative variables.
    integer(I1P)::                           ss               !< Counter.
    integer(I4P)::                           s                !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    call Ud%alloc(Ns=Ns)
    do ss=1_I1P,s1-1_I1P
      KSa(:,ss) = KS(ss)%cons2array()
    enddo
    do s=1_I4P,Ns
      Ud%rs(s) = rk_stage(s1=s1,Dt=Dt,Un=U%rs(s),KS=KSa(s,1:s1-1_I1P))
    enddo
    Ud%rv%x = rk_stage(s1=s1,Dt=Dt,Un=U%rv%x,KS=KSa(Ns+1_I4P,1:s1-1_I1P))
    Ud%rv%y = rk_stage(s1=s1,Dt=Dt,Un=U%rv%y,KS=KSa(Ns+2_I4P,1:s1-1_I1P))
    Ud%rv%z = rk_stage(s1=s1,Dt=Dt,Un=U%rv%z,KS=KSa(Ns+3_I4P,1:s1-1_I1P))
    Ud%re   = rk_stage(s1=s1,Dt=Dt,Un=U%re  ,KS=KSa(Ns+4_I4P,1:s1-1_I1P))
    call cons2prim(cons=Ud,prim=P,species0=species0)
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine rk_stages_sum_backend
  endsubroutine rk_stages_sum

  !> @brief Procedure for computing Runge-Kutta one time step integration.
  pure subroutine rk_time_integration(block,RU)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: block               !< Block-level data.
  real(R8P),          intent(OUT)::   RU(1:)              !< NormL2 of residuals of conservative variables.
  real(R8P)                           R(1:size(RU,dim=1)) !< Residuals of conservative variables.
  integer(I4P)::                      i,j,k               !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  RU = 0._R8P
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k,R)       &
  !$OMP SHARED(block)          &
  !$OMP REDUCTION(+: RU)
  !$OMP DO
  do k=1,block%dims%Nk
    do j=1,block%dims%Nj
      do i=1,block%dims%Ni
        call rk_time_integration_backend(Ns=block%dims%Ns,rk_ord=block%dims%Nrk,&
                                         Dt=block%C(i,j,k)%Dt,KS=block%C(i,j,k)%KS,U=block%C(i,j,k)%U,R=R)
        RU = RU + (R*R)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  RU = sqrt(RU)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure subroutine rk_time_integration_backend(Ns,rk_ord,Dt,KS,U,R)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I4P),            intent(IN)::    Ns                            !< Number of species.
    integer(I1P),            intent(IN)::    rk_ord                        !< Number of Runge-Kutta stages.
    real(R8P),               intent(IN)::    Dt                            !< Current time step.
    type(Type_Conservative), intent(IN)::    KS(1:)                        !< Runge-Kutta stages.
    type(Type_Conservative), intent(INOUT):: U                             !< Conservative variables.
    real(R8P),               intent(OUT)::   R(1:)                         !< Residuals of conservative variables.
    real(R8P)::                              KSa(1:size(R,dim=1),1:rk_ord) !< Dummy Runge-Kutta stages.
    integer(I1P)::                           s1                            !< Counter.
    integer(I4P)::                           s                             !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    do s1=1_I1P,rk_ord
      KSa(:,s1) = KS(s1)%cons2array()
    enddo
    do s=1_I4P,Ns
      call rk_time_integ(Dt=Dt,KS=KSa(s,:),R=R(s),Unp1=U%rs(s))
    enddo
    call rk_time_integ(Dt=Dt,KS=KSa(Ns+1_I4P,:),R=R(Ns+1_I4P),Unp1=U%rv%x)
    call rk_time_integ(Dt=Dt,KS=KSa(Ns+2_I4P,:),R=R(Ns+2_I4P),Unp1=U%rv%y)
    call rk_time_integ(Dt=Dt,KS=KSa(Ns+3_I4P,:),R=R(Ns+3_I4P),Unp1=U%rv%z)
    call rk_time_integ(Dt=Dt,KS=KSa(Ns+4_I4P,:),R=R(Ns+4_I4P),Unp1=U%re  )
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine rk_time_integration_backend
  endsubroutine rk_time_integration

  !> @brief Procedure for loading block.
  subroutine load_block(block,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),     intent(INOUT):: block   !< Block data.
  integer(I8P), optional, intent(IN)::    pos     !< Position specifier.
  integer(I4P), optional, intent(OUT)::   iostat  !< IO error.
  character(*), optional, intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),           intent(IN)::    unit    !< Logic unit.
  integer(I4P)::                          iostatd !< IO error.
  character(500)::                        iomsgd  !< Temporary variable for IO error message.
  integer(I4P)::                          i,j,k   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    call block%dims%load(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
  else
    call block%dims%load(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  endif
  call block%alloc ; call block%alloc(members=.true.)
  associate(gc => block%dims%gc(1:6),Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
    do k=0-gc(5),Nk+gc(6) ; do j=0-gc(3),Nj+gc(4) ; do i=0-gc(1),Ni+gc(2)
      call block%node(i,j,k)%load(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
    do k=1-gc(5),Nk+gc(6) ; do j=1-gc(3),Nj+gc(4) ; do i=0-gc(1),Ni+gc(2)
      call block%Fi(i,j,k)%load(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
    do k=1-gc(5),Nk+gc(6) ; do j=0-gc(3),Nj+gc(4) ; do i=1-gc(1),Ni+gc(2)
      call block%Fj(i,j,k)%load(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
    do k=0-gc(5),Nk+gc(6) ; do j=1-gc(3),Nj+gc(4) ; do i=1-gc(1),Ni+gc(2)
      call block%Fk(i,j,k)%load(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
    do k=1-gc(5),Nk+gc(6) ; do j=1-gc(3),Nj+gc(4) ; do i=1-gc(1),Ni+gc(2)
      call block%C(i,j,k)%load(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
  endassociate
  call block%BC%load(unit=unit,iostat=iostatd,iomsg=iomsgd)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_block

  !> @brief Procedure for saving block.
  subroutine save_block(block,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),     intent(IN)::  block   !< Block data.
  integer(I8P), optional, intent(IN)::  pos     !< Position specifier.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  integer(I4P)::                        i,j,k   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    call block%dims%save(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)
  else
    call block%dims%save(unit=unit,        iostat=iostatd,iomsg=iomsgd)
  endif
  associate(gc => block%dims%gc(1:6),Ni => block%dims%Ni,Nj => block%dims%Nj,Nk => block%dims%Nk)
    do k=0-gc(5),Nk+gc(6) ; do j=0-gc(3),Nj+gc(4) ; do i=0-gc(1),Ni+gc(2)
      call block%node(i,j,k)%save(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
    do k=1-gc(5),Nk+gc(6) ; do j=1-gc(3),Nj+gc(4) ; do i=0-gc(1),Ni+gc(2)
      call block%Fi(i,j,k)%save(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
    do k=1-gc(5),Nk+gc(6) ; do j=0-gc(3),Nj+gc(4) ; do i=1-gc(1),Ni+gc(2)
      call block%Fj(i,j,k)%save(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
    do k=0-gc(5),Nk+gc(6) ; do j=1-gc(3),Nj+gc(4) ; do i=1-gc(1),Ni+gc(2)
      call block%Fk(i,j,k)%save(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
    do k=1-gc(5),Nk+gc(6) ; do j=1-gc(3),Nj+gc(4) ; do i=1-gc(1),Ni+gc(2)
      call block%C(i,j,k)%save(unit=unit,iostat=iostatd,iomsg=iomsgd)
    enddo ; enddo ; enddo
  endassociate
  call block%BC%save(unit=unit,iostat=iostatd,iomsg=iomsgd)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_block

  !> @brief Procedure for printing mesh infos with a pretty format.
  subroutine print_info_mesh(block,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),     intent(IN)::  block       !< Block-level data.
  character(*), optional, intent(IN)::  pref        !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat      !< IO error.
  character(*), optional, intent(OUT):: iomsg       !< IO error message.
  integer(I4P),           intent(IN)::  unit        !< Logic unit.
  character(len=:), allocatable::       prefd       !< Prefixing string.
  integer(I4P)::                        iostatd     !< IO error.
  character(500)::                      iomsgd      !< Temporary variable for IO error message.
  type(Type_Vector), allocatable::      NFi(:,:,:)  !< |
  type(Type_Vector), allocatable::      NFj(:,:,:)  !< |
  type(Type_Vector), allocatable::      NFk(:,:,:)  !< |
  real(R_P),         allocatable::      Si (:,:,:)  !< | Dummy variables for printing only internal cells info.
  real(R8P),         allocatable::      Sj (:,:,:)  !< |
  real(R8P),         allocatable::      Sk (:,:,:)  !< |
  real(R8P),         allocatable::      V  (:,:,:)  !< |
  character(len=:),  allocatable::      vmax,vmin   !< String for printing max and min of variables.
  character(len=:),  allocatable::      xmm,ymm,zmm !< String for printing max and min of variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Block dimensions:'
  call block%dims%print(unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd//'  ')
  associate(Ni => block%dims%Ni, Nj => block%dims%Nj, Nk => block%dims%Nk)
    NFi = block%Fi(0:Ni,1:Nj,1:Nk)%N
    NFj = block%Fj(1:Ni,0:Nj,1:Nk)%N
    NFk = block%Fk(1:Ni,1:Nj,0:Nk)%N
    Si  = block%Fi(0:Ni,1:Nj,1:Nk)%S
    Sj  = block%Fj(1:Ni,0:Nj,1:Nk)%S
    Sk  = block%Fk(1:Ni,1:Nj,0:Nk)%S
    V   = block%C (1:Ni,1:Nj,1:Nk)%V
  endassociate
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Block normals:'
  xmm=trim(str(n=maxval(NFi(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFi(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFi(:,:,:)%z)))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   NFi-x-max '//xmm//' NFi-y-max '//ymm//' NFi-z-max '//zmm
  xmm=trim(str(n=minval(NFi(:,:,:)%x))) ; ymm=trim(str(n=minval(NFi(:,:,:)%y))) ; zmm=trim(str(n=minval(NFi(:,:,:)%z)))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   NFi-x-min '//xmm//' NFi-y-min '//ymm//' NFi-z-min '//zmm
  xmm=trim(str(n=maxval(NFj(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFj(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFj(:,:,:)%z)))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   NFj-x-max '//xmm//' NFj-y-max '//ymm//' NFj-z-max '//zmm
  xmm=trim(str(n=minval(NFj(:,:,:)%x))) ; ymm=trim(str(n=minval(NFj(:,:,:)%y))) ; zmm=trim(str(n=minval(NFj(:,:,:)%z)))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   NFj-x-min '//xmm//' NFj-y-min '//ymm//' NFj-z-min '//zmm
  xmm=trim(str(n=maxval(NFk(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFk(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFk(:,:,:)%z)))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   NFk-x-max '//xmm//' NFk-y-max '//ymm//' NFk-z-max '//zmm
  xmm=trim(str(n=minval(NFk(:,:,:)%x))) ; ymm=trim(str(n=minval(NFk(:,:,:)%y))) ; zmm=trim(str(n=minval(NFk(:,:,:)%z)))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   NFk-x-min '//xmm//' NFk-y-min '//ymm//' NFk-z-min '//zmm
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Block faces area:'
  vmax=trim(str(n=maxval(Si(:,:,:)))) ; vmin=trim(str(n=minval(Si(:,:,:))))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Si-max '//vmax//' Si-min '//vmin
  vmax=trim(str(n=maxval(Sj(:,:,:)))) ; vmin=trim(str(n=minval(Sj(:,:,:))))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Sj-max '//vmax//' Sj-min '//vmin
  vmax=trim(str(n=maxval(Sk(:,:,:)))) ; vmin=trim(str(n=minval(Sk(:,:,:))))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   Sk-max '//vmax//' Sk-min '//vmin
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Block volumes:'
  vmax=trim(str(n=maxval(V(:,:,:)))) ; vmin=trim(str(n=minval(V(:,:,:))))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'   V-max '//vmax//' V-min '//vmin
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_info_mesh

  !> @brief Procedure for printing mesh infos with a pretty format.
  subroutine print_info_bc(block,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),     intent(IN)::  block       !< Block-level data.
  character(*), optional, intent(IN)::  pref        !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat      !< IO error.
  character(*), optional, intent(OUT):: iomsg       !< IO error message.
  integer(I4P),           intent(IN)::  unit        !< Logic unit.
  character(len=:), allocatable::       prefd       !< Prefixing string.
  integer(I4P)::                        iostatd     !< IO error.
  character(500)::                      iomsgd      !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  call block%BC%print(unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_info_bc

  !> @brief Procedure for printing fluid infos with a pretty format.
  subroutine print_info_fluid(block,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),     intent(IN)::  block          !< Block-level data.
  character(*), optional, intent(IN)::  pref           !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat         !< IO error.
  character(*), optional, intent(OUT):: iomsg          !< IO error message.
  integer(I4P),           intent(IN)::  unit           !< Logic unit.
  character(len=:), allocatable::       prefd          !< Prefixing string.
  integer(I4P)::                        iostatd        !< IO error.
  character(500)::                      iomsgd         !< Temporary variable for IO error message.
  character(len=:), allocatable::       vmax,vmin      !< String for printing max and min of variables.
  real(R8P),        allocatable::       dummy(:,:,:,:) !< Dummy variables for printing only internal cells info.
  integer(I4P)::                        i,j,k,s        !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  associate(Ni => block%dims%Ni, Nj => block%dims%Nj, Nk => block%dims%Nk, Nc => block%dims%Nc, Np => block%dims%Np)
    if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Nc,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R8P
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          dummy(:,i,j,k) = block%C(i,j,k)%U%cons2array()
        enddo
      enddo
    enddo
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'  Conservative variables'
    do s=1,Nc
      vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
      write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'    U('//trim(str(.true.,s))//'): max '//vmax//' min '//vmin
    enddo
    if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Np,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R8P
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          dummy(:,i,j,k) = block%C(i,j,k)%P%prim2array()
        enddo
      enddo
    enddo
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'  Primitive variables'
    do s=1,Np
      vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
      write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'    P('//trim(str(.true.,s))//'): max '//vmax//' min '//vmin
    enddo
    if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          dummy(1,i,j,k) = block%C(i,j,k)%Dt
        enddo
      enddo
    enddo
    vmax=trim(str(n=maxval(dummy(1,:,:,:)))) ; vmin=trim(str(n=minval(dummy(1,:,:,:))))
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//'  Time steps Dt: max '//vmax//' min '//vmin
  endassociate
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_info_fluid

  !> @brief Procedure for printing block with a pretty format.
  subroutine print_block(block,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock),     intent(IN)::  block          !< Block-level data.
  character(*), optional, intent(IN)::  pref           !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat         !< IO error.
  character(*), optional, intent(OUT):: iomsg          !< IO error message.
  integer(I4P),           intent(IN)::  unit           !< Logic unit.
  character(len=:), allocatable::       prefd          !< Prefixing string.
  integer(I4P)::                        iostatd        !< IO error.
  character(500)::                      iomsgd         !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Mesh infos:'
  call block%print_info_mesh(unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd//'  ')
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' BC infos:'
  call block%print_info_bc(unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd//'  ')
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Fluid infos:'
  call block%print_info_fluid(unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd//'  ')
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_block

  !> @brief Procedure for generating a mirrored block.
  elemental function mirror(block,mirrorX,mirrorY,mirrorZ) result(blockmir)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(IN):: block    !< Block-level data.
  logical, optional,  intent(IN):: mirrorX  !< Flag for activating mirroring on X axis versor.
  logical, optional,  intent(IN):: mirrorY  !< Flag for activating mirroring on Y axis versor.
  logical, optional,  intent(IN):: mirrorZ  !< Flag for activating mirroring on Z axis versor.
  type(Type_SBlock)::              blockmir !< Mirrored block.
  integer(I4P)::                   i,j,k    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  blockmir = block
  ! mirroring mesh
  if (present(mirrorX)) then
    if (mirrorX) then
    endif
  endif
  if (present(mirrorY)) then
    if (mirrorY) then
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP PRIVATE(i,j,k)         &
      !$OMP SHARED(blockmir)
      !$OMP DO
      do k=1-blockmir%dims%gc(5),blockmir%dims%Nk+blockmir%dims%gc(6)
        do j=1-blockmir%dims%gc(3),blockmir%dims%Nj+blockmir%dims%gc(4)
          do i=1-blockmir%dims%gc(1),blockmir%dims%Ni+blockmir%dims%gc(2)
            blockmir%node(i,j,k)%y = -blockmir%node(i,j,k)%y
          enddo
        enddo
      enddo
      !$OMP END PARALLEL
    endif
  endif
  if (present(mirrorZ)) then
    if (mirrorZ) then
    endif
  endif
  call blockmir%metrics
  call blockmir%metrics_correction
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction mirror

  ! Assignment (=)
  !> @brief Procedure for assignment between two blocks variables.
  elemental subroutine assign_blk(blk1,blk2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_SBlock), intent(INOUT):: blk1
  type(Type_SBlock),  intent(IN)::    blk2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  blk1%exts = blk2%exts
  blk1%BC   = blk2%BC
  blk1%IC   = blk2%IC
  blk1%dims = blk2%dims
  call blk1%alloc
  call blk1%alloc(members=.true.)
  if (allocated(blk2%node)) blk1%node = blk2%node
  if (allocated(blk2%Fi  )) blk1%Fi   = blk2%Fi
  if (allocated(blk2%Fj  )) blk1%Fj   = blk2%Fj
  if (allocated(blk2%Fk  )) blk1%Fk   = blk2%Fk
  if (allocated(blk2%C   )) blk1%C    = blk2%C
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_blk
  !> @}
endmodule Data_Type_SBlock
