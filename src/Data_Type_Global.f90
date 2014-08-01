!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_GlobalDerivedType Data_Type_Global
!> Module definition of Type_Global
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_GlobalInterface Data_Type_Global
!> Module definition of Type_Global
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_GlobalPublicProcedure Data_Type_Global
!> Module definition of Type_Global
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_GlobalPrivateProcedure Data_Type_Global
!> Module definition of Type_Global
!> @}

!> @brief Module Data_Type_Global contains the definition of Type_Global and useful procedures for its handling.
!> Global-level data are referred to those informations of global interest.
module Data_Type_Global
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                            ! Integers and reals precision definition.
USE Data_Type_BC                                            ! Definition of Type_BC.
USE Data_Type_Adimensional,     only: Type_Adimensional     ! Definition of Type_Adimensional.
USE Data_Type_BC_in1,           only: Type_BC_in1           ! Definition of Type_BC_in1.
USE Data_Type_Block_Dimensions, only: Type_Block_Dimensions ! Definition of Type_Block_Dimensions.
USE Data_Type_Cell,             only: Type_Cell             ! Definition of Type_Cell.
USE Data_Type_Cell,             only: cells_bc_set_ext      ! Procedure for setting etrapolation BC of cells array.
USE Data_Type_Cell,             only: cells_bc_set_ref      ! Procedure for setting reflective BC of cells array.
USE Data_Type_Cell,             only: cells_bc_set_per      ! Procedure for setting periodic BC of cells array.
USE Data_Type_CompiledCode,     only: Type_CompiledCode     ! Definition of Type_CompiledCode.
USE Data_Type_Face,             only: Type_Face             ! Definition of Type_Face.
USE Data_Type_File_Profile,     only: Type_File_Profile     ! Definition of Type_File_Profile.
USE Data_Type_OS,               only: Type_OS               ! Definition of Type_OS.
USE Data_Type_Parallel,         only: Type_Parallel         ! Definition of Type_Parallel.
USE Data_Type_SBlock,           only: Type_SBlock           ! Definition of Type_SBlock.
USE Data_Type_Space_Step,       only: Type_Space_Step       ! Definition of Type_Space_Step.
USE Data_Type_Species,          only: Type_Species          ! Definition of Type_Species.
USE Data_Type_Time_Step,        only: Type_Time_Step        ! Definition of Type_Time_Step.
USE Data_Type_Tree,             only: Type_Tree             ! Definition of Type_Tree.
USE Data_Type_Vector,           only: Type_Vector           ! Definition of Type_Vector.
#ifdef _MPI
USE Lib_Parallel,               only: prim_sendrecv         ! Library for send/receive data for parallel (MPI) operations.
USE MPI                                                     ! MPI runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the global-level data.
!> @ingroup Data_Type_GlobalDerivedType
type, public:: Type_Global
  type(Type_CompiledCode):: cco        !< Compiled code used options.
  type(Type_OS)::           OS         !< Running architecture.
  type(Type_Parallel)::     parallel   !< Parallelization data.
  type(Type_Space_Step)::   space_step !< Space step data.
  type(Type_Time_Step)::    time_step  !< Time step data.
  type(Type_Species)::      species0   !< Initial species.
  type(Type_Adimensional):: adim       !< Non-dimensionalization data.
  type(Type_BC_in1)::       bc_in1     !< Inflow 1 boundary conditions.
  type(Type_Tree)::         block_dims !< Mesh global dimensions (blocks of all processes), tree of Type_Block_Dimensions.
  type(Type_Tree)::         block      !< Block-level data (blocks of myrank), tree of Type_SBlock.
  contains
    procedure:: free                 ! Procedure for freeing dynamic memory.
    procedure:: set_Ns_from_species0 ! Procedure for setting the number of species accordingly species0.
    procedure:: set_Nrk_from_rk_ord  ! Procedure for setting the number of Runge-Kutta stages from time order accuracy, rk_ord.
    procedure:: blocks_init          ! Procedure for initializing blocks memory.
    procedure:: blocks_dims_update   ! Procedure for updating blocks dimensions form mesh ones.
    procedure:: solve                ! Procedure for solving the conservation equations for one grid level.
    procedure:: boundary_conditions  ! Procedure for imposing the boundary conditions for one grid level.
    final::     finalize             ! Procedure for freeing dynamic memory when finalizing.
endtype Type_Global
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_GlobalPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic memory.
  subroutine free(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call global%cco%free
  call global%parallel%free
  call global%species0%free
  call global%bc_in1%free
  call global%block_dims%free
  call global%block%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Subroutine for freeing dynamic memory when finalizing.
  subroutine finalize(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(INOUT):: global !< Global data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call global%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> @brief Procedure for setting the number of species accordingly species0.
  subroutine set_Ns_from_species0(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT)::    global  !< Global data.
  type(Type_Block_Dimensions), pointer:: blkdims !< Pointer for scanning global%block_dims tree.
  integer(I8P)::                         ID      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do while(global%block_dims%loopID(ID=ID))
  blkdims => global%block_dims%dat(ID=ID)
    blkdims%Np = global%species0%Np
    blkdims%Nc = global%species0%Nc
    blkdims%Ns = global%species0%Ns
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_Ns_from_species0

  !> @brief Procedure for setting the number of Runge-Kutta stages from time order accuracy, rk_ord.
  elemental subroutine set_Nrk_from_rk_ord(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT)::    global  !< Global data.
  type(Type_Block_Dimensions), pointer:: blkdims !< Pointer for scanning global%block_dims tree.
  integer(I8P)::                         ID      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do while(global%block_dims%loopID(ID=ID))
  blkdims => global%block_dims%dat(ID=ID)
    blkdims%Nrk = global%time_step%rk_ord
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_Nrk_from_rk_ord

  !> @brief Procedure for initializing blocks memory.
  !> @note Note that global%parallel and global%block_dims must be already initialized when invoking this procedure.
  subroutine blocks_init(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT)::    global  !< Global data.
  type(Type_Block_Dimensions), pointer:: blkdims !< Pointer for scanning global%block_dims tree.
  integer(I4P),                pointer:: proc    !< Pointer for scanning global%parallel%BPmap      tree.
  type(Type_SBlock)::                    block   !< Block prototype.
  integer(I8P)::                         ID      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call global%block%free
  call global%block%init(source=global%block_dims)
  do while(global%block_dims%loopID(ID=ID))
    proc => global%parallel%BPmap%dat(ID=ID)
    if (global%parallel%is_process(p=proc)) then
      call block%free
      blkdims => global%block_dims%dat(ID=ID)
      block%dims = blkdims
      call block%alloc
      call global%block%put(ID=ID,d=block)
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine blocks_init

  !> @brief Procedure for updating blocks dimensions form mesh ones.
  elemental subroutine blocks_dims_update(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT)::    global  !< Global data.
  type(Type_SBlock),           pointer:: block   !< Pointer for scanning global%block tree.
  type(Type_Block_Dimensions), pointer:: blkdims !< Pointer for scanning global%block_dims tree.
  integer(I8P)::                         ID      !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do while(global%block%loopID(ID=ID))
    block   => global%block%dat(ID=ID)
    blkdims => global%block_dims%dat(ID=ID)
    block%dims = blkdims
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine blocks_dims_update

  !> @brief Procedure for imposing bc over a generic direction.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) arrays are declared as assumed-shape ones.
  pure subroutine set_bc(global,gc,ic,N,F,C)
  !-------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN)::    global       !< Global-level data.
  integer(I1P),      intent(IN)::    gc(1:)       !< Number of ghost cells [1:2].
  integer(I4P),      intent(IN)::    ic           !< Number of internal cells used for extrapolation (1 or gc).
  integer(I4P),      intent(IN)::    N            !< Number of internal cells.
  type(Type_Face),   intent(IN)::    F(0-gc(1):)  !< Faces data [0-gc(1):N+gc(2)].
  type(Type_Cell),   intent(INOUT):: C(1-gc(1):)  !< Cells data [1-gc(1):N+gc(2)].
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  ! left
  select case(F(0)%BC%tp)
  case(bc_ext)
    call cells_bc_set_ext(gc=gc,ic=ic,N=N,          boundary='l',C=C(1-gc(1):N+gc(2)))
  case(bc_ref)
    call cells_bc_set_ref(gc=gc,ic=ic,N=N,NF=F(0)%N,boundary='l',C=C(1-gc(1):N+gc(2)))
  case(bc_per)
    call cells_bc_set_per(gc=gc,ic=ic,N=N,          boundary='l',C=C(1-gc(1):N+gc(2)))
  case(bc_adj)
    call cells_bc_set_adj(gc=gc(1),F=F(1-gc(1):0),               C=C(1-gc(1):0))
  case(bc_in1)
    call cells_bc_set_in1(gc=gc(1),F=F(1-gc(1):0),               C=C(1-gc(1):0))
  endselect
  ! right
  select case(F(N)%BC%tp)
  case(bc_ext)
    call cells_bc_set_ext(gc=gc,ic=ic,N=N,          boundary='r',C=C(1-gc(1):N+gc(2)))
  case(bc_ref)
    call cells_bc_set_ref(gc=gc,ic=ic,N=N,NF=F(N)%N,boundary='r',C=C(1-gc(1):N+gc(2)))
  case(bc_per)
    call cells_bc_set_per(gc=gc,ic=ic,N=N,          boundary='r',C=C(1-gc(1):N+gc(2)))
  case(bc_adj)
    call cells_bc_set_adj(gc=gc(2),F=F(N:N+gc(2)-1),             C=C(N+1:N+gc(2)))
  case(bc_in1)
    call cells_bc_set_in1(gc=gc(2),F=F(N:N+gc(2)-1),             C=C(N+1:N+gc(2)))
  endselect
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Procedure for imposing adjacent boundary conditions.
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b F and \b C are declared as
    !> assumed-shape with only the lower bound defined. Their extentions are: F [1-gc:0], C [1-gc:0].
    !> When this procedure is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays F and C must be
    !> properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1] for F and [N+1:N+gc] for C.
    pure subroutine cells_bc_set_adj(gc,F,C)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I1P),      intent(IN)::    gc       !< Number of ghost cells.
    type(Type_Face),   intent(IN)::    F(1-gc:) !< Faces data [1-gc:0].
    type(Type_Cell),   intent(INOUT):: C(1-gc:) !< Cells data [1-gc:0].
    type(Type_SBlock), pointer::       block    !< Pointer for scanning global%block tree.
    integer(I4P)::                     i        !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
#ifdef _MPI
    ! for multi-processes simulation it is possible that the adjacent cell is in other processes than the actual and in case the
    ! data have been already exchanged by the MPI procedure prim_sendrecv
    if (global%parallel%procmap(F(0)%BC%adj%b)/=global%parallel%myrank) return
#endif
    do i=1-gc,0
      block => global%block%dat(ID=F(i)%BC%adj%ID)
      C(i)%P = block%C(F(i)%BC%adj%i,F(i)%BC%adj%j,F(i)%BC%adj%k)%P
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine cells_bc_set_adj

    !> @brief Procedure for imposing inflow 1 boundary conditions.
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b F and \b C are declared as
    !> assumed-shape with only the lower bound defined. Their extentions are: F [1-gc:0], C [1-gc:0].
    !> When this procedure is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays F and C must be
    !> properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1] for F and [N+1:N+gc] for C.
    pure subroutine cells_bc_set_in1(gc,F,C)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I1P),    intent(IN)::    gc       !< Number of ghost cells.
    type(Type_Face), intent(IN)::    F(1-gc:) !< Faces data [1-gc:0].
    type(Type_Cell), intent(INOUT):: C(1-gc:) !< Cells data [1-gc:0].
    integer(I_P)::                   i        !< Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    do i=1-gc,0
      C(i)%P = global%BC_in1%P(F(i)%BC%inf)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine cells_bc_set_in1
  endsubroutine set_bc

  !> @brief Procedure for imposing the boundary conditions of blocks of grid level "l" updating the ghost cells.
  !> @note Considering the ghost cell \f$ c^0 \f$ and the inner one \f$ c^1 \f$ (being \f$ c^N \f$ the other bound)
  !> the available boundary conditions are:
  !> - \b REF: reflective boundary condition \f$ \begin{array}{*{20}{c}} P^0 = P^1 \\ P_{\vec v\cdot \vec n}^0 =
  !>           - P_{\vec v\cdot \vec n}^1\end{array} \f$;
  !> - \b EXT: extrapolation boundary condition \f$ P^0 = P^1 \f$;
  !> - \b PER: periodic boundary condition \f$ P^0 = P^N \f$;
  !> - \b ADJ: adjacent (cell) boundary condition \f$ P^0 = P^a \f$ where "a" is the adjacent cell (b,i,j,k indexes must be
  !>           specified);
  !> - \b IN1: supersonic inflow steady boundary condition \f$ P^0 = P^{in1} \f$. \n
  !> where \f$P\f$ are the primitive variables.
  subroutine boundary_conditions(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global-level data.
  type(Type_SBlock), pointer::        block  !< Pointer for scanning global%block tree.
  integer(I8P)::                      ID     !< Counter.
  integer(I4P)::                      i,j,k  !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef _MPI
  ! doing the multi-processes communications if necessary
  do b=1,global%mesh_dims%Nb
    call prim_sendrecv(l=l,parallel=global%parallel,mesh_dims=global%mesh_dims,block=global%block(:,l))
  enddo
#endif
  do while(global%block%loopID(ID=ID))
    block => global%block%dat(ID=ID)
    associate(gc=>block%dims%gc,Ni=>block%dims%Ni,Nj=>block%dims%Nj,Nk=>block%dims%Nk)
    !!$OMP PARALLEL DEFAULT(NONE) &
    !!$OMP PRIVATE(i,j,k)         &
    !!$OMP SHARED(Ni,Nj,Nk,gc,global)
    !!$OMP DO
    do k=1,Nk
      do j=1,Nj
        call set_bc(global=global,gc=gc(1:2),ic=0,N=Ni, &
                    F=block%Fi(0-gc(1):Ni+gc(2),j,k),&
                    C=block%C( 1-gc(1):Ni+gc(2),j,k))
      enddo
    enddo
    !!$OMP DO
    do k=1,Nk
      do i=1,Ni
        call set_bc(global=global,gc=gc(3:4),ic=0,N=Nj, &
                    F=block%Fj(i,0-gc(3):Nj+gc(4),k),&
                    C=block%C( i,1-gc(3):Nj+gc(4),k))
      enddo
    enddo
    !!$OMP DO
    do j=1,Nj
      do i=1,Ni
        call set_bc(global=global,gc=gc(5:6),ic=0,N=Nk, &
                    F=block%Fk(i,j,0-gc(5):Nk+gc(6)),&
                    C=block%C( i,j,1-gc(5):Nk+gc(6)))
      enddo
    enddo
    !!$OMP END PARALLEL
    endassociate
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine boundary_conditions

  !> @brief Procedure for solving (performing one time step integration) the conservation equations.
  subroutine solve(global,prof)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global),      intent(INOUT):: global    !< Global data.
  type(Type_File_Profile), intent(INOUT):: prof      !< Profiling file.
  real(R8P), allocatable::                 Dtmin(:)  !< Min t step of actual process for each blk.
  real(R8P)::                              DtminL    !< Min t step of actual process over all blks.
  real(R8P)::                              gDtmin    !< Global (all processes/all blks) min t step.
  real(R8P), allocatable::                 RU  (:,:) !< NormL2 of conservartive residuals.
  real(R8P), allocatable::                 mRU (:)   !< Maximum of RU of actual process.
  real(R8P), allocatable::                 gmRU(:)   !< Global (all processes) maximum of RU.
  integer(I4P)::                           b,Nb      !< Blocks counters.
  integer(I1P)::                           s1        !< Runge-Kutta stages counters.
  integer(I4P)::                           i,j,k     !< Counters.
  type(Type_SBlock), pointer::             block     !< Pointer for scanning global%block tree.
  integer(I8P)::                           ID        !< Counter.
#ifdef _MPI
  integer(I4P)::                           err       !< Error trapping flag: 0 no errors, >0 errors.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(n=>global%time_step%n,myrank=>global%parallel%myrank,Nthreads=>global%parallel%Nthreads,Nproc=>global%parallel%Nproc)
    ! allocating dynamic arrays
    Nb = global%block%length()
    allocate(Dtmin(1:Nb),RU(1:global%species0%Nc,1:Nb),mRU(1:global%species0%Nc),gmRU(1:global%species0%Nc))
    ! converting conservative variables to primitive ones
#ifdef PROFILING
    call prof%profile(p=2,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    do while(global%block%loopID(ID=ID))
      block => global%block%dat(ID=ID)
      call block%conservative2primitive(species0=global%species0)
    enddo
#ifdef PROFILING
    call prof%profile(p=2,pstop=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    ! imposing the boundary conditions
#ifdef PROFILING
    call prof%profile(p=3,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    call global%boundary_conditions
#ifdef PROFILING
    call prof%profile(p=3,pstop=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    ! updating time varying variables: Dt,Dtmin
    n = n + 1_I8P
#ifdef PROFILING
    call prof%profile(p=4,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    do while(global%block%loopID(ID=ID))
      block => global%block%dat(ID=ID)
      call block%compute_time(time_step=global%time_step,Dtmin=Dtmin(b))
    enddo
#ifdef PROFILING
    call prof%profile(p=4,pstop=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    DtminL = minval(Dtmin)
#ifdef _MPI
    ! for multi-processes simulation all processes must exchange their DtminL for computing the global variables
    call MPI_ALLREDUCE(DtminL,gDtmin,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,err)
#else
    ! for single processes DtminL are already global variables
    gDtmin = DtminL
#endif
    if (global%time_step%unsteady) then
      call global%time_step%update(gDtmin=gDtmin)
      do while(global%block%loopID(ID=ID))
        block => global%block%dat(ID=ID)
        block%C%Dt = gDtmin
      enddo
    endif
    ! evaluating the Runge-Kutta stages
    ! Runge-Kutta stages initialization
    do while(global%block%loopID(ID=ID))
      block => global%block%dat(ID=ID)
      associate(gc=>block%dims%gc,Ni=>block%dims%Ni,Nj=>block%dims%Nj,Nk=>block%dims%Nk)
        do k=1-gc(5),Nk+gc(6)
          do j=1-gc(3),Nj+gc(4)
            do i=1-gc(1),Ni+gc(2)
              block%C(i,j,k)%KS = 0._R8P
            enddo
          enddo
        enddo
      endassociate
    enddo
    do s1=1,global%time_step%rk_ord
      if (s1>1) then
        ! summing the stages up to s1-1: $K_{s1}=R\left(U^n+\sum_{s2=1}^{s1-1}{Dt \cdot rk_{c2}(s1,s2) \cdot K_{s2}}\right)$
        ! updating primitive variables: $P=conservative2primitive(K_{s1})$
        do while(global%block%loopID(ID=ID))
          block => global%block%dat(ID=ID)
          call block%rk_stages_sum(s1=s1,species0=global%species0)
        enddo
        ! imposing the boundary conditions
        call global%boundary_conditions
      endif
      ! computing the s1-th Runge-Kutta stage
#ifdef PROFILING
      call prof%profile(p=5,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
      do while(global%block%loopID(ID=ID))
        block => global%block%dat(ID=ID)
        call block%residuals(s1=s1,space_step=global%space_step,species0=global%species0)
      enddo
#ifdef PROFILING
      call prof%profile(p=5,pstop=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    enddo
    ! Runge-Kutta time integration
#ifdef PROFILING
    call prof%profile(p=6,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    do while(global%block%loopID(ID=ID))
      block => global%block%dat(ID=ID)
      call block%rk_time_integration(RU=RU(:,b))
    enddo
#ifdef PROFILING
    call prof%profile(p=6,pstop=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    ! finding the maximum value of residuals of actual process
    mRU = maxval(RU,dim=2)
#ifdef _MPI
    ! for multi-processes simulation all processes must exchange their mRU for computing the global gmRU
    call MPI_ALLREDUCE(mRU,gmRU,global%species0%Nc,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,err)
#else
    ! for single processes mRU is already global gmRU
    gmRU = mRU
#endif
    ! deallocating dynamic arrays
    deallocate(Dtmin,RU,mRU,gmRU)
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine solve
  !> @}
endmodule Data_Type_Global
