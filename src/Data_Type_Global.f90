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
USE IR_Precision                                          ! Integers and reals precision definition.
USE Data_Type_BC                                          ! Definition of Type_BC.
USE Data_Type_Adimensional,    only: Type_Adimensional    ! Definition of Type_Adimensional.
USE Data_Type_BC_in1,          only: Type_BC_in1          ! Definition of Type_BC_in1.
USE Data_Type_Cell,            only: Type_Cell            ! Definition of Type_Cell.
USE Data_Type_Cell,            only: cells_bc_set_ext     ! Procedure for setting etrapolation BC of cells array.
USE Data_Type_Cell,            only: cells_bc_set_ref     ! Procedure for setting reflective BC of cells array.
USE Data_Type_Cell,            only: cells_bc_set_per     ! Procedure for setting periodic BC of cells array.
USE Data_Type_CompiledCode,    only: Type_CompiledCode    ! Definition of Type_CompiledCode.
USE Data_Type_Face,            only: Type_Face            ! Definition of Type_Face.
USE Data_Type_File_Profile,    only: Type_File_Profile    ! Definition of Type_File_Profile.
USE Data_Type_Mesh_Dimensions, only: Type_Mesh_Dimensions ! Definition of Type_Mesh_Dimensions.
USE Data_Type_OS,              only: Type_OS              ! Definition of Type_OS.
USE Data_Type_Parallel,        only: Type_Parallel        ! Definition of Type_Parallel.
USE Data_Type_SBlock,          only: Type_SBlock          ! Definition of Type_SBlock.
USE Data_Type_Space_Step,      only: Type_Space_Step      ! Definition of Type_Space_Step.
USE Data_Type_Species,         only: Type_Species         ! Definition of Type_Species.
USE Data_Type_Time_Step,       only: Type_Time_Step       ! Definition of Type_Time_Step.
USE Data_Type_Vector,          only: Type_Vector          ! Definition of Type_Vector.
#ifdef _MPI
USE Lib_Parallel,              only: prim_sendrecv        ! Library for send/receive data for parallel (MPI) operations.
USE MPI                                                   ! MPI runtime library.
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
  type(Type_CompiledCode)::        cco        !< Compiled code used options.
  type(Type_OS)::                  OS         !< Running architecture.
  type(Type_Parallel)::            parallel   !< Parallelization data.
  type(Type_Space_Step)::          space_step !< Space step data.
  type(Type_Time_Step)::           time_step  !< Time step data.
  type(Type_Species)::             species0   !< Initial species.
  type(Type_Adimensional)::        adim       !< Non-dimensionalization data.
  type(Type_BC_in1)::              bc_in1     !< Inflow 1 boundary conditions.
  type(Type_Mesh_Dimensions)::     mesh_dims  !< Mesh global dimensions.
  type(Type_SBlock), allocatable:: block(:,:) !< Block-level data [1:mesh_dims%Nb,1:mesh_dims%Nl].
  contains
    procedure:: free                 ! Procedure for freeing memory.
    procedure:: alloc                ! Procedure for allocating memory.
    procedure:: set_Ns_from_species0 ! Procedure for setting the number of species accordingly species0.
    procedure:: set_Nrk_from_rk_ord  ! Procedure for setting the number of Runge-Kutta stages from time order accuracy, rk_ord.
    procedure:: blocks_dims_update   ! Procedure for updating blocks dimensions form mesh ones.
    procedure:: solve_grl            ! Procedure for solving the conservation equations for one grid level.
    procedure:: boundary_conditions  ! Procedure for imposing the boundary conditions for one grid level.
    final::     finalize             ! Procedure for freeing dynamic memory when finalizing.
endtype Type_Global
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_GlobalPrivateProcedure
  !> @{
  !> @brief Procedure for freeing dynamic data of Type_Global variables.
  elemental subroutine free(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call global%bc_in1%free
  call global%mesh_dims%free
  if (allocated(global%block)) then
    call global%block%free ; deallocate(global%block)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free

  !> @brief Procedure for allocating dynamic data of Type_Global variables.
  elemental subroutine alloc(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global data.
  integer(I4P)::                      b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(global%block)) then
    call global%block%free ; deallocate(global%block)
  endif
  allocate(global%block(1:global%mesh_dims%Nb,1:global%mesh_dims%Nl))
  do l=1,global%mesh_dims%Nl
    do b=1,global%mesh_dims%Nb
      global%block(b,l)%dims = global%mesh_dims%block_dims(global%parallel%blockmap(b),l)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc

  !> @brief Subroutine for freeing dynamic memory when finalizing.
  elemental subroutine finalize(global)
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
  elemental subroutine set_Ns_from_species0(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global data.
  integer(I4P)::                      b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,global%mesh_dims%Nl
    do b=1,global%mesh_dims%Nb_tot
      global%mesh_dims%block_dims(b,l)%Np = global%species0%Np
      global%mesh_dims%block_dims(b,l)%Nc = global%species0%Nc
      global%mesh_dims%block_dims(b,l)%Ns = global%species0%Ns
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_Ns_from_species0

  !> @brief Procedure for setting the number of Runge-Kutta stages from time order accuracy, rk_ord.
  elemental subroutine set_Nrk_from_rk_ord(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global data.
  integer(I4P)::                      b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do l=1,global%mesh_dims%Nl
    do b=1,global%mesh_dims%Nb_tot
      global%mesh_dims%block_dims(b,l)%Nrk = global%time_step%rk_ord
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_Nrk_from_rk_ord

  !> @brief Procedure for updating blocks dimensions form mesh ones.
  elemental subroutine blocks_dims_update(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global data.
  integer(I4P)::                      b,l    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(mesh_dims => global%mesh_dims, blockmap => global%parallel%blockmap)
    do l=1,global%mesh_dims%Nl
      do b=1,global%mesh_dims%Nb
        global%block(b,l)%dims = mesh_dims%block_dims(blockmap(b),l)
      enddo
    enddo
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine blocks_dims_update

  !> @brief Procedure for imposing bc over a generic direction.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) arrays are declared as assumed-shape ones.
  pure subroutine set_bc(global,l,gc,ic,N,F,C)
  !-------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN)::    global       !< Global-level data.
  integer(I4P),      intent(IN)::    l            !< Current grid level.
  integer(I1P),      intent(IN)::    gc(1: )      !< Number of ghost cells [1:2].
  integer(I4P),      intent(IN)::    ic           !< Number of internal cells used for extrapolation (1 or gc).
  integer(I4P),      intent(IN)::    N            !< Number of internal cells.
  type(Type_Face),   intent(IN)::    F(0-gc(1):)  !< Faces data [0-gc(1):N+gc(2)].
  type(Type_Cell),   intent(INOUT):: C(1-gc(1):)  !< Cells data [1-gc(1):N+gc(2)].
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  associate(block=>global%block(:,l))
    ! left
    select case(F(0)%BC%tp)
    case(bc_ext)
      call cells_bc_set_ext(gc=gc,ic=ic,N=N,          boundary='l',C=C(1-gc(1):N+gc(2)))
    case(bc_ref)
      call cells_bc_set_ref(gc=gc,ic=ic,N=N,NF=F(0)%N,boundary='l',C=C(1-gc(1):N+gc(2)))
    case(bc_per)
      call cells_bc_set_per(gc=gc,ic=ic,N=N,          boundary='l',C=C(1-gc(1):N+gc(2)))
    case(bc_adj)
      call set_adj(global=global,block=block,gc=gc(1),F=F(1-gc(1):0),C=C(1-gc(1):0))
    case(bc_in1)
      call set_in1(global=global,gc=gc(1),F=F(1-gc(1):0),C=C(1-gc(1):0))
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
      call set_adj(global=global,block=block,gc=gc(2),F=F(N:N+gc(2)-1),C=C(N+1:N+gc(2)))
    case(bc_in1)
      call set_in1(global=global,gc=gc(2),F=F(N:N+gc(2)-1),C=C(N+1:N+gc(2)))
    endselect
  endassociate
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Subroutine for imposing adjacent boundary conditions.
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b bc and \b P are declared as
    !> assumed-shape with only the lower bound defined. Their extentions are: bc [1-gc:0], P [1-gc:0].
    !> @note When this subroutine is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays BC and P must be
    !> properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1] for BC and [N+1:N+gc] for P.
    pure subroutine set_adj(global,block,gc,F,C)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    type(Type_Global), intent(IN)::    global       !< Global-level data.
    type(Type_SBlock), intent(IN)::    block(1:)    !< Block-level data.
    integer(I1P),      intent(IN)::    gc           !< Number of ghost cells.
    type(Type_Face),   intent(IN)::    F(1-gc:)     !< Faces data [1-gc:0].
    type(Type_Cell),   intent(INOUT):: C(1-gc:)     !< Cells data [1-gc:0].
    integer(I_P)::                     il,bg        !< Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
#ifdef _MPI
    ! for multi-processes simulation it is possible that the adjacent cell is in other processes than the actual and in case the
    ! data have been already exchanged by the MPI procedure prim_sendrecv
    if (global%parallel%procmap(F(0)%BC%adj%b)/=global%parallel%myrank) return
#endif
    do il=1-gc,0
      bg = minloc(array=global%parallel%blockmap,dim=1,mask=global%parallel%blockmap==F(il)%BC%adj%b)
      C(il)%P = block(bg)%C(F(il)%BC%adj%i,F(il)%BC%adj%j,F(il)%BC%adj%k)%P
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_adj

    !> @brief Subroutine for imposing inflow 1 boundary conditions.
    !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b bc and \b P are declared as
    !> assumed-shape with only the lower bound defined. Their extentions are: bc [1-gc:0], P [1-gc:0].
    !> @note When this subroutine is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays BC and P must be
    !> properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1] for BC and [N+1:N+gc] for P.
    pure subroutine set_in1(global,gc,F,C)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    type(Type_Global), intent(IN)::    global   !< Global-level data.
    integer(I1P),      intent(IN)::    gc       !< Number of ghost cells.
    type(Type_Face),   intent(IN)::    F(1-gc:) !< Faces data [1-gc:0].
    type(Type_Cell),   intent(INOUT):: C(1-gc:) !< Cells data [1-gc:0].
    integer(I_P)::                     i        !< Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    do i=1-gc,0
      C(i)%P = global%BC_in1%P(F(i)%BC%inf)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_in1
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
  subroutine boundary_conditions(global,l)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global    !< Global-level data.
  integer(I4P),       intent(IN)::    l         !< Current grid level.
  integer(I4P)::                      b,i,j,k   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef _MPI
  ! doing the multi-processes communications if necessary
  do b=1,global%mesh_dims%Nb
    call prim_sendrecv(l=l,parallel=global%parallel,mesh_dims=global%mesh_dims,block=global%block(:,l))
  enddo
#endif
  associate(Nb=>global%mesh_dims%Nb,block=>global%block(:,l))
  do b=1,Nb
    associate(gc=>block(b)%dims%gc,Ni=>block(b)%dims%Ni,Nj=>block(b)%dims%Nj,Nk=>block(b)%dims%Nk)
    !!$OMP PARALLEL DEFAULT(NONE) &
    !!$OMP FIRSTPRIVATE(b)        &
    !!$OMP PRIVATE(i,j,k)         &
    !!$OMP SHARED(Ni,Nj,Nk,gc,global)
    !!$OMP DO
    do k=1,Nk
      do j=1,Nj
        call set_bc(global=global,l=l,gc=gc(1:2),ic=0,N=Ni,&
                    F=block(b)%Fi(0-gc(1):Ni+gc(2),j,k),   &
                    C=block(b)%C( 1-gc(1):Ni+gc(2),j,k))
      enddo
    enddo
    !!$OMP DO
    do k=1,Nk
      do i=1,Ni
        call set_bc(global=global,l=l,gc=gc(3:4),ic=0,N=Nj,&
                    F=block(b)%Fj(i,0-gc(3):Nj+gc(4),k),   &
                    C=block(b)%C( i,1-gc(3):Nj+gc(4),k))
      enddo
    enddo
    !!$OMP DO
    do j=1,Nj
      do i=1,Ni
        call set_bc(global=global,l=l,gc=gc(5:6),ic=0,N=Nk,&
                    F=block(b)%Fk(i,j,0-gc(5):Nk+gc(6)),   &
                    C=block(b)%C( i,j,1-gc(5):Nk+gc(6)))
      enddo
    enddo
    !!$OMP END PARALLEL
    endassociate
  enddo
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine boundary_conditions

  !> @brief Procedure for solving (performing one time step integration) the conservation equations for grid level "l".
  subroutine solve_grl(global,l,prof)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global),      intent(INOUT):: global    !< Global data.
  integer(I4P),            intent(IN)::    l         !< Current grid level.
  type(Type_File_Profile), intent(INOUT):: prof      !< Profiling file.
  real(R8P), allocatable::                 Dtmin(:)  !< Min t step of actual process for each blk.
  real(R8P)::                              DtminL    !< Min t step of actual process over all blks.
  real(R8P)::                              gDtmin    !< Global (all processes/all blks) min t step.
  real(R8P), allocatable::                 RU  (:,:) !< NormL2 of conservartive residuals.
  real(R8P), allocatable::                 mRU (:)   !< Maximum of RU of actual process.
  real(R8P), allocatable::                 gmRU(:)   !< Global (all processes) maximum of RU.
  integer(I4P)::                           b         !< Blocks counter.
  integer(I1P)::                           s1        !< Runge-Kutta stages counters.
  integer(I4P)::                           i,j,k     !< Counters.
#ifdef _MPI
  integer(I4P)::                           err       !< Error trapping flag: 0 no errors, >0 errors.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(Nb=>global%mesh_dims%Nb,block=>global%block(:,l),n=>global%time_step%n,&
            myrank=>global%parallel%myrank,Nthreads=>global%parallel%Nthreads,Nproc=>global%parallel%Nproc)
    ! allocating dynamic arrays
    allocate(Dtmin(1:Nb),RU(1:global%species0%Nc,1:Nb),mRU(1:global%species0%Nc),gmRU(1:global%species0%Nc))
    ! converting conservative variables to primitive ones
#ifdef PROFILING
    call prof%profile(p=2,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    do b=1,Nb
      call block(b)%conservative2primitive(species0=global%species0)
    enddo
#ifdef PROFILING
    call prof%profile(p=2,pstop=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    ! imposing the boundary conditions
#ifdef PROFILING
    call prof%profile(p=3,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    call global%boundary_conditions(l=l)
#ifdef PROFILING
    call prof%profile(p=3,pstop=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    ! updating time varying variables: Dt,Dtmin
    n = n + 1_I8P
#ifdef PROFILING
    call prof%profile(p=4,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    do b=1,Nb
      call block(b)%compute_time(time_step=global%time_step,Dtmin=Dtmin(b))
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
      do b=1,Nb
        block(b)%C%Dt = gDtmin
      enddo
    endif
    ! evaluating the Runge-Kutta stages
    ! Runge-Kutta stages initialization
    do b=1,Nb
      associate(gc=>block(b)%dims%gc,Ni=>block(b)%dims%Ni,Nj=>block(b)%dims%Nj,Nk=>block(b)%dims%Nk)
        do k=1-gc(5),Nk+gc(6)
          do j=1-gc(3),Nj+gc(4)
            do i=1-gc(1),Ni+gc(2)
              block(b)%C(i,j,k)%KS = 0._R8P
            enddo
          enddo
        enddo
      endassociate
    enddo
    do s1=1,global%time_step%rk_ord
      if (s1>1) then
        ! summing the stages up to s1-1: $K_{s1}=R\left(U^n+\sum_{s2=1}^{s1-1}{Dt \cdot rk_{c2}(s1,s2) \cdot K_{s2}}\right)$
        ! updating primitive variables: $P=conservative2primitive(K_{s1})$
        do b=1,Nb
            call block(b)%rk_stages_sum(s1=s1,species0=global%species0)
        enddo
        ! imposing the boundary conditions
        call global%boundary_conditions(l=l)
      endif
      ! computing the s1-th Runge-Kutta stage
#ifdef PROFILING
      call prof%profile(p=5,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
      do b=1,Nb
        call block(b)%residuals(s1=s1,space_step=global%space_step,species0=global%species0)
      enddo
#ifdef PROFILING
      call prof%profile(p=5,pstop=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    enddo
    ! Runge-Kutta time integration
#ifdef PROFILING
    call prof%profile(p=6,pstart=.true.,myrank=myrank,Nthreads=Nthreads,Nproc=Nproc)
#endif
    do b=1,Nb
      call block(b)%rk_time_integration(RU=RU(:,b))
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
  endsubroutine solve_grl
  !> @}
endmodule Data_Type_Global
