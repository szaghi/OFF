module Lib_Fluidynamic
!-----------------------------------------------------------------------------------------------------------------------------------
! The module Lib_Fluidynamic contains the definition of fluidynamic functions and subroutines.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                 ! Integers and reals precision definition.
USE Data_Type_BC, &                                                              ! Definition of Type_BC.
                                        init_bc=>init, set_bc=>set               ! Function for initializing Type_BC.
USE Data_Type_Conservative, &                                                    ! Definition of Type_Conservative.
                                        init_cons=>init, set_cons=>set           ! Function for initializing Type_Conservative.
USE Data_Type_Globals                                                            ! Definition of Type_Global and Type_Block.
USE Data_Type_Primitive, &                                                       ! Definition of Type_Primitive.
                                        init_prim=>init, set_prim=>set           ! Function for initializing Type_Primitive.
USE Data_Type_Time                                                               ! Definition of Type_Time.
USE Data_Type_Vector, &                                                          ! Definition of Type_Vector.
                                        init_vec=>init, set_vec=>set             ! Function for initializing Type_Vector.
USE Lib_IO_Misc                                                                  ! Procedures for IO and strings operations.
USE Lib_Math,                     only: interpolate1                             ! Function for computing linear interpolation.
USE Lib_Parallel,                 only: blockmap, &                              ! Local/global blocks map.
                                        Psendrecv                                ! Subroutine for send/recive Primitive variables.
USE Lib_Riemann,                  only: &
#ifdef SMSWliu
                                        chk_smooth_liu,        &                 ! Liu's smoothness function.
#elif defined SMSWvanleer || defined SMSWvanalbada || defined SMSWharten
                                        chk_smooth_limiter,    &                 ! Slope limiter-based smoothness function.
#else
                                        chk_smooth_z,          &                 ! Riemann-Solver-like smoothness function.
#endif
#if defined RSHLLCb || defined RSHLLCc || defined RSHLLCp || defined RSHLLCt || defined RSHLLCz
                                        Riem_Solver=>Riem_Solver_HLLC            ! HLLC Riemann solver.
#elif RSEXA
                                        Riem_Solver=>Riem_Solver_Exact_U         ! Exact Riemann solver.
#elif RSPVL
                                        Riem_Solver=>Riem_Solver_PVL             ! PVL Riemann solver.
#elif RSTR
                                        Riem_Solver=>Riem_Solver_TR              ! TR Riemann solver.
#elif RSTS
                                        Riem_Solver=>Riem_Solver_TS              ! TS Riemann solver.
#elif RSAPRS
                                        Riem_Solver=>Riem_Solver_APRS            ! APRS Riemann solver.
#elif RSALFR
                                        Riem_Solver=>Riem_Solver_ALFR            ! ALFR Riemann solver.
#elif defined RSLFp || defined RSLFz
                                        Riem_Solver=>Riem_Solver_LaxFriedrichs   ! Lax-Friedrichs Riemann solver.
#elif RSROE
                                        Riem_Solver=>Riem_Solver_Roe             ! Roe Riemann solver.
#else
                                        Riem_Solver=>Riem_Solver_HLLC            ! HLLC Riemann solver.
#endif
USE Lib_Thermodynamic_Laws_Ideal, only: a                                        ! Function for computing speed of sound.
USE Lib_WENO,                     only: &
#ifdef HYBRIDC
                                        noweno_central, &                        ! Function for comp. central diff. reconctruction.
#elif HYBRID
                                        weno_optimal,   &                        ! Function for comp. optimal WENO reconctruction.
#endif
                                        weno                                     ! Function for computing WENO reconctruction.
#ifdef MPI2
USE Lib_Parallel,                 only: procmap, &                               ! Proc/blocks map.
                                        blockmap                                 ! Local/global blocks map.
USE MPI                                                                          ! MPI runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: primitive2conservative
public:: conservative2primitive
public:: residuals
public:: boundary_conditions
public:: solve_grl
public:: rk_init
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Runge-Kutta coefficients
real(R_P), allocatable:: rk_c1(:)   ! c1 coefficients of Runge-Kutta's integration [1:rk_ord]
real(R_P), allocatable:: rk_c2(:,:) ! c2 coefficients of Runge-Kutta's integration [1:rk_ord,1:rk_ord]

integer(I1P):: flip = 0_I_P ! Flip-Flop flag for restart solution file.
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  pure subroutine prim2cons(prim,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for converting primitive variables to conservative variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive),    intent(IN)::    prim ! Primitive variables.
  type(Type_Conservative), intent(INOUT):: cons ! Conservative variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = prim%r
  cons%rv = prim%d*prim%v
  cons%re = prim%p/(prim%g-1._R_P) + 0.5_R_P*prim%d*sq_norm(prim%v)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim2cons

  pure subroutine cons2prim(cp0,cv0,cons,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for converting conservative variables to primitive variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN)::    cp0(:)                  ! Specific heat at constant p of initial species.
  real(R_P),               intent(IN)::    cv0(:)                  ! Specific heat at constant v of initial species.
  type(Type_Conservative), intent(IN)::    cons                    ! Conservative variables.
  type(Type_Primitive),    intent(INOUT):: prim                    ! Primitive variables.
  real(R_P)::                              c(1:size(cp0(:),dim=1)) ! Species concentration.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r     = cons%rs
  prim%d     = sum(cons%rs(:))
  c          = prim%r/prim%d
  prim%g     = dot_product(c,cp0)/dot_product(c,cv0)
  prim%v     = cons%rv/prim%d
  prim%p     = (prim%g - 1._R_P)*(cons%re - 0.5_R_P*prim%d*(sq_norm(prim%v)))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine cons2prim

  subroutine primitive2conservative(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for converting primitive variables to conservative variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Block),  intent(INOUT):: block  ! Block-level data.
  integer(I_P)::                     i,j,k  ! Space counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(block)
  !$OMP DO
  do k=1,block%mesh%Nk
    do j=1,block%mesh%Nj
      do i=1,block%mesh%Ni
        call prim2cons(prim = block%fluid%P(i,j,k), cons = block%fluid%U(i,j,k))
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------

  endsubroutine primitive2conservative

  subroutine conservative2primitive(global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for converting conservative variables to primitive variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN)::    global ! Global-level data.
  type(Type_Block),  intent(INOUT):: block  ! Block-level data.
  integer(I_P)::                     i,j,k  ! Space counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(global,block)
  !$OMP DO
  do k=1,block%mesh%Nk
    do j=1,block%mesh%Nj
      do i=1,block%mesh%Ni
        call cons2prim(cp0 = global%fluid%cp0, cv0 = global%fluid%cv0, cons = block%fluid%U(i,j,k), prim = block%fluid%P(i,j,k))
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine conservative2primitive

  subroutine compute_time(global,block,Dtmin)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for evaluating the local and global time step value by CFL conditions.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN)::    global                        ! Global-level data.
  type(Type_Block),  intent(INOUT):: block                         ! Block-level data.
  real(R_P),         intent(OUT)::   Dtmin                         ! Minimum Dt.
  real(R_P)::                        vmax                          ! Maximum speed of waves.
  real(R_P)::                        a                             ! Speed of sound.
  real(R_P)::                        vmiL,vmiR,vmjL,vmjR,vmkL,vmkR ! Dummy velocities.
  type(Type_Vector)::                vm                            ! Dummy vectorial velocities.
  integer(I_P)::                     Ni,Nj,Nk,gc(1:6)              ! Temp var for storing block dimensions.
  integer(I_P)::                     i,j,k                         ! Space counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = block%mesh%Ni
  Nj = block%mesh%Nj
  Nk = block%mesh%Nk
  gc = block%mesh%gc
  ! computing the minimum Dt into the inner cells
  !$OMP PARALLEL DEFAULT(NONE)                                  &
  !$OMP PRIVATE(i,j,k,vmax,a,vmiL,vmiR,vmjL,vmjR,vmkL,vmkR,vm)  &
  !$OMP SHARED(Ni,Nj,Nk,global,block)
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        ! computing the local speed of sound
        a = sqrt((block%fluid%P(i,j,k)%g*block%fluid%P(i,j,k)%p)/(block%fluid%P(i,j,k)%d))
        ! evaluating the maximum propagation speed of acoustic segnals multiplied for face area
        ! left i
        vm   = 0.5_R_P*(block%fluid%P(i-1,j,k)%v+block%fluid%P(i,j,k)%v)
        vmiL = (vm.dot.block%mesh%NFi(i-1,j,k))*block%mesh%Si(i-1,j,k)
        vmiL = abs(vmiL) + a
        ! right i
        vm   = 0.5_R_P*(block%fluid%P(i,j,k)%v+block%fluid%P(i+1,j,k)%v)
        vmiR = (vm.dot.block%mesh%NFi(i,j,k))*block%mesh%Si(i,j,k)
        vmiR = abs(vmiR) + a
        ! left j
        vm   = 0.5_R_P*(block%fluid%P(i,j-1,k)%v+block%fluid%P(i,j,k)%v)
        vmjL = (vm.dot.block%mesh%NFj(i,j-1,k))*block%mesh%Sj(i,j-1,k)
        vmjL = abs(vmjL) + a
        ! right j
        vm   = 0.5_R_P*(block%fluid%P(i,j,k)%v+block%fluid%P(i,j+1,k)%v)
        vmjR = (vm.dot.block%mesh%NFj(i,j,k))*block%mesh%Sj(i,j,k)
        vmjR = abs(vmjR) + a
        ! left k
        vm   = 0.5_R_P*(block%fluid%P(i,j,k-1)%v+block%fluid%P(i,j,k)%v)
        vmkL = (vm.dot.block%mesh%NFk(i,j,k-1))*block%mesh%Sk(i,j,k-1)
        vmkL = abs(vmkL) + a
        ! right k
        vm   = 0.5_R_P*(block%fluid%P(i,j,k)%v+block%fluid%P(i,j,k+1)%v)
        vmkR = (vm.dot.block%mesh%NFk(i,j,k))*block%mesh%Sk(i,j,k)
        vmkR = abs(vmkR) + a
        ! vmax
        vmax = max(vmiL,vmiR,vmjL,vmjR,vmkL,vmkR)
        block%fluid%Dt(i,j,k) = block%mesh%V(i,j,k)/vmax*global%fluid%CFL
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  ! computing minimum Dt
  Dtmin = minval(block%fluid%Dt(1:Ni,1:Nj,1:Nk))
  ! ghost cells estrapolation: imposing the minum value of Dt
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(gc,Ni,Nj,Nk,Dtmin,block)
  ! left i frame
  !$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),0
        block%fluid%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! right i frame
  !$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=Ni+1,Ni+gc(2)
        block%fluid%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! left j frame
  !$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),0
      do i=1-gc(1),Ni+gc(2)
        block%fluid%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! right j frame
  !$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=Nj+1,Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        block%fluid%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! left k frame
  !$OMP DO
  do k=1-gc(5),0
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        block%fluid%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  ! right k frame
  !$OMP DO
  do k=Nk+1,Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        block%fluid%Dt(i,j,k) = Dtmin
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_time

  subroutine residuals(s1,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing the residuals. This the space operator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    s1                                                  ! Current Runge-kutta stage.
  type(Type_Global), intent(IN)::    global                                              ! Global-level data.
  type(Type_Block),  intent(INOUT):: block                                               ! Block-level data.
  type(Type_Conservative)::          Fi(0:block%mesh%Ni,1:block%mesh%Nj,1:block%mesh%Nk) ! I fluxes.
  type(Type_Conservative)::          Fj(1:block%mesh%Ni,0:block%mesh%Nj,1:block%mesh%Nk) ! J fluxes.
  type(Type_Conservative)::          Fk(1:block%mesh%Ni,1:block%mesh%Nj,0:block%mesh%Nk) ! K fluxes.
  integer(I1P)::                     gcu                                                 ! Number of ghost cells used.
  integer(I_P)::                     Ni,Nj,Nk,Ns                                         ! Temp var for storing block dims.
  integer(I1P)::                     gc(1:6)                                             ! Temp var for storing ghost cells num.
  integer(I_P)::                     i,j,k                                               ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = block%mesh%Ni
  Nj = block%mesh%Nj
  Nk = block%mesh%Nk
  gc = block%mesh%gc
  Ns = global%fluid%Ns
#ifndef NULi
  ! computing the fluxes in i direction
  gcu = min(global%mesh%gco,gc(1),gc(2))
  do k=1,Nk
    do j=1,Nj
      call fluxes(gc  = gcu,                              &
                  N   = Ni,                               &
                  Ns  = Ns,                               &
                  cp0 = global%fluid%cp0,                 &
                  cv0 = global%fluid%cv0,                 &
                  NF  = block%mesh%NFi(0-gcu:Ni+gcu,j,k), &
                  P   = block%fluid%P (1-gcu:Ni+gcu,j,k), &
                  F   = Fi            (    0:Ni    ,j,k))
    enddo
  enddo
#endif
#ifndef NULj
  ! computing the fluxes in j direction
  gcu = min(global%mesh%gco,gc(3),gc(4))
  do k=1,Nk
    do i=1,Ni
      call fluxes(gc  = gcu,                              &
                  N   = Nj,                               &
                  Ns  = Ns,                               &
                  cp0 = global%fluid%cp0,                 &
                  cv0 = global%fluid%cv0,                 &
                  NF  = block%mesh%NFj(i,0-gcu:Nj+gcu,k), &
                  P   = block%fluid%P (i,1-gcu:Nj+gcu,k), &
                  F   = Fj            (i,    0:Nj    ,k))
    enddo
  enddo
#endif
#ifndef NULk
  ! computing the fluxes in k direction
  gcu = min(global%mesh%gco,gc(5),gc(6))
  do j=1,Nj
    do i=1,Ni
      call fluxes(gc  = gcu,                              &
                  N   = Nk,                               &
                  Ns  = Ns,                               &
                  cp0 = global%fluid%cp0,                 &
                  cv0 = global%fluid%cv0,                 &
                  NF  = block%mesh%NFk(i,j,0-gcu:Nk+gcu), &
                  P   = block%fluid%P (i,j,1-gcu:Nk+gcu), &
                  F   = Fk            (i,j,    0:Nk    ))
    enddo
  enddo
#endif
#if defined NULi || NULj || NULk
  ! nullify fluxes and velocity for 2D and 1D simulations
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(Ni,Nj,Nk,Fi,Fj,Fk)
#ifdef NULi
  ! nullify i direction
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=0,Ni
        Fi(i,j,k) = 0._R_P
      enddo
    enddo
  enddo
#endif
#ifdef NULj
  ! nullify j direction
  !$OMP DO
  do k=1,Nk
    do j=0,Nj
      do i=1,Ni
        Fj(i,j,k) = 0._R_P
      enddo
    enddo
  enddo
#endif
#ifdef NULk
  ! nullify k direction
  !$OMP DO
  do k=0,Nk
    do j=1,Nj
      do i=1,Ni
        Fk(i,j,k) = 0._R_P
      enddo
    enddo
  enddo
#endif
  !$OMP END PARALLEL
#endif

  ! computing the residuals
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(s1,Ni,Nj,Nk,Ns,block,Fi,Fj,Fk)
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        block%fluid%KS(i,j,k,s1) = (                                                                              &
#ifndef NULi
                                    + block%mesh%Si(i-1,j,  k  )*Fi(i-1,j,  k  ) - block%mesh%Si(i,j,k)*Fi(i,j,k) &
#endif
#ifndef NULj
                                    + block%mesh%Sj(i,  j-1,k  )*Fj(i,  j-1,k  ) - block%mesh%Sj(i,j,k)*Fj(i,j,k) &
#endif
#ifndef NULk
                                    + block%mesh%Sk(i,  j,  k-1)*Fk(i,  j,  k-1) - block%mesh%Sk(i,j,k)*Fk(i,j,k) &
#endif
                                    )/block%mesh%V(i,j,k)
      enddo
    enddo
  enddo
#ifdef NULi
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        block%fluid%KS(i,j,k,s1)%rv%x = 0._R_P
      enddo
    enddo
  enddo
#endif
#ifdef NULj
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        block%fluid%KS(i,j,k,s1)%rv%y = 0._R_P
      enddo
    enddo
  enddo
#endif
#ifdef NULk
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        block%fluid%KS(i,j,k,s1)%rv%z = 0._R_P
      enddo
    enddo
  enddo
#endif
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    subroutine fluxes(gc,N,Ns,cp0,cv0,NF,P,F)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for computing interfaces fluxes. The 3D primitive variables are first projected in the 1D normal of the interface;
    ! then the 1D Riemann problem is solved and the 3D fluxes are computed.
    ! 1D remapping:
    ! 1)    density of species 1    (r1)
    ! 2)    density of species 2    (r2)
    ! ...
    ! s)    density of species s-th (rs)
    ! ...
    ! Ns)   density of species Ns   (rNs)
    ! Ns+1) velocity                (u normal component)
    ! Ns+2) pressure                (p)
    ! Ns+3) density                 (r=sum(rs))
    ! Ns+4) specific heats ratio    (g)
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I1P),            intent(IN)::    gc                         ! Number of ghost cells used.
    integer(I_P),            intent(IN)::    N                          ! Number of cells.
    integer(I_P),            intent(IN)::    Ns                         ! Number of species.
    real(R_P),               intent(IN)::    cp0(1:Ns)                  ! Initial specific heat cp.
    real(R_P),               intent(IN)::    cv0(1:Ns)                  ! Initial specific heat cv.
    type(Type_Vector),       intent(IN)::    NF (           0-gc:N+gc ) ! Face normal.
    type(Type_Primitive),    intent(IN)::    P  (           1-gc:N+gc ) ! Primitive variables (3D).
    type(Type_Conservative), intent(INOUT):: F  (              0:N    ) ! Flux of mass (3D).
    type(Type_Vector)::                      ut (       1:2,1-gc:N+gc ) ! Tangential velocity: left (1) and right (2) interfaces.
    real(R_P)::                              P1D(1:Ns+2,1:2,1-gc:-1+gc) ! Primitive variables 1D fixed to the current interface.
#ifdef RECVC
    real(R_P)::                              Pm (1:Ns+4,       0:N    ) ! Mean of primitive variables across interfaces.
    real(R_P)::                              LPm(1:Ns+2,1:Ns+2,0:N    ) ! Mean left eigenvectors matrix.
    real(R_P)::                              RPm(1:Ns+2,1:Ns+2,0:N    ) ! Mean right eigenvectors matrix.
    real(R_P)::                              C  (1:Ns+2,1:2,1-gc:-1+gc) ! Interface value of characteristic variables.
    real(R_P)::                              CR (1:Ns+2,1:2           ) ! Reconstructed interface value of characteristic variables.
#endif
    real(R_P)::                              PR (1:Ns+4,1:2,   0:N+1  ) ! Reconstructed interface values of primitive variables 1D.
    real(R_P)::                              F_r                        ! Flux of mass (1D).
    real(R_P)::                              F_u                        ! Flux of momentum (1D).
    real(R_P)::                              F_e                        ! Flux of energy (1D).
    integer(I_P)::                           i,j,k,v                    ! Counters.
    logical::                                ROR(1:2)                   ! Logical flag for testing the result of WENO reconst.
    integer(I_P)::                           or                         ! Counter of order for ROR algorithm.
#if defined HYBRID || HYBRIDC
    logical::                                smooth  (1-gc:N+gc  )      ! Logical flag for testing the smoothness of the stencils.
    logical::                                smooth_f(1-gc:N+gc-1)      ! Logical flag for testing the smoothness of the stencils.
#endif
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    call init_cons(Ns=Ns,cons=F)
    !$OMP PARALLEL DEFAULT(NONE)                 &
    !$OMP PRIVATE(i,j,k,v,or,ROR,P1D,F_r,F_u,F_e &
#ifdef RECVC
    !$OMP         ,C,CR)                         &
#else
    !$OMP         )                              &
#endif
    !$OMP SHARED(gc,N,Ns,NF,ut,cp0,cv0,P,PR,F    &
#ifdef RECVC
    !$OMP        ,Pm,LPm,RPm                     &
#endif
#if defined HYBRID || HYBRIDC
    !$OMP        ,smooth,smooth_f)
#else
    !$OMP        )
#endif

    ! computing velocity vectors and its tangential component
    !$OMP DO
    do i=1-gc,N+gc
      ! computing tangential velocity component
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==1-gc.AND.j==1) cycle
        if (i==N+gc.AND.j==2) cycle
        ut(j,i) = P(i)%v - (P(i)%v.paral.NF(i+j-2)) ! ut
      enddo
    enddo

    ! computing the left and right states or Riemann Problems
    select case(gc)
    case(2_I1P,3_I1P,4_I1P) ! 3rd, 5th or 7th order WENO reconstruction
#if defined HYBRID || HYBRIDC
      ! checking the smoothness of the interface if high-order method is used
#ifdef SMSWliu
      !$OMP DO
      do i=0,N+1 ! loop over cells
        smooth(i) = chk_smooth_liu(pm1  =             P(i-1)%p,               &
                                   rm1  =             P(i-1)%d,               &
                                   um1  =            (P(i-1)%v.dot.NF(i-1)),  &
                                   gm1  =             P(i-1)%g,               &
                                   cpm1 = dot_product(P(i-1)%r/P(i-1)%d,cp0), &
                                   cvm1 = dot_product(P(i-1)%r/P(i-1)%d,cv0), &
                                   p1   =             P(i  )%p,               &
                                   r1   =             P(i  )%d,               &
                                   u1m  =            (P(i  )%v.dot.NF(i-1)),  &
                                   u1p  =            (P(i  )%v.dot.NF(i  )),  &
                                   g1   =             P(i  )%g,               &
                                   cp1  = dot_product(P(i  )%r/P(i  )%d,cp0), &
                                   cv1  = dot_product(P(i  )%r/P(i  )%d,cv0), &
                                   pp1  =             P(i+1)%p,               &
                                   rp1  =             P(i+1)%d,               &
                                   up1  =            (P(i+1)%v.dot.NF(i  )),  &
                                   gp1  =             P(i+1)%g,               &
                                   cpp1 = dot_product(P(i+1)%r/P(i+1)%d,cp0), &
                                   cvp1 = dot_product(P(i+1)%r/P(i+1)%d,cv0))
      enddo
#elif defined SMSWvanalbada || defined SMSWvanleer || defined SMSWharten
      !$OMP DO
      do i=0,N+1 ! loop over cells
        smooth(i) = chk_smooth_limiter(pm1 =  P(i  )%p              -  P(i-1)%p,              &
                                       pp1 =  P(i+1)%p              -  P(i  )%p,              &
                                       rm1 =  P(i  )%d              -  P(i-1)%d,              &
                                       rp1 =  P(i+1)%d              -  P(i  )%d,              &
                                       um1 = (P(i  )%v.dot.NF(i-1)) - (P(i-1)%v.dot.NF(i-1)), &
                                       up1 = (P(i+1)%v.dot.NF(i  )) - (P(i  )%v.dot.NF(i  )), &
                                       gm1 =  P(i  )%g              -  P(i-1)%g,              &
                                       gp1 =  P(i+1)%g              -  P(i  )%g)
      enddo
#elif defined SMSWz
      !$OMP DO
      do i=1-gc,N+gc-1 ! loop over interfaces
        smooth_f(i) = chk_smooth_z(p1 =  P(i  )%p,            &
                                   r1 =  P(i  )%d,            &
                                   u1 = (P(i  )%v.dot.NF(i)), &
                                   g1 =  P(i  )%g,            &
                                   p4 =  P(i+1)%p,            &
                                   r4 =  P(i+1)%d,            &
                                   u4 = (P(i+1)%v.dot.NF(i)), &
                                   g4 =  P(i+1)%g)
      enddo
      !$OMP DO
      do i=1-gc+1,N+gc-1 ! loop over cells
        smooth(i) = (smooth_f(i-1).AND.smooth_f(i))
      enddo
      !$OMP SINGLE
      smooth(1-gc)=smooth(1-gc+1)
      smooth(N+gc)=smooth(N+gc-1)
      !$OMP END SINGLE
#endif
      !$OMP SINGLE
      do i=1-gc,N+gc ! loop over cells
      enddo
      !$OMP END SINGLE
#endif

#ifdef RECVC
      ! computing Pm, LPm and RPm across the interfaces  for the local characteristic projection
      !$OMP DO
      do i=0,N ! loop over interfaces
        ! computing mean of primitive variables across the interface
        Pm(1:Ns,i) = 0.5_R_P*( P(i)%r(1:NS)      +  P(i+1)%r(1:NS)     ) ! rs
        Pm(Ns+1,i) = 0.5_R_P*((P(i)%v.dot.NF(i)) + (P(i+1)%v.dot.NF(i))) ! un
        Pm(Ns+2,i) = 0.5_R_P*( P(i)%p            +  P(i+1)%p           ) ! p
        Pm(Ns+3,i) = 0.5_R_P*( P(i)%d            +  P(i+1)%d           ) ! r
        Pm(Ns+4,i) = 0.5_R_P*( P(i)%g            +  P(i+1)%g           ) ! g
        ! computing mean left and right eigenvectors matrix across the interface
        LPm(1:Ns+2,1:Ns+2,i) = LP(Ns,Pm(1:Ns+4,i))
        RPm(1:Ns+2,1:Ns+2,i) = RP(Ns,Pm(1:Ns+4,i))
      enddo
#endif

      ! computing the reconstruction
      !$OMP DO
      do i=0,N+1 ! loop over cells
        ! computing 1D primitive variables for the stencil [i+1-gc,i-1+gc] fixing the velocity across the interfaces i+-1/2
        do k=i+1-gc,i-1+gc
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            P1D(1:Ns,j,k-i) =  P(k)%r(1:Ns)
            P1D(Ns+1,j,k-i) = (P(k)%v.dot.NF(i+j-2))
            P1D(Ns+2,j,k-i) =  P(k)%p
          enddo
        enddo

#ifdef RECVC
        ! transforming variables into local characteristic fields for the stencil [1-gc,-1+gc]
        do k=i+1-gc,i-1+gc
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            do v=1,Ns+2
              C(v,j,k-i) = dot_product(LPm(v,1:Ns+2,i+j-2),P1D(1:Ns+2,j,k-i))
            enddo
          enddo
        enddo
#endif

#ifdef HYBRID
        ! if stencil is smooth the WENO algorithm is computed directly with optimal weights
        if (all(smooth(i+1-gc:i-1+gc))) then
          ! computing WENO reconstruction with ideal weights (without the computation of smoothness indicator)
#ifdef RECVC
          do v=1,Ns+2
            CR(v,1:2) = weno_optimal(S=gc,V=C(v,1:2,1-gc:-1+gc))
          enddo
          ! trasforming local reconstructed characteristic variables to primitive ones
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            do v=1,Ns+2
              PR(v,j,i) = dot_product(RPm(v,1:Ns+2,i+j-2),CR(1:Ns+2,j))
            enddo
            PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
            PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
          enddo
#else
          do v=1,Ns+2
            PR(v,1:2,i) = weno_optimal(S=gc,V=P1D(v,1:2,1-gc:-1+gc))
          enddo
          ! computing the last reconstructed variables
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
            PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
          enddo
#endif
          cycle ! the full WENO reconstructions are not necessary
        endif
#elif HYBRIDC
        ! if stencil is smooth the WENO algorithm is switched to the central difference reconstruction
        if (all(smooth(i+1-gc:i-1+gc))) then
          ! computing central difference reconstruction
#ifdef RECVC
          do v=1,Ns+2
            CR(v,1:2) = noweno_central(S=gc,V=C(v,1:2,1-gc:-1+gc))
          enddo
          ! trasforming local reconstructed characteristic variables to primitive ones
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            do v=1,Ns+2
              PR(v,j,i) = dot_product(RPm(v,1:Ns+2,i+j-2),CR(1:Ns+2,j))
            enddo
            PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
            PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
          enddo
#else
          do v=1,Ns+2
            PR(v,1:2,i) = noweno_central(S=gc,V=P1D(v,1:2,1-gc:-1+gc))
          enddo
          ! computing the last reconstructed variables
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
            PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
          enddo
#endif
          cycle ! the full WENO reconstructions are not necessary
        endif
#endif

        ! computing WENO reconstruction with ROR (Recursive Order Reduction) algorithm
        ROR_check: do or=gc,2,-1
          ! computing WENO reconstruction
#ifdef RECVC
          do v=1,Ns+2
            CR(v,1:2) = weno(S=or,V=C(v,1:2,1-or:-1+or))
          enddo
          ! trasforming local reconstructed characteristic variables to primitive ones
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            do v=1,Ns+2
              PR(v,j,i) = dot_product(RPm(v,1:Ns+2,i+j-2),CR(1:Ns+2,j))
            enddo
          enddo
#else
          do v=1,Ns+2
            PR(v,1:2,i) = weno(S=gc,V=P1D(v,1:2,1-gc:-1+gc))
          enddo
#endif
          ! extrapolation of the reconstructed values for the non computed boundaries
          if (i==0  ) PR(1:Ns+2,1,i) = PR(1:Ns+2,2,i)
          if (i==N+1) PR(1:Ns+2,2,i) = PR(1:Ns+2,1,i)

#ifdef PPL
          ! applaying the maximum-principle-satisfying limiter to the reconstructed densities and pressure
          do v=1,Ns
            call positivity_preserving_limiter(o = gc, vmean = P(i)%r(v), vr = PR(v,1:2,i)) ! densities
          enddo
          call positivity_preserving_limiter(o = gc, vmean = P(i)%p, vr = PR(Ns+2,1:2,i))   ! pressure
#endif

          ! computing the reconstructed variables of total density and specific heats ratio
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
            PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
          enddo

          ! verifing the WENO reconstruction at order or-th
          ROR = .true.
          do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            if (i==0  .AND.j==1) cycle
            if (i==N+1.AND.j==2) cycle
            ROR(j) = ((PR(Ns+2,j,i)>0._R_P).AND.(PR(Ns+3,j,i)>0._R_P).AND.(PR(Ns+4,j,i)>0._R_P))
          enddo
          if (ROR(1).AND.ROR(2)) exit ROR_check
          if (or==2) then
            ! the high order reconstruction is falied; used 1st order reconstruction
            do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
              if (i==0  .AND.j==1) cycle
              if (i==N+1.AND.j==2) cycle
              PR(1:Ns,j,i) =  P(i)%r(1:Ns)          ! rs
              PR(Ns+1,j,i) = (P(i)%v.dot.NF(i+j-2)) ! un
              PR(Ns+2,j,i) =  P(i)%p                ! p
              PR(Ns+3,j,i) =  P(i)%d                ! r
              PR(Ns+4,j,i) =  P(i)%g                ! g
            enddo
          endif
        enddo ROR_check
      enddo
    case(1_I1P) ! 1st order piecewise constant reconstruction
      !$OMP DO
      do i=0,N+1
        do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          if (i==0  .AND.j==1) cycle
          if (i==N+1.AND.j==2) cycle
          PR(1:Ns,j,i) =  P(i)%r(1:Ns)          ! rs
          PR(Ns+1,j,i) = (P(i)%v.dot.NF(i+j-2)) ! un
          PR(Ns+2,j,i) =  P(i)%p                ! p
          PR(Ns+3,j,i) =  P(i)%d                ! r
          PR(Ns+4,j,i) =  P(i)%g                ! g
        enddo
      enddo
    endselect

    ! solving Rimeann Problems
    !$OMP DO
    do i=0,N
      ! face normal fluxes
      call Riem_Solver(p1  = PR(Ns+2,2,i  ), & !| right value of cell i   (letf  of interface i+1/2)
                       r1  = PR(Ns+3,2,i  ), & !|
                       u1  = PR(Ns+1,2,i  ), & !|
                       g1  = PR(Ns+4,2,i  ), & !|
#if defined RSHLLCb || defined RSROE
                       cp1 = dot_product(PR(1:Ns,2,i)/PR(Ns+3,2,i),cp0), &
                       cv1 = dot_product(PR(1:Ns,2,i)/PR(Ns+3,2,i),cv0), &
#endif
                       p4  = PR(Ns+2,1,i+1), & !| left  value of cell i+1 (right of interface i+1/2)
                       r4  = PR(Ns+3,1,i+1), & !|
                       u4  = PR(Ns+1,1,i+1), & !|
                       g4  = PR(Ns+4,1,i+1), & !|
#if defined RSHLLCb || defined RSROE
                       cp4 = dot_product(PR(1:Ns,1,i+1)/PR(Ns+3,1,i+1),cp0), &
                       cv4 = dot_product(PR(1:Ns,1,i+1)/PR(Ns+3,1,i+1),cv0), &
#endif
                       F_r = F_r           , &
                       F_u = F_u           , &
                       F_E = F_e)
      ! uptdating fluxes with tangential components
      if (F_r>0._R_P) then
        F(i)%rs = F_r*PR(1:Ns,2,i  )/PR(Ns+3,2,i  )
        F(i)%rv = F_u*NF(i) + F_r*ut(2,i)
        F(i)%re = F_e + F_r*sq_norm(ut(2,i))*0.5_R_P
      else
        F(i)%rs = F_r*PR(1:Ns,1,i+1)/PR(Ns+3,1,i+1)
        F(i)%rv = F_u*NF(i) + F_r*ut(1,i+1)
        F(i)%re = F_e + F_r*sq_norm(ut(1,i+1))*0.5_R_P
      endif
    enddo

    !$OMP END PARALLEL
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine fluxes

    pure function LP(Ns,primitive) result(L)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Function for computing left eigenvectors from primitive variables.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN):: Ns                ! Number of species.
    real(R_P),    intent(IN):: primitive(1:Ns+4) ! Primitive variables.
    real(R_P)::                L(1:Ns+2,1:Ns+2)  ! Left eigenvectors matrix.
    real(R_P)::                a                 ! Speed of sound.
    real(R_P)::                a2                ! a^2.
    real(R_P)::                r_inv             ! 1/r.
    real(R_P)::                ar_inv            ! 1/(a*r).
    integer(I_P)::             i                 ! Counter.
    !real(R_P)::                gp                ! g*p.
    !real(R_P)::                gp_a              ! g*p/a.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ! old not correct
    !a    = sqrt(primitive(Ns+4)*primitive(Ns+2)/primitive(Ns+3))
    !gp   = primitive(Ns+4)*primitive(Ns+2)
    !gp_a = gp/a

    !L = 0._R_P
    !do i=1,Ns
      !L(i,i) = gp
      !L(i,Ns+2) =-primitive(i)
    !enddo
    !L(Ns+1,Ns+1) =-gp_a ; L(Ns+1,Ns+2) = 1._R_P
    !L(Ns+2,Ns+1) = gp_a ; L(Ns+2,Ns+2) = 1._R_P
    ! old not correct

    a      = sqrt(primitive(Ns+4)*primitive(Ns+2)/primitive(Ns+3))
    a2     = a*a
    r_inv  = 1._R_P/(primitive(Ns+3))
    ar_inv = r_inv*1._R_P/(a)

    L = 0._R_P
    do i=1,Ns
      L(i,i   ) = a2
      L(i,Ns+2) = -primitive(i)*r_inv
    enddo
    L(Ns+1,Ns+1) = -1._R_P ; L(Ns+1,Ns+2) = ar_inv
    L(Ns+2,Ns+1) =  1._R_P ; L(Ns+2,Ns+2) = ar_inv
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction LP

    pure function RP(Ns,primitive) result(R)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Function for computing right eigenvectors from primitive variables.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN):: Ns                ! Number of species.
    real(R_P),    intent(IN):: primitive(1:Ns+4) ! Primitive variables.
    real(R_P)::                R(1:Ns+2,1:Ns+2)  ! Right eigenvectors matrix.
    real(R_P)::                a                 ! Speed of sound.
    real(R_P)::                a_inv             ! 1/a.
    real(R_P)::                a2_inv            ! 1/(a^2).
    real(R_P)::                a_2inv            ! 1/(2*a).
    real(R_P)::                ar_2              ! a*r/2.
    integer(I_P)::             i                 ! Counter.
    !real(R_P)::                gp_inv            ! 1/(g*p).
    !real(R_P)::                a_2gp             ! a/(2*g*p).
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ! old not correct
    !a      = sqrt(primitive(Ns+4)*primitive(Ns+2)/primitive(Ns+3))
    !gp_inv = 1._R_P/(primitive(Ns+4)*primitive(Ns+2))
    !a_2gp  = a*gp_inv*0.5_R_P

    !R = 0._R_P
    !do i=1,Ns
      !R(i,i) = gp_inv
      !R(i,Ns+1) = primitive(i)*gp_inv*0.5_R_P
      !R(i,Ns+2) = R(i,Ns+1)
    !enddo
    !R(Ns+1,Ns+1) =-a_2gp   ; R(Ns+1,Ns+2) = a_2gp
    !R(Ns+2,Ns+1) = 0.5_R_P ; R(Ns+2,Ns+2) = 0.5_R_P
    ! old not correct

    a      = sqrt(primitive(Ns+4)*primitive(Ns+2)/primitive(Ns+3))
    a_inv  = 1._R_P/a
    a2_inv = a_inv*a_inv
    a_2inv = 0.5_R_P*a_inv
    ar_2   = 0.5_R_P*a*primitive(Ns+3)

    R = 0._R_P
    do i=1,Ns
      R(i,i   ) = a2_inv
      R(i,Ns+1) = primitive(i)*a_2inv ; R(i,Ns+2) = R(i,Ns+1)
    enddo
    R(Ns+1,Ns+1) = -0.5_R_P ; R(Ns+1,Ns+2) = -R(Ns+1,Ns+1)
    R(Ns+2,Ns+1) =  ar_2    ; R(Ns+2,Ns+2) =  R(Ns+2,Ns+1)
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction RP

    pure subroutine positivity_preserving_limiter(o,vmean,vr)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Positivity preserving limiter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN)::    o                                        ! Order of space reconstruction.
    real(R_P),    intent(IN)::    vmean                                    ! Mean value.
    real(R_P),    intent(INOUT):: vr(1:2)                                  ! Reconstructed values.
    real(R_P), parameter::        w1(2:4)=[1._R_P/6._R_P,  &               ! Weights of Gauss-Lobatto's 3 points quadrature.
                                           1._R_P/12._R_P, &               ! Weights of Gauss-Lobatto's 4 points quadrature.
                                           1._R_P/20._R_P]                 ! Weights of Gauss-Lobatto's 5 points quadrature.
    real(R_P), parameter::        w2(2:4)=[1._R_P/(1._R_P-2._R_P*w1(2)), & ! 1/(1-2*w1)
                                           1._R_P/(1._R_P-2._R_P*w1(3)), & ! 1/(1-2*w1)
                                           1._R_P/(1._R_P-2._R_P*w1(4))]   ! 1/(1-2*w1)
    real(R_P), parameter::        e  = tiny(1._R_P)                        ! Small number for avoid division for zero.
    real(R_P)::                   vstar                                    ! Star value.
    real(R_P)::                   theta                                    ! Limiter coefficient.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    vstar = (vmean - w1(o)*(vr(1) + vr(2)))*w2(o)
    theta = min((vmean-0.1_R_P*vmean)/(vmean - min(vstar,vr(1),vr(2)) + e),1._R_P)
    vr(1) = theta*(vr(1) - vmean) + vmean
    vr(2) = theta*(vr(2) - vmean) + vmean
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine positivity_preserving_limiter
  endsubroutine residuals

  subroutine boundary_conditions(myrank,l,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! The subroutine boundary_conditions imposes the boundary conditions of blocks of grid level "l" updating the ghost cells.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank                  ! Actual rank process.
  integer(I_P),      intent(IN)::    l                       ! Actual grid level.
  type(Type_Global), intent(IN)::    global                  ! Global-level data.
  type(Type_Block),  intent(INOUT):: block(1:global%mesh%Nb) ! Block-level data.
  integer(I_P)::                     Ni,Nj,Nk,gc(1:6)        ! Temp var for storing block dimensions.
  integer(I_P)::                     b,i,j,k                 ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef MPI2
  ! doing the multi-processes comunications if necessary
  call Psendrecv(myrank=myrank,l=l,global=global,block=block)
#endif
  !$OMP PARALLEL DEFAULT(NONE)       &
  !$OMP PRIVATE(b,i,j,k,Ni,Nj,Nk,gc) &
  !$OMP SHARED(l,global,block)
  !$OMP DO
  do b=1,global%mesh%Nb
    Ni = block(b)%mesh%Ni
    Nj = block(b)%mesh%Nj
    Nk = block(b)%mesh%Nk
    gc = block(b)%mesh%gc
    do k=1,Nk
      do j=1,Nj
        ! left i
        select case(block(b)%bc%BCi(0,j,k)%tp)
        case(bc_ext)
          call set_ext_l(gc=gc(1),ic=0,N=Ni,P=block(b)%fluid%P(1-gc(1):0+gc(1),j,k))
        case(bc_ref)
          call set_ref_l(gc=gc(1),ic=0,N=Ni,NF=block(b)%mesh%NFi(0,j,k),P=block(b)%fluid%P(1-gc(1):0+gc(1),j,k))
        case(bc_per)
          call set_per(gc=gc(1),ic=0,N=Ni,boundary='l',P=block(b)%fluid%P(1-gc(1):Ni+gc(1),j,k))
        case(bc_adj)
          call set_adj(gc=gc(1),bc=block(b)%bc%BCi(1-gc(1):0,j,k),P=block(b)%fluid%P(1-gc(1):0,j,k))
        case(bc_in1)
          call set_in1(gc=gc(1),bc=block(b)%bc%BCi(1-gc(1):0,j,k),P=block(b)%fluid%P(1-gc(1):0,j,k))
        case(bc_in2)
          !call set_in2(gc=gc(1),Np=Np,bc=block(b)%bc%BCi(1-gc(1):0,j,k),t=global%fluid%t,P=block(b)%fluid%P(1-gc(1):0,j,k))
        endselect
        ! right i
        select case(block(b)%bc%BCi(Ni,j,k)%tp)
        case(bc_ext)
          call set_ext_r(gc=gc(2),ic=0,N=Ni,P=block(b)%fluid%P(Ni-gc(2):Ni+gc(2),j,k))
        case(bc_ref)
          call set_ref_r(gc=gc(2),ic=0,N=Ni,NF=block(b)%mesh%NFi(Ni,j,k),P=block(b)%fluid%P(Ni-gc(2):Ni+gc(2),j,k))
        case(bc_per)
          call set_per(gc=gc(2),ic=0,N=Ni,boundary='r',P=block(b)%fluid%P(1-gc(2):Ni+gc(2),j,k))
        case(bc_adj)
          call set_adj(gc=gc(2),bc=block(b)%bc%BCi(Ni:Ni+gc(2)-1,j,k),P=block(b)%fluid%P(Ni+1:Ni+gc(2),j,k))
        case(bc_in1)
          call set_in1(gc=gc(2),bc=block(b)%bc%BCi(Ni:Ni+gc(2)-1,j,k),P=block(b)%fluid%P(Ni+1:Ni+gc(2),j,k))
        case(bc_in2)
          !call set_in2(gc=gc(2),Np=Np,bc=block(b)%bc%BCi(Ni:Ni+gc(2)-1,j,k),t=global%fluid%t,P=block(b)%fluid%P(Ni+1:Ni+gc(2),j,k))
        endselect
      enddo
    enddo
    do k=1,Nk
      do i=1,Ni
        ! left j
        select case(block(b)%bc%BCj(i,0,k)%tp)
        case(bc_ext)
          call set_ext_l(gc=gc(3),ic=0,N=Nj,P=block(b)%fluid%P(i,1-gc(3):0+gc(3),k))
        case(bc_ref)
          call set_ref_l(gc=gc(3),ic=0,N=Nj,NF=block(b)%mesh%NFj(i,0,k),P=block(b)%fluid%P(i,1-gc(3):0+gc(3),k))
        case(bc_per)
          call set_per(gc=gc(3),ic=0,N=Nj,boundary='l',P=block(b)%fluid%P(i,1-gc(3):Nj+gc(3),k))
        case(bc_adj)
          call set_adj(gc=gc(3),bc=block(b)%bc%BCj(i,1-gc(3):0,k),P=block(b)%fluid%P(i,1-gc(3):0,k))
        case(bc_in1)
          call set_in1(gc=gc(3),bc=block(b)%bc%BCj(i,1-gc(3):0,k),P=block(b)%fluid%P(i,1-gc(3):0,k))
        case(bc_in2)
          !call set_in2(gc=gc(3),Np=Np,bc=block(b)%bc%BCj(i,1-gc(3):0,k),t=global%fluid%t,P=block(b)%fluid%P(i,1-gc(3):0,k))
        endselect
        ! right j
        select case(block(b)%bc%BCj(i,Nj,k)%tp)
        case(bc_ext)
          call set_ext_r(gc=gc(4),ic=0,N=Nj,P=block(b)%fluid%P(i,Nj-gc(4):Nj+gc(4),k))
        case(bc_ref)
          call set_ref_r(gc=gc(4),ic=0,N=Nj,NF=block(b)%mesh%NFj(i,Nj,k),P=block(b)%fluid%P(i,Nj-gc(4):Nj+gc(4),k))
        case(bc_per)
          call set_per(gc=gc(4),ic=0,N=Nj,boundary='r',P=block(b)%fluid%P(i,1-gc(4):Nj+gc(4),k))
        case(bc_adj)
          call set_adj(gc=gc(4),bc=block(b)%bc%BCj(i,Nj:Nj+gc(4)-1,k),P=block(b)%fluid%P(i,Nj+1:Nj+gc(4),k))
        case(bc_in1)
          call set_in1(gc=gc(4),bc=block(b)%bc%BCj(i,Nj:Nj+gc(4)-1,k),P=block(b)%fluid%P(i,Nj+1:Nj+gc(4),k))
        case(bc_in2)
          !call set_in2(gc=gc(4),Np=Np,bc=block(b)%bc%BCj(i,Nj:Nj+gc(4)-1,k),t=global%fluid%t,P=block(b)%fluid%P(i,Nj+1:Nj+gc(4),k))
        endselect
      enddo
    enddo
    do j=1,Nj
      do i=1,Ni
        ! left k
        select case(block(b)%bc%BCk(i,j,0)%tp)
        case(bc_ext)
          call set_ext_l(gc=gc(5),ic=0,N=Nk,P=block(b)%fluid%P(i,j,1-gc(5):0+gc(5)))
        case(bc_ref)
          call set_ref_l(gc=gc(5),ic=0,N=Nk,NF=block(b)%mesh%NFk(i,j,0),P=block(b)%fluid%P(i,j,1-gc(5):0+gc(5)))
        case(bc_per)
          call set_per(gc=gc(5),ic=0,N=Nk,boundary='l',P=block(b)%fluid%P(i,j,1-gc(5):Nk+gc(5)))
        case(bc_adj)
          call set_adj(gc=gc(5),bc=block(b)%bc%BCk(i,j,1-gc(5):0),P=block(b)%fluid%P(i,j,1-gc(5):0))
        case(bc_in1)
          call set_in1(gc=gc(5),bc=block(b)%bc%BCk(i,j,1-gc(5):0),P=block(b)%fluid%P(i,j,1-gc(5):0))
        case(bc_in2)
          !call set_in2(gc=gc(5),Np=Np,bc=block(b)%bc%BCk(i,j,1-gc(5):0),t=global%fluid%t,P=block(b)%fluid%P(i,j,1-gc(5):0))
        endselect
        ! right k
        select case(block(b)%bc%BCk(i,j,Nk)%tp)
        case(bc_ext)
          call set_ext_r(gc=gc(6),ic=0,N=Nk,P=block(b)%fluid%P(i,j,Nk-gc(6):Nk+gc(6)))
        case(bc_ref)
          call set_ref_r(gc=gc(6),ic=0,N=Nk,NF=block(b)%mesh%NFk(i,j,Nk),P=block(b)%fluid%P(i,j,Nk-gc(6):Nk+gc(6)))
        case(bc_per)
          call set_per(gc=gc(6),ic=0,N=Nk,boundary='r',P=block(b)%fluid%P(i,j,1-gc(6):Nk+gc(6)))
        case(bc_adj)
          call set_adj(gc=gc(6),bc=block(b)%bc%BCk(i,j,Nk:Nk+gc(6)-1),P=block(b)%fluid%P(i,j,Nk+1:Nk+gc(6)))
        case(bc_in1)
          call set_in1(gc=gc(6),bc=block(b)%bc%BCk(i,j,Nk:Nk+gc(6)-1),P=block(b)%fluid%P(i,j,Nk+1:Nk+gc(6)))
        case(bc_in2)
          !call set_in2(gc=gc(6),Np=Np,bc=block(b)%bc%BCk(i,j,Nk:Nk+gc(6)-1),t=global%fluid%t,P=block(b)%fluid%P(i,j,Nk+1:Nk+gc(6)))
        endselect
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure subroutine set_ext_l(gc,ic,N,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for imposing extrapolation of ghost cells from internal ones (left boundary).
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc           ! Number of ghost cells.
    integer(I_P),         intent(IN)::    ic           ! Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N            ! Number of internal cells.
    type(Type_Primitive), intent(INOUT):: P(1-gc:0+gc) ! Primitive variables.
    integer(I_P)::                        i            ! Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (ic==1.OR.N<gc) then
      ! extrapolation using only the cell 1
      do i=1-gc,0
        P(i) = P(1)
      enddo
    else
      ! extrapolation using the cells 1,2,...,gc
      do i=1-gc,0
        P(i) = P(-i+1)
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_ext_l

    pure subroutine set_ext_r(gc,ic,N,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for imposing extrapolation of ghost cells from internal ones (right boundary).
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc           ! Number of ghost cells.
    integer(I_P),         intent(IN)::    ic           ! Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N            ! Number of internal cells.
    type(Type_Primitive), intent(INOUT):: P(N-gc:N+gc) ! Primitive variables.
    integer(I_P)::                        i            ! Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (ic==1.OR.N<gc) then
      ! extrapolation using only the cell N
      do i=N+1,N+gc
        P(i) = P(N)
      enddo
    else
      ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
      do i=N+1,N+gc
        P(i) = P(N+1-(i-N))
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_ext_r

    pure subroutine set_ref_l(gc,ic,N,NF,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for imposing reflective boundary conditions (left boundary).
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc           ! Number of ghost cells.
    integer(I_P),         intent(IN)::    ic           ! Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N            ! Number of internal cells.
    type(Type_Vector),    intent(IN)::    NF           ! Left face normal.
    type(Type_Primitive), intent(INOUT):: P(1-gc:0+gc) ! Left section of primitive variables.
    integer(I_P)::                        i            ! Cell counter.
    type(Type_Vector)::                   vr           ! Reflected velocity vector.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (ic==1.OR.N<gc) then
      ! reflection using only the cell 1
      vr = P(1)%v - (2._R_P*(P(1)%v.paral.NF)) ! reflected velocity
      do i=1-gc,0
        P(i)%r = P(1)%r
        P(i)%v = vr
        P(i)%p = P(1)%p
        P(i)%d = P(1)%d
        P(i)%g = P(1)%g
      enddo
    else
      ! reflection using the cells 1,2,...,gc
      do i=1-gc,0
        vr = P(-i+1)%v - (2._R_P*(P(-i+1)%v.paral.NF)) ! reflected velocity
        P(i)%r = P(-i+1)%r
        P(i)%v = vr
        P(i)%p = P(-i+1)%p
        P(i)%d = P(-i+1)%d
        P(i)%g = P(-i+1)%g
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_ref_l

    pure subroutine set_ref_r(gc,ic,N,NF,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for imposing reflective boundary conditions (right boundary).
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc           ! Number of ghost cells.
    integer(I_P),         intent(IN)::    ic           ! Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N            ! Number of internal cells.
    type(Type_Vector),    intent(IN)::    NF           ! Right face normal.
    type(Type_Primitive), intent(INOUT):: P(N-gc:N+gc) ! Right section of primitive variables.
    integer(I_P)::                        i            ! Cell counter.
    type(Type_Vector)::                   vr           ! Reflected velocity vector.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (ic==1.OR.N<gc) then
      ! reflection using only the cell N
      vr = P(N)%v - (2._R_P*(P(N)%v.paral.NF)) ! reflected velocity
      do i=N+1,N+gc
        P(i)%r = P(N)%r
        P(i)%v = vr
        P(i)%p = P(N)%p
        P(i)%d = P(N)%d
        P(i)%g = P(N)%g
      enddo
    else
      ! reflection using the cells N-gc,N-gc+1,N-gc+2,...,N
      do i=N+1,N+gc
        vr = P(N+1-(i-N))%v - (2._R_P*(P(N+1-(i-N))%v.paral.NF)) ! reflected velocity
        P(i)%r = P(N+1-(i-N))%r
        P(i)%v = vr
        P(i)%p = P(N+1-(i-N))%p
        P(i)%d = P(N+1-(i-N))%d
        P(i)%g = P(N+1-(i-N))%g
      enddo
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_ref_r

    pure subroutine set_per(gc,ic,N,boundary,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for imposing periodic boundary conditions.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc           ! Number of ghost cells.
    integer(I_P),         intent(IN)::    ic           ! Number of internal cells used for extrapolation (1 or gc).
    integer(I_P),         intent(IN)::    N            ! Number of internal cells.
    character(1),         intent(IN)::    boundary     ! Boundary left ('l') or right ('r').
    type(Type_Primitive), intent(INOUT):: P(1-gc:N+gc) ! Left section of primitive variables.
    integer(I_P)::                        i            ! Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if (boundary=='l') then
      if (ic==1.OR.N<gc) then
        ! extrapolation using only the cell N
        do i=1-gc,0
          P(i) = P(N)
        enddo
      else
        ! extrapolation using the cells N-gc,N-gc+1,N-gc+2,...,N
        do i=1-gc,0
          P(i) = P(i+N)
        enddo
      endif
    endif
    if (boundary=='r') then
      if (ic==1.OR.N<gc) then
        ! extrapolation using only the cell 1
        do i=N+1,N+gc
          P(i) = P(1)
        enddo
      else
        ! extrapolation using the cells 1,2,...,gc
        do i=N+1,N+gc
          P(i) = P(i-N)
        enddo
      endif
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_per

    pure subroutine set_adj(gc,bc,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for imposing adjacent boundary conditions.
    ! Note that when this subroutine is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays BC and P must be
    ! properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1] for BC and [N+1:N+gc] for P.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc         ! Number of ghost cells.
    type(Type_BC),        intent(IN)::    bc(1-gc:0) ! Boundary conditions infos.
    type(Type_Primitive), intent(INOUT):: P (1-gc:0) ! Left section of primitive variables.
    integer(I_P)::                        i,b        ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
#ifdef MPI2
    ! for multi-processes simulation it is possible that the adjacent cell is in other processes than the actual and in case the
    ! data have been already exchanged by the MPI subroutine Psendrecv
    if (procmap(bc(0)%adj%b)/=myrank) return
#endif
    do i=1-gc,0
      b = minloc(array=blockmap,dim=1,mask=blockmap==bc(i)%adj%b)
      P(i) = block(b)%fluid%P(bc(i)%adj%i,bc(i)%adj%j,bc(i)%adj%k)
      !call get_cell_lbijk(P = P(i),        &
      !                    l = l,           &
      !                    b = b,           &
      !                    i = BC(i)%adj%i, &
      !                    j = BC(i)%adj%j, &
      !                    k = BC(i)%adj%k)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_adj

    pure subroutine set_in1(gc,bc,P)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for imposing inflow 1 boundary conditions.
    ! Note that when this subroutine is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays BC(i,j,k) and P must be
    ! properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1].
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P),         intent(IN)::    gc         ! Number of ghost cells.
    type(Type_BC),        intent(IN)::    bc(1-gc:0) ! Boundary conditions.
    type(Type_Primitive), intent(INOUT):: P (1-gc:0) ! Left section of primitive variables.
    integer(I_P)::                        i          ! Cell counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    do i=1-gc,0
      P(i) = global%bc%in1(bc(i)%inf)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine set_in1

   !subroutine set_in2(gc,bc,t,P)
   !!-------------------------------------------------------------------------------------------------------------------------------
   !! Subroutine for imposing inflow 2 boundary conditions.
   !! Note that when this subroutine is called for a 'right' (Ni,Nj,Nk) boundary the section of arrays BC(i,j,k) and P must be
   !! properly remapped: in the section [1-gc:0] must be passed the actual section [N:N+gc-1].
   !!-------------------------------------------------------------------------------------------------------------------------------
   !
   !!-------------------------------------------------------------------------------------------------------------------------------
   !implicit none
   !integer(I_P),  intent(IN)::  gc         ! Number of ghost cells.
   !type(Type_BC), intent(IN)::  bc(1-gc:0) ! Boundary conditions.
   !real(R_P),     intent(IN)::  t          ! Time.
   !real(R_P),     intent(OUT):: P (1-gc:0) ! Left section of primitive variables.
   !integer(I_P)::               i          ! Cell counter.
   !integer(I_P)::               n          ! Time counter.
   !integer(I_P)::               v          ! Variable counter.
   !!-------------------------------------------------------------------------------------------------------------------------------
   !
   !!-------------------------------------------------------------------------------------------------------------------------------
   !do i=1-gc,0
   !  time_search: do n=bc_global%bb_in2(1,bc(i)%inf)+1,bc_global%bb_in2(2,bc(i)%inf)
   !    if (t<bc_global%in2(Np+1,n))          exit time_search
   !    if (n==bc_global%bb_in2(2,bc(i)%inf)) exit time_search
   !  enddo time_search
   !  do v=1,size()
   !    P(i)%r(v) = interpolate1(x1 = bc_global%in2(Np+1,n-1),x2 = bc_global%in2(Np+1,n), &
   !                             V1 = bc_global%in2(v   ,n-1),V2 = bc_global%in2(v   ,n), &
   !                             x  = t)
   !  enddo
   !enddo
   !return
   !!-------------------------------------------------------------------------------------------------------------------------------
   !endsubroutine set_in2
  endsubroutine boundary_conditions

  subroutine solve_grl(myrank,l,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for solving the conservation equations for a given grid level
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank                                   ! Actual rank process.
  integer(I_P),      intent(IN)::    l                                        ! Actual grid level.
  type(Type_Global), intent(INOUT):: global                                   ! Global-level data.
  type(Type_Block),  intent(INOUT):: block(1:global%mesh%Nb)                  ! Block-level data.
  real(R_P)::                        Dtmin(1:global%mesh%Nb)                  ! Min t step of actual process for each blk.
  real(R_P)::                        DtminL                                   ! Min t step of actual process over all blks.
  real(R_P)::                        gDtmin                                   ! Global (all processes/all blks) min t step.
  real(R_P)::                        RU  (1:global%fluid%Nc,1:global%mesh%Nb) ! NormL2 of conservartive residuals.
  real(R_P)::                        mRU (1:global%fluid%Nc)                  ! Maximum of RU of actual process.
  real(R_P)::                        gmRU(1:global%fluid%Nc)                  ! Global (all processes) maximum of RU.
  integer(I_P)::                     err                                      ! Error traping flag: 0 no errors, >0 errors.
  integer(I_P)::                     b                                        ! Blocks counter.
  integer(I_P)::                     s1                                       ! Runge-Kutta stages counters.
  real(R_P)::                        sec_elp                                  ! Seconds elapsed from the simulation start.
  real(R_P)::                        sec_res                                  ! Seconds residual from the simulation end.
  type(Type_Time)::                  time_elp                                 ! Time elapsed (in days,hours,min... format).
  type(Type_Time)::                  time_res                                 ! Time residual (in days,hours,min... format).
  real(R_P)::                        progress                                 ! Status (%) of simulation progress.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! converting conservative variables to primitive ones
  do b=1,global%mesh%Nb
    call conservative2primitive(global = global, block = block(b))
  enddo

  ! imposing the boundary conditions
  call boundary_conditions(myrank = myrank, l = l, global = global, block = block)

  ! saving the actual solution
  if (global%file%sol_out>0) then
    if ((mod(global%fluid%n,global%file%sol_out)==0).OR. &
        (global%fluid%t==global%fluid%Tmax).OR.          &
        (global%fluid%n==global%fluid%Nmax)) then
      do b=1,global%mesh%Nb
        err=save_bfluid(filename=file_name(basename=trim(global%file%Path_OutPut)//global%file%File_Sol,&
                                           suffix='.sol',blk=blockmap(b),grl=l,n=global%fluid%n),       &
                        global=global,block=block(b))
      enddo
    endif
  endif

  ! saving the restart solution file
  if ((mod(global%fluid%n,global%file%restart_out)==0).OR. &
      (global%fluid%t==global%fluid%Tmax).OR.              &
      (global%fluid%n==global%fluid%Nmax)) then
    flip = 1_I1P - flip
    do b=1,global%mesh%Nb
      err=save_bfluid(filename=file_name(basename=trim(global%file%Path_OutPut)//global%file%File_Sol,&
                                         suffix='.sol',blk=blockmap(b),grl=l,flip=flip),              &
                      global=global,block=block(b))
    enddo
  endif

  ! updating time varying variables: Dt,Dtmin
  global%fluid%n = global%fluid%n + 1_I8P
  do b=1,global%mesh%Nb
    call compute_time(global=global,block=block(b),Dtmin=Dtmin(b))
  enddo
  DtminL = minval(Dtmin)
#ifdef MPI2
  ! for multi-processes simulation all processes must exchange their DtminL for computing the global variables
  call MPI_ALLREDUCE(DtminL,gDtmin,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,err)
#else
  ! for single processes DtminL are already global variables
  gDtmin = DtminL
#endif
  if (global%fluid%unsteady) then ! for an unsteady accurate simulation each cell is updated by means of global minimum time step
    ! control for the last iterate
    if (global%fluid%Nmax<=0) then
      if ((global%fluid%t+gDtmin)>global%fluid%Tmax) then
        ! the global minimum time step is so high that the last iteration will go over Tmax
        ! it is decreased both for in order to achieve exactly Tmax
        gDtmin=abs(global%fluid%Tmax-global%fluid%t)
      endif
    endif
    global%fluid%t = global%fluid%t + gDtmin
    do b=1,global%mesh%Nb
      block(b)%fluid%Dt = gDtmin
    enddo
  endif

  ! updating console
  if (myrank==0) then
    if ((mod(global%fluid%n,global%file%screen_out)==0).OR. &
        (global%fluid%t==global%fluid%Tmax).OR. &
        (global%fluid%n==global%fluid%Nmax).OR. &
        (global%fluid%n==1)) then
      sec_elp=Crono()
      if (global%fluid%Nmax>0) then
        progress = global%fluid%n*100/(global%fluid%Nmax*1._R_P)
        sec_res  = global%fluid%Nmax*sec_elp/global%fluid%n - sec_elp
      elseif (global%fluid%Tmax>0._R_P) then
        progress = 100*global%fluid%t/global%fluid%Tmax
        sec_res  = global%fluid%Tmax*sec_elp/global%fluid%t - sec_elp
      else
        progress = 0._R_P
        sec_res  = 0._R_P
      endif
      time_elp = Seconds_To_Time(sec_elp)
      time_res = Seconds_To_Time(sec_res)
      write(stdout,'(A)',                      iostat=err)'----------------------------------------------------------------------'
      write(stdout,'(A39,I30)',                iostat=err)' Actual grid level                  l: ',l
      write(stdout,'(A39,23X,F6.2,A)',         iostat=err)' Simulation progress                p: ',progress,'%'
      write(stdout,'(A39,I30)',                iostat=err)' Actual step number                 n: ',global%fluid%n
      write(stdout,'(A39,ES30.12)',            iostat=err)' Actual simulation time             t: ',global%fluid%t
      write(stdout,'(A39,ES30.12)',            iostat=err)' Actual time step value           gDt: ',gDtmin
      write(stdout,*)
      write(stdout,'(10X,A15,20X,A15)',        iostat=err)' Elapsed Time','Residual Time'
      write(stdout,'(10X,A9,I6,20X,A9,I6)',    iostat=err)' Days    =',time_elp%Days,   'Days    =',time_res%Days
      write(stdout,'(10X,A9,I6,20X,A9,I6)',    iostat=err)' Hours   =',time_elp%Hours,  'Hours   =',time_res%Hours
      write(stdout,'(10X,A9,I6,20X,A9,I6)',    iostat=err)' Minutes =',time_elp%Minutes,'Minutes =',time_res%Minutes
      write(stdout,'(10X,A9,F6.1,20X,A9,F6.1)',iostat=err)' Seconds =',time_elp%Seconds,'Seconds =',time_res%Seconds
      write(stdout,'(A)',                      iostat=err)'----------------------------------------------------------------------'
      write(stdout,*)
    endif
  endif

  ! evaluating the Runge-Kutta stages
  ! Runge-Kutta stages initialization
  do b=1,global%mesh%Nb
    call init_cons(Ns=global%fluid%Ns,cons=block(b)%fluid%KS)
  enddo
  do s1=1,global%fluid%rk_ord
    if (s1>1) then
      ! summing the stages up to s1-1 for update P
      do b=1,global%mesh%Nb
          call rk_stages_sum(s1=s1,global=global,block=block(b))
      enddo
      ! imposing the boundary conditions
      call boundary_conditions(myrank = myrank, l = l, global = global, block = block)
    endif
    ! computing the s1-th Runge-Kutta stage: K_s1=R(Un+sum_s2=1-s1-1(Dt*rk_c2(s1,s2)*K_s2))
    do b=1,global%mesh%Nb
      call residuals(s1=s1,global=global,block=block(b))
    enddo
  enddo
  ! Runge-Kutta time integration
  do b=1,global%mesh%Nb
    call rk_time_integration(global=global,block=block(b),RU=RU(:,b))
  enddo
  ! finding the maximum value of residuals of actual process
  mRU = maxval(RU,dim=2)
#ifdef MPI2
  ! for multi-processes simulation all processes must exchange their mRU for computing the global gmRU
  call MPI_ALLREDUCE(mRU,gmRU,global%fluid%Nc,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,err)
#else
  ! for single processes mRU is already global gmRU
  gmRU = mRU
#endif

  ! updating the log file of residuals
  if (myrank==0) then
    if ((mod(global%fluid%n,global%file%probe_out)==0).OR. &
       (global%fluid%t==global%fluid%Tmax).OR. &
       (global%fluid%n==global%fluid%Nmax).OR. &
       (global%fluid%n==1)) then
      write(global%file%unit_res,trim(global%file%varform_res))global%fluid%n,global%fluid%t,(gmRU(b),b=1,global%fluid%Nc)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine solve_grl

  subroutine rk_init(S)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for initializing Runge-Kutta coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: S ! Number of stages used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! allocating variables
  if (allocated(rk_c1)) deallocate(rk_c1) ; allocate(rk_c1(    1:S)) ; rk_c1 = 0._R_P
  if (allocated(rk_c2)) deallocate(rk_c2) ; allocate(rk_c2(1:S,1:S)) ; rk_c2 = 0._R_P
  ! inizializing the coefficients
  select case(S)
  case(1_I1P)
    ! RK(1,1) Forward-Euler
    rk_c1(1) = 1._R_P
  case(2_I1P)
    ! SSPRK(2,2)
    rk_c1(1) = 0.5_R_P
    rk_c1(2) = 0.5_R_P

    rk_c2(2,1) = 1._R_P
  case(3_I1P)
    ! SSPRK(3,3)
    rk_c1(1) = 1._R_P/6._R_P
    rk_c1(2) = 1._R_P/6._R_P
    rk_c1(3) = 2._R_P/3._R_P

    rk_c2(2,1) = 1._R_P
    rk_c2(3,1) = 0.25_R_P ; rk_c2(3,2) = 0.25_R_P
  case(4_I1P)
    ! NSSPRK(4,4)
    rk_c1(1) = 55._R_P/108._R_P
    rk_c1(2) = 1._R_P/3._R_P
    rk_c1(3) = 1._R_P/3._R_P
    rk_c1(4) = 1._R_P/6._R_P

    rk_c2(2,1) = 1._R_P
    rk_c2(3,1) = 0.25_R_P       ; rk_c2(3,2) = 0.5_R_P
    rk_c2(4,1) = 5._R_P/18._R_P ; rk_c2(4,2) = 0._R_P  ; rk_c2(4,3) = 1._R_P
  case(5_I1P)
    ! SSPRK(5,4)
    rk_c1(1) = 0.14681187618661_R_P
    rk_c1(2) = 0.24848290924556_R_P
    rk_c1(3) = 0.10425883036650_R_P
    rk_c1(4) = 0.27443890091960_R_P
    rk_c1(5) = 0.22600748319395_R_P

    rk_c2(2,1)=0.39175222700392_R_P
    rk_c2(3,1)=0.21766909633821_R_P;rk_c2(3,2)=0.36841059262959_R_P
    rk_c2(4,1)=0.08269208670950_R_P;rk_c2(4,2)=0.13995850206999_R_P;rk_c2(4,3)=0.25189177424738_R_P
    rk_c2(5,1)=0.06796628370320_R_P;rk_c2(5,2)=0.11503469844438_R_P;rk_c2(5,3)=0.20703489864929_R_P;rk_c2(5,4)=0.54497475021237_R_P

    ! NSSPRK(5,3)
    !rk_c1(1) = 0.25_R_P
    !rk_c1(5) = 0.75_R_P

    !rk_c2(2,1)=1._R_P/7._R_P
    !rk_c2(3,1)=0._R_P        ; rk_c2(3,2)=3._R_P/16._R_P
    !rk_c2(4,1)=0._R_P        ; rk_c2(4,2)=0._R_P         ; rk_c2(4,3)=1._R_P/3._R_P
    !rk_c2(5,1)=0._R_P        ; rk_c2(5,2)=0._R_P         ; rk_c2(5,3)=0._R_P        ; rk_c2(5,4)=2._R_P/3._R_P
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rk_init

  subroutine rk_stages_sum(s1,global,block)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for summing Runge-Kutta stages for updating primitive variables (block%fluid%P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    s1      ! Current Runge-kutta stage.
  type(Type_Global), intent(IN)::    global  ! Global-level data.
  type(Type_Block),  intent(INOUT):: block   ! Block-level data.
  type(Type_Conservative)::          U,rksum ! Dummy conservative varibales.
  integer(I_P)::                     i,j,k,s ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE)   &
  !$OMP PRIVATE(i,j,k,s,rksum,U) &
  !$OMP SHARED(s1,global,block,rk_c2)
  call init_cons(Ns = global%fluid%Ns, cons = rksum)
  call init_cons(Ns = global%fluid%Ns, cons = U    )
  !$OMP DO
  do k=1,block%mesh%Nk
    do j=1,block%mesh%Nj
      do i=1,block%mesh%Ni
        rksum = 0._R_P
        do s=1,s1-1
          rksum = rksum + rk_c2(s1,s)*block%fluid%KS(i,j,k,s)
        enddo
        U = block%fluid%U(i,j,k) + block%fluid%Dt(i,j,k)*rksum
        call cons2prim(cp0 = global%fluid%cp0, cv0 = global%fluid%cv0, cons = U, prim = block%fluid%P(i,j,k))
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rk_stages_sum

  subroutine rk_time_integration(global,block,RU)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing Runge-Kutta integration.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(IN)::    global                ! Global-level data.
  type(Type_Block),  intent(INOUT):: block                 ! Block-level data.
  real(R_P),         intent(OUT)::   RU(1:global%fluid%Nc) ! NormL2 of residuals of conservative variables.
  type(Type_Conservative)::          rksum                 ! Dummy conservative varibale.
  integer(I_P)::                     i,j,k,s               ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  RU = 0._R_P
  !$OMP PARALLEL DEFAULT(NONE)     &
  !$OMP PRIVATE(i,j,k,s,rksum)     &
  !$OMP SHARED(global,block,rk_c1) &
  !$OMP REDUCTION(+: RU)
  call init_cons(Ns = global%fluid%Ns, cons = rksum)
  !$OMP DO
  do k=1,block%mesh%Nk
    do j=1,block%mesh%Nj
      do i=1,block%mesh%Ni
        rksum = 0._R_P
        do s=1,global%fluid%rk_ord
          rksum = rksum + rk_c1(s)*block%fluid%KS(i,j,k,s)
        enddo
        block%fluid%U(i,j,k) = block%fluid%U(i,j,k) + block%fluid%Dt(i,j,k)*rksum
        RU = RU + cons2array(rksum*rksum)
      enddo
    enddo
  enddo
  !$OMP DO
  do s=1,global%fluid%Nc
    RU(s) = sqrt(RU(s))
  enddo
  !$OMP END PARALLEL
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rk_time_integration
endmodule Lib_Fluidynamic
