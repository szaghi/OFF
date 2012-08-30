!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Fluxes_ConvectivePublicProcedure Lib_Fluxes_Convective
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Fluxes_ConvectivePrivateProcedure Lib_Fluxes_Convective
!> @}

!> This module contains the definition of procedures for computing convective fluxes (hyperbolic operator).
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @todo \b RotatedRieman: Implement Rotated Rieman Solver technique
!> @ingroup Library
module Lib_Fluxes_Convective
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                              ! Integers and reals precision definition.
USE Data_Type_Conservative, init_cons=>init, set_cons=>set    ! Definition of Type_Conservative.
USE Data_Type_Primitive, init_prim=>init, set_prim=>set       ! Definition of Type_Primitive.
USE Data_Type_Vector, set_vec=>set                            ! Definition of Type_Vector.
USE Lib_Riemann, only: &
#ifdef SMSWliu
                       chk_smooth_liu,        &               ! Liu's smoothness function.
#elif defined SMSWvanleer || defined SMSWvanalbada || defined SMSWharten
                       chk_smooth_limiter,    &               ! Slope limiter-based smoothness function.
#else
                       chk_smooth_z,          &               ! Riemann-Solver-like smoothness function.
#endif
#if defined RSHLLCb || defined RSHLLCc || defined RSHLLCp || defined RSHLLCt || defined RSHLLCz
                       Riem_Solver=>Riem_Solver_HLLC          ! HLLC Riemann solver.
#elif RSEXA
                       Riem_Solver=>Riem_Solver_Exact_U       ! Exact Riemann solver.
#elif RSPVL
                       Riem_Solver=>Riem_Solver_PVL           ! PVL Riemann solver.
#elif RSTR
                       Riem_Solver=>Riem_Solver_TR            ! TR Riemann solver.
#elif RSTS
                       Riem_Solver=>Riem_Solver_TS            ! TS Riemann solver.
#elif RSAPRS
                       Riem_Solver=>Riem_Solver_APRS          ! APRS Riemann solver.
#elif RSALFR
                       Riem_Solver=>Riem_Solver_ALFR          ! ALFR Riemann solver.
#elif defined RSLFp || defined RSLFz
                       Riem_Solver=>Riem_Solver_LaxFriedrichs ! Lax-Friedrichs Riemann solver.
#elif RSROE
                       Riem_Solver=>Riem_Solver_Roe           ! Roe Riemann solver.
#else
                       Riem_Solver=>Riem_Solver_HLLC          ! HLLC Riemann solver.
#endif
USE Lib_WENO,    only: &
#ifdef HYBRIDC
                       noweno_central, &                      ! Function for comp. central diff. reconstruction.
#elif HYBRID
                       weno_optimal,   &                      ! Function for comp. optimal WENO reconstruction.
#endif
                       weno                                   ! Function for computing WENO reconstruction.
USE Lib_Thermodynamic_Laws_Ideal, only: a                     ! Function for computing speed of sound.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: fluxes_convective
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> Subroutine for computing interfaces convective fluxes (hyperbolic operator). Taking in input the primitive variables along a
  !> 1D slice of N+2*gc cells the subroutine returns in ouput the N+1 convective fluxes at cells interfaces.
  !> @note The 3D primitive variables are projected in the 1D normal direction to the interface.
  !> The 1D remapping is as following:
  !>  - 1)    density of species 1    (r1)
  !>  - 2)    density of species 2    (r2)
  !>  - ...
  !>  - s)    density of species s-th (rs)
  !>  - ...
  !>  - Ns)   density of species Ns   (rNs)
  !>  - Ns+1) velocity                (interface's normal component)
  !>  - Ns+2) pressure                (p)
  !>  - Ns+3) density                 (r=sum(rs))
  !>  - Ns+4) specific heats ratio    (g)
  subroutine fluxes_convective(gc,N,Ns,cp0,cv0,NF,P,F)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN)::    gc                         !< Number of ghost cells used.
  integer(I_P),            intent(IN)::    N                          !< Number of cells.
  integer(I_P),            intent(IN)::    Ns                         !< Number of species.
  real(R_P),               intent(IN)::    cp0(1:Ns)                  !< Initial specific heat cp.
  real(R_P),               intent(IN)::    cv0(1:Ns)                  !< Initial specific heat cv.
  type(Type_Vector),       intent(IN)::    NF (           0-gc:N+gc ) !< Interface normal.
  type(Type_Primitive),    intent(IN)::    P  (           1-gc:N+gc ) !< Primitive variables (3D format).
  type(Type_Conservative), intent(INOUT):: F  (              0:N    ) !< Convective fluxes (3D format).
  type(Type_Vector)::                      ut (       1:2,1-gc:N+gc ) !< Tangential velocity: left (1) and right (2) interfaces.
  real(R_P)::                              P1D(1:Ns+2,1:2,1-gc:-1+gc) !< Primitive variables 1D fixed to the current interface.
#ifdef RECVC
  real(R_P)::                              Pm (1:Ns+4,       0:N    ) !< Mean of primitive variables across interfaces.
  real(R_P)::                              LPm(1:Ns+2,1:Ns+2,0:N    ) !< Mean left eigenvectors matrix.
  real(R_P)::                              RPm(1:Ns+2,1:Ns+2,0:N    ) !< Mean right eigenvectors matrix.
  real(R_P)::                              C  (1:Ns+2,1:2,1-gc:-1+gc) !< Interface value of characteristic variables.
  real(R_P)::                              CR (1:Ns+2,1:2           ) !< Reconstructed interface value of charac. variables.
#endif
#ifdef LMA
  real(R_P)::                              ulma(1:2)                  !< Left (1) and right (2) reconstructed normal speed adjusted
                                                                      !< for low Mach number flows.
  real(R_P)::                              z                          !< Scaling coefficient for Low Mach number Adjustment.
#endif
  real(R_P)::                              PR (1:Ns+4,1:2,   0:N+1  ) !< Reconstructed interface values of primitive variables 1D.
  real(R_P)::                              F_r                        !< Flux of mass (1D).
  real(R_P)::                              F_u                        !< Flux of momentum (1D).
  real(R_P)::                              F_E                        !< Flux of energy (1D).
  integer(I_P)::                           i,j,k,v                    !< Counters.
  logical::                                ROR(1:2)                   !< Logical flag for testing the result of WENO reconst.
  integer(I_P)::                           or                         !< Counter of order for ROR algorithm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
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
          PR(v,1:2,i) = weno(S=or,V=P1D(v,1:2,1-or:-1+or))
        enddo
#endif
        ! extrapolation of the reconstructed values for the non computed boundaries
        if (i==0  ) PR(1:Ns+2,1,i) = PR(1:Ns+2,2,i)
        if (i==N+1) PR(1:Ns+2,2,i) = PR(1:Ns+2,1,i)
#ifdef PPL
        ! applaying the maximum-principle-satisfying limiter to the reconstructed densities and pressure
        do v=1,Ns
          call positivity_preserving_limiter(o = or, vmean = P(i)%r(v), vr = PR(v,1:2,i)) ! densities
        enddo
        call positivity_preserving_limiter(o = or, vmean = P(i)%p, vr = PR(Ns+2,1:2,i))   ! pressure
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
#ifdef LMA
  ! applying Low Mach number Adjustment for decreasing the dissipation in local low Mach number region
  do i=0,N
    z = min(1._R_P,max(sqrt(PR(Ns+1,2,i  )*PR(Ns+1,2,i  )+sq_norm(ut(2,i  )))/a(p=P(i  )%p,r=P(i  )%d,g=P(i  )%g),&
                       sqrt(PR(Ns+1,1,i+1)*PR(Ns+1,1,i+1)+sq_norm(ut(1,i+1)))/a(p=P(i+1)%p,r=P(i+1)%d,g=P(i+1)%g)))
    ulma(1) = 0.5_R_P*(PR(Ns+1,2,i) + PR(Ns+1,1,i+1)) + z*0.5_R_P*(PR(Ns+1,2,i) - PR(Ns+1,1,i+1))
    ulma(2) = 0.5_R_P*(PR(Ns+1,2,i) + PR(Ns+1,1,i+1)) - z*0.5_R_P*(PR(Ns+1,2,i) - PR(Ns+1,1,i+1))
    PR(Ns+1,2,i  ) = ulma(1)
    PR(Ns+1,1,i+1) = ulma(2)
  enddo
#endif
  ! solving Riemann Problems
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
                     F_E = F_E)
    ! uptdating fluxes with tangential components
    if (F_r>0._R_P) then
      F(i)%rs = F_r*PR(1:Ns,2,i  )/PR(Ns+3,2,i  )
      F(i)%rv = F_u*NF(i) + F_r*ut(2,i)
      F(i)%re = F_E + F_r*sq_norm(ut(2,i))*0.5_R_P
    else
      F(i)%rs = F_r*PR(1:Ns,1,i+1)/PR(Ns+3,1,i+1)
      F(i)%rv = F_u*NF(i) + F_r*ut(1,i+1)
      F(i)%re = F_E + F_r*sq_norm(ut(1,i+1))*0.5_R_P
    endif
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> Function for computing left eigenvectors from primitive variables.
    pure function LP(Ns,primitive) result(L)
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

    !> Function for computing right eigenvectors from primitive variables.
    pure function RP(Ns,primitive) result(R)
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

    !> Positivity preserving limiter.
    pure subroutine positivity_preserving_limiter(o,vmean,vr)
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
  endsubroutine fluxes_convective
endmodule Lib_Fluxes_Convective

 !subroutine convective_fluxes(gc,N,Ns,cp0,cv0,NF,P,F)
  !!-------------------------------------------------------------------------------------------------------------------------------
  !implicit none
  !integer(I1P),            intent(IN)::    gc                         !< Number of ghost cells used.
  !integer(I_P),            intent(IN)::    N                          !< Number of cells.
  !integer(I_P),            intent(IN)::    Ns                         !< Number of species.
  !real(R_P),               intent(IN)::    cp0(1:Ns)                  !< Initial specific heat cp.
  !real(R_P),               intent(IN)::    cv0(1:Ns)                  !< Initial specific heat cv.
  !type(Type_Vector),       intent(IN)::    NF (           0-gc:N+gc ) !< Face normal.
  !type(Type_Primitive),    intent(IN)::    P  (           1-gc:N+gc ) !< Primitive variables (3D).
  !type(Type_Conservative), intent(INOUT):: F  (              0:N    ) !< Flux of mass (3D).
  !type(Type_Vector)::                      ut (       1:2,1-gc:N+gc ) !< Tangential velocity: left (1) and right (2) interfaces.
  !real(R_P)::                              P1D(1:Ns+2,1:2,1-gc:-1+gc) !< Primitive variables 1D fixed to the current interface.
!#ifdef RECVC
  !real(R_P)::                              Pm (1:Ns+4,       0:N    ) !< Mean of primitive variables across interfaces.
  !real(R_P)::                              LPm(1:Ns+2,1:Ns+2,0:N    ) !< Mean left eigenvectors matrix.
  !real(R_P)::                              RPm(1:Ns+2,1:Ns+2,0:N    ) !< Mean right eigenvectors matrix.
  !real(R_P)::                              C  (1:Ns+2,1:2,1-gc:-1+gc) !< Interface value of characteristic variables.
  !real(R_P)::                              CR (1:Ns+2,1:2           ) !< Reconstructed interface value of charac. variables.
!#endif
  !real(R_P)::                              PR (1:Ns+4,1:2,   0:N+1  ) !< Reconstructed interface values of primitive variables 1D.
  !real(R_P)::                              F_r                        !< Flux of mass (1D).
  !real(R_P)::                              F_u                        !< Flux of momentum (1D).
  !real(R_P)::                              F_e                        !< Flux of energy (1D).
  !integer(I_P)::                           i,j,k,v                    !< Counters.
  !logical::                                ROR(1:2)                   !< Logical flag for testing the result of WENO reconst.
  !integer(I_P)::                           or                         !< Counter of order for ROR algorithm.
!#if defined HYBRID || HYBRIDC
  !logical::                                smooth  (1-gc:N+gc  )      !< Logical flag for testing the smoothness of the stencils.
  !logical::                                smooth_f(1-gc:N+gc-1)      !< Logical flag for testing the smoothness of the stencils.
!#endif
  !!-------------------------------------------------------------------------------------------------------------------------------

  !!-------------------------------------------------------------------------------------------------------------------------------
  !!$OMP PARALLEL DEFAULT(NONE)                 &
  !!$OMP PRIVATE(i,j,k,v,or,ROR,P1D,F_r,F_u,F_e &
!#ifdef RECVC
  !!$OMP         ,C,CR)                         &
!#else
  !!$OMP         )                              &
!#endif
  !!$OMP SHARED(gc,N,Ns,NF,ut,cp0,cv0,P,PR,F    &
!#ifdef RECVC
  !!$OMP        ,Pm,LPm,RPm                     &
!#endif
!#if defined HYBRID || HYBRIDC
  !!$OMP        ,smooth,smooth_f)
!#else
  !!$OMP        )
!#endif

  !! computing velocity vectors and its tangential component
  !!$OMP DO
  !do i=1-gc,N+gc
    !! computing tangential velocity component
    !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      !if (i==1-gc.AND.j==1) cycle
      !if (i==N+gc.AND.j==2) cycle
      !ut(j,i) = P(i)%v - (P(i)%v.paral.NF(i+j-2)) ! ut
    !enddo
  !enddo

  !! computing the left and right states or Riemann Problems
  !select case(gc)
  !case(2_I1P,3_I1P,4_I1P) ! 3rd, 5th or 7th order WENO reconstruction
!#if defined HYBRID || HYBRIDC
    !! checking the smoothness of the interface if high-order method is used
!#ifdef SMSWliu
    !!$OMP DO
    !do i=0,N+1 ! loop over cells
      !smooth(i) = chk_smooth_liu(pm1  =             P(i-1)%p,               &
                                 !rm1  =             P(i-1)%d,               &
                                 !um1  =            (P(i-1)%v.dot.NF(i-1)),  &
                                 !gm1  =             P(i-1)%g,               &
                                 !cpm1 = dot_product(P(i-1)%r/P(i-1)%d,cp0), &
                                 !cvm1 = dot_product(P(i-1)%r/P(i-1)%d,cv0), &
                                 !p1   =             P(i  )%p,               &
                                 !r1   =             P(i  )%d,               &
                                 !u1m  =            (P(i  )%v.dot.NF(i-1)),  &
                                 !u1p  =            (P(i  )%v.dot.NF(i  )),  &
                                 !g1   =             P(i  )%g,               &
                                 !cp1  = dot_product(P(i  )%r/P(i  )%d,cp0), &
                                 !cv1  = dot_product(P(i  )%r/P(i  )%d,cv0), &
                                 !pp1  =             P(i+1)%p,               &
                                 !rp1  =             P(i+1)%d,               &
                                 !up1  =            (P(i+1)%v.dot.NF(i  )),  &
                                 !gp1  =             P(i+1)%g,               &
                                 !cpp1 = dot_product(P(i+1)%r/P(i+1)%d,cp0), &
                                 !cvp1 = dot_product(P(i+1)%r/P(i+1)%d,cv0))
    !enddo
!#elif defined SMSWvanalbada || defined SMSWvanleer || defined SMSWharten
    !!$OMP DO
    !do i=0,N+1 ! loop over cells
      !smooth(i) = chk_smooth_limiter(pm1 =  P(i  )%p              -  P(i-1)%p,              &
                                     !pp1 =  P(i+1)%p              -  P(i  )%p,              &
                                     !rm1 =  P(i  )%d              -  P(i-1)%d,              &
                                     !rp1 =  P(i+1)%d              -  P(i  )%d,              &
                                     !um1 = (P(i  )%v.dot.NF(i-1)) - (P(i-1)%v.dot.NF(i-1)), &
                                     !up1 = (P(i+1)%v.dot.NF(i  )) - (P(i  )%v.dot.NF(i  )), &
                                     !gm1 =  P(i  )%g              -  P(i-1)%g,              &
                                     !gp1 =  P(i+1)%g              -  P(i  )%g)
    !enddo
!#elif defined SMSWz
    !!$OMP DO
    !do i=1-gc,N+gc-1 ! loop over interfaces
      !smooth_f(i) = chk_smooth_z(p1 =  P(i  )%p,            &
                                 !r1 =  P(i  )%d,            &
                                 !u1 = (P(i  )%v.dot.NF(i)), &
                                 !g1 =  P(i  )%g,            &
                                 !p4 =  P(i+1)%p,            &
                                 !r4 =  P(i+1)%d,            &
                                 !u4 = (P(i+1)%v.dot.NF(i)), &
                                 !g4 =  P(i+1)%g)
    !enddo
    !!$OMP DO
    !do i=1-gc+1,N+gc-1 ! loop over cells
      !smooth(i) = (smooth_f(i-1).AND.smooth_f(i))
    !enddo
    !!$OMP SINGLE
    !smooth(1-gc)=smooth(1-gc+1)
    !smooth(N+gc)=smooth(N+gc-1)
    !!$OMP END SINGLE
!#endif
    !!$OMP SINGLE
    !do i=1-gc,N+gc ! loop over cells
    !enddo
    !!$OMP END SINGLE
!#endif

!#ifdef RECVC
    !! computing Pm, LPm and RPm across the interfaces  for the local characteristic projection
    !!$OMP DO
    !do i=0,N ! loop over interfaces
      !! computing mean of primitive variables across the interface
      !Pm(1:Ns,i) = 0.5_R_P*( P(i)%r(1:NS)      +  P(i+1)%r(1:NS)     ) ! rs
      !Pm(Ns+1,i) = 0.5_R_P*((P(i)%v.dot.NF(i)) + (P(i+1)%v.dot.NF(i))) ! un
      !Pm(Ns+2,i) = 0.5_R_P*( P(i)%p            +  P(i+1)%p           ) ! p
      !Pm(Ns+3,i) = 0.5_R_P*( P(i)%d            +  P(i+1)%d           ) ! r
      !Pm(Ns+4,i) = 0.5_R_P*( P(i)%g            +  P(i+1)%g           ) ! g
      !! computing mean left and right eigenvectors matrix across the interface
      !LPm(1:Ns+2,1:Ns+2,i) = LP(Ns,Pm(1:Ns+4,i))
      !RPm(1:Ns+2,1:Ns+2,i) = RP(Ns,Pm(1:Ns+4,i))
    !enddo
!#endif

    !! computing the reconstruction
    !!$OMP DO
    !do i=0,N+1 ! loop over cells
      !! computing 1D primitive variables for the stencil [i+1-gc,i-1+gc] fixing the velocity across the interfaces i+-1/2
      !do k=i+1-gc,i-1+gc
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !P1D(1:Ns,j,k-i) =  P(k)%r(1:Ns)
          !P1D(Ns+1,j,k-i) = (P(k)%v.dot.NF(i+j-2))
          !P1D(Ns+2,j,k-i) =  P(k)%p
        !enddo
      !enddo

!#ifdef RECVC
      !! transforming variables into local characteristic fields for the stencil [1-gc,-1+gc]
      !do k=i+1-gc,i-1+gc
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !do v=1,Ns+2
            !C(v,j,k-i) = dot_product(LPm(v,1:Ns+2,i+j-2),P1D(1:Ns+2,j,k-i))
          !enddo
        !enddo
      !enddo
!#endif

!#ifdef HYBRID
      !! if stencil is smooth the WENO algorithm is computed directly with optimal weights
      !if (all(smooth(i+1-gc:i-1+gc))) then
        !! computing WENO reconstruction with ideal weights (without the computation of smoothness indicator)
!#ifdef RECVC
        !do v=1,Ns+2
          !CR(v,1:2) = weno_optimal(S=gc,V=C(v,1:2,1-gc:-1+gc))
        !enddo
        !! trasforming local reconstructed characteristic variables to primitive ones
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !do v=1,Ns+2
            !PR(v,j,i) = dot_product(RPm(v,1:Ns+2,i+j-2),CR(1:Ns+2,j))
          !enddo
          !PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
          !PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
        !enddo
!#else
        !do v=1,Ns+2
          !PR(v,1:2,i) = weno_optimal(S=gc,V=P1D(v,1:2,1-gc:-1+gc))
        !enddo
        !! computing the last reconstructed variables
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
          !PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
        !enddo
!#endif
        !cycle ! the full WENO reconstructions are not necessary
      !endif
!#elif HYBRIDC
      !! if stencil is smooth the WENO algorithm is switched to the central difference reconstruction
      !if (all(smooth(i+1-gc:i-1+gc))) then
        !! computing central difference reconstruction
!#ifdef RECVC
        !do v=1,Ns+2
          !CR(v,1:2) = noweno_central(S=gc,V=C(v,1:2,1-gc:-1+gc))
        !enddo
        !! trasforming local reconstructed characteristic variables to primitive ones
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !do v=1,Ns+2
            !PR(v,j,i) = dot_product(RPm(v,1:Ns+2,i+j-2),CR(1:Ns+2,j))
          !enddo
          !PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
          !PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
        !enddo
!#else
        !do v=1,Ns+2
          !PR(v,1:2,i) = noweno_central(S=gc,V=P1D(v,1:2,1-gc:-1+gc))
        !enddo
        !! computing the last reconstructed variables
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
          !PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
        !enddo
!#endif
        !cycle ! the full WENO reconstructions are not necessary
      !endif
!#endif

      !! computing WENO reconstruction with ROR (Recursive Order Reduction) algorithm
      !ROR_check: do or=gc,2,-1
        !! computing WENO reconstruction
!#ifdef RECVC
        !do v=1,Ns+2
          !CR(v,1:2) = weno(S=or,V=C(v,1:2,1-or:-1+or))
        !enddo
        !! trasforming local reconstructed characteristic variables to primitive ones
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !do v=1,Ns+2
            !PR(v,j,i) = dot_product(RPm(v,1:Ns+2,i+j-2),CR(1:Ns+2,j))
          !enddo
        !enddo
!#else
        !do v=1,Ns+2
          !PR(v,1:2,i) = weno(S=gc,V=P1D(v,1:2,1-gc:-1+gc))
        !enddo
!#endif
        !! extrapolation of the reconstructed values for the non computed boundaries
        !if (i==0  ) PR(1:Ns+2,1,i) = PR(1:Ns+2,2,i)
        !if (i==N+1) PR(1:Ns+2,2,i) = PR(1:Ns+2,1,i)

!#ifdef PPL
        !! applaying the maximum-principle-satisfying limiter to the reconstructed densities and pressure
        !do v=1,Ns
          !call positivity_preserving_limiter(o = gc, vmean = P(i)%r(v), vr = PR(v,1:2,i)) ! densities
        !enddo
        !call positivity_preserving_limiter(o = gc, vmean = P(i)%p, vr = PR(Ns+2,1:2,i))   ! pressure
!#endif

        !! computing the reconstructed variables of total density and specific heats ratio
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !PR(Ns+3,j,i) = sum(PR(1:Ns,j,i))
          !PR(Ns+4,j,i) = dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cp0(1:Ns))/dot_product(PR(1:Ns,j,i)/PR(Ns+3,j,i),cv0(1:Ns))
        !enddo

        !! verifing the WENO reconstruction at order or-th
        !ROR = .true.
        !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          !if (i==0  .AND.j==1) cycle
          !if (i==N+1.AND.j==2) cycle
          !ROR(j) = ((PR(Ns+2,j,i)>0._R_P).AND.(PR(Ns+3,j,i)>0._R_P).AND.(PR(Ns+4,j,i)>0._R_P))
        !enddo
        !if (ROR(1).AND.ROR(2)) exit ROR_check
        !if (or==2) then
          !! the high order reconstruction is falied; used 1st order reconstruction
          !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            !if (i==0  .AND.j==1) cycle
            !if (i==N+1.AND.j==2) cycle
            !PR(1:Ns,j,i) =  P(i)%r(1:Ns)          ! rs
            !PR(Ns+1,j,i) = (P(i)%v.dot.NF(i+j-2)) ! un
            !PR(Ns+2,j,i) =  P(i)%p                ! p
            !PR(Ns+3,j,i) =  P(i)%d                ! r
            !PR(Ns+4,j,i) =  P(i)%g                ! g
          !enddo
        !endif
      !enddo ROR_check
    !enddo
  !case(1_I1P) ! 1st order piecewise constant reconstruction
    !!$OMP DO
    !do i=0,N+1
      !do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        !if (i==0  .AND.j==1) cycle
        !if (i==N+1.AND.j==2) cycle
        !PR(1:Ns,j,i) =  P(i)%r(1:Ns)          ! rs
        !PR(Ns+1,j,i) = (P(i)%v.dot.NF(i+j-2)) ! un
        !PR(Ns+2,j,i) =  P(i)%p                ! p
        !PR(Ns+3,j,i) =  P(i)%d                ! r
        !PR(Ns+4,j,i) =  P(i)%g                ! g
      !enddo
    !enddo
  !endselect

  !! solving Rimeann Problems
  !!$OMP DO
  !do i=0,N
    !! face normal fluxes
    !call Riem_Solver(p1  = PR(Ns+2,2,i  ), & !| right value of cell i   (letf  of interface i+1/2)
                     !r1  = PR(Ns+3,2,i  ), & !|
                     !u1  = PR(Ns+1,2,i  ), & !|
                     !g1  = PR(Ns+4,2,i  ), & !|
!#if defined RSHLLCb || defined RSROE
                     !cp1 = dot_product(PR(1:Ns,2,i)/PR(Ns+3,2,i),cp0), &
                     !cv1 = dot_product(PR(1:Ns,2,i)/PR(Ns+3,2,i),cv0), &
!#endif
                     !p4  = PR(Ns+2,1,i+1), & !| left  value of cell i+1 (right of interface i+1/2)
                     !r4  = PR(Ns+3,1,i+1), & !|
                     !u4  = PR(Ns+1,1,i+1), & !|
                     !g4  = PR(Ns+4,1,i+1), & !|
!#if defined RSHLLCb || defined RSROE
                     !cp4 = dot_product(PR(1:Ns,1,i+1)/PR(Ns+3,1,i+1),cp0), &
                     !cv4 = dot_product(PR(1:Ns,1,i+1)/PR(Ns+3,1,i+1),cv0), &
!#endif
                     !F_r = F_r           , &
                     !F_u = F_u           , &
                     !F_E = F_e)
    !! uptdating fluxes with tangential components
    !if (F_r>0._R_P) then
      !F(i)%rs = F_r*PR(1:Ns,2,i  )/PR(Ns+3,2,i  )
      !F(i)%rv = F_u*NF(i) + F_r*ut(2,i)
      !F(i)%re = F_e + F_r*sq_norm(ut(2,i))*0.5_R_P
    !else
      !F(i)%rs = F_r*PR(1:Ns,1,i+1)/PR(Ns+3,1,i+1)
      !F(i)%rv = F_u*NF(i) + F_r*ut(1,i+1)
      !F(i)%re = F_e + F_r*sq_norm(ut(1,i+1))*0.5_R_P
    !endif
  !enddo

  !!$OMP END PARALLEL
  !return
  !!-------------------------------------------------------------------------------------------------------------------------------
  !endsubroutine convective_fluxes
