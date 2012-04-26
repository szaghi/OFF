module Lib_Multigrid
!-----------------------------------------------------------------------------------------------------------------------------------
!! The module Lib_Multigrid contains the definition of multigrid functions and subroutines.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                       ! Integers and reals precision definition.
USE Data_Type_BC, init_bc => init, set_bc => set       ! Definition of Type_BC.
USE Data_Type_Vector, init_vec => init, set_vec => set ! Definition of Type_Vector.
USE Lib_Fluidynamic, only: conservative2primitive, &   ! Function for converting conservative variables to primitive ones.
                           residuals,              &   ! Subroutine for computing conservative variables residuals.
                           boundary_conditions         ! Subroutine for setting boundary conditions.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: UR
public:: TE
public:: alloc_multigrid
public:: restrict
public:: correct
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Multigrid variables
real(R_P), allocatable:: UR(:,:) ! Restriction of conservative variables      [1:Nc,bb1:bb2].
real(R_P), allocatable:: TE(:,:) ! Truncation error of conservative variables [1:Nc,bb1:bb2].
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine alloc_multigrid(Nc,bb1,bb2)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for safe allocation of multigrid algorithm variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: Nc      ! Number of conservative variables.
  integer(I_P), intent(IN):: bb1,bb2 ! Bounds of allocation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(UR)) deallocate(UR) ; allocate(UR(1:Nc,bb1:bb2)) ; UR = 0._R_P
  if (allocated(TE)) deallocate(TE) ; allocate(TE(1:Nc,bb1:bb2)) ; TE = 0._R_P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_multigrid

  subroutine restrict(myrank,l,Nb,gco,gc,Ni_f,Nj_f,Nk_f,Ni_c,Nj_c,Nk_c,bb_f,bb_c,bb_Fi_c,bb_Fj_c,bb_Fk_c,Np,Nc,Ns, &
                      NFi_c,NFj_c,NFk_c,BCi_c,BCj_c,BCk_c,Nin1,in1,Nin2,bb_in2,in2,Si_c,Sj_c,Sk_c,                 &
                      V_f,V_c,                                                                                     &
                      U_f,R_f,                                                                                     &
                      t,U_c,                                                                                       &
                      UR,TE)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for the restriction of the fine solution on the coarse one and for the estimation of the truncation error.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::  myrank                                  ! Actual rank process.
  integer(I_P),      intent(IN)::  l                                       ! Actual grid level.
  integer(I_P),      intent(IN)::  Nb                                      ! Number of blocks.
  integer(I_P),      intent(IN)::  gco                                     ! Number of ghost cells necessary for achieving the
                                                                           ! space order.
  integer(I_P),      intent(IN)::  gc(1:6,1:Nb)                            ! Number of ghost cells.
  integer(I_P),      intent(IN)::  Ni_f(1:Nb)                              ! Number of cells in i direction.
  integer(I_P),      intent(IN)::  Nj_f(1:Nb)                              ! Number of cells in j direction.
  integer(I_P),      intent(IN)::  Nk_f(1:Nb)                              ! Number of cells in k direction.
  integer(I_P),      intent(IN)::  Ni_c(1:Nb)                              ! Number of cells in i direction.
  integer(I_P),      intent(IN)::  Nj_c(1:Nb)                              ! Number of cells in j direction.
  integer(I_P),      intent(IN)::  Nk_c(1:Nb)                              ! Number of cells in k direction.
  integer(I_P),      intent(IN)::  bb_f   (1:2,1:Nb)                       ! Bounds of cell centers.
  integer(I_P),      intent(IN)::  bb_c   (1:2,1:Nb)                       ! Bounds of cell centers.
  integer(I_P),      intent(IN)::  bb_Fi_c(1:2,1:Nb)                       ! Bounds of interfaces in i direction.
  integer(I_P),      intent(IN)::  bb_Fj_c(1:2,1:Nb)                       ! Bounds of interfaces in j direction.
  integer(I_P),      intent(IN)::  bb_Fk_c(1:2,1:Nb)                       ! Bounds of interfaces in k direction.
  integer(I_P),      intent(IN)::  Np                                      ! Number of primitive variables.
  integer(I_P),      intent(IN)::  Nc                                      ! Number of conservative variables.
  integer(I_P),      intent(IN)::  Ns                                      ! Number of species.
  type(Type_Vector), intent(IN)::  NFi_c (     bb_Fi_c(1,1):bb_Fi_c(2,Nb)) ! Face i normal.
  type(Type_Vector), intent(IN)::  NFj_c (     bb_Fj_c(1,1):bb_Fj_c(2,Nb)) ! Face j normal.
  type(Type_Vector), intent(IN)::  NFk_c (     bb_Fk_c(1,1):bb_Fk_c(2,Nb)) ! Face k normal.
  type(Type_BC),     intent(IN)::  BCi_c (     bb_Fi_c(1,1):bb_Fi_c(2,Nb)) ! Bc of i interfaces.
  type(Type_BC),     intent(IN)::  BCj_c (     bb_Fj_c(1,1):bb_Fj_c(2,Nb)) ! Bc of j interfaces.
  type(Type_BC),     intent(IN)::  BCk_c (     bb_Fk_c(1,1):bb_Fk_c(2,Nb)) ! Bc of k interfaces.
  integer(I_P),      intent(IN)::  Nin1                                    ! Number of inflow 1 boundary conditions.
  real(R_P),         intent(IN)::  in1(1:Np,1:Nin1)                        ! Primitive variables of inflow 1 bc.
  integer(I_P),      intent(IN)::  Nin2                                    ! Number of inflow 2 boundary conditions.
  integer(I_P),      intent(IN)::  bb_in2(1:2,1:Nin2)                      ! Bounds of inflow 2 boundary conditions.
  real(R_P),         intent(IN)::  in2(1:Np+1,bb_in2(1,1):bb_in2(2,Nin2))  ! Primitive variables (plus time) of inflow 2 bc.
  real(R_P),         intent(IN)::  Si_c  (     bb_Fi_c(1,1):bb_Fi_c(2,Nb)) ! I surface area.
  real(R_P),         intent(IN)::  Sj_c  (     bb_Fj_c(1,1):bb_Fj_c(2,Nb)) ! J surface area.
  real(R_P),         intent(IN)::  Sk_c  (     bb_Fk_c(1,1):bb_Fk_c(2,Nb)) ! K surface area.
  real(R_P),         intent(IN)::  V_f   (     bb_f   (1,1):bb_f   (2,Nb)) ! Cell volume.
  real(R_P),         intent(IN)::  V_c   (     bb_c   (1,1):bb_c   (2,Nb)) ! Cell volume.
  real(R_P),         intent(IN)::  U_f   (1:Nc,bb_f   (1,1):bb_f   (2,Nb)) ! Conservative variables.
  real(R_P),         intent(IN)::  R_f   (1:Nc,bb_f   (1,1):bb_f   (2,Nb)) ! Residuals of conservative variables.
  real(R_P),         intent(IN)::  t                                       ! Time.
  real(R_P),         intent(OUT):: U_c   (1:Nc,bb_c   (1,1):bb_c   (2,Nb)) ! Conservative variables.
  real(R_P),         intent(OUT):: UR    (1:Nc,bb_c   (1,1):bb_c   (2,Nb)) ! Restriction of conservative variables.
  real(R_P),         intent(OUT):: TE    (1:Nc,bb_c   (1,1):bb_c   (2,Nb)) ! Truncation error of conservative variables.
  real(R_P)::                      P_c   (1:Np,bb_c   (1,1):bb_c   (2,Nb)) ! Primitive variables.
  real(R_P)::                      R_c   (1:Nc,bb_c   (1,1):bb_c   (2,Nb)) ! Residuals of conservative variables.
  integer(I_P)::                   b                                       ! Block counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do b=1,Nb
    ! restric the solution U
    call fine2coarse_U(gc = gc(1:6,b),                                 &
                       Ni_f = Ni_f(b), Nj_f = Nj_f(b), Nk_f = Nk_f(b), &
                       Ni_c = Ni_c(b), Nj_c = Nj_c(b), Nk_c = Nk_c(b), &
                       Nc = Nc,                                        &
                       V_f = V_f(     bb_f(1,b):bb_f(2,b)),            &
                       V_c = V_c(     bb_c(1,b):bb_c(2,b)),            &
                       U_f = U_f(1:Nc,bb_f(1,b):bb_f(2,b)),            &
                       U_c = U_c(1:Nc,bb_c(1,b):bb_c(2,b)))
    ! restric the residual R
    call fine2coarse_R(gc = gc(1:6,b),                                 &
                       Ni_f = Ni_f(b), Nj_f = Nj_f(b), Nk_f = Nk_f(b), &
                       Ni_c = Ni_c(b), Nj_c = Nj_c(b), Nk_c = Nk_c(b), &
                       Nc = Nc,                                        &
                       R_f = R_f(1:Nc,bb_f(1,b):bb_f(2,b)),            &
                       R_c = R_c(1:Nc,bb_c(1,b):bb_c(2,b)))
    ! converting conservative variables to primitive ones
    call conservative2primitive(gc   = gc(1:6,b),                     &
                                Ni   = Ni_c(  b),                     &
                                Nj   = Nj_c(  b),                     &
                                Nk   = Nk_c(  b),                     &
                                Np   = Np,                            &
                                Nc   = Nc,                            &
                                Ns   = Ns,                            &
                                cons = U_c(1:Nc,bb_c(1,b):bb_c(2,b)), &
                                prim = P_c(1:Np,bb_c(1,b):bb_c(2,b)))
  enddo
  ! imposing the boundary conditions
  do b=1,Nb
    call boundary_conditions(myrank = myrank,                                 &
                             l      = l,                                      &
                             b      = b,                                      &
                             gc     = gc(1:6,b),                              &
                             Ni     = Ni_c(  b),                              &
                             Nj     = Nj_c(  b),                              &
                             Nk     = Nk_c(  b),                              &
                             Np     = Np,                                     &
                             Ns     = Ns,                                     &
                             NFi    = NFi_c(     bb_Fi_c(1,b):bb_Fi_c(2,b)),  &
                             NFj    = NFj_c(     bb_Fj_c(1,b):bb_Fj_c(2,b)),  &
                             NFk    = NFk_c(     bb_Fk_c(1,b):bb_Fk_c(2,b)),  &
                             BCi    = BCi_c(     bb_Fi_c(1,b):bb_Fi_c(2,b)),  &
                             BCj    = BCj_c(     bb_Fj_c(1,b):bb_Fj_c(2,b)),  &
                             BCk    = BCk_c(     bb_Fk_c(1,b):bb_Fk_c(2,b)),  &
                             Nin1   = Nin1,                                   &
                             in1    = in1(1:Np,1:Nin1),                       &
                             Nin2   = Nin2,                                   &
                             bb_in2 = bb_in2(1:2,1:Nin2),                     &
                             in2    = in2(1:Np+1,bb_in2(1,1):bb_in2(2,Nin2)), &
                             t      = t,                                      &
                             P      = P_c  (1:Np,bb_c   (1,b):bb_c   (2,b)))
  enddo
  ! computing the residuals (saved temporary in TE) of restricted solution
  do b=1,Nb
    !call residuals(gco    = gco,                                   &
                   !gc     = gc(1:6,b),                             &
                   !Ni     = Ni_c(  b),                             &
                   !Nj     = Nj_c(  b),                             &
                   !Nk     = Nk_c(  b),                             &
                   !Np     = Np,                                    &
                   !Nc     = Nc,                                    &
                   !Ns     = Ns,                                    &
                   !NFi    = NFi_c(     bb_Fi_c(1,b):bb_Fi_c(2,b)), &
                   !NFj    = NFj_c(     bb_Fj_c(1,b):bb_Fj_c(2,b)), &
                   !NFk    = NFk_c(     bb_Fk_c(1,b):bb_Fk_c(2,b)), &
                   !Si     = Si_c (     bb_Fi_c(1,b):bb_Fi_c(2,b)), &
                   !Sj     = Sj_c (     bb_Fj_c(1,b):bb_Fj_c(2,b)), &
                   !Sk     = Sk_c (     bb_Fk_c(1,b):bb_Fk_c(2,b)), &
                   !V      = V_c  (     bb_c   (1,b):bb_c   (2,b)), &
                   !BCi    = BCi_c(     bb_Fi_c(1,b):bb_Fi_c(2,b)), &
                   !BCj    = BCj_c(     bb_Fj_c(1,b):bb_Fj_c(2,b)), &
                   !BCk    = BCk_c(     bb_Fk_c(1,b):bb_Fk_c(2,b)), &
                   !P      = P_c  (1:Np,bb_c   (1,b):bb_c   (2,b)), &
                   !Pmax   = Pmax (1:Np                          ), &
                   !Pmin   = Pmin (1:Np                          ), &
                   !R      = TE   (1:Nc,bb_c   (1,b):bb_c   (2,b)))
  enddo
  ! storing the restricted solution in UR and estimate the truncation error TE
  do b=1,Nb
    call te_estimation(gc = gc(1:6,b),                     &
                       Ni = Ni_c(  b),                     &
                       Nj = Nj_c(  b),                     &
                       Nk = Nk_c(  b),                     &
                       Nc = Nc,                            &
                       U  = U_c(1:Nc,bb_c(1,b):bb_c(2,b)), &
                       R  = R_c(1:Nc,bb_c(1,b):bb_c(2,b)), &
                       UR = UR (1:Nc,bb_c(1,b):bb_c(2,b)), &
                       TE = TE (1:Nc,bb_c(1,b):bb_c(2,b)))

  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    subroutine fine2coarse_U(gc,Ni_f,Nj_f,Nk_f,Ni_c,Nj_c,Nk_c,Nc,V_f,V_c,U_f,U_c)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for restricting the conservative variables on the coarse grid from the fine one.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN)::  gc(1:6)                                                            ! Ghost cells number.
    integer(I_P), intent(IN)::  Ni_f,Nj_f,Nk_f                                                     ! Cells number of fine grid.
    integer(I_P), intent(IN)::  Ni_c,Nj_c,Nk_c                                                     ! Cells number of coarse grid.
    integer(I_P), intent(IN)::  Nc                                                                 ! Number of cons. variables.
    real(R_P),    intent(IN)::  V_f(     1-gc(1):Ni_f+gc(2),1-gc(3):Nj_f+gc(4),1-gc(5):Nk_f+gc(6)) ! Fine grid finite volume.
    real(R_P),    intent(IN)::  V_c(     1-gc(1):Ni_c+gc(2),1-gc(3):Nj_c+gc(4),1-gc(5):Nk_c+gc(6)) ! Coarse grid finite volume.
    real(R_P),    intent(IN)::  U_f(1:Nc,1-gc(1):Ni_f+gc(2),1-gc(3):Nj_f+gc(4),1-gc(5):Nk_f+gc(6)) ! Fine grid cons. variables.
    real(R_P),    intent(OUT):: U_c(1:Nc,1-gc(1):Ni_c+gc(2),1-gc(3):Nj_c+gc(4),1-gc(5):Nk_c+gc(6)) ! Coarse grid cons. variables.
    integer(I_P)::              i,j,k,c                                                            ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    do k=1,Nk_c
      do j=1,Nj_c
        do i=1,Ni_c
          do c=1,Nc
            U_c(c,i,j,k) = (U_f(c,i*2-1,j*2-1,k*2-1)*V_f(i*2-1,j*2-1,k*2-1) + &
                            U_f(c,i*2-1,j*2  ,k*2-1)*V_f(i*2-1,j*2  ,k*2-1) + &
                            U_f(c,i*2-1,j*2-1,k*2  )*V_f(i*2-1,j*2-1,k*2  ) + &
                            U_f(c,i*2-1,j*2  ,k*2  )*V_f(i*2-1,j*2  ,k*2  ) + &
                            U_f(c,i*2  ,j*2-1,k*2-1)*V_f(i*2  ,j*2-1,k*2-1) + &
                            U_f(c,i*2  ,j*2  ,k*2-1)*V_f(i*2  ,j*2  ,k*2-1) + &
                            U_f(c,i*2  ,j*2-1,k*2  )*V_f(i*2  ,j*2-1,k*2  ) + &
                            U_f(c,i*2  ,j*2  ,k*2  )*V_f(i*2  ,j*2  ,k*2  ))/V_c(i,j,k)
          enddo
        enddo
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine fine2coarse_U

    subroutine fine2coarse_R(gc,Ni_f,Nj_f,Nk_f,Ni_c,Nj_c,Nk_c,Nc,R_f,R_c)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine for restricting the residuals of conservative variables on the coarse grid from the fine one.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN)::  gc(1:6)                                                            ! Ghost cells number.
    integer(I_P), intent(IN)::  Ni_f,Nj_f,Nk_f                                                     ! Cells number of fine grid.
    integer(I_P), intent(IN)::  Ni_c,Nj_c,Nk_c                                                     ! Cells number of coarse grid.
    integer(I_P), intent(IN)::  Nc                                                                 ! Number of cons. variables.
    real(R_P),    intent(IN)::  R_f(1:Nc,1-gc(1):Ni_f+gc(2),1-gc(3):Nj_f+gc(4),1-gc(5):Nk_f+gc(6)) ! Fine grid res. of cons. vars.
    real(R_P),    intent(OUT):: R_c(1:Nc,1-gc(1):Ni_c+gc(2),1-gc(3):Nj_c+gc(4),1-gc(5):Nk_c+gc(6)) ! Coarse grid res. of cons. vars.
    integer(I_P)::              i,j,k,c                                          ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    R_c = 0._R_P
    do k=1,Nk_c
      do j=1,Nj_c
        do i=1,Ni_c
          do c=1,Nc
            R_c(c,i,j,k) = R_f(c,i*2-1,j*2-1,k*2-1) + &
                           R_f(c,i*2-1,j*2  ,k*2-1) + &
                           R_f(c,i*2-1,j*2-1,k*2  ) + &
                           R_f(c,i*2-1,j*2  ,k*2  ) + &
                           R_f(c,i*2  ,j*2-1,k*2-1) + &
                           R_f(c,i*2  ,j*2  ,k*2-1) + &
                           R_f(c,i*2  ,j*2-1,k*2  ) + &
                           R_f(c,i*2  ,j*2  ,k*2  )
          enddo
        enddo
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine fine2coarse_R

    subroutine te_estimation(gc,Ni,Nj,Nk,Nc,U,R,UR,TE)
    !-------------------------------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN)::    gc(1:6)                                                     ! Number of ghost cells.
    integer(I_P), intent(IN)::    Ni                                                          ! Number of cells in i direction.
    integer(I_P), intent(IN)::    Nj                                                          ! Number of cells in j direction.
    integer(I_P), intent(IN)::    Nk                                                          ! Number of cells in k direction.
    integer(I_P), intent(IN)::    Nc                                                          ! Number of conservative variables.
    real(R_P),    intent(IN)::    U (1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Conservative variables.
    real(R_P),    intent(IN)::    R (1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Residuals of conservative variables.
    real(R_P),    intent(OUT)::   UR(1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Conservative variables.
    real(R_P),    intent(INOUT):: TE(1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Truncation error of cons. variables.
    integer(I_P)::                c,i,j,k                                                     ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(c,i,j,k)       &
    !$OMP SHARED(Nc,Ni,Nj,Nk,U,R,UR,TE)
    !$OMP DO
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          do c=1,Nc
            UR(c,i,j,k) = U(c,i,j,k)
            TE(c,i,j,k) = R(c,i,j,k) - TE(c,i,j,k)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine te_estimation
  endsubroutine restrict

  subroutine correct(myrank,l,Nb,gco,gc,Ni_f,Nj_f,Nk_f,Ni_c,Nj_c,Nk_c,bb_f,bb_c,bb_Fi_f,bb_Fj_f,bb_Fk_f,Np,Nc,Ns, &
                     NFi_f,NFj_f,NFk_f,BCi_f,BCj_f,BCk_f,Nin1,in1,Nin2,bb_in2,in2,Si_f,Sj_f,Sk_f,V_f,UR,U_c,      &
                     t,U_f,P_f)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for the correction of the fine solution.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank                                  ! Actual rank process.
  integer(I_P),      intent(IN)::    l                                       ! Actual grid level.
  integer(I_P),      intent(IN)::    Nb                                      ! Number of blocks.
  integer(I_P),      intent(IN)::    gco                                     ! Number of ghost cells necessary for achieving the
                                                                             ! space order.
  integer(I_P),      intent(IN)::    gc(1:6,1:Nb)                            ! Number of ghost cells.
  integer(I_P),      intent(IN)::    Ni_f(  1:Nb)                            ! Number of cells in i direction (fine grid).
  integer(I_P),      intent(IN)::    Nj_f(  1:Nb)                            ! Number of cells in j direction (fine grid).
  integer(I_P),      intent(IN)::    Nk_f(  1:Nb)                            ! Number of cells in k direction (fine grid).
  integer(I_P),      intent(IN)::    Ni_c(  1:Nb)                            ! Number of cells in i direction (coarse grid).
  integer(I_P),      intent(IN)::    Nj_c(  1:Nb)                            ! Number of cells in j direction (coarse grid).
  integer(I_P),      intent(IN)::    Nk_c(  1:Nb)                            ! Number of cells in k direction (coarse grid).
  integer(I_P),      intent(IN)::    bb_f(1:2,1:Nb)                          ! Bounds of cell centers.
  integer(I_P),      intent(IN)::    bb_c(1:2,1:Nb)                          ! Bounds of cell centers.
  integer(I_P),      intent(IN)::    bb_Fi_f(1:2,1:Nb)                       ! Bounds of interfaces in i direction.
  integer(I_P),      intent(IN)::    bb_Fj_f(1:2,1:Nb)                       ! Bounds of interfaces in j direction.
  integer(I_P),      intent(IN)::    bb_Fk_f(1:2,1:Nb)                       ! Bounds of interfaces in k direction.
  integer(I_P),      intent(IN)::    Np                                      ! Number of primitive variables.
  integer(I_P),      intent(IN)::    Nc                                      ! Number of conservative variables.
  integer(I_P),      intent(IN)::    Ns                                      ! Number of species.
  type(Type_Vector), intent(IN)::    NFi_f (     bb_Fi_f(1,1):bb_Fi_f(2,Nb)) ! Face i normal.
  type(Type_Vector), intent(IN)::    NFj_f (     bb_Fj_f(1,1):bb_Fj_f(2,Nb)) ! Face j normal.
  type(Type_Vector), intent(IN)::    NFk_f (     bb_Fk_f(1,1):bb_Fk_f(2,Nb)) ! Face k normal.
  type(Type_BC),     intent(IN)::    BCi_f (     bb_Fi_f(1,1):bb_Fi_f(2,Nb)) ! Bc of i interfaces.
  type(Type_BC),     intent(IN)::    BCj_f (     bb_Fj_f(1,1):bb_Fj_f(2,Nb)) ! Bc of j interfaces.
  type(Type_BC),     intent(IN)::    BCk_f (     bb_Fk_f(1,1):bb_Fk_f(2,Nb)) ! Bc of k interfaces.
  integer(I_P),      intent(IN)::    Nin1                                    ! Number of inflow 1 boundary conditions.
  real(R_P),         intent(IN)::    in1(1:Np,1:Nin1)                        ! Primitive variables of inflow 1 bc.
  integer(I_P),      intent(IN)::    Nin2                                    ! Number of inflow 2 boundary conditions.
  integer(I_P),      intent(IN)::    bb_in2(1:2,1:Nin2)                      ! Bounds of inflow 2 boundary conditions.
  real(R_P),         intent(IN)::    in2(1:Np+1,bb_in2(1,1):bb_in2(2,Nin2))  ! Primitive variables (plus time) of inflow 2 bc.
  real(R_P),         intent(IN)::    Si_f  (     bb_Fi_f(1,1):bb_Fi_f(2,Nb)) ! I surface area.
  real(R_P),         intent(IN)::    Sj_f  (     bb_Fj_f(1,1):bb_Fj_f(2,Nb)) ! J surface area.
  real(R_P),         intent(IN)::    Sk_f  (     bb_Fk_f(1,1):bb_Fk_f(2,Nb)) ! K surface area.
  real(R_P),         intent(IN)::    V_f   (     bb_f   (1,1):bb_f   (2,Nb)) ! Cell volume.
  real(R_P),         intent(IN)::    UR    (1:Nc,bb_c   (1,1):bb_c   (2,Nb)) ! Restriction of conservative variables.
  real(R_P),         intent(IN)::    U_c   (1:Nc,bb_c   (1,1):bb_c   (2,Nb)) ! Conservative variables (coarse grid).
  real(R_P),         intent(IN)::    t                                       ! Time.
  real(R_P),         intent(INOUT):: U_f   (1:Nc,bb_f   (1,1):bb_f   (2,Nb)) ! Conservative variables (fine grid).
  real(R_P),         intent(OUT)::   P_f   (1:Np,bb_f   (1,1):bb_f   (2,Nb)) ! Primitive variables (fine grid).
  real(R_P)::                        R_f   (1:Nc,bb_f   (1,1):bb_f   (2,Nb)) ! Residuals of conservative variables (fine grid).
  real(R_P)::                        CU    (1:Nc,bb_c   (1,1):bb_c   (2,Nb)) ! Correction of conservative variables (coarse grid).
  real(R_P)::                        CUp   (1:Nc,bb_f   (1,1):bb_f   (2,Nb)) ! Correction prolongation of conservative variables.
  integer(I_P)::                     b                                       ! Block counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do b=1,Nb
    ! compute the coarse correction
    call get_correction(gc = gc(1:6,b),                     &
                        Ni = Ni_c(  b),                     &
                        Nj = Nj_c(  b),                     &
                        Nk = Nk_c(  b),                     &
                        Nc = Nc,                            &
                        UR = UR (1:Nc,bb_c(1,b):bb_c(2,b)), &
                        U  = U_c(1:Nc,bb_c(1,b):bb_c(2,b)), &
                        CU = CU (1:Nc,bb_c(1,b):bb_c(2,b)))
    ! prolongate the coarse correction solution
    call coarse2fine_U(gc = gc(1:6,b), Nc = Nc,                        &
                       Ni_c = Ni_c(b), Nj_c = Nj_c(b), Nk_c = Nk_c(b), &
                       Ni_f = Ni_f(b), Nj_f = Nj_f(b), Nk_f = Nk_f(b), &
                       U_c = CU (1:Nc,bb_c(1,b):bb_c(2,b)),            &
                       U_f = CUp(1:Nc,bb_f(1,b):bb_f(2,b)))
    ! correct the fine solution
    call update_fine(gc = gc(1:6,b),                     &
                     Ni = Ni_f(  b),                     &
                     Nj = Nj_f(  b),                     &
                     Nk = Nk_f(  b),                     &
                     Nc = Nc,                            &
                     CU = CUp(1:Nc,bb_f(1,b):bb_f(2,b)), &
                     U  = U_f(1:Nc,bb_f(1,b):bb_f(2,b)))
    ! converting conservative variables to primitive ones
    call conservative2primitive(gc   = gc(1:6,b),                     &
                                Ni   = Ni_f(  b),                     &
                                Nj   = Nj_f(  b),                     &
                                Nk   = Nk_f(  b),                     &
                                Np   = Np,                            &
                                Nc   = Nc,                            &
                                Ns   = Ns,                            &
                                cons = U_f(1:Nc,bb_f(1,b):bb_f(2,b)), &
                                prim = P_f(1:Np,bb_f(1,b):bb_f(2,b)))
  enddo
  ! imposing the boundary conditions
  do b=1,Nb
    call boundary_conditions(myrank = myrank,                                 &
                             l      = l,                                      &
                             b      = b,                                      &
                             gc     = gc(1:6,b),                              &
                             Ni     = Ni_f(  b),                              &
                             Nj     = Nj_f(  b),                              &
                             Nk     = Nk_f(  b),                              &
                             Np     = Np,                                     &
                             Ns     = Ns,                                     &
                             NFi    = NFi_f(     bb_Fi_f(1,b):bb_Fi_f(2,b)),  &
                             NFj    = NFj_f(     bb_Fj_f(1,b):bb_Fj_f(2,b)),  &
                             NFk    = NFk_f(     bb_Fk_f(1,b):bb_Fk_f(2,b)),  &
                             BCi    = BCi_f(     bb_Fi_f(1,b):bb_Fi_f(2,b)),  &
                             BCj    = BCj_f(     bb_Fj_f(1,b):bb_Fj_f(2,b)),  &
                             BCk    = BCk_f(     bb_Fk_f(1,b):bb_Fk_f(2,b)),  &
                             Nin1   = Nin1,                                   &
                             in1    = in1(1:Np,1:Nin1),                       &
                             Nin2   = Nin2,                                   &
                             bb_in2 = bb_in2(1:2,1:Nin2),                     &
                             in2    = in2(1:Np+1,bb_in2(1,1):bb_in2(2,Nin2)), &
                             t      = t,                                      &
                             P      = P_f  (1:Np,bb_f   (1,b):bb_f   (2,b)))
  enddo
  ! computing the residuals for uptdating normal velocity at interface (unF) doing a residuals computation with the updated fine
  ! solution
  do b=1,Nb
    !call residuals(gco    = gco,                                   &
                   !gc     = gc(1:6,b),                             &
                   !Ni     = Ni_f(  b),                             &
                   !Nj     = Nj_f(  b),                             &
                   !Nk     = Nk_f(  b),                             &
                   !Np     = Np,                                    &
                   !Nc     = Nc,                                    &
                   !Ns     = Ns,                                    &
                   !NFi    = NFi_f(     bb_Fi_f(1,b):bb_Fi_f(2,b)), &
                   !NFj    = NFj_f(     bb_Fj_f(1,b):bb_Fj_f(2,b)), &
                   !NFk    = NFk_f(     bb_Fk_f(1,b):bb_Fk_f(2,b)), &
                   !Si     = Si_f (     bb_Fi_f(1,b):bb_Fi_f(2,b)), &
                   !Sj     = Sj_f (     bb_Fj_f(1,b):bb_Fj_f(2,b)), &
                   !Sk     = Sk_f (     bb_Fk_f(1,b):bb_Fk_f(2,b)), &
                   !V      = V_f  (     bb_f   (1,b):bb_f   (2,b)), &
                   !BCi    = BCi_f(     bb_Fi_f(1,b):bb_Fi_f(2,b)), &
                   !BCj    = BCj_f(     bb_Fj_f(1,b):bb_Fj_f(2,b)), &
                   !BCk    = BCk_f(     bb_Fk_f(1,b):bb_Fk_f(2,b)), &
                   !P      = P_f  (1:Np,bb_f   (1,b):bb_f   (2,b)), &
                   !Pmax   = Pmax (1:Np                          ), &
                   !Pmin   = Pmin (1:Np                          ), &
                   !R      = R_f  (1:Nc,bb_f   (1,b):bb_f   (2,b)))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    subroutine get_correction(gc,Ni,Nj,Nk,Nc,UR,U,CU)
    !-------------------------------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN)::  gc(1:6)                                                     ! Ghost cells number.
    integer(I_P), intent(IN)::  Ni,Nj,Nk                                                    ! Cells number (i,j,k).
    integer(I_P), intent(IN)::  Nc                                                          ! Number of conservative variables.
    real(R_P),    intent(IN)::  UR(1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Restriction of conservative variables.
    real(R_P),    intent(IN)::  U (1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Conservative variables.
    real(R_P),    intent(OUT):: CU(1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Correction of conservative variables.
    integer(I_P)::              i,j,k                                                       ! Spaces counters.
    integer(I_P)::              c                                                           ! Conservative variables counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    CU = 0._R_P
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(c,i,j,k)       &
    !$OMP SHARED(Nc,Ni,Nj,Nk,CU,U,UR)
    !$OMP DO
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          do c=1,Nc
            CU(c,i,j,k) = U(c,i,j,k) - UR(c,i,j,k)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine get_correction

    subroutine coarse2fine_U(gc,Ni_c,Nj_c,Nk_c,Ni_f,Nj_f,Nk_f,Nc,U_c,U_f)
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Subroutine coarse2fine prolongates the conservative variables on the fine grid from the coarse one.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN)::  gc(1:6)                                                            ! Ghost cells number.
    integer(I_P), intent(IN)::  Ni_c,Nj_c,Nk_c                                                     ! Cells number of coarse grid.
    integer(I_P), intent(IN)::  Ni_f,Nj_f,Nk_f                                                     ! Cells number of fine grid.
    integer(I_P), intent(IN)::  Nc                                                                 ! Number of cons. variables.
    real(R_P),    intent(IN)::  U_c(1:Nc,1-gc(1):Ni_c+gc(2),1-gc(3):Nj_c+gc(4),1-gc(5):Nk_c+gc(6)) ! Coarse grid cons. variables.
    real(R_P),    intent(OUT):: U_f(1:Nc,1-gc(1):Ni_f+gc(2),1-gc(3):Nj_f+gc(4),1-gc(5):Nk_f+gc(6)) ! Fine grid cons. variables.
    real(R_P)::                 wi,wi_m,wj,wj_m,wk,wk_m                                            ! Weights of trilinear interp.
    integer(I_P)::              i,j,k,i_c,j_c,k_c,c                                                ! Counters.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    !do k=1,Nk_c
      !do j=1,Nj_c
        !do i=1,Ni_c
          !do c=1,Nc
            !U_f(c,i*2-1,j*2-1,k*2-1) = U_c(c,i,j,k)
            !U_f(c,i*2-1,j*2  ,k*2-1) = U_c(c,i,j,k)
            !U_f(c,i*2-1,j*2  ,k*2  ) = U_c(c,i,j,k)
            !U_f(c,i*2-1,j*2-1,k*2  ) = U_c(c,i,j,k)
            !U_f(c,i*2  ,j*2-1,k*2-1) = U_c(c,i,j,k)
            !U_f(c,i*2  ,j*2  ,k*2-1) = U_c(c,i,j,k)
            !U_f(c,i*2  ,j*2  ,k*2  ) = U_c(c,i,j,k)
            !U_f(c,i*2  ,j*2-1,k*2  ) = U_c(c,i,j,k)
          !enddo
        !enddo
      !enddo
    !enddo
    U_f = 0._R_P
    do k=1,Nk_f
      k_c  = k/2 + 1
      wk   = real(1 + 2*mod(k,2),R_P)/4._R_P
      wk_m = 1._R_P - wk
      do j=1,Nj_f
        j_c  = j/2 + 1
        wj   = real(1 + 2*mod(j,2),R_P)/4._R_P
        wj_m = 1._R_P - wj
        do i=1,Ni_f
          i_c  = i/2 + 1
          wi   = real(1 + 2*mod(i,2),R_P)/4._R_P
          wi_m = 1._R_P - wi
          do c=1,Nc
            U_f(c,i,j,k) = wi_m*((U_c(c,i_c-1,j_c-1,k_c-1)*wj_m + U_c(c,i_c-1,j_c,k_c-1)*wj)*wk_m  + &
                                 (U_c(c,i_c-1,j_c-1,k_c  )*wj_m + U_c(c,i_c-1,j_c,k_c  )*wj)*wk  ) + &
                             wi*((U_c(c,i_c  ,j_c-1,k_c-1)*wj_m + U_c(c,i_c  ,j_c,k_c-1)*wj)*wk_m  + &
                                 (U_c(c,i_c  ,j_c-1,k_c  )*wj_m + U_c(c,i_c  ,j_c,k_c  )*wj)*wk  )

          enddo
        enddo
      enddo
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine coarse2fine_U

    subroutine update_fine(gc,Ni,Nj,Nk,Nc,CU,U)
    !-------------------------------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I_P), intent(IN)::    gc(1:6)                                                     ! Ghost cells number.
    integer(I_P), intent(IN)::    Ni,Nj,Nk                                                    ! Cells number (i,j,k).
    integer(I_P), intent(IN)::    Nc                                                          ! Number of conservative variables.
    real(R_P),    intent(IN)::    CU(1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Correction for fine grid.
    real(R_P),    intent(INOUT):: U (1:Nc,1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)) ! Conservative variables.
    integer(I_P)::                i,j,k                                                       ! Spaces counters.
    integer(I_P)::                c                                                           ! Conservative variables counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(c,i,j,k)       &
    !$OMP SHARED(Nc,Ni,Nj,Nk,CU,U)
    !$OMP DO
    do k=1,Nk
      do j=1,Nj
        do i=1,Ni
          do c=1,Nc
            U(c,i,j,k) = U(c,i,j,k) + CU(c,i,j,k)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine update_fine
  endsubroutine correct
endmodule Lib_Multigrid
