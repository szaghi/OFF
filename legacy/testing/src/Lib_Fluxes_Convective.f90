!> @ingroup Library
!> @{
!> @defgroup Lib_Fluxes_ConvectiveLibrary Lib_Fluxes_Convective
!> @}

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
!> @ingroup Lib_Fluxes_ConvectiveLibrary
module Lib_Fluxes_Convective
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                        ! Integers and reals precision definition.
USE Data_Type_Cell,                   only: Type_Cell                   ! Definition of Type_Cell.
USE Data_Type_Conservative,           only: Type_Conservative           ! Definition of Type_Conservative.
USE Data_Type_Face,                   only: Type_Face                   ! Definition of Type_Face.
USE Data_Type_Primitive,              only: Type_Primitive              ! Definition of Type_Primitive.
USE Data_Type_Primitive1D,            only: Type_Primitive1D            ! Definition of Type_Primitive1D.
USE Data_Type_Riemann_Conservative1D, only: Type_Riemann_Conservative1D ! Definition of Type_Riemann_Conservative1D.
USE Data_Type_Riemann_Primitive1D,    only: Type_Riemann_Primitive1D    ! Definition of Type_Riemann_Primitive1D.
USE Data_Type_Species,                only: Type_Species                ! Definition of Type_Species.
USE Data_Type_Vector,                 only: Type_Vector                 ! Definition of Type_Vector.
USE Lib_Riemann_Solvers,              only: riemann_solver              ! Procedure for solving a Riemann problem.
USE Lib_WENO,                         only: weno                        ! Procedure for computing WENO reconstruction.
USE Lib_Variables_Conversions                                           ! Procedure for converting primitives to 1D primitives.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: set_interface_reconstruction
public:: fluxes_convective
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
procedure(IR_interface), pointer:: interface_reconstruction => null()
!> Abstract interfaces for pointer procedures.
abstract interface
  pure subroutine IR_interface(gc,N,Ns,species0,F,C,primR)
  use IR_Precision
  use Data_Type_Cell
  use Data_Type_Face
  use Data_Type_Primitive1D
  use Data_Type_Species
  integer(I1P),           intent(IN)::  gc
  integer(I4P),           intent(IN)::  N
  integer(I4P),           intent(IN)::  Ns
  type(Type_Species),     intent(IN)::  species0
  type(Type_Face),        intent(IN)::  F(0-gc:N+gc)
  type(Type_Cell),        intent(IN)::  C(1-gc:N+gc)
  type(Type_Primitive1D), intent(OUT):: primR(1:2,0:N+1)
  endsubroutine IR_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Fluxes_ConvectivePublicProcedure
  !> @{
  !> @brief Procedure for setting the interface values reconstruction algorithm accordingly to the user options.
  subroutine set_interface_reconstruction(reconstruction_type)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: reconstruction_type !< Type of variables to be reconstructed:
                                                 !< - "primitive_reconstruction"  => reconstruction on primitive vars;
                                                 !< - "characteristic_reconstruction" => reconstruction on characteristic vars.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(trim(adjustl(reconstruction_type)))
  case("primitive_reconstruction")
    interface_reconstruction => interface_reconstruction_prim
  case("characteristic_reconstruction")
    interface_reconstruction => interface_reconstruction_charac
  case default
    interface_reconstruction => interface_reconstruction_charac
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_interface_reconstruction

  !> @brief Procedure for computing interfaces convective fluxes (hyperbolic operator). Taking in input the primitive variables
  !> along a 1D slice of N+2*gc cells the procedure returns in ouput the N+1 convective fluxes at cells interfaces.
  !> @note For avoiding the creation of temporary arrays (improving the efficiency) the arrays \b NF, \b P and \b F are declared as
  !> assumed-shape with only the lower bound defined. Their extentions are: NF [0-gc:N+gc], P [1-gc:N+gc], F [0:N].
  pure subroutine fluxes_convective(gc,N,species0,F,C,Fl)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN)::    gc                !< Number of ghost cells used.
  integer(I4P),            intent(IN)::    N                 !< Number of cells.
  type(Type_Species),      intent(IN)::    species0          !< Initial species.
  type(Type_Face),         intent(IN)::    F(0-gc:)          !< Faces data [0-gc:N+gc].
  type(Type_Cell),         intent(IN)::    C(1-gc:)          !< Cells data [1-gc:N+gc].
  type(Type_Conservative), intent(INOUT):: Fl(0:)            !< Convective fluxes of 3D conservative variables [0:N].
  type(Type_Primitive1D)::                 primR(1:2,0:N+1)  !< Reconstructed interface values of 1D primitive variables.
  type(Type_Riemann_Primitive1D)::         state1            !< State 1 (left) of interface Riemann Problem.
  type(Type_Riemann_Primitive1D)::         state4            !< State 4 (right) of interface Riemann Problem.
  type(Type_Riemann_Conservative1D)::      flux              !< Convective fluxes of 1D conservative variables.
  type(Type_Vector)::                      ut(1:2,1-gc:N+gc) !< Tangential velocity: left (1) and right (2) interfaces.
  integer(I4P)::                           i,j               !< Counters.
#ifdef LMA
  real(R8P)::                              ulma(1:2)         !< Left (1) and right (2) reconstructed normal speed adjusted
                                                             !< for low Mach number flows.
  real(R8P)::                              z                 !< Scaling coefficient for Low Mach number Adjustment.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing velocity vectors and its tangential component
  do i=1-gc,N+gc
    ! computing tangential velocity component
    do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      if (i==1-gc.AND.j==1) cycle
      if (i==N+gc.AND.j==2) cycle
      ut(j,i) = C(i)%P%v - (C(i)%P%v.paral.F(i+j-2)%N) ! ut
    enddo
  enddo
  ! computing the left and right states or Riemann Problems
  select case(gc)
  case(2_I1P,3_I1P,4_I1P) ! 3rd, 5th or 7th order WENO reconstruction
    ! computing the reconstruction
    call interface_reconstruction(gc=gc,N=N,Ns=species0%Ns,species0=species0,F=F,C=C,primR=primR)
  case(1_I1P) ! 1st order piecewise constant reconstruction
    do i=0,N+1
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        call prim3D2prim1D(prim3D=C(i)%P,N=F(i+j-2)%N,prim1D=primR(j,i))
      enddo
    enddo
  endselect
#ifdef LMA
  ! applying Low Mach number Adjustment for decreasing the dissipation in local low Mach number region
  do i=0,N
    z = min(1._R8P,max(sqrt(primR(2,i  )%u*primR(2,i  )%u+ut(2,i  )%sq_norm())/C(i  )%P%a(),&
                       sqrt(primR(1,i+1)%u*primR(1,i+1)%u+ut(1,i+1)%sq_norm())/C(i+1)%P%a()))
    ulma(1) = 0.5_R8P*(primR(2,i)%u + primR(1,i+1)%u) + z*0.5_R8P*(primR(2,i)%u - primR(1,i+1)%u)
    ulma(2) = 0.5_R8P*(primR(2,i)%u + primR(1,i+1)%u) - z*0.5_R8P*(primR(2,i)%u - primR(1,i+1)%u)
    primR(2,i  )%u = ulma(1)
    primR(1,i+1)%u = ulma(2)
  enddo
#endif
  ! solving Riemann Problems
  do i=0,N
    ! face normal fluxes
    call prim1D2riemprim1D(prim1D=primR(2,i  ),riemprim1D=state1)
    call prim1D2riemprim1D(prim1D=primR(1,i+1),riemprim1D=state4)
    call riemann_solver(state1=state1,state4=state4,flux=flux)
    ! uptdating fluxes with tangential components
    if (flux%r>0._R8P) then
      Fl(i)%rs = flux%r*primR(2,i  )%r/primR(2,i  )%d
      Fl(i)%rv = flux%u*F(i)%N + flux%r*ut(2,i)
      Fl(i)%re = flux%E + flux%r*ut(2,i)%sq_norm()*0.5_R8P
    else
      Fl(i)%rs = flux%r*primR(1,i+1)%r/primR(1,i+1)%d
      Fl(i)%rv = flux%u*F(i)%N + flux%r*ut(1,i+1)
      Fl(i)%re = flux%E + flux%r*ut(1,i+1)%sq_norm()*0.5_R8P
    endif
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine fluxes_convective
  !> @}

  !> @ingroup Lib_Fluxes_ConvectivePrivateProcedure
  !> @{
  !> @brief Procedure for imposing the positivity preserving limiter.
  pure subroutine positivity_preserving_limiter(o,vmean,vrL,vrR)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN)::    o                                        ! Order of space reconstruction.
  real(R8P),    intent(IN)::    vmean                                    ! Mean value.
  real(R8P),    intent(INOUT):: vrL                                      ! Left reconstructed value.
  real(R8P),    intent(INOUT):: vrR                                      ! Right reconstructed value.
  real(R8P), parameter::        w1(2:4)=[1._R8P/6._R8P,  &               ! Weights of Gauss-Lobatto's 3 points quadrature.
                                         1._R8P/12._R8P, &               ! Weights of Gauss-Lobatto's 4 points quadrature.
                                         1._R8P/20._R8P]                 ! Weights of Gauss-Lobatto's 5 points quadrature.
  real(R8P), parameter::        w2(2:4)=[1._R8P/(1._R8P-2._R8P*w1(2)), & ! 1/(1-2*w1)
                                         1._R8P/(1._R8P-2._R8P*w1(3)), & ! 1/(1-2*w1)
                                         1._R8P/(1._R8P-2._R8P*w1(4))]   ! 1/(1-2*w1)
  real(R8P), parameter::        e  = tiny(1._R8P)                        ! Small number for avoid division for zero.
  real(R8P)::                   vstar                                    ! Star value.
  real(R8P)::                   theta                                    ! Limiter coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vstar = (vmean - w1(o)*(vrL + vrR))*w2(o)
  theta = min((vmean-0.1_R8P*vmean)/(vmean - min(vstar,vrL,vrR) + e),1._R8P)
  vrL   = theta*(vrL - vmean) + vmean
  vrR   = theta*(vrR - vmean) + vmean
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine positivity_preserving_limiter

  !> @brief Procedure for reconstructing interface values of primitive variables along the interface F using WENO interpolation on
  !> (local) chracacteristic projection.
  pure subroutine interface_reconstruction_charac(gc,N,Ns,species0,F,C,primR)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),           intent(IN)::  gc                         !< Number of ghost cells used.
  integer(I4P),           intent(IN)::  N                          !< Number of cells.
  integer(I4P),           intent(IN)::  Ns                         !< Number of initial species.
  type(Type_Species),     intent(IN)::  species0                   !< Initial species.
  type(Type_Face),        intent(IN)::  F(0-gc:N+gc)               !< Faces data [0-gc:N+gc].
  type(Type_Cell),        intent(IN)::  C(1-gc:N+gc)               !< Cells data [1-gc:N+gc].
  type(Type_Primitive1D), intent(OUT):: primR(1:2,0:N+1)           !< Reconstructed interface values of 1D primitives.
  type(Type_Primitive1D)::              prim1D(1:2,1-gc:-1+gc)     !< 1D primitive variables into the stencil used.
  real(R8P)::                           WI(1:Ns+4,1:2,1-gc:-1+gc)  !< WENO Input: can be the interface value of either
                                                                   !< the characteristic variables or the primitive ones.
  real(R8P)::                           WO(1:Ns+4,1:2)             !< WENO Ouput: can be the reconstructed interface value of either
                                                                   !< the characteristic variables or the primitive ones.
  type(Type_Primitive1D)::              prim1Dm                    !< Mean of 1D primitive.
  real(R8P)::                           eigenLm(1:Ns+4,1:Ns+4,1:2) !< Left eigenvectors matrix of prim1Dm.
  real(R8P)::                           eigenRm(1:Ns+4,1:Ns+4,1:2) !< Left eigenvectors matrix of prim1Dm.
  logical::                             ROR(1:2)                   !< Logical flag for testing the result of WENO reconst.
  integer(I1P)::                        or                         !< Counter of order for ROR algorithm.
  integer(I4P)::                        i,j,k,v                    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=0,N+1 ! loop over cells
    ! computing 1D primitive variables for the stencil [i+1-gc,i-1+gc] fixing the velocity across the interfaces i+-1/2
    do k=i+1-gc,i-1+gc
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        call prim3D2prim1D(prim3D=C(k)%P,N=F(i+j-2)%N,prim1D=prim1D(j,k-i))
      enddo
    enddo
    ! computing mean left and right eigenvectors matrix across the interface
    do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      if (i==0  .AND.j==1) cycle
      if (i==N+1.AND.j==2) cycle
      ! computing mean of primitive variables across the interface
      prim1Dm = 0.5_R8P*(prim1D(j,j-1)+prim1D(j,j-2))
      ! computing mean left and right eigenvectors matrix across the interface
      eigenLm(1:Ns+4,1:Ns+4,j) = prim1Dm%eigenvectL()
      eigenRm(1:Ns+4,1:Ns+4,j) = prim1Dm%eigenvectR()
    enddo
    ! preparing the WENO input
    do k=i+1-gc,i-1+gc
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        ! transforming variables into local characteristic fields for the stencil [1-gc,-1+gc]
        call prim1D2charac(prim=prim1D(j,k-i),eigenvectL=eigenLm(:,:,j),charac=WI(:,j,k-i))
      enddo
    enddo
    ! computing WENO reconstruction with ROR (Recursive Order Reduction) algorithm
    ROR_check: do or=gc,2_I1P,-1_I1P
      ! computing WENO reconstruction
      do v=1,Ns+2
        call weno(S=or,V=WI(v,1:2,1-or:-1+or),VR=WO(v,1:2))
      enddo
      ! converting WENO output
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        ! trasforming local reconstructed characteristic variables to primitive ones
        call charac2prim1D(charac=WO(:,j),eigenvectR=eigenRm(:,:,j),prim=primR(j,i))
      enddo
      ! extrapolation of the reconstructed values for the non computed boundaries
      if (i==0  ) primR(1,i) = primR(2,i)
      if (i==N+1) primR(2,i) = primR(1,i)
#ifdef PPL
      ! applaying the maximum-principle-satisfying limiter to the reconstructed densities and pressure
      do v=1,Ns
        call positivity_preserving_limiter(o=or,vmean=C(i)%P%r(v),vrL=primR(1,i)%r(v),vrR=primR(2,i)%r(v)) ! densities
      enddo
        call positivity_preserving_limiter(o=or,vmean=C(i)%P%p   ,vrL=primR(1,i)%p   ,vrR=primR(2,i)%p   ) ! pressure
#endif
      ! computing the reconstructed variables of total density and specific heats ratio
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        call primR(j,i)%compute_d
        call primR(j,i)%compute_g(species0=species0)
      enddo
      ! verifing the WENO reconstruction at order or-th
      ROR = .true.
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        ROR(j) = ((primR(j,i)%p>0._R8P).AND.(primR(j,i)%d>0._R8P).AND.(primR(j,i)%g>0._R8P))
      enddo
      if (ROR(1).AND.ROR(2)) exit ROR_check
      if (or==2) then
        ! the high order reconstruction is falied; used 1st order reconstruction
        do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          if (i==0  .AND.j==1) cycle
          if (i==N+1.AND.j==2) cycle
          call prim3D2prim1D(prim3D=C(i)%P,N=F(i+j-2)%N,prim1D=primR(j,i))
        enddo
      endif
    enddo ROR_check
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interface_reconstruction_charac

  !> @brief Procedure for reconstructing interface values of primitive variables along the interface F using WENO interpolation on
  !> primitive variables directly.
  pure subroutine interface_reconstruction_prim(gc,N,Ns,species0,F,C,primR)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),           intent(IN)::  gc                         !< Number of ghost cells used.
  integer(I4P),           intent(IN)::  N                          !< Number of cells.
  integer(I4P),           intent(IN)::  Ns                         !< Number of initial species.
  type(Type_Species),     intent(IN)::  species0                   !< Initial species.
  type(Type_Face),        intent(IN)::  F(0-gc:N+gc)               !< Faces data [0-gc:N+gc].
  type(Type_Cell),        intent(IN)::  C(1-gc:N+gc)               !< Cells data [1-gc:N+gc].
  type(Type_Primitive1D), intent(OUT):: primR(1:2,0:N+1)           !< Reconstructed interface values of 1D primitives.
  type(Type_Primitive1D)::              prim1D(1:2,1-gc:-1+gc)     !< 1D primitive variables into the stencil used.
  real(R8P)::                           WI(1:Ns+4,1:2,1-gc:-1+gc)  !< WENO Input: can be the interface value of either
                                                                   !< the characteristic variables or the primitive ones.
  real(R8P)::                           WO(1:Ns+4,1:2)             !< WENO Ouput: can be the reconstructed interface value
                                                                   !< of either the characteristic variables or the
                                                                   !< primitive ones.
  logical::                             ROR(1:2)                   !< Logical flag for testing the result of WENO reconst.
  integer(I1P)::                        or                         !< Counter of order for ROR algorithm.
  integer(I4P)::                        i,j,k,v                    !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=0,N+1 ! loop over cells
    ! computing 1D primitive variables for the stencil [i+1-gc,i-1+gc] fixing the velocity across the interfaces i+-1/2
    do k=i+1-gc,i-1+gc
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        call prim3D2prim1D(prim3D=C(k)%P,N=F(i+j-2)%N,prim1D=prim1D(j,k-i))
      enddo
    enddo
    ! preparing the WENO input
    do k=i+1-gc,i-1+gc
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        WI(:,j,k-i) = prim1D(j,k-i)%prim2array()
      enddo
    enddo
    ! computing WENO reconstruction with ROR (Recursive Order Reduction) algorithm
    ROR_check: do or=gc,2_I1P,-1_I1P
      ! computing WENO reconstruction
      do v=1,Ns+2
        call weno(S=or,V=WI(v,1:2,1-or:-1+or),VR=WO(v,1:2))
      enddo
      ! converting WENO output
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        call primR(j,i)%array2prim(WO(:,j))
      enddo
      ! extrapolation of the reconstructed values for the non computed boundaries
      if (i==0  ) primR(1,i) = primR(2,i)
      if (i==N+1) primR(2,i) = primR(1,i)
#ifdef PPL
      ! applaying the maximum-principle-satisfying limiter to the reconstructed densities and pressure
      do v=1,Ns
        call positivity_preserving_limiter(o=or,vmean=C(i)%P%r(v),vrL=primR(1,i)%r(v),vrR=primR(2,i)%r(v)) ! densities
      enddo
        call positivity_preserving_limiter(o=or,vmean=C(i)%P%p   ,vrL=primR(1,i)%p   ,vrR=primR(2,i)%p   ) ! pressure
#endif
      ! computing the reconstructed variables of total density and specific heats ratio
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        call primR(j,i)%compute_d
        call primR(j,i)%compute_g(species0=species0)
      enddo
      ! verifing the WENO reconstruction at order or-th
      ROR = .true.
      do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        if (i==0  .AND.j==1) cycle
        if (i==N+1.AND.j==2) cycle
        ROR(j) = ((primR(j,i)%p>0._R8P).AND.(primR(j,i)%d>0._R8P).AND.(primR(j,i)%g>0._R8P))
      enddo
      if (ROR(1).AND.ROR(2)) exit ROR_check
      if (or==2) then
        ! the high order reconstruction is falied; used 1st order reconstruction
        do j=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
          if (i==0  .AND.j==1) cycle
          if (i==N+1.AND.j==2) cycle
          call prim3D2prim1D(prim3D=C(i)%P,N=F(i+j-2)%N,prim1D=primR(j,i))
        enddo
      endif
    enddo ROR_check
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interface_reconstruction_prim
  !> @}
endmodule Lib_Fluxes_Convective
