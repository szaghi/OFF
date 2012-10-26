!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws.
!> This is a library module.
!> The Riemann Problem solvers provide as solution the convective fluxes
!> \f$(\overline{\overline F} - \overline{\overline {{Q_S}}})_{conv} \f$ normal to
!> the interface direction. The solvers contained into this library have a unique API. They take as input the primitive variables
!> in the left (state 1) and right (state 4) cells with respect the interface and provide the convective fluxes as output:
!> @code
!> ...
!> call Riem_Solver(p1,r1,u1,g1,p4,r4,u4,g4,F_r,F_u,F_E)
!> ...
!> @endcode
!> where \f$p_1,\rho_1,u_1,\gamma_1\f$ define the left state, \f$p_4,\rho_4,u_4,\gamma_4\f$ define the right state and \f$F_r\f$,
!> \f$F_u\f$ and \f$F_E\f$ are the convective fluxes of mass, momentum and energy conservation, respectively.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Library
module Lib_Riemann
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                             ! precision of integers and reals
USE Lib_Thermodynamic_Laws_Ideal, only: E, & ! Function for computing total energy.
                                        H, & ! Function for computing total entalpy.
                                        a, & ! Function for computing speed of sound.
                                        r    ! Function for computing density.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: chk_smooth_z
public:: chk_smooth_limiter
public:: chk_smooth_liu
public:: Riem_Solver_LaxFriedrichs
public:: Riem_Solver_PVL
public:: Riem_Solver_TR
public:: Riem_Solver_TS
public:: Riem_Solver_APRS
public:: Riem_Solver_ALFR
public:: Riem_Solver_HLLC
public:: Riem_Solver_Roe
public:: Riem_Solver_Exact_U
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
real(R_P), parameter:: toll = 1.D-10 ! Tollerance.
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! auxiliary subroutines
  elemental subroutine fluxes(p,r,u,g,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing the conservative fluxes from primitive variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p   ! Pressure.
  real(R_P), intent(IN)::  r   ! Density.
  real(R_P), intent(IN)::  u   ! Velocity.
  real(R_P), intent(IN)::  g   ! Specific heats ratio.
  real(R_P), intent(OUT):: F_r ! Flux of mass conservation.
  real(R_P), intent(OUT):: F_u ! Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E ! Flux of energy conservation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  F_r = r*u
  F_u = F_r*u + p
  F_E = F_r*H(p=p,r=r,u=u,g=g)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine fluxes

  elemental subroutine Rarefaction(sgn,g,delta,eta,u0,p0,a0,ux,rx,px,ax,s0,sx)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! subroutine for computing an unknown state "x" from a known state "0" when the two states are separated by a rarefaction; it is
  ! assumed that the velocity of the unknown state "ux" is known. There is an input variable that indicates if the rarefaction
  ! propagates on the "u-a" direction (left) or on the "u+a" one (right).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),    intent(IN)::  sgn      ! sign for distinguishing "left" (-1) from "right" (1) wave
  real(R_P),    intent(IN)::  g        ! specific heats ratio
  real(R_P),    intent(IN)::  delta    ! (g-1)/2
  real(R_P),    intent(IN)::  eta      ! 2*g/(g-1)
  real(R_P),    intent(IN)::  u0,p0,a0 ! known state (speed, pressure and speed of sound)
  real(R_P),    intent(IN)::  ux       ! known speed of unknown state
  real(R_P),    intent(OUT):: rx,px,ax ! unknown pressure and density
  real(R_P),    intent(OUT):: s0,sx    ! wave speeds (head and back fronts)
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ax = a0 + sgn*delta*(ux - u0) ! unknown speed of sound
  px = p0*((ax/a0)**(eta))      ! unknown pressure
  rx = r(p=px,a=ax,g=g)         ! unknown density
  s0 = u0 + sgn*a0              ! left wave speed
  sx = ux + sgn*ax              ! right wave speed
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Rarefaction

  elemental subroutine Shock(sgn,g,gm1,gp1,u0,p0,a0,ux,rx,px,ax,ss)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! subroutine for computing an unknown state "x" from a known state "0" when the two states are separated by a shock; it is
  ! assumed that the velocity of the unknown state "ux" is known. There is an input variable that indicates if the shock propagates
  ! on the "u-a" direction (left) or on the "u+a" one (right).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),    intent(IN)::  sgn      ! sign for distinguishing "left" (-1) from "right" (1) wave
  real(R_P),    intent(IN)::  g        ! specific heats ratio
  real(R_P),    intent(IN)::  gm1,gp1  ! gm1 = g - 1 , gp1 = g + 1
  real(R_P),    intent(IN)::  u0,p0,a0 ! known state (speed, pressure and speed of sound)
  real(R_P),    intent(IN)::  ux       ! unknown speed
  real(R_P),    intent(OUT):: rx,px,ax ! unknown state (density, pressure and speed of sound).
  real(R_P),    intent(OUT):: ss       ! shock wave speed
  real(R_P)::                 M0       ! relative Mach number of known state
  real(R_P)::                 x        ! dummy variable
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x   = 0.25_R_P*gp1*(ux - u0)/a0              ! dummy variable
  M0  = x + sgn*sqrt(1.0_R_P + x*x)            ! relative Mach number of known state
  x   = 1._R_P + 2._R_P*g*(M0*M0 - 1._R_P)/gp1 ! dummy variable (pressure ratio px/p0)

  ax = a0*sqrt((gp1 + gm1*x)/(gp1 + gm1/x)) ! unknown speed of sound
  px = p0*x                                 ! unknown pressure
  rx = r(p=px,a=ax,g=g)                     ! unknown density
  ss = u0 + a0*M0                           ! shock wave speed
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Shock

  elemental subroutine WavesSpeed14u(u1,a1,g1,u4,a4,g4,u23,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing waves speed 1 and 4 knowing the value of speed (u23) of intermediates states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  u1   ! Velocity of state 1.
  real(R_P), intent(IN)::  a1   ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1   ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  u4   ! Velocity of state 4.
  real(R_P), intent(IN)::  a4   ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4   ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  u23  ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1   ! Left signal velocities.
  real(R_P), intent(OUT):: S4   ! Right signal velocities.
  real(R_P)::              x    ! Dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (u23<u1) then
    ! shock
    x  = 0.25_R_P*(g1 + 1._R_P)*(u23-u1)/a1
    S1 = u1 + a1*(x - sqrt(1.0_R_P+x*x))
  else
    ! rarefaction
    S1 = u1 - a1
  endif
  ! computing right state
  if (u23>u4) then
    ! shock
    x  = 0.25_R_P*(g4 + 1._R_P)*(u23-u4)/a4
    S4 = u4 + a4*(x + sqrt(1.0_R_P+x*x))
  else
    ! rarefaction
    S4 = u4  + a4
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine WavesSpeed14u

  elemental subroutine WavesSpeed14up(p1,u1,a1,g1,p4,u4,a4,g4,u23,p23,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing waves speed 1 and 4 knowing the value of speed and pressure (u23, p23) of intermediates states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1   ! Pressure of state 1.
  real(R_P), intent(IN)::  u1   ! Velocity of state 1.
  real(R_P), intent(IN)::  a1   ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1   ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4   ! Pressure of state 4.
  real(R_P), intent(IN)::  u4   ! Velocity of state 4.
  real(R_P), intent(IN)::  a4   ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4   ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  u23  ! Velocity of intermediate states.
  real(R_P), intent(IN)::  p23  ! Pressure of intermediate states.
  real(R_P), intent(OUT):: S1   ! Left signal velocities.
  real(R_P), intent(OUT):: S4   ! Right signal velocities.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (u23<u1) then
    ! shock
    S1 = u1 - a1*sqrt(1._R_P + 0.5_R_P*(g1 + 1._R_P)/g1*(p23/p1-1._R_P))
  else
    ! rarefaction
    S1 = u1 - a1
  endif
  ! computing right state
  if (u23>u4) then
    ! shock
    S4 = u4 + a4*sqrt(1._R_P + 0.5_R_P*(g4 + 1._R_P)/g4*(p23/p4-1._R_P))
  else
    ! rarefaction
    S4 = u4  + a4
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine WavesSpeed14up

  elemental subroutine WavesSpeed1234u(u1,a1,g1,u4,a4,g4,u23,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing waves speed knowing the value of speed (u23) of intermediates states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  u1    ! Velocity of state 1.
  real(R_P), intent(IN)::  a1    ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1    ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  u4    ! Velocity of state 4.
  real(R_P), intent(IN)::  a4    ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4    ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  u23   ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1,S2 ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4 ! Right signal velocities.
  real(R_P)::              x     ! Dummy variable.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (u23<u1) then
    ! shock
    x  = 0.25_R_P*(g1 + 1._R_P)*(u23-u1)/a1
    S1 = u1 + a1*(x - sqrt(1.0_R_P+x*x))
    S2 = S1
  else
    ! rarefaction
    S1 = u1  - a1
    S2 = u23 - a1 - 0.5_R_P*(g1 - 1._R_P)*(u23 - u1)
  endif
  ! computing right state
  if (u23>u4) then
    ! shock
    x  = 0.25_R_P*(g4 + 1._R_P)*(u23-u4)/a4
    S4 = u4 + a4*(x + sqrt(1.0_R_P+x*x))
    S3 = S4
  else
    ! rarefaction
    S4 = u4  + a4
    S3 = u23 + a4 + 0.5_R_P*(g4 - 1._R_P)*(u23 - u4)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine WavesSpeed1234u

  elemental subroutine WavesSpeed1234up(p1,u1,a1,g1,p4,u4,a4,g4,u23,p23,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing waves speed knowing the value of speed and pressure (u23, p23) of intermediates states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    ! Pressure of state 1.
  real(R_P), intent(IN)::  u1    ! Velocity of state 1.
  real(R_P), intent(IN)::  a1    ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1    ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4    ! Pressure of state 4.
  real(R_P), intent(IN)::  u4    ! Velocity of state 4.
  real(R_P), intent(IN)::  a4    ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4    ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  u23   ! Velocity of intermediate states.
  real(R_P), intent(IN)::  p23   ! Pressure of intermediate states.
  real(R_P), intent(OUT):: S1,S2 ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4 ! Right signal velocities.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (u23<u1) then
    ! shock
    S1 = u1 - a1*sqrt(1._R_P + 0.5_R_P*(g1 + 1._R_P)/g1*(p23/p1-1._R_P))
    S2 = S1
  else
    ! rarefaction
    S1 = u1  - a1
    S2 = u23 - a1 - 0.5_R_P*(g1 - 1._R_P)*(u23 - u1)
  endif
  ! computing right state
  if (u23>u4) then
    ! shock
    S4 = u4 + a4*sqrt(1._R_P + 0.5_R_P*(g4 + 1._R_P)/g4*(p23/p4-1._R_P))
    S3 = S4
  else
    ! rarefaction
    S4 = u4  + a4
    S3 = u23 + a4 + 0.5_R_P*(g4 - 1._R_P)*(u23 - u4)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine WavesSpeed1234up

  elemental subroutine InterStates23u(p1,u1,a1,g1,gm1_1,gp1_1,delta1,eta1, &
                                      p4,u4,a4,g4,gm1_4,gp1_4,delta4,eta4, &
                                      u23,                                 &
                                      r2,p2,a2,r3,p3,a3,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing intermediates states knowing the value of speed (u23) of intermediates states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1          ! Pressure of state 1.
  real(R_P), intent(IN)::  u1          ! Velocity of state 1.
  real(R_P), intent(IN)::  a1          ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1          ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  gm1_1,gp1_1 ! g-1, g+1 of state 1.
  real(R_P), intent(IN)::  delta1,eta1 ! (g-1)/2, 2*g/(g-1) of state 1.
  real(R_P), intent(IN)::  p4          ! Pressure of state 4.
  real(R_P), intent(IN)::  u4          ! Velocity of state 4.
  real(R_P), intent(IN)::  a4          ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4          ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  gm1_4,gp1_4 ! g-1, g+1 of state 4.
  real(R_P), intent(IN)::  delta4,eta4 ! (g-1)/2, 2*g/(g-1) of state 4.
  real(R_P), intent(IN)::  u23         ! Velocity of intermediate states.
  real(R_P), intent(OUT):: r2,p2,a2    ! Density, pressure and speed of sound of state 2.
  real(R_P), intent(OUT):: r3,p3,a3    ! Density, pressure and speed of sound of state 3.
  real(R_P), intent(OUT):: S1,S2       ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4       ! Right signal velocities.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left wave
  if (abs(u23-u1)<=Zero) then
    call Rarefaction(sgn   = -1._R_P, &
                     g     = g1,      &
                     delta = delta1,  &
                     eta   = eta1,    &
                     u0    = u1,      &
                     p0    = p1,      &
                     a0    = a1,      &
                     ux    = u23,     &
                     rx    = r2,      &
                     px    = p2,      &
                     ax    = a2,      &
                     s0    = S1,      &
                     sx    = S2)
  else
    if (u23<u1) then
      call Shock(sgn = -1._R_P, &
                 g   = g1,      &
                 gm1 = gm1_1,   &
                 gp1 = gp1_1,   &
                 u0  = u1,      &
                 p0  = p1,      &
                 a0  = a1,      &
                 ux  = u23,     &
                 rx  = r2,      &
                 px  = p2,      &
                 ax  = a2,      &
                 ss  = S1)
      S2 = S1
    else
      call Rarefaction(sgn   = -1._R_P, &
                       g     = g1,      &
                       delta = delta1,  &
                       eta   = eta1,    &
                       u0    = u1,      &
                       p0    = p1,      &
                       a0    = a1,      &
                       ux    = u23,     &
                       rx    = r2,      &
                       px    = p2,      &
                       ax    = a2,      &
                       s0    = S1,      &
                       sx    = S2)
    endif
  endif
  ! computing right wave
  if (abs(u23-u4)<=Zero) then
    call Rarefaction(sgn   = 1._R_P, &
                     g     = g4,     &
                     delta = delta4, &
                     eta   = eta4,   &
                     u0    = u4,     &
                     p0    = p4,     &
                     a0    = a4,     &
                     ux    = u23,    &
                     rx    = r3,     &
                     px    = p3,     &
                     ax    = a3,     &
                     s0    = S4,     &
                     sx    = S3)
  else
    if (u23>u4) then
      call Shock(sgn = 1._R_P, &
                 g   = g4,     &
                 gm1 = gm1_4,  &
                 gp1 = gp1_4,  &
                 u0  = u4,     &
                 p0  = p4,     &
                 a0  = a4,     &
                 ux  = u23,    &
                 rx  = r3,     &
                 px  = p3,     &
                 ax  = a3,     &
                 ss  = S4)
      S3 = S4
    else
      call Rarefaction(sgn   = 1._R_P, &
                       g     = g4,     &
                       delta = delta4, &
                       eta   = eta4,   &
                       u0    = u4,     &
                       p0    = p4,     &
                       a0    = a4,     &
                       ux    = u23,    &
                       rx    = r3,     &
                       px    = p3,     &
                       ax    = a3,     &
                       s0    = S4,     &
                       sx    = S3)
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine InterStates23u

  elemental subroutine InterStates23up(p1,u1,a1,g1,p4,u4,a4,g4,u23,p23,r2,r3,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing intermediates states knowing the value of speed and pressure (u23, p23) of intermediates states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1      ! Pressure of state 1.
  real(R_P), intent(IN)::  u1      ! Velocity of state 1.
  real(R_P), intent(IN)::  a1      ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1      ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4      ! Pressure of state 4.
  real(R_P), intent(IN)::  u4      ! Velocity of state 4.
  real(R_P), intent(IN)::  a4      ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4      ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  u23     ! Velocity of intermediate states.
  real(R_P), intent(IN)::  p23     ! Pressure of intermediate states.
  real(R_P), intent(OUT):: r2      ! Density of state 2.
  real(R_P), intent(OUT):: r3      ! Density of state 3.
  real(R_P), intent(OUT):: S1,S2   ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4   ! Right signal velocities.
  real(R_P)::              a2,a3   ! Speed of sound of state 2 and 3.
  real(R_P)::              gm1,gp1 ! g - 1, g + 1.
  real(R_P)::              p_p     ! Pressure ratio.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing left state
  if (u23<u1) then
    ! shock
    gm1 = g1 - 1._R_P
    gp1 = g1 + 1._R_P
    p_p = p23/p1
    a2  = a1*sqrt((gp1 + gm1*p_p)/(gp1 + gm1/p_p))
    r2  = r(p=p23,a=a2,g=g1)
    S1  = u1 - a1*sqrt(1._R_P + 0.5_R_P*gp1/g1*(p_p-1._R_P))
    S2  = S1
  else
    ! rarefaction
    a2 = a1 - 0.5_R_P*(g1-1._R_P)*(u23-u1)
    r2 = r(p=p23,a=a2,g=g1)
    S1 = u1  - a1
    S2 = u23 - a2
  endif
  ! computing right state
  if (u23>u4) then
    ! shock
    gm1 = g4 - 1._R_P
    gp1 = g4 + 1._R_P
    p_p = p23/p4
    a3  = a4*sqrt((gp1 + gm1*p_p)/(gp1 + gm1/p_p))
    r3  = r(p=p23,a=a3,g=g4)
    S4  = u4 + a4*sqrt(1._R_P + 0.5_R_P*gp1/g4*(p_p-1._R_P))
    S3  = S4
  else
    ! rarefaction
    a3 = a4 + 0.5_R_P*(g4-1._R_P)*(u23-u4)
    r3 = r(p=p23,a=a3,g=g4)
    S4 = u4  + a4
    S3 = u23 + a3
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine InterStates23up

  subroutine sampling_fluxes(p1,r1,u1,a1,g1,delta1,eta1, &
                             p4,r4,u4,a4,g4,delta4,eta4, &
                             p23,r2,r3,S,S1,S2,S3,S4,    &
                             F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for sampling interface solution and computing fluxes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1       ! Pressure of state 1.
  real(R_P), intent(IN)::  r1       ! Density of state 1.
  real(R_P), intent(IN)::  u1       ! Velocity of state 1.
  real(R_P), intent(IN)::  a1       ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1       ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  delta1   ! (g-1)/2 of state 1.
  real(R_P), intent(IN)::  eta1     ! 2*g/(g-1) of state 1.
  real(R_P), intent(IN)::  p4       ! Pressure of state 4.
  real(R_P), intent(IN)::  r4       ! Density of state 4.
  real(R_P), intent(IN)::  u4       ! Velocity of state 4.
  real(R_P), intent(IN)::  a4       ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4       ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  delta4   ! (g-1)/2 of state 4.
  real(R_P), intent(IN)::  eta4     ! 2*g/(g-1) of state 4.
  real(R_P), intent(IN)::  p23      ! Pressure of intermediate states.
  real(R_P), intent(IN)::  r2       ! Density of state 2.
  real(R_P), intent(IN)::  r3       ! Density of state 3.
  real(R_P), intent(IN)::  S        ! Velocity of intermediate states.
  real(R_P), intent(IN)::  S1,S2    ! Left signal velocities.
  real(R_P), intent(IN)::  S3,S4    ! Right signal velocities.
  real(R_P), intent(OUT):: F_r      ! Flux of mass conservation.
  real(R_P), intent(OUT):: F_u      ! Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E      ! Flux of energy conservation.
  real(R_P)::              pF,rF,aF ! Primitive variables at interface.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(minloc([-S1,S1*S2,S2*S,S*S3,S3*S4,S4],dim=1))
  case(1) ! left supersonic
    call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(2) ! left transonic
    aF = (a1 + u1*delta1)/(1._R_P + delta1)
    pF = p1*(aF/a1)**eta1
    rF = r(p=pF,a=aF,g=g1)
    call fluxes(p = pF, r = rF, u = aF, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(3) ! left subsonic
    call fluxes(p = p23, r = r2, u = S, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(4) ! right subsonic
    call fluxes(p = p23, r = r3, u = S, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  case(5) ! right transonic
    aF = (a4 - u4*delta4)/(1._R_P + delta4)
    pF = p4*(aF/a4)**eta4
    rF = r(p=pF,a=aF,g=g4)
    call fluxes(p = pF, r = rF, u = -aF, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  case(6) ! right supersonic
    call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine sampling_fluxes

  elemental subroutine sampling_fluxes_old(p1,r1,u1,a1,g1,delta1,eta1, &
                                           p4,r4,u4,a4,g4,delta4,eta4, &
                                           p23,r2,r3,S,S1,S2,S3,S4,    &
                                           F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for sampling interface solution and computing fluxes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1       ! Pressure of state 1.
  real(R_P), intent(IN)::  r1       ! Density of state 1.
  real(R_P), intent(IN)::  u1       ! Velocity of state 1.
  real(R_P), intent(IN)::  a1       ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1       ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  delta1   ! (g-1)/2 of state 1.
  real(R_P), intent(IN)::  eta1     ! 2*g/(g-1) of state 1.
  real(R_P), intent(IN)::  p4       ! Pressure of state 4.
  real(R_P), intent(IN)::  r4       ! Density of state 4.
  real(R_P), intent(IN)::  u4       ! Velocity of state 4.
  real(R_P), intent(IN)::  a4       ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4       ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  delta4   ! (g-1)/2 of state 4.
  real(R_P), intent(IN)::  eta4     ! 2*g/(g-1) of state 4.
  real(R_P), intent(IN)::  p23      ! Pressure of intermediate states.
  real(R_P), intent(IN)::  r2       ! Density of state 2.
  real(R_P), intent(IN)::  r3       ! Density of state 3.
  real(R_P), intent(IN)::  S        ! Velocity of intermediate states.
  real(R_P), intent(IN)::  S1,S2    ! Left signal velocities.
  real(R_P), intent(IN)::  S3,S4    ! Right signal velocities.
  real(R_P), intent(OUT):: F_r      ! Flux of mass conservation.
  real(R_P), intent(OUT):: F_u      ! Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E      ! Flux of energy conservation.
  real(R_P)::              pF,rF,aF ! Primitive variables at interface.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (S>toll) then
    ! left subsonic
    call fluxes(p = p23, r = r2, u = S, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
    if (S<=u1) then
      if (S1>0._R_P) then
        ! left supersonic
        call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
      endif
    else
      if (S1>0._R_P) then
        ! left supersonic
        call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
      else
        if (S2>0._R_P) then
          ! left transonic
          aF = (a1 + u1*delta1)/(1._R_P + delta1)
          pF = p1*(aF/a1)**eta1
          rF = r(p=pF,a=aF,g=g1)
          call fluxes(p = pF, r = rF, u = aF, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
        endif
      endif
    endif
  else
    ! right subsonic
    call fluxes(p = p23, r = r3, u = S, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
    if (S>=u4) then
      if (S4<0._R_P) then
        ! right supersonic
        call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
      endif
    else
      if (S4<0._R_P) then
        ! right supersonic
        call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
      else
        if (S3<0._R_P) then
          ! right transonic
          aF = (a4 - u4*delta4)/(1._R_P + delta4)
          pF = p4*(aF/a4)**eta4
          rF = r(p=pF,a=aF,g=g4)
          call fluxes(p = pF, r = rF, u = -aF, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
        endif
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine sampling_fluxes_old

  elemental subroutine sampling_fluxes_vacuum(p1,r1,u1,a1,g1,delta1,eta1, &
                                              p4,r4,u4,a4,g4,delta4,eta4, &
                                              S1,S2,S3,S4,                &
                                              F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for sampling interface solution and computing fluxes.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1       ! Pressure of state 1.
  real(R_P), intent(IN)::  r1       ! Density of state 1.
  real(R_P), intent(IN)::  u1       ! Velocity of state 1.
  real(R_P), intent(IN)::  a1       ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1       ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  delta1   ! (g-1)/2 of state 1.
  real(R_P), intent(IN)::  eta1     ! 2*g/(g-1) of state 1.
  real(R_P), intent(IN)::  p4       ! Pressure of state 4.
  real(R_P), intent(IN)::  r4       ! Density of state 4.
  real(R_P), intent(IN)::  u4       ! Velocity of state 4.
  real(R_P), intent(IN)::  a4       ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4       ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  delta4   ! (g-1)/2 of state 4.
  real(R_P), intent(IN)::  eta4     ! 2*g/(g-1) of state 4.
  real(R_P), intent(IN)::  S1,S2    ! Left signal velocities.
  real(R_P), intent(IN)::  S3,S4    ! Right signal velocities.
  real(R_P), intent(OUT):: F_r      ! Flux of mass conservation.
  real(R_P), intent(OUT):: F_u      ! Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E      ! Flux of energy conservation.
  real(R_P)::              pF,rF,aF ! Primitive variables at interface.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(minloc([-S1,S1*S2,S2*S3,S3*S4,S4],dim=1))
  case(1) ! left supersonic
    call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(2) ! left transonic
    aF = (a1 + u1*delta1)/(1._R_P + delta1)
    pF = p1*(aF/a1)**eta1
    rF = r(p=pF,a=aF,g=g1)
    call fluxes(p = pF, r = rF, u = aF, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(3) ! vacuum
    F_r = toll
    F_u = toll
    F_E = toll
  case(4) ! right transonic
    aF = (a4 - u4*delta4)/(1._R_P + delta4)
    pF = p4*(aF/a4)**eta4
    rF = r(p=pF,a=aF,g=g4)
    call fluxes(p = pF, r = rF, u = -aF, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  case(5) ! right supersonic
    call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine sampling_fluxes_vacuum

  ! subroutines for evaluating slope limiters
  elemental function VanLeer_limiter(Sm1,Sp1) result(limit)
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Symmetric (van Leer 1974) limiter based on van Leer function.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: Sm1,Sp1 ! Upwind (S[i-1]-S[i]) and downwind (S[i]-S[i+1]) differences.
  real(R_P)::             limit   ! Value of limit.
  real(R_P)::             r       ! Sm1/Sp1.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  if (Sm1*Sp1>=0._R_P) then
    r     = (Sm1+Zero)/(Sp1+Zero)
    limit = 2._R_P*r/(1._R_P + r)
  else
    limit = 0._R_P
  endif
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VanLeer_limiter

  elemental function VanAlbada_limiter(Sm1,Sp1) result(limit)
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Symmetric (van Albada 1982) limiter based on van Albada function.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: Sm1,Sp1 ! Upwind (S[i-1]-S[i]) and downwind (S[i]-S[i+1]) differences.
  real(R_P)::             limit   ! Value of limit.
  real(R_P)::             r       ! Sm1/Sp1.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  if (Sm1*Sp1>=0._R_P) then
    r     = (Sm1+Zero)/(Sp1+Zero)
    limit = (r*r+r)/(1._R_P+r*r)
  else
    limit = 0._R_P
  endif
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VanAlbada_limiter

  elemental function Harten_limiter(Sm1,Sp1) result(limit)
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Symmetric limiter based on Harten's shock switch.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: Sm1,Sp1 ! Upwind (S[i-1]-S[i]) and downwind (S[i]-S[i+1]) differences.
  real(R_P)::             limit   ! Value of limit.
  real(R_P)::             r       ! Sm1/Sp1.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  if (Sm1*Sp1>=0._R_P) then
    r     = (Sm1+Zero)/(Sp1+Zero)
    limit = 1._R_P - abs(1._R_P - r)/(1._R_P + r)
  else
    limit = 0._R_P
  endif
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction Harten_limiter

  ! subroutines for evaluating the smoothness of cell or interface
  !> Function for evaluating the smoothness of an interface using a Riemann solver-like algorithm.
  function chk_smooth_z(p1,r1,u1,g1,p4,r4,u4,g4) result(smooth)
  !-------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: p1     !< Pressure of state 1.
  real(R_P), intent(IN):: r1     !< Density of state 1.
  real(R_P), intent(IN):: u1     !< Velocity of state 1.
  real(R_P), intent(IN):: g1     !< Specific heats ratio of state 1.
  real(R_P), intent(IN):: p4     !< Pressure of state 4.
  real(R_P), intent(IN):: r4     !< Density of state 4.
  real(R_P), intent(IN):: u4     !< Velocity of state 4.
  real(R_P), intent(IN):: g4     !< Specific heats ratio of state 4.
  logical::               smooth !< Logical flag for testing the smoothness of the stencil.
  real(R_P)::             a1,a4  !< Speeds of sound of states 1 and 4.
  real(R_P)::             pmin   !< min(p1,p4).
  real(R_P)::             pmax   !< max(p1,p4).
  real(R_P)::             u23    !< Velocity of intermediate states.
  real(R_P)::             p23    !< Pressure of intermediate states.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound of state 1 and 4
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  ! evaluating intermediate states pressure and velocity
  pmin = min(p1,p4)
  pmax = max(p1,p4)
  call PVL_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = u23, p23 = p23)
  if (pmax/pmin>2._R_P.OR.pmin>p23.OR.p23>pmax) then
    if (p23<pmin) then
      call TR_u23(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = u23)
    else
      call TS_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, u23 = u23, p23 = p23)
    endif
  endif
  smooth = ((u23>=u1).AND.(u23<=u4))
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endfunction chk_smooth_z

  !> Function for evaluating the smoothness of a cell using a discontinuity detector based on slope limiter function.
  function chk_smooth_limiter(pm1,pp1,rm1,rp1,um1,up1,gm1,gp1) result(smooth)
  !-------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: pm1     !< Upwind difference of pressure.
  real(R_P), intent(IN):: pp1     !< Downwind difference of pressure.
  real(R_P), intent(IN):: rm1     !< Upwind difference of density.
  real(R_P), intent(IN):: rp1     !< Downwind difference of density.
  real(R_P), intent(IN):: um1     !< Upwind difference of speed.
  real(R_P), intent(IN):: up1     !< Downwind difference of speed.
  real(R_P), intent(IN):: gm1     !< Upwind difference of specific heats ratio.
  real(R_P), intent(IN):: gp1     !< Downwind difference of specific heats ratio.
  logical::               smooth  !< Logical flag for testing the smoothness of the stencil.
  real(R_P)::             p_limit !< Smoothness limit (1 => smooth) of stencil according to pressure values.
  real(R_P)::             r_limit !< Smoothness limit (1 => smooth) of stencil according to density values.
  real(R_P)::             u_limit !< Smoothness limit (1 => smooth) of stencil according to speed values.
  real(R_P)::             g_limit !< Smoothness limit (1 => smooth) of stencil according to g values.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
#ifdef SMSWvanleer
  p_limit = VanLeer_limiter(Sm1 = pm1, Sp1 = pp1)
  r_limit = VanLeer_limiter(Sm1 = rm1, Sp1 = rp1)
  u_limit = VanLeer_limiter(Sm1 = um1, Sp1 = up1)
  g_limit = VanLeer_limiter(Sm1 = gm1, Sp1 = gp1)
#elif SMSWharten
  p_limit = Harten_limiter(Sm1 = pm1, Sp1 = pp1)
  r_limit = Harten_limiter(Sm1 = rm1, Sp1 = rp1)
  u_limit = Harten_limiter(Sm1 = um1, Sp1 = up1)
  g_limit = Harten_limiter(Sm1 = gm1, Sp1 = gp1)
#else
  p_limit = VanAlbada_limiter(Sm1 = pm1, Sp1 = pp1)
  r_limit = VanAlbada_limiter(Sm1 = rm1, Sp1 = rp1)
  u_limit = VanAlbada_limiter(Sm1 = um1, Sp1 = up1)
  g_limit = VanAlbada_limiter(Sm1 = gm1, Sp1 = gp1)
#endif
  smooth  = ((int(p_limit + 0.1_R_P)==1).AND. & ! degenerate wave (flow     discontinuity)
             (int(r_limit + 0.1_R_P)==1).AND. & ! degenerate wave (species  discontinuity)
             (int(u_limit + 0.1_R_P)==1).AND. & ! non linear wave (momentum discontinuity)
             (int(g_limit + 0.1_R_P)==1))       ! non linear wave (energy   discontinuity)
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endfunction chk_smooth_limiter

  !> Function for evaluating the smoothness of a cell using the Liu's entropy conditions.
  function chk_smooth_liu(pm1,rm1,um1,gm1,cpm1,cvm1,     &
                          p1, r1, u1m,u1p, g1, cp1, cv1, &
                          pp1,rp1,up1,gp1,cpp1,cvp1) result(smooth)
  !-------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: pm1      !< Pressure of state i-1.
  real(R_P), intent(IN):: rm1      !< Density of state i-1.
  real(R_P), intent(IN):: um1      !< Velocity of state i-1.
  real(R_P), intent(IN):: gm1      !< Specific heats ratio of state i-1.
  real(R_P), intent(IN):: cpm1     !< Specific heats cp of state i-1.
  real(R_P), intent(IN):: cvm1     !< Specific heats cv of state i-1.
  real(R_P), intent(IN):: p1       !< Pressure of state i.
  real(R_P), intent(IN):: r1       !< Density of state i.
  real(R_P), intent(IN):: u1m      !< Velocity of state i at left interface.
  real(R_P), intent(IN):: u1p      !< Velocity of state i at left interface.
  real(R_P), intent(IN):: g1       !< Specific heats ratio of state i.
  real(R_P), intent(IN):: cp1      !< Specific heats cp of state i.
  real(R_P), intent(IN):: cv1      !< Specific heats cv of state i.
  real(R_P), intent(IN):: pp1      !< Pressure of state i+1.
  real(R_P), intent(IN):: rp1      !< Density of state i+1.
  real(R_P), intent(IN):: up1      !< Velocity of state i+1.
  real(R_P), intent(IN):: gp1      !< Specific heats ratio of state i+1.
  real(R_P), intent(IN):: cpp1     !< Specific heats cp of state i+1.
  real(R_P), intent(IN):: cvp1     !< Specific heats cv of state i+1.
  logical::               smooth   !< Logical flag for testing the smoothness of the stencil.
  real(R_P)::             a1,a4    !< Speeds of sound of states 1 and 4.
  real(R_P)::             u_roe    !< Velocity of Roe average state.
  real(R_P)::             a_roe    !< Speed of sound of Roe average state.
  logical::               smooth_l !< Logical flag for testing the smoothness across the left interface.
  logical::               smooth_r !< Logical flag for testing the smoothness across the right interface.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  ! verifing the presence of a shock across the left interface
  call Roe_ua(p1 = pm1, r1 = rm1, u1 = um1, g1 = gm1, cp1 = cpm1, cv1 = cvm1, &
              p4 = p1,  r4 = r1,  u4 = u1m, g4 = g1,  cp4 = cp1,  cv4 = cv1,  &
              u_roe = u_roe, a_roe = a_roe)
  a1 = a(p=pm1,r=rm1,g=gm1)
  a4 = a(p=p1, r=r1, g=g1 )
  smooth_l = ((abs(u1m + a4) < abs(u_roe + a_roe)).AND.(abs(u_roe + a_roe) < abs(um1 + a1)).AND.&
              (abs(u1m - a4) < abs(u_roe - a_roe)).AND.(abs(u_roe - a_roe) < abs(um1 - a1)))
  ! verifing the presence of a shock across the right interface
  call Roe_ua(p1 = p1,  r1 = r1,  u1 = u1p, g1 = g1,  cp1 = cp1,  cv1 = cv1,  &
              p4 = pp1, r4 = rp1, u4 = up1, g4 = gp1, cp4 = cpp1, cv4 = cvp1, &
              u_roe = u_roe, a_roe = a_roe)
  a1 = a(p=p1, r=r1, g=g1 )
  a4 = a(p=pp1,r=rp1,g=gp1)
  smooth_r = ((abs(up1 + a4) < abs(u_roe + a_roe)).AND.(abs(u_roe + a_roe) < abs(u1p + a1)).AND.&
              (abs(up1 - a4) < abs(u_roe - a_roe)).AND.(abs(u_roe - a_roe) < abs(u1p - a1)))

  smooth = (smooth_l.AND.smooth_r)

  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endfunction chk_smooth_liu

  ! subroutines for evaluating intermediate states
  elemental subroutine PVL_u23(p1,r1,u1,a1,p4,r4,u4,a4,u23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating speed of intermediate states using PVL (Primitive Variables Linearization) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1  ! Pressure of state 1.
  real(R_P), intent(IN)::  r1  ! Density of state 1.
  real(R_P), intent(IN)::  u1  ! Velocity of state 1.
  real(R_P), intent(IN)::  a1  ! Speed of sound of state 1.
  real(R_P), intent(IN)::  p4  ! Pressure of state 4.
  real(R_P), intent(IN)::  r4  ! Density of state 4.
  real(R_P), intent(IN)::  u4  ! Velocity of state 4.
  real(R_P), intent(IN)::  a4  ! Speed of sound of state 4.
  real(R_P), intent(OUT):: u23 ! Velocity of intermediate states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  u23 = 0.5_R_P*(u1 + u4) - 2.0_R_P*(p4 - p1)/((r1 + r4)*(a1 + a4))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_u23

  elemental subroutine PVL_p23(p1,r1,u1,a1,p4,r4,u4,a4,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating pressure of intermediate states using PVL (Primitive Variables Linearization) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1  ! Pressure of state 1.
  real(R_P), intent(IN)::  r1  ! Density of state 1.
  real(R_P), intent(IN)::  u1  ! Velocity of state 1.
  real(R_P), intent(IN)::  a1  ! Speed of sound of state 1.
  real(R_P), intent(IN)::  p4  ! Pressure of state 4.
  real(R_P), intent(IN)::  r4  ! Density of state 4.
  real(R_P), intent(IN)::  u4  ! Velocity of state 4.
  real(R_P), intent(IN)::  a4  ! Speed of sound of state 4.
  real(R_P), intent(OUT):: p23 ! Pressure of intermediate states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  p23 = 0.5_R_P*((p1 + p4) - 0.25_R_P*(u4 - u1)*(r1 + r4)*(a1 + a4))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_p23

  elemental subroutine PVL_up23(p1,r1,u1,a1,p4,r4,u4,a4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating speed and pressure of intermediate states using PVL (Primitive Variables Linearization) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1  ! Pressure of state 1.
  real(R_P), intent(IN)::  r1  ! Density of state 1.
  real(R_P), intent(IN)::  u1  ! Velocity of state 1.
  real(R_P), intent(IN)::  a1  ! Speed of sound of state 1.
  real(R_P), intent(IN)::  p4  ! Pressure of state 4.
  real(R_P), intent(IN)::  r4  ! Density of state 4.
  real(R_P), intent(IN)::  u4  ! Velocity of state 4.
  real(R_P), intent(IN)::  a4  ! Speed of sound of state 4.
  real(R_P), intent(OUT):: u23 ! Velocity of intermediate states.
  real(R_P), intent(OUT):: p23 ! Pressure of intermediate states.
  real(R_P)::              ram ! Mean value of r*a.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ram = 0.25_R_P*(r1 + r4)*(a1 + a4)
  u23 = 0.5_R_P*((u1 + u4) - (p4 - p1)/ram)
  p23 = 0.5_R_P*((p1 + p4) - (u4 - u1)*ram)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_up23

  elemental subroutine PVL_WavesSpeed14(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,S,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the fastest (1,4) waves speed using PVL (Primitive Variables Linearization) approximation.
  ! Two variants are provided: the first uses both the approximations of u and p provided by the PVL approximation for computing
  ! the waves speed (WSup algorithm), while the second uses only the approximation of u (WSu algorithm).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1  ! Pressure of state 1.
  real(R_P), intent(IN)::  r1  ! Density of state 1.
  real(R_P), intent(IN)::  u1  ! Velocity of state 1.
  real(R_P), intent(IN)::  a1  ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1  ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4  ! Pressure of state 4.
  real(R_P), intent(IN)::  r4  ! Density of state 4.
  real(R_P), intent(IN)::  u4  ! Velocity of state 4.
  real(R_P), intent(IN)::  a4  ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4  ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: S   ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1  ! Fastest left signal velocities.
  real(R_P), intent(OUT):: S4  ! Fastest right signal velocities.
#ifdef WSu
#else
  real(R_P)::              p23 ! Pressure of intermediate states 2 and 3.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef WSu
  ! evaluating intermediate states velocity
  call PVL_u23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = S)
  ! computing waves speed 1 and 4
  call WavesSpeed14u(u1 = u1, a1 = a1, g1 = g1, u4 = u4, a4 = a4, g4 = g4, u23 = S, S1 = S1, S4 = S4)
#else
  ! evaluating intermediate states pressure and velocity
  call PVL_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = S, p23 = p23)
  ! computing waves speed 1 and 4
  call WavesSpeed14up(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23, S1 = S1, S4 = S4)
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_WavesSpeed14

  elemental subroutine PVL_WavesSpeed1234(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,S,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating waves speed using PVL (Primitive Variables Linearization) approximation.
  ! Two variants are provided: the first uses both the approximations of u and p provided by the PVL approximation for computing
  ! the waves speed (WSup algorithm), while the second uses only the approximation of u (WSu algorithm).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    ! Pressure of state 1.
  real(R_P), intent(IN)::  r1    ! Density of state 1.
  real(R_P), intent(IN)::  u1    ! Velocity of state 1.
  real(R_P), intent(IN)::  a1    ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1    ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4    ! Pressure of state 4.
  real(R_P), intent(IN)::  r4    ! Density of state 4.
  real(R_P), intent(IN)::  u4    ! Velocity of state 4.
  real(R_P), intent(IN)::  a4    ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4    ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: S     ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1,S2 ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4 ! Right signal velocities.
#ifdef WSu
#else
  real(R_P)::              p23   ! Pressure of intermediate states 2 and 3.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef WSu
  ! evaluating intermediate states pressure and velocity
  call PVL_u23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = S)
  ! computing waves speed
  call WavesSpeed1234u(u1 = u1, a1 = a1, g1 = g1, &
                       u4 = u4, a4 = a4, g4 = g4, &
                       u23 = S,                   &
                       S1 = S1, S2 = S2, S3 = S3, S4 = S4)
#else
  ! evaluating intermediate states pressure and velocity
  call PVL_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = S, p23 = p23)
  ! computing waves speed
  call WavesSpeed1234up(p1 = p1, u1 = u1, a1 = a1, g1 = g1, &
                        p4 = p4, u4 = u4, a4 = a4, g4 = g4, &
                        u23 = S, p23 = p23,                 &
                        S1 = S1, S2 = S2, S3 = S3, S4 = S4)
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_WavesSpeed1234

  elemental subroutine PVL_States(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,u23,p23,r2,r3,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the intermediate states using PVL (Primitive Variables Linearization) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    ! Pressure of state 1.
  real(R_P), intent(IN)::  r1    ! Density of state 1.
  real(R_P), intent(IN)::  u1    ! Velocity of state 1.
  real(R_P), intent(IN)::  a1    ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1    ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4    ! Pressure of state 4.
  real(R_P), intent(IN)::  r4    ! Density of state 4.
  real(R_P), intent(IN)::  u4    ! Velocity of state 4.
  real(R_P), intent(IN)::  a4    ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4    ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: u23   ! Velocity of intermediate states.
  real(R_P), intent(OUT):: p23   ! Pressure of intermediate states.
  real(R_P), intent(OUT):: r2    ! Density of state 2.
  real(R_P), intent(OUT):: r3    ! Density of state 3.
  real(R_P), intent(OUT):: S1,S2 ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4 ! Right signal velocities.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call PVL_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = u23, p23 = p23)
  ! computing intermediate states
  call InterStates23up(p1 = p1, u1 = u1, a1 = a1, g1 = g1,     &
                       p4 = p4, u4 = u4, a4 = a4, g4 = g4,     &
                       u23 = u23, p23 = p23, r2 = r2, r3 = r3, &
                       S1 = S1, S2 = S2, S3 = S3, S4 = S4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_States

  elemental subroutine CVL_u23(p1,r1,u1,a1,p4,r4,u4,a4,u23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating speed of intermediate states using CVL (Charactheristic Variables Linearization) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1  ! Pressure of state 1.
  real(R_P), intent(IN)::  r1  ! Density of state 1.
  real(R_P), intent(IN)::  u1  ! Velocity of state 1.
  real(R_P), intent(IN)::  a1  ! Speed of sound of state 1.
  real(R_P), intent(IN)::  p4  ! Pressure of state 4.
  real(R_P), intent(IN)::  r4  ! Density of state 4.
  real(R_P), intent(IN)::  u4  ! Velocity of state 4.
  real(R_P), intent(IN)::  a4  ! Speed of sound of state 4.
  real(R_P), intent(OUT):: u23 ! Velocity of intermediate states.
  real(R_P)::              c1  ! Caracteristic of state 1.
  real(R_P)::              c4  ! Caracteristic of state 4.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  c1  = r1*a1
  c4  = r4*a4
  u23 = (c1*u1 + c4*u4 + (p1 - p4))/(c1 + c4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine CVL_u23

  elemental subroutine CVL_up23(p1,r1,u1,a1,p4,r4,u4,a4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating speed and pressure of intermediate states using CVL (Charactheristic Variables Linearization)
  ! approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1      ! Pressure of state 1.
  real(R_P), intent(IN)::  r1      ! Density of state 1.
  real(R_P), intent(IN)::  u1      ! Velocity of state 1.
  real(R_P), intent(IN)::  a1      ! Speed of sound of state 1.
  real(R_P), intent(IN)::  p4      ! Pressure of state 4.
  real(R_P), intent(IN)::  r4      ! Density of state 4.
  real(R_P), intent(IN)::  u4      ! Velocity of state 4.
  real(R_P), intent(IN)::  a4      ! Speed of sound of state 4.
  real(R_P), intent(OUT):: u23     ! Velocity of intermediate states.
  real(R_P), intent(OUT):: p23     ! Pressure of intermediate states.
  real(R_P)::              c1      ! Caracteristic of state 1.
  real(R_P)::              c4      ! Caracteristic of state 4.
  real(R_P)::              c14_inv ! 1/(c1+c4).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  c1      = r1*a1
  c4      = r4*a4
  c14_inv = 1._R_P/(c1 + c4)
  u23     = c14_inv*(c1*u1 + c4*u4 +       (p1 - p4))
  p23     = c14_inv*(c4*p1 + c1*p4 + c1*c4*(u1 - u4))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine CVL_up23

  elemental subroutine CVL_WavesSpeed14(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,S,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the fastest (1,4) waves speed using CVL (Charactheristic Variables Linearization) approximation.
  ! Two variants are provided: the first uses both the approximations of u and p provided by the CVL approximation for computing
  ! the waves speed (WSup algorithm), while the second uses only the approximation of u (WSu algorithm).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1  ! Pressure of state 1.
  real(R_P), intent(IN)::  r1  ! Density of state 1.
  real(R_P), intent(IN)::  u1  ! Velocity of state 1.
  real(R_P), intent(IN)::  a1  ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1  ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4  ! Pressure of state 4.
  real(R_P), intent(IN)::  r4  ! Density of state 4.
  real(R_P), intent(IN)::  u4  ! Velocity of state 4.
  real(R_P), intent(IN)::  a4  ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4  ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: S   ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1  ! Fastest left signal velocities.
  real(R_P), intent(OUT):: S4  ! Fastest right signal velocities.
#ifdef WSu
#else
  real(R_P)::              p23 ! Pressure of intermediate states 2 and 3.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef WSu
  ! evaluating intermediate states pressure and velocity
  call CVL_u23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = S)
  ! computing waves speed 1 and 4
  call WavesSpeed14u(u1 = u1, a1 = a1, g1 = g1, u4 = u4, a4 = a4, g4 = g4, u23 = S, S1 = S1, S4 = S4)
#else
  ! evaluating intermediate states pressure and velocity
  call CVL_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = S, p23 = p23)
  ! computing waves speed 1 and 4
  call WavesSpeed14up(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23, S1 = S1, S4 = S4)
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine CVL_WavesSpeed14

  elemental subroutine CVL_States(p1,r1,u1,g1,p4,r4,u4,g4,u23,p23,r2,r3,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the intermediate states using CVL (Charactheristic Variables Linearization) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    ! Pressure of state 1.
  real(R_P), intent(IN)::  r1    ! Density of state 1.
  real(R_P), intent(IN)::  u1    ! Velocity of state 1.
  real(R_P), intent(IN)::  g1    ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4    ! Pressure of state 4.
  real(R_P), intent(IN)::  r4    ! Density of state 4.
  real(R_P), intent(IN)::  u4    ! Velocity of state 4.
  real(R_P), intent(IN)::  g4    ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: u23   ! Velocity of intermediate states.
  real(R_P), intent(OUT):: p23   ! Pressure of intermediate states.
  real(R_P), intent(OUT):: r2    ! Density of state 2.
  real(R_P), intent(OUT):: r3    ! Density of state 3.
  real(R_P), intent(OUT):: S1,S2 ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4 ! Right signal velocities.
  real(R_P)::              a1,a4 ! Speed of sound of state 1 and 4.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  call CVL_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = u23, p23 = p23)
  ! computing intermediate states
  call InterStates23up(p1 = p1, u1 = u1, a1 = a1, g1 = g1,     &
                       p4 = p4, u4 = u4, a4 = a4, g4 = g4,     &
                       u23 = u23, p23 = p23, r2 = r2, r3 = r3, &
                       S1 = S1, S2 = S2, S3 = S3, S4 = S4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine CVL_States

  elemental subroutine TR_u23(p1,u1,a1,g1,p4,u4,a4,g4,u23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating speed of intermediate states using TR (Two Rarefactions) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1        ! Pressure of state 1.
  real(R_P), intent(IN)::  u1        ! Velocity of state 1.
  real(R_P), intent(IN)::  a1        ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1        ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4        ! Pressure of state 4.
  real(R_P), intent(IN)::  u4        ! Velocity of state 4.
  real(R_P), intent(IN)::  a4        ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4        ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: u23       ! Velocity of intermediate states.
  real(R_P)::              gm,gm1    ! Mean value of specific heats ratios and g-1.
  real(R_P)::              dm        ! Mean value of d = (g-1)/2.
  real(R_P)::              zm,zm_inv ! Mean value of z = d/g and 1/z.
  real(R_P)::              p14       ! p14 = (p1/p4)**zm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  gm     = 0.5_R_P*(g1 + g4)
  gm1    = gm - 1._R_P
  dm     = 0.5_R_P*gm1
  zm     = dm/gm
  zm_inv = 1._R_P/zm
  p14    = (p1/p4)**zm
  u23    = (p14*u1/a1 + u4/a4 + 2._R_P*(p14 - 1._R_P)/gm1)/(p14/a1 + 1._R_P/a4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TR_u23

  elemental subroutine TR_up23(p1,u1,a1,g1,p4,u4,a4,g4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating speed and pressure of intermediate states using TR (Two Rarefactions) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1        ! Pressure of state 1.
  real(R_P), intent(IN)::  u1        ! Velocity of state 1.
  real(R_P), intent(IN)::  a1        ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1        ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4        ! Pressure of state 4.
  real(R_P), intent(IN)::  u4        ! Velocity of state 4.
  real(R_P), intent(IN)::  a4        ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4        ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: u23       ! Velocity of intermediate states.
  real(R_P), intent(OUT):: p23       ! Pressure of intermediate states.
  real(R_P)::              gm,gm1    ! Mean value of specific heats ratios and g-1.
  real(R_P)::              dm        ! Mean value of d = (g-1)/2.
  real(R_P)::              zm,zm_inv ! Mean value of z = d/g and 1/z.
  real(R_P)::              p14       ! p14 = (p1/p4)**zm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  gm     = 0.5_R_P*(g1 + g4)
  gm1    = gm - 1._R_P
  dm     = 0.5_R_P*gm1
  zm     = dm/gm
  zm_inv = 1._R_P/zm
  p14    = (p1/p4)**zm
  u23    = (p14*u1/a1 + u4/a4 + 2._R_P*(p14 - 1._R_P)/gm1)/(p14/a1 + 1._R_P/a4)
  p23    = 0.5_R_P*(p1*(1._R_P + dm*(u1 - u23)/a1)**zm_inv + p4*(1._R_P + dm*(u23 - u4)/a4)**zm_inv)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TR_up23

  elemental subroutine TR_WavesSpeed14(p1,u1,a1,g1,p4,u4,a4,g4,S,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the fastest (1,4) waves speed using TR (Two Rarefactions) approximation.
  ! Two variants are provided: the first uses both the approximations of u and p provided by the TR approximation for computing
  ! the waves speed (WSup algorithm), while the second uses only the approximation of u (WSu algorithm).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1  ! Pressure of state 1.
  real(R_P), intent(IN)::  u1  ! Velocity of state 1.
  real(R_P), intent(IN)::  a1  ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1  ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4  ! Pressure of state 4.
  real(R_P), intent(IN)::  u4  ! Velocity of state 4.
  real(R_P), intent(IN)::  a4  ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4  ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: S   ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1  ! Fastest left signal velocities.
  real(R_P), intent(OUT):: S4  ! Fastest right signal velocities.
#ifdef WSu
#else
  real(R_P)::              p23 ! Pressure of intermediate states 2 and 3.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef WSu
  ! evaluating intermediate states pressure and velocity
  call TR_u23(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S)
  ! computing waves speed 1 and 4
  call WavesSpeed14u(u1 = u1, a1 = a1, g1 = g1, u4 = u4, a4 = a4, g4 = g4, u23 = S, S1 = S1, S4 = S4)
#else
  ! evaluating intermediate states pressure and velocity
  call TR_up23(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23)
  ! computing waves speed 1 and 4
  call WavesSpeed14up(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23, S1 = S1, S4 = S4)
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TR_WavesSpeed14

  elemental subroutine TR_States(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,u23,p23,r2,r3,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating intermediate states using TR (Two Rarefactions) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    ! Pressure of state 1.
  real(R_P), intent(IN)::  r1    ! Density of state 1.
  real(R_P), intent(IN)::  u1    ! Velocity of state 1.
  real(R_P), intent(IN)::  a1    ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1    ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4    ! Pressure of state 4.
  real(R_P), intent(IN)::  r4    ! Density of state 4.
  real(R_P), intent(IN)::  u4    ! Velocity of state 4.
  real(R_P), intent(IN)::  a4    ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g4    ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: u23   ! Velocity of intermediate states.
  real(R_P), intent(OUT):: p23   ! Pressure of intermediate states.
  real(R_P), intent(OUT):: r2    ! Density of state 2.
  real(R_P), intent(OUT):: r3    ! Density of state 3.
  real(R_P), intent(OUT):: S1,S2 ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4 ! Right signal velocities.
  real(R_P)::              a2,a3 ! Speed of sound of state 2 and 3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  call TR_up23(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = u23, p23 = p23)
  ! evaluating left state
  r2 = r1*(p23/p1)**(1._R_P/g1)
  a2 = a(p=p23,r=r2,g=g1)
  S1 = u1  - a1
  S2 = u23 - a2
  ! evaluating right state
  r3 = r4*(p23/p4)**(1._R_P/g4)
  a3 = a(p=p23,r=r3,g=g4)
  S4 = u4  + a4
  S3 = u23 + a3
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TR_States

  elemental subroutine TS_up23(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating speed and pressure of intermediate states using TS (Two Shocks) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1          ! Pressure of state 1.
  real(R_P), intent(IN)::  r1          ! Density of state 1.
  real(R_P), intent(IN)::  u1          ! Velocity of state 1.
  real(R_P), intent(IN)::  a1          ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1          ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4          ! Pressure of state 4.
  real(R_P), intent(IN)::  r4          ! Density of state 4.
  real(R_P), intent(IN)::  u4          ! Velocity of state 4.
  real(R_P), intent(IN)::  a4          ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4          ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: u23         ! Velocity of intermediate states.
  real(R_P), intent(OUT):: p23         ! Pressure of intermediate states.
  real(R_P)::              g1p,g4p     ! Dummy variables for computing TS pressure of intermediate states.
  real(R_P)::              gp1,gm1_gp1 ! g+1, (g-1)/(g+1).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call PVL_p23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, p23 = p23)
  gp1 = g1 + 1._R_P ; gm1_gp1 = (g1 - 1._R_P)/gp1 ; g1p = sqrt((2._R_P/(gp1*r1))/(p23 + gm1_gp1*p1))
  gp1 = g4 + 1._R_P ; gm1_gp1 = (g4 - 1._R_P)/gp1 ; g4p = sqrt((2._R_P/(gp1*r4))/(p23 + gm1_gp1*p4))
  p23 = (g1p*p1 + g4p*p4 + u1 - u4)/(g1p + g4p)
  u23 = 0.5_R_P*(u1 + u4 + (p23 - p4)*g4p - (p23 - p1)*g1p)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TS_up23

  elemental subroutine TS_WavesSpeed14(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,S,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the fastest (1,4) waves speed using TS (Two Shocks) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1  ! Pressure of state 1.
  real(R_P), intent(IN)::  r1  ! Density of state 1.
  real(R_P), intent(IN)::  u1  ! Velocity of state 1.
  real(R_P), intent(IN)::  a1  ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1  ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4  ! Pressure of state 4.
  real(R_P), intent(IN)::  r4  ! Density of state 1.
  real(R_P), intent(IN)::  u4  ! Velocity of state 4.
  real(R_P), intent(IN)::  a4  ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4  ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: S   ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1  ! Fastest left signal velocities.
  real(R_P), intent(OUT):: S4  ! Fastest right signal velocities.
  real(R_P)::              p23 ! Pressure of intermediate states 2 and 3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  call TS_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23)
  ! computing waves speed 1 and 4
  call WavesSpeed14up(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23, S1 = S1, S4 = S4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TS_WavesSpeed14

  elemental subroutine TS_States(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,u23,p23,r2,r3,S1,S2,S3,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating intermediate states using TS (Two Shocks) approximation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1      ! Pressure of state 1.
  real(R_P), intent(IN)::  r1      ! Density of state 1.
  real(R_P), intent(IN)::  u1      ! Velocity of state 1.
  real(R_P), intent(IN)::  a1      ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1      ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4      ! Pressure of state 4.
  real(R_P), intent(IN)::  r4      ! Density of state 4.
  real(R_P), intent(IN)::  u4      ! Velocity of state 4.
  real(R_P), intent(IN)::  a4      ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g4      ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: u23     ! Velocity of intermediate states.
  real(R_P), intent(OUT):: p23     ! Pressure of intermediate states.
  real(R_P), intent(OUT):: r2      ! Density of state 2.
  real(R_P), intent(OUT):: r3      ! Density of state 3.
  real(R_P), intent(OUT):: S1,S2   ! Left signal velocities.
  real(R_P), intent(OUT):: S3,S4   ! Right signal velocities.
  real(R_P)::              gp1,gm1 ! g+1, g-1.
  real(R_P)::              p_p     ! Pressure ratio.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  call TS_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, u23 = u23, p23 = p23)
  ! evaluating left state
  p_p = p23/p1
  gp1 = g1 + 1._R_P
  gm1 = g1 - 1._R_P
  r2  = r1*((gm1 + gp1*p_p)/(gp1 + gm1*p_p))
  S1  = u1  - a1*sqrt(1._R_P + 0.5_R_P*gp1/g1*(p_p-1._R_P))
  S2  = S1
  ! evaluating right state
  p_p = p23/p4
  gp1 = g4 + 1._R_P
  gm1 = g4 - 1._R_P
  r3  = r4*((gm1 + gp1*p_p)/(gp1 + gm1*p_p))
  S4  = u4  + a4*sqrt(1._R_P + 0.5_R_P*gp1/g4*(p_p-1._R_P))
  S3  = S4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TS_States

  elemental subroutine Roe_up23(p1,r1,u1,g1,cp1,cv1,p4,r4,u4,g4,cp4,cv4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating speed and pressure of intermediate states using Roe averages.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1     ! Pressure of state 1.
  real(R_P), intent(IN)::  r1     ! Density of state 1.
  real(R_P), intent(IN)::  u1     ! Velocity of state 1.
  real(R_P), intent(IN)::  g1     ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  cp1    ! Specific heats cp of state 1.
  real(R_P), intent(IN)::  cv1    ! Specific heats cv of state 1.
  real(R_P), intent(IN)::  p4     ! Pressure of state 4.
  real(R_P), intent(IN)::  r4     ! Density of state 4.
  real(R_P), intent(IN)::  u4     ! Velocity of state 4.
  real(R_P), intent(IN)::  g4     ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  cp4    ! Specific heats cp of state 4.
  real(R_P), intent(IN)::  cv4    ! Specific heats cv of state 4.
  real(R_P), intent(OUT):: u23    ! Velocity of Roe average state.
  real(R_P), intent(OUT):: p23    ! Pressure of Roe average state.
  real(R_P)::              x,omx  ! x = sqrt(r1)/(sqrt(r1)+sqrt(r4)),  omx = 1-x = sqrt(r4)/(sqrt(r1)+sqrt(r4)).
  real(R_P)::              r_roe  ! Density of Roe average state.
  real(R_P)::              cp_roe ! Specific heats cp of Roe average state
  real(R_P)::              cv_roe ! Specific heats cv of Roe average state
  real(R_P)::              g_roe  ! Specific heats ratio of Roe average state.
  real(R_P)::              H_roe  ! Total specific entalpy of Roe average state
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x = sqrt(r1)/(sqrt(r1)+sqrt(r4)) ; omx = 1._R_P-x

  r_roe  = sqrt(r1*r4)
  u23    =                     u1*x +                     u4*omx
  H_roe  = H(p=p1,r=r1,u=u1,g=g1)*x + H(p=p4,r=r4,u=u4,g=g4)*omx
  cp_roe =                    cp1*x +                    cp4*omx
  cv_roe =                    cv1*x +                    cv4*omx
  g_roe  = cp_roe/cv_roe
  p23    = (H_roe-0.5_R_P*u23*u23)*(g_roe-1._R_P)/g_roe*r_roe
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Roe_up23

  elemental subroutine Roe_WavesSpeed14(p1,r1,u1,g1,cp1,cv1,p4,r4,u4,g4,cp4,cv4,S,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the fastest (1,4) waves speed using Roe averages.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    ! Pressure of state 1.
  real(R_P), intent(IN)::  r1    ! Density of state 1.
  real(R_P), intent(IN)::  u1    ! Velocity of state 1.
  real(R_P), intent(IN)::  g1    ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  cp1   ! Specific heats cp of state 1.
  real(R_P), intent(IN)::  cv1   ! Specific heats cv of state 1.
  real(R_P), intent(IN)::  p4    ! Pressure of state 4.
  real(R_P), intent(IN)::  r4    ! Density of state 4.
  real(R_P), intent(IN)::  u4    ! Velocity of state 4.
  real(R_P), intent(IN)::  g4    ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  cp4   ! Specific heats cp of state 4.
  real(R_P), intent(IN)::  cv4   ! Specific heats cv of state 4.
  real(R_P), intent(OUT):: S     ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1    ! Fastest left signal velocities.
  real(R_P), intent(OUT):: S4    ! Fastest right signal velocities.
  real(R_P)::              a1,a4 ! Speed of sound of state 1 and 4.
  real(R_P)::              p23   ! Pressure of intermediate states 2 and 3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  call Roe_up23(p1 = p1, r1 = r1, u1 = u1, g1 = g1, cp1 = cp1, cv1 = cv1, &
                p4 = p4, r4 = r4, u4 = u4, g4 = g4, cp4 = cp4, cv4 = cv4, &
                u23 = S, p23 = p23)
  ! computing waves speed 1 and 4
  call WavesSpeed14up(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23, S1 = S1, S4 = S4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Roe_WavesSpeed14

  elemental subroutine Roe_ua(p1,r1,u1,g1,cp1,cv1,p4,r4,u4,g4,cp4,cv4,u_roe,a_roe)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing speed and speed of sound of Roe average state.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1     ! Pressure of state 1.
  real(R_P), intent(IN)::  r1     ! Density of state 1.
  real(R_P), intent(IN)::  u1     ! Velocity of state 1.
  real(R_P), intent(IN)::  g1     ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  cp1    ! Specific heats cp of state 1.
  real(R_P), intent(IN)::  cv1    ! Specific heats cv of state 1.
  real(R_P), intent(IN)::  p4     ! Pressure of state 4.
  real(R_P), intent(IN)::  r4     ! Density of state 4.
  real(R_P), intent(IN)::  u4     ! Velocity of state 4.
  real(R_P), intent(IN)::  g4     ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  cp4    ! Specific heats cp of state 4.
  real(R_P), intent(IN)::  cv4    ! Specific heats cv of state 4.
  real(R_P), intent(OUT):: u_roe  ! Velocity of Roe average state.
  real(R_P), intent(OUT):: a_roe  ! Speed of sound of Roe average state.
  real(R_P)::              x,omx  ! x = sqrt(r1)/(sqrt(r1)+sqrt(r4)),  omx = 1-x = sqrt(r4)/(sqrt(r1)+sqrt(r4)).
  real(R_P)::              cp_roe ! Specific heats cp of Roe average state
  real(R_P)::              cv_roe ! Specific heats cv of Roe average state
  real(R_P)::              H_roe  ! Total specific entalpy of Roe average state
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x = sqrt(r1)/(sqrt(r1)+sqrt(r4)) ; omx = 1._R_P-x

  u_roe  =                     u1*x +                     u4*omx
  H_roe  = H(p=p1,r=r1,u=u1,g=g1)*x + H(p=p4,r=r4,u=u4,g=g4)*omx
  cp_roe =                    cp1*x +                    cp4*omx
  cv_roe =                    cv1*x +                    cv4*omx
  a_roe  = sqrt((cp_roe/cv_roe - 1._R_P)*(H_roe - 0.5_R_P*u_roe*u_roe))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Roe_ua

  elemental subroutine Roe_State(p1,r1,u1,g1,cp1,cv1,p4,r4,u4,g4,cp4,cv4,r_roe,u_roe,a_roe,H_roe)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the intermediate state Uroe from the known states U1,U4 using the Roe linearization. This subroutine
  ! is optimized for using in Roe Riemann solver. There is another subroutine (Roe_Mean) that computes the same variables that is
  ! optimized for using in WENO reconstruction.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1     ! Pressure of state 1.
  real(R_P), intent(IN)::  r1     ! Density of state 1.
  real(R_P), intent(IN)::  u1     ! Velocity of state 1.
  real(R_P), intent(IN)::  g1     ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  cp1    ! Specific heats cp of state 1.
  real(R_P), intent(IN)::  cv1    ! Specific heats cv of state 1.
  real(R_P), intent(IN)::  p4     ! Pressure of state 4.
  real(R_P), intent(IN)::  r4     ! Density of state 4.
  real(R_P), intent(IN)::  u4     ! Velocity of state 4.
  real(R_P), intent(IN)::  g4     ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  cp4    ! Specific heats cp of state 4.
  real(R_P), intent(IN)::  cv4    ! Specific heats cv of state 4.
  real(R_P), intent(OUT):: r_roe  ! Density of Roe average state.
  real(R_P), intent(OUT):: u_roe  ! Velocity of Roe average state.
  real(R_P), intent(OUT):: a_roe  ! Speed of sound of Roe average state.
  real(R_P), intent(OUT):: H_roe  ! Total specific entalpy of Roe average state.
  real(R_P)::              x,omx  ! x = sqrt(r1)/(sqrt(r1)+sqrt(r4)),  omx = 1-x = sqrt(r4)/(sqrt(r1)+sqrt(r4)).
  real(R_P)::              cp_roe ! Specific heats cp of Roe average state (dummy variable).
  real(R_P)::              cv_roe ! Specific heats cv of Roe average state (dummy variable).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x = sqrt(r1)/(sqrt(r1)+sqrt(r4)) ; omx = 1._R_P - x

  r_roe  = sqrt(r1*r4)
  u_roe  =                     u1*x +                     u4*omx
  H_roe  = H(p=p1,r=r1,u=u1,g=g1)*x + H(p=p4,r=r4,u=u4,g=g4)*omx
  cp_roe =                    cp1*x +                    cp4*omx
  cv_roe =                    cv1*x +                    cv4*omx
  a_roe  = sqrt((cp_roe/cv_roe - 1._R_P)*(H_roe - 0.5_R_P*u_roe*u_roe))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Roe_State

  elemental subroutine BCLC_WavesSpeed14(p1,r1,u1,a1,g1,cp1,cv1,p4,r4,u4,a4,g4,cp4,cv4,S,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the fastest (1,4) waves speed using BCLC (Batten, Clarke, Lambert, Causon) approximation. It is based
  ! on the Roe averages (see "On the Choice of Wavespeeds for the HLLC Riemann Solver", J. Sci. Comput., 1997).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    ! Pressure of state 1.
  real(R_P), intent(IN)::  r1    ! Density of state 1.
  real(R_P), intent(IN)::  u1    ! Velocity of state 1.
  real(R_P), intent(IN)::  a1    ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1    ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  cp1   ! Specific heats cp of state 1.
  real(R_P), intent(IN)::  cv1   ! Specific heats cv of state 1.
  real(R_P), intent(IN)::  p4    ! Pressure of state 4.
  real(R_P), intent(IN)::  r4    ! Density of state 4.
  real(R_P), intent(IN)::  u4    ! Velocity of state 4.
  real(R_P), intent(IN)::  a4    ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4    ! Specific heats ratio of state 4.
  real(R_P), intent(IN)::  cp4   ! Specific heats cp of state 4.
  real(R_P), intent(IN)::  cv4   ! Specific heats cv of state 4.
  real(R_P), intent(OUT):: S     ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1    ! Fastest left signal velocities.
  real(R_P), intent(OUT):: S4    ! Fastest right signal velocities.
  real(R_P)::              a_roe ! Roe average of speed of sound.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing the Roe average of speed and speed of sound
  call Roe_ua(p1 = p1, r1 = r1, u1 = u1, g1 = g1, cp1 = cp1, cv1 = cv1, &
              p4 = p4, r4 = r4, u4 = u4, g4 = g4, cp4 = cp4, cv4 = cv4, &
              u_roe = S, a_roe = a_roe)
  ! evaluating waves speed with BCLC algorithm
  S1 = min((u1 - a1),(S - a_roe))
  S4 = max((u4 + a4),(S + a_roe))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine BCLC_WavesSpeed14

  elemental subroutine Z_WavesSpeed14(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,S,S1,S4)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for evaluating the fastest (1,4) waves speed using Z algorithm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1   ! Pressure of state 1.
  real(R_P), intent(IN)::  r1   ! Density of state 1.
  real(R_P), intent(IN)::  u1   ! Velocity of state 1.
  real(R_P), intent(IN)::  a1   ! Speed of sound of state 1.
  real(R_P), intent(IN)::  g1   ! Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4   ! Pressure of state 4.
  real(R_P), intent(IN)::  r4   ! Density of state 4.
  real(R_P), intent(IN)::  u4   ! Velocity of state 4.
  real(R_P), intent(IN)::  a4   ! Speed of sound of state 4.
  real(R_P), intent(IN)::  g4   ! Specific heats ratio of state 4.
  real(R_P), intent(OUT):: S    ! Velocity of intermediate states.
  real(R_P), intent(OUT):: S1   ! Fastest left signal velocities.
  real(R_P), intent(OUT):: S4   ! Fastest right signal velocities.
  real(R_P)::              p23  ! Pressure of intermediate states 2 and 3.
  real(R_P)::              pmin ! min(p1,p4).
  real(R_P)::              pmax ! max(p1,p4).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pmin = min(p1,p4)
  pmax = max(p1,p4)
  ! evaluating intermediate states pressure and velocity
  call PVL_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = S, p23 = p23)
  if (pmax/pmin>2._R_P.OR.pmin>p23.OR.p23>pmax) then
    if (p23<pmin) then
      call TR_up23(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23)
    else
      call TS_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23)
    endif
  endif
  ! computing waves speed 1 and 4
  call WavesSpeed14up(p1 = p1, u1 = u1, a1 = a1, g1 = g1, p4 = p4, u4 = u4, a4 = a4, g4 = g4, u23 = S, p23 = p23, S1 = S1, S4 = S4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Z_WavesSpeed14

  ! Riemann Solvers
  !> Approximate Riemann solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm.
  elemental subroutine Riem_Solver_LaxFriedrichs(p1,r1,u1,g1,p4,r4,u4,g4,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1             !< Pressure of state 1.
  real(R_P), intent(IN)::  r1             !< Density of state 1.
  real(R_P), intent(IN)::  u1             !< Velocity of state 1.
  real(R_P), intent(IN)::  g1             !< Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4             !< Pressure of state 4.
  real(R_P), intent(IN)::  r4             !< Density of state 4.
  real(R_P), intent(IN)::  u4             !< Velocity of state 4.
  real(R_P), intent(IN)::  g4             !< Specific heats ratio of state 4.
  real(R_P), intent(OUT):: F_r            !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u            !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E            !< Flux of energy conservation.
  real(R_P)::              a1,a4          !< Speed of sound of state 1 and 4.
  real(R_P)::              S1,S4          !< Maximum wave speed of state 1 and 4.
  real(R_P)::              S              !< Contact discontinuity wave speed.
  real(R_P)::              lmax           !< Maximum wave speed estimation.
  real(R_P)::              F_r1,F_u1,F_E1 !< Fluxes of state 1.
  real(R_P)::              F_r4,F_u4,F_E4 !< Fluxes of state 4.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound of state 1 and 4
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  ! evaluating the waves speed S1, S and S4
#ifdef RSLFz
  call Z_WavesSpeed14(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                      p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                      S = S, S1 = S1, S4 = S4)
#else
  call PVL_WavesSpeed14(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                        p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                        S = S, S1 = S1, S4 = S4)
#endif
  ! evalutaing the maximum waves speed
  lmax = max(abs(S1),abs(S4))
  ! computing the fluxes of state 1 and 4
  call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r1, F_u = F_u1, F_E = F_E1)
  call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r4, F_u = F_u4, F_E = F_E4)
  ! computing fluxes
  F_r = 0.5_R_P*(F_r1 + F_r4 - lmax*(r4                        - r1                       ))
  F_u = 0.5_R_P*(F_u1 + F_u4 - lmax*(r4*u4                     - r1*u1                    ))
  F_E = 0.5_R_P*(F_E1 + F_E4 - lmax*(r4*E(p=p4,r=r4,u=u4,g=g4) - r1*E(p=p1,r=r1,u=u1,g=g1)))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_LaxFriedrichs

  !> Approximate Riemann solver based on Primitive Variables Linearization algorithm.
  subroutine Riem_Solver_PVL(p1,r1,u1,g1,p4,r4,u4,g4,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    !< Pressure of state 1.
  real(R_P), intent(IN)::  r1    !< Density of state 1.
  real(R_P), intent(IN)::  u1    !< Velocity of state 1.
  real(R_P), intent(IN)::  g1    !< Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4    !< Pressure of state 4.
  real(R_P), intent(IN)::  r4    !< Density of state 4.
  real(R_P), intent(IN)::  u4    !< Velocity of state 4.
  real(R_P), intent(IN)::  g4    !< Specific heats ratio of state 4.
  real(R_P), intent(OUT):: F_r   !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u   !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E   !< Flux of energy conservation.
  real(R_P)::              a1,a4 !< Speed of sound of state 1 and 4.
  real(R_P)::              S1,S2 !< Left signal velocities.
  real(R_P)::              S4,S3 !< Right signal velocities.
  real(R_P)::              S     !< Velocity of intermediate states.
  real(R_P)::              p23   !< Pressure of intermediate states.
  real(R_P)::              r2    !< Density of state 2.
  real(R_P)::              r3    !< Density of state 3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound of state 1 and 4
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  ! evaluating intermediate states
  call PVL_States(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                  p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                  p23 = p23, u23 = S, r2 = r2, r3 = r3,        &
                  S1 = S1, S2 = S2, S3 = S3, S4 = S4)
  ! sampling interface solution and computing fluxes
  call sampling_fluxes(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, delta1 = 0.5_R_P*(g1-1.0_R_P), eta1 = 2._R_P*g1/(g1-1.0_R_P), &
                       p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, delta4 = 0.5_R_P*(g4-1.0_R_P), eta4 = 2._R_P*g4/(g4-1.0_R_P), &
                       p23 = p23, r2 = r2, r3 = r3, S = S, S1 = S1, S2 = S2, S3 = S3, S4 = S4,                                    &
                       F_r = F_r, F_u = F_u, F_E = F_E)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_PVL

  !> Approximate Riemann solver based on Two Rarefactions algorithm.
  subroutine Riem_Solver_TR(p1,r1,u1,g1,p4,r4,u4,g4,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    !< Pressure of state 1.
  real(R_P), intent(IN)::  r1    !< Density of state 1.
  real(R_P), intent(IN)::  u1    !< Velocity of state 1.
  real(R_P), intent(IN)::  g1    !< Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4    !< Pressure of state 4.
  real(R_P), intent(IN)::  r4    !< Density of state 4.
  real(R_P), intent(IN)::  u4    !< Velocity of state 4.
  real(R_P), intent(IN)::  g4    !< Specific heats ratio of state 4.
  real(R_P), intent(OUT):: F_r   !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u   !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E   !< Flux of energy conservation.
  real(R_P)::              a1,a4 !< Speed of sound of state 1 and 4.
  real(R_P)::              S1,S2 !< Left signal velocities.
  real(R_P)::              S4,S3 !< Right signal velocities.
  real(R_P)::              S     !< Velocity of intermediate states.
  real(R_P)::              p23   !< Pressure of intermediate states.
  real(R_P)::              r2    !< Density of state 2.
  real(R_P)::              r3    !< Density of state 3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound of state 1 and 4
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  ! evaluating intermediate states
  call TR_States(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                 p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                 p23 = p23, u23 = S, r2 = r2, r3 = r3,        &
                 S1 = S1, S2 = S2, S3 = S3, S4 = S4)
  ! sampling interface solution and computing fluxes
  call sampling_fluxes(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, delta1 = 0.5_R_P*(g1-1.0_R_P), eta1 = 2._R_P*g1/(g1-1.0_R_P), &
                       p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, delta4 = 0.5_R_P*(g4-1.0_R_P), eta4 = 2._R_P*g4/(g4-1.0_R_P), &
                       p23 = p23, r2 = r2, r3 = r3, S = S, S1 = S1, S2 = S2, S3 = S3, S4 = S4,                                    &
                       F_r = F_r, F_u = F_u, F_E = F_E)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_TR

  !> Approximate Riemann solver based on Two Shocks algorithm.
  elemental subroutine Riem_Solver_TS(p1,r1,u1,g1,p4,r4,u4,g4,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    !< Pressure of state 1.
  real(R_P), intent(IN)::  r1    !< Density of state 1.
  real(R_P), intent(IN)::  u1    !< Velocity of state 1.
  real(R_P), intent(IN)::  g1    !< Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4    !< Pressure of state 4.
  real(R_P), intent(IN)::  r4    !< Density of state 4.
  real(R_P), intent(IN)::  u4    !< Velocity of state 4.
  real(R_P), intent(IN)::  g4    !< Specific heats ratio of state 4.
  real(R_P), intent(OUT):: F_r   !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u   !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E   !< Flux of energy conservation.
  real(R_P)::              a1,a4 !< Speed of sound of state 1 and 4.
  real(R_P)::              S1,S2 !< Left signal velocities.
  real(R_P)::              S4,S3 !< Right signal velocities.
  real(R_P)::              S     !< Velocity of intermediate states.
  real(R_P)::              p23   !< Pressure of intermediate states.
  real(R_P)::              r2    !< Density of state 2.
  real(R_P)::              r3    !< Density of state 3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound of state 1 and 4
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  ! evaluating intermediate states
  call TS_States(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                 p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                 p23 = p23, u23 = S, r2 = r2, r3 = r3,        &
                 S1 = S1, S2 = S2, S3 = S3, S4 = S4)
  ! computing fluxes
  select case(minloc([-S1,S1*S,S*S4,S4],dim=1))
  case(1)
    call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(2)
    call fluxes(p = p23, r = r2, u = S, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(3)
    call fluxes(p = p23, r = r3, u = S, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  case(4)
    call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_TS

  !> Approximate Riemann solver based on Adaptive (non iterative) PVL-TR-TS algorithm.
  subroutine Riem_Solver_APRS(p1,r1,u1,g1,p4,r4,u4,g4,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1          !< Pressure of state 1.
  real(R_P), intent(IN)::  r1          !< Density of state 1.
  real(R_P), intent(IN)::  u1          !< Velocity of state 1.
  real(R_P), intent(IN)::  g1          !< Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4          !< Pressure of state 4.
  real(R_P), intent(IN)::  r4          !< Density of state 4.
  real(R_P), intent(IN)::  u4          !< Velocity of state 4.
  real(R_P), intent(IN)::  g4          !< Specific heats ratio of state 4.
  real(R_P), intent(OUT):: F_r         !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u         !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E         !< Flux of energy conservation.
  real(R_P)::              a1,a4       !< Speed of sound of state 1 and 4.
  real(R_P)::              gm1_1,gm1_4 !< g-1 of state 1 and 4.
  real(R_P)::              pmin        !< min(p1,p4).
  real(R_P)::              pmax        !< max(p1,p4).
  real(R_P)::              S1,S2       !< Left signal velocities.
  real(R_P)::              S4,S3       !< Right signal velocities.
  real(R_P)::              S           !< Velocity of intermediate states.
  real(R_P)::              p23         !< Pressure of intermediate states.
  real(R_P)::              r2          !< Density and speed of sound of state 2.
  real(R_P)::              r3          !< Density and speed of sound of state 3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound of state 1 and 4
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  ! evaluating intermediate states pressure and velocity
  pmin = min(p1,p4)
  pmax = max(p1,p4)
  call PVL_up23(p1 = p1, r1 = r1, u1 = u1, a1 = a1, p4 = p4, r4 = r4, u4 = u4, a4 = a4, u23 = S, p23 = p23)
  if (pmax/pmin>2._R_P.OR.pmin>p23.OR.p23>pmax) then
    if (p23<pmin) then
      call TR_States(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                     p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                     p23 = p23, u23 = S, r2 = r2, r3 = r3,        &
                     S1 = S1, S2 = S2, S3 = S3, S4 = S4)
    else
      call TS_States(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                     p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                     p23 = p23, u23 = S, r2 = r2, r3 = r3,        &
                     S1 = S1, S2 = S2, S3 = S3, S4 = S4)
    endif
  else
    call PVL_States(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                    p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                    p23 = p23, u23 = S, r2 = r2, r3 = r3,        &
                    S1 = S1, S2 = S2, S3 = S3, S4 = S4)
  endif
  ! sampling interface solution and computing fluxes
  gm1_1 = g1 - 1.0_R_P ; gm1_4 = g4 - 1.0_R_P
  call sampling_fluxes(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, delta1 = 0.5_R_P*gm1_1, eta1 = 2._R_P*g1/gm1_1, &
                       p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, delta4 = 0.5_R_P*gm1_4, eta4 = 2._R_P*g4/gm1_4, &
                       p23 = p23, r2 = r2, r3 = r3, S = S, S1 = S1, S2 = S2, S3 = S3, S4 = S4,                      &
                       F_r = F_r, F_u = F_u, F_E = F_E)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_APRS

  !> Approximate Riemann solver based on Adaptive (non iterative) LF-TR algorithm.
  subroutine Riem_Solver_ALFR(p1,r1,u1,g1,p4,r4,u4,g4,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1                !< Pressure of state 1.
  real(R_P), intent(IN)::  r1                !< Density of state 1.
  real(R_P), intent(IN)::  u1                !< Velocity of state 1.
  real(R_P), intent(IN)::  g1                !< Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4                !< Pressure of state 4.
  real(R_P), intent(IN)::  r4                !< Density of state 4.
  real(R_P), intent(IN)::  u4                !< Velocity of state 4.
  real(R_P), intent(IN)::  g4                !< Specific heats ratio of state 4.
  real(R_P), intent(OUT):: F_r               !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u               !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E               !< Flux of energy conservation.
  real(R_P)::              a1,a2,a3,a4       !< Speed of sound of state 1, 2, 3 and 4.
  real(R_P)::              gm1_1,delta1,eta1 !< g-1, (g-1)/2, 2*g/(g-1)
  real(R_P)::              gm1_4,delta4,eta4 !< g-1, (g-1)/2, 2*g/(g-1)
  real(R_P)::              p2,r2             !< Primitive variables of state 2.
  real(R_P)::              p3,r3             !< Primitive variables of state 3.
  real(R_P)::              S1,S2             !< Left signal velocities.
  real(R_P)::              S4,S3             !< Right signal velocities.
  logical::                RCR               !< Flag for identifying Two Rarefactions condition.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound and gas constants of state 1 and 4
  a1     = a(p=p1,r=r1,g=g1)
  a4     = a(p=p4,r=r4,g=g4)
  gm1_1  = g1 - 1.0_R_P
  gm1_4  = g4 - 1.0_R_P
  delta1 = 0.5_R_P*gm1_1
  delta4 = 0.5_R_P*gm1_4
  eta1   = g1/delta1
  eta4   = g4/delta4
  ! checking RCR condition
  RCR = .false.
  if (p4>=p1) then
    if (u4>=u1) then
      ! NCR: computing the right Rarefaction
      call Rarefaction(sgn   = 1._R_P, &
                       g     = g4,     &
                       delta = delta4, &
                       eta   = eta4,   &
                       u0    = u4,     &
                       p0    = p4,     &
                       a0    = a4,     &
                       ux    = u1,     &
                       rx    = r3,     &
                       px    = p3,     &
                       ax    = a3,     &
                       s0    = S4,     &
                       sx    = S3)
      if (p3<=p1) RCR = .true.
    endif
  else
    if (u4>=u1) then
      ! RCN: computing the left Rarefaction
      call Rarefaction(sgn   = -1._R_P, &
                       g     = g1,      &
                       delta = delta1,  &
                       eta   = eta1,    &
                       u0    = u1,      &
                       p0    = p1,      &
                       a0    = a1,      &
                       ux    = u4,      &
                       rx    = r2,      &
                       px    = p2,      &
                       ax    = a2,      &
                       s0    = S1,      &
                       sx    = S2)
      if (p2<=p4) RCR = .true.
    endif
  endif
  if (RCR) then
    ! strong rarefactions may cause LF failure avoided using TR approximation
    call Riem_Solver_TR(p1=p1,r1=r1,u1=u1,g1=g1,p4=p4,r4=r4,u4=u4,g4=g4,F_r=F_r,F_u=F_u,F_E=F_e)
  else
    ! using LF approximation
    call Riem_Solver_LaxFriedrichs(p1=p1,r1=r1,u1=u1,g1=g1,p4=p4,r4=r4,u4=u4,g4=g4,F_r=F_r,F_u=F_u,F_E=F_e)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_ALFR

  !> Approximate Riemann solver based on HLLC algorithm.
  !> @note Five variants are provided: HLLCb, using BCLC estimation of waves speed, HLLCc using CVL algorithm, HLLCp using PVL
  !> algorithm, HLLCt using TR algorithm and HLLCz using Z one.
  elemental subroutine Riem_Solver_HLLC(p1,r1,u1,g1, &
#if defined RSHLLCb
                                        cp1,cv1,     &
#endif
                                        p4,r4,u4,g4, &
#if defined RSHLLCb
                                        cp4,cv4,     &
#endif
                                        F_r,F_u,F_E)
  !--------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1    !< Pressure of state 1.
  real(R_P), intent(IN)::  r1    !< Density of state 1.
  real(R_P), intent(IN)::  u1    !< Velocity of state 1.
  real(R_P), intent(IN)::  g1    !< Specific heats ratio of state 1.
#if defined RSHLLCb
  real(R_P), intent(IN)::  cp1   !< Specific heats cp of state 1.
  real(R_P), intent(IN)::  cv1   !< Specific heats cv of state 1.
#endif
  real(R_P), intent(IN)::  p4    !< Pressure of state 4.
  real(R_P), intent(IN)::  r4    !< Density of state 4.
  real(R_P), intent(IN)::  u4    !< Velocity of state 4.
  real(R_P), intent(IN)::  g4    !< Specific heats ratio of state 4
#if defined RSHLLCb
  real(R_P), intent(IN)::  cp4   !< Specific heats cp of state 4.
  real(R_P), intent(IN)::  cv4   !< Specific heats cv of state 4.
#endif
  real(R_P), intent(OUT):: F_r   !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u   !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E   !< Flux of energy conservation.
  real(R_P)::              a1,a4 !< Speed of sound of state 1 and 4.
  real(R_P)::              S1,S4 !< Left and right wave speeds.
  real(R_P)::              S     !< Contact discontinuity wave speed.
  real(R_P)::              U1S   !< Mass conservation.
  real(R_P)::              U2S   !< Momentum conservation.
  real(R_P)::              U3S   !< Energy conservation.
  real(R_P)::              E1    !< Entalpy and internal energy of left state.
  real(R_P)::              E4    !< Entalpy and internal energy of right state.
#ifdef VACUUM
  real(R_P)::              d1,d4 !< 2/(g-1) of state 1 and 4.
  real(R_P)::              S2,S3 !< Left and right speeds of vacuum fronts.
#endif
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound of state 1 and 4
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
#ifdef VACUUM
  d1 = 2._R_P/(g1-1._R_P)
  d4 = 2._R_P/(g4-1._R_P)
  ! check the vacuum condition
  if ((d1*a1+d4*a4)<=(u4-u1)) then
    ! there is vacuum!
    ! computing the speed of rarefactions bounding the vacuum
    S1 = u1 - a1
    S2 = u1 + d1*a1
    S3 = u4 - d4*a4
    S4 = u4 + a4
    ! sampling interface solution and computing fluxes
    call sampling_fluxes_vacuum(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, delta1 = 1._R_P/d1, eta1 = g1*d1, &
                                p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, delta4 = 1._R_P/d4, eta4 = g4*d4, &
                                S1 = S1, S2 = S2, S3 = S3, S4 = S4,                                        &
                                F_r = F_r, F_u = F_u, F_E = F_E)
  endif
#endif
  ! evaluating the waves speed S1, S and S4
#if defined RSHLLCb
  call BCLC_WavesSpeed14(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, cp1 = cp1, cv1 = cv1, &
                         p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, cp4 = cp4, cv4 = cv4, &
                         S = S, S1 = S1, S4 = S4)
#elif defined RSHLLCc
  call CVL_WavesSpeed14(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                        p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                        S = S, S1 = S1, S4 = S4)
#elif defined RSHLLCt
  call TR_WavesSpeed14(p1 = p1, u1 = u1, a1 = a1, g1 = g1, &
                       p4 = p4, u4 = u4, a4 = a4, g4 = g4, &
                       S = S, S1 = S1, S4 = S4)
#elif defined RSHLLCz
  call Z_WavesSpeed14(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                      p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                      S = S, S1 = S1, S4 = S4)
#else
  call PVL_WavesSpeed14(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                        p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                        S = S, S1 = S1, S4 = S4)
#endif
  S = (r4*u4*(S4-u4) - r1*u1*(S1-u1) + p1 - p4)/(r4*(S4-u4) - r1*(S1-u1))
  ! computing fluxes
  select case(minloc([-S1,S1*S,S*S4,S4],dim=1))
  case(1)
    call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(2)
    call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
    E1  = E(p=p1,r=r1,u=u1,g=g1)
    U1S = r1*(S1-u1)/(S1-S)
    U2S = U1S*S
    U3S = U1S*(E1+(S-u1)*(S+p1/(r1*(S1-u1))))

    F_r = F_r + S1*(U1S - r1)
    F_u = F_u + S1*(U2S - r1*u1)
    F_E = F_E + S1*(U3S - r1*E1)
  case(3)
    call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
    E4  = E(p=p4,r=r4,u=u4,g=g4)
    U1S = r4*(S4-u4)/(S4-S)
    U2S = U1S*S
    U3S = U1S*(E4+(S-u4)*(S+p4/(r4*(S4-u4))))

    F_r = F_r + S4*(U1S - r4)
    F_u = F_u + S4*(U2S - r4*u4)
    F_E = F_E + S4*(U3S - r4*E4)
  case(4)
    call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_HLLC

  !> Approximate Riemann solver based on Roe linearization.
  !> @note The Harten-Hyman entropy fix is used to avoid shock-expansions (entropy violation).
  elemental subroutine Riem_Solver_Roe(p1,r1,u1,g1,cp1,cv1,p4,r4,u4,g4,cp4,cv4,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1              !< Pressure of state 1.
  real(R_P), intent(IN)::  r1              !< Density of state 1.
  real(R_P), intent(IN)::  u1              !< Velocity of state 1.
  real(R_P), intent(IN)::  g1              !< Specific heats ratio of state 1.
  real(R_P), intent(IN)::  cp1             !< Specific heats cp of state 1.
  real(R_P), intent(IN)::  cv1             !< Specific heats cv of state 1.
  real(R_P), intent(IN)::  p4              !< pressure of state 4.
  real(R_P), intent(IN)::  r4              !< Density of state 4.
  real(R_P), intent(IN)::  u4              !< Velocity of state 4.
  real(R_P), intent(IN)::  g4              !< Specific heats ratio of state 4.
  real(R_P), intent(IN)::  cp4             !< Specific heats cp of state 4.
  real(R_P), intent(IN)::  cv4             !< Specific heats cv of state 4.
  real(R_P), intent(OUT):: F_r             !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u             !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E             !< Flux of energy conservation.
  real(R_P)::              a1,a4           !< Speed of sound of state 1 and 4.
  real(R_P)::              S1,S2           !< Left signal velocities.
  real(R_P)::              S4,S3           !< Right signal velocities.
  real(R_P)::              S               !< Velocity of intermediate states.
  real(R_P)::              F_r1,F_u1,F_E1  !< Fluxes of state 1.
  real(R_P)::              F_r4,F_u4,F_E4  !< Fluxes of state 4.
  real(R_P)::              Dr              !< Density difference  Dr = r4-r1.
  real(R_P)::              Du              !< Velocity difference Du = u4-u1.
  real(R_P)::              Dp              !< Pressure difference Dp = p4-p1.
  real(R_P)::              r_r,u_r,a_r,H_r !< Roe's state.
  real(R_P)::              aa1,aa2,aa3     !< Wawes amplitudes Roe's estimation.
  real(R_P)::              ll1,ll2,ll3     !< Wawes speeds Roe's estimation.
  real(R_P)::              ls1,    ls3     !< Wawes speeds Roe's estimation with entropy fix of Harten-Hyman.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing speed of sound of state 1 and 4
  a1 = a(p=p1,r=r1,g=g1)
  a4 = a(p=p4,r=r4,g=g4)
  ! evaluating the intermediates states for providing the sampling information necessary for the entropy fix
  call PVL_WavesSpeed1234(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, &
                          p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, &
                          S = S, S1 = S1, S2 = S2, S3 = S3, S4 = S4)
  ! fluxes computation
  select case(minloc([-S1,S1*S2,S2*S,S*S3,S3*S4,S4],dim=1))
  case(1) ! left supersonic
    call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r, F_u = F_u, F_E = F_E)
  case(2) ! left transonic
    call Roe_State(r1 = r1, u1 = u1, p1 = p1, g1 = g1, cp1 = cp1, cv1 = cv1, &
                   r4 = r4, u4 = u4, p4 = p4, g4 = g4, cp4 = cp4, cv4 = cv4, &
                   r_roe = r_r, u_roe = u_r, a_roe = a_r, H_roe = H_r)
    Du  = u4 - u1
    Dp  = p4 - p1
    aa1 = 0.5_R_P*(Dp - r_r*a_r*Du)/(a_r*a_r)
    ll1 = u_r - a_r
    ls1 = S1*(S2-ll1)/(S2-S1)
    call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r1, F_u = F_u1, F_E = F_E1)

    F_r = F_r1 + aa1*ls1
    F_u = F_u1 + aa1*ls1*ll1
    F_E = F_E1 + aa1*ls1*(H_r-u_r*a_r)
  case(3,4) ! subsonic
    call Roe_State(r1 = r1, u1 = u1, p1 = p1, g1 = g1, cp1 = cp1, cv1 = cv1, &
                   r4 = r4, u4 = u4, p4 = p4, g4 = g4, cp4 = cp4, cv4 = cv4, &
                   r_roe = r_r, u_roe = u_r, a_roe = a_r, H_roe = H_r)
    Dr  = r4 - r1
    Du  = u4 - u1
    Dp  = p4 - p1
    aa1 = 0.5_R_P*(Dp - r_r*a_r*Du)/(a_r*a_r)
    aa2 = Dr - Dp/(a_r*a_r)
    aa3 = 0.5_R_P*(Dp + r_r*a_r*Du)/(a_r*a_r)
    ll1 = u_r - a_r
    ll2 = u_r
    ll3 = u_r + a_r
    call fluxes(p = p1, r = r1, u = u1, g = g1, F_r = F_r1, F_u = F_u1, F_E = F_E1)
    call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r4, F_u = F_u4, F_E = F_E4)

    F_r = 0.5_R_P*(F_r1 + F_r4 - (aa1*abs(ll1)              +aa2*abs(ll2)                +aa3*abs(ll3)              ))
    F_u = 0.5_R_P*(F_u1 + F_u4 - (aa1*abs(ll1)*ll1          +aa2*abs(ll2)*ll2            +aa3*abs(ll3)*ll3          ))
    F_E = 0.5_R_P*(F_E1 + F_E4 - (aa1*abs(ll1)*(H_r-u_r*a_r)+aa2*abs(ll2)*ll2*ll2*0.5_R_P+aa3*abs(ll3)*(H_r+u_r*a_r)))
  case(5) ! right transonic
    call Roe_State(r1 = r1, u1 = u1, p1 = p1, g1 = g1, cp1 = cp1, cv1 = cv1, &
                   r4 = r4, u4 = u4, p4 = p4, g4 = g4, cp4 = cp4, cv4 = cv4, &
                   r_roe = r_r, u_roe = u_r, a_roe = a_r, H_roe = H_r)
    Du  = u4 - u1
    Dp  = p4 - p1
    aa3 = 0.5_R_P*(Dp + r_r*a_r*Du)/(a_r*a_r)
    ll3 = u_r + a_r
    ls3 = S4*(ll3-S3)/(S4-S3)
    call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r4, F_u = F_u4, F_E = F_E4)

    F_r = F_r4 - aa3*ls3
    F_u = F_u4 - aa3*ls3*ll3
    F_E = F_E4 - aa3*ls3*(H_r+u_r*a_r)
  case(6) ! right supersonic
    call fluxes(p = p4, r = r4, u = u4, g = g4, F_r = F_r, F_u = F_u, F_E = F_E)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_Roe

  !> Exact Riemann solver based on iterative solution of u-function.
  subroutine Riem_Solver_Exact_U(p1,r1,u1,g1,p4,r4,u4,g4,F_r,F_u,F_E)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::  p1                      !< Pressure of state 1.
  real(R_P), intent(IN)::  r1                      !< Density of state 1.
  real(R_P), intent(IN)::  u1                      !< U speed of state 1.
  real(R_P), intent(IN)::  g1                      !< Specific heats ratio of state 1.
  real(R_P), intent(IN)::  p4                      !< Pressure of state 4.
  real(R_P), intent(IN)::  r4                      !< Density of state 4.
  real(R_P), intent(IN)::  u4                      !< U speed of state 4.
  real(R_P), intent(IN)::  g4                      !< Specific heats ratio of state 4.
  real(R_P), intent(OUT):: F_r                     !< Flux of mass conservation.
  real(R_P), intent(OUT):: F_u                     !< Flux of momentum conservation.
  real(R_P), intent(OUT):: F_E                     !< Flux of energy conservation.
  real(R_P)::              a1,a2,a3,a4             !< Speed of sound of state 1, 2, 3 and 4.
  real(R_P)::              gm1_1,gp1_1,delta1,eta1 !< g-1, g+1, (g-1)/2, 2*g/(g-1)
  real(R_P)::              gm1_4,gp1_4,delta4,eta4 !< g-1, g+1, (g-1)/2, 2*g/(g-1)
  real(R_P)::              p2,r2                   !< Primitive variables of state 2.
  real(R_P)::              p3,r3                   !< Primitive variables of state 3.
  real(R_P)::              dp2,dp3                 !< Derivate of pessure (dp/du) of state 2 and 3.
  real(R_P)::              u23                     !< Speed of state 2 and 3.
  real(R_P)::              S1                      !< Downstream front of C1 wave.
  real(R_P)::              S2                      !< Upstream front of C1 wave.
  real(R_P)::              S3                      !< Downstream front of C2 wave.
  real(R_P)::              S4                      !< Upstream front of C1 wave.
  real(R_P)::              dum,alfa,beta           !< Dummies coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  a1     = a(p=p1,r=r1,g=g1)        ! computing speed of sound of state 1
  a4     = a(p=p4,r=r4,g=g4)        ! computing speed of sound of state 4
  gm1_1  = g1 - 1.0_R_P             ! computing g-1 for state 1
  gm1_4  = g4 - 1.0_R_P             ! computing g-1 for state 4
  gp1_1  = g1 + 1.0_R_P             ! computing g+1 for state 1
  gp1_4  = g4 + 1.0_R_P             ! computing g+1 for state 4
  delta1 = 0.5_R_P*gm1_1            ! computing (g-1)/2 for state 1
  delta4 = 0.5_R_P*gm1_4            ! computing (g-1)/2 for state 4
  eta1   = g1/delta1                ! computing 2*g/(g-1) for state 1
  eta4   = g4/delta4                ! computing 2*g/(g-1) for state 4
  if (p1<p4) then
    dum  = 0.5_R_P*gm1_4/g4         ! computing (g-1)/(g*2)
  else
    dum  = 0.5_R_P*gm1_1/g1         ! computing (g-1)/(g*2)
  endif
  alfa   = (p1/p4)**dum             ! computing alfa coefficient
  beta   = alfa*delta1/a1+delta4/a4 ! computing beta coefficient
  ! computing first approximation of intemediate speed u23
  u23    = (alfa-1.0_R_P)/beta + &
           0.5_R_P*(u1+u4)     + &
           0.5_R_P*(u1-u4)*(alfa*delta1/a1-delta4/a4)/beta ! initializing u speed
  ! computing Newton-Rapson iterations
  Newton: do
    ! computing intermediate states using the current approximation of u23
    call InterStates23u(p1 = p1, u1 = u1, a1 = a1, g1 = g1, gm1_1 = gm1_1, gp1_1 = gp1_1, delta1 = delta1, eta1 = eta1, &
                        p4 = p4, u4 = u4, a4 = a4, g4 = g4, gm1_4 = gm1_4, gp1_4 = gp1_4, delta4 = delta4, eta4 = eta4, &
                        u23 = u23,                                                                                      &
                        r2 = r2, p2 = p2, a2 = a2, r3 = r3, p3 = p3, a3 = a3 ,S1 = S1, S2 = S2, S3 = S3, S4 =S4)
    ! computing dp/du for the computed intermediates states
    dp2 = -1._R_P*g1*p2/a2
    dp3 =  1._R_P*g4*p3/a3
    ! evaluating the Newton-Rapson convergence
    if (abs(1.0_R_P - (p2/p3))>=toll) then ! Newton iterative step
      u23  = u23 - ((p2-p3)/(dp2-dp3))     ! updating u value
    else                                   ! Newton iterations have been converged
      exit Newton
    endif
  enddo Newton
  ! sampling interface solution and computing fluxes
  call sampling_fluxes(p1 = p1, r1 = r1, u1 = u1, a1 = a1, g1 = g1, delta1 = delta1, eta1 = eta1, &
                       p4 = p4, r4 = r4, u4 = u4, a4 = a4, g4 = g4, delta4 = delta4, eta4 = eta4, &
                       p23 = p2, r2 = r2, r3 = r3, S = u23, S1 = S1, S2 = S2, S3 = S3, S4 = S4,   &
                       F_r = F_r, F_u = F_u, F_E = F_E)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Riem_Solver_Exact_U
endmodule Lib_Riemann
