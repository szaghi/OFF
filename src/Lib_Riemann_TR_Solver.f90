!> @ingroup Library
!> @{
!> @defgroup Lib_Riemann_TR_SolverLibrary Lib_Riemann_TR_Solver
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Riemann_TR_SolverPublicProcedure Lib_Riemann_TR_Solver
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Riemann_TR_SolverPrivateProcedure Lib_Riemann_TR_Solver
!> @}

!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws my means of
!> TR (Two Rarefactions) approximation.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Riemann_TR_SolverLibrary
module Lib_Riemann_TR_Solver
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                    ! precision of integers and reals
USE Data_Type_Riemann_Primitive1D,  only: Type_Riemann_Primitive1D  ! Definition of Type_Riemann_Primitive1D.
USE Data_Type_Riemann_InterState1D, only: Type_Riemann_InterState1D ! Definition of Type_Riemann_InterState1D.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: TR_u23
public:: TR_up23
public:: TR_compute_waves14_u
public:: TR_compute_waves14_up
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Riemann_TR_SolverPublicProcedure
  !> @{
  !> @brief Procedure for evaluating speed of intermediate states using TR (Two Rarefactions) approximation.
  elemental function TR_u23(state1,state4) result(u23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D), intent(IN)::  state1    !< State 1 (left).
  type(Type_Riemann_Primitive1D), intent(IN)::  state4    !< State 4 (right).
  real(R8P)::                                   u23       !< Velocity of intermediate states.
  real(R8P)::                                   gm,gm1    !< Mean value of specific heats ratios and g-1.
  real(R8P)::                                   dm        !< Mean value of d = (g-1)/2.
  real(R8P)::                                   zm,zm_inv !< Mean value of z = d/g and 1/z.
  real(R8P)::                                   p14       !< p14 = (p1/p4)**zm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  gm     = 0.5_R8P*(state1%g + state4%g)
  gm1    = gm - 1._R8P
  dm     = 0.5_R8P*gm1
  zm     = dm/gm
  zm_inv = 1._R8P/zm
  p14    = (state1%p/state4%p)**zm
  u23    = (p14*state1%u/state1%a + state4%u/state4%a + 2._R8P*(p14 - 1._R8P)/gm1)/(p14/state1%a + 1._R8P/state4%a)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction TR_u23

  !> @brief Procedure for evaluating speed and pressure of intermediate states using TR (Two Rarefactions) approximation.
  elemental subroutine TR_up23(state1,state4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D), intent(IN)::  state1    !< State 1 (left).
  type(Type_Riemann_Primitive1D), intent(IN)::  state4    !< State 4 (right).
  real(R8P),                      intent(OUT):: u23       !< Velocity of intermediate states.
  real(R8P),                      intent(OUT):: p23       !< Pressure of intermediate states.
  real(R8P)::                                   gm,gm1    !< Mean value of specific heats ratios and g-1.
  real(R8P)::                                   dm        !< Mean value of d = (g-1)/2.
  real(R8P)::                                   zm,zm_inv !< Mean value of z = d/g and 1/z.
  real(R8P)::                                   p14       !< p14 = (p1/p4)**zm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  gm     = 0.5_R8P*(state1%g + state4%g)
  gm1    = gm - 1._R8P
  dm     = 0.5_R8P*gm1
  zm     = dm/gm
  zm_inv = 1._R8P/zm
  p14    = (state1%p/state4%p)**zm
  u23    = (p14*state1%u/state1%a + state4%u/state4%a + 2._R8P*(p14 - 1._R8P)/gm1)/(p14/state1%a + 1._R8P/state4%a)
  p23    = 0.5_R8P*(state1%p*(1._R8P+dm*(state1%u-u23)/state1%a)**zm_inv+state4%p*(1._R8P+dm*(u23-state4%u)/state4%a)**zm_inv)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TR_up23

  !> @brief Procedure for evaluating the fastest (1,4) waves speed using TR (Two Rarefactions) approximation.
  !> This variant uses only the approximation of u provided by the TR approximation for computing the waves speed (WSup algorithm).
  elemental subroutine TR_compute_waves14_u(state1,state4,state23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),  intent(IN)::    state4  !< State 4 (right).
  type(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state 2,3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  state23%u = TR_u23(state1=state1,state4=state4)
  ! computing waves speed 1 and 4
  call state23%compute_w14_from_u(state1=state1,state4=state4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TR_compute_waves14_u

  !> @brief Procedure for evaluating the fastest (1,4) waves speed using TR (Two Rarefactions) approximation.
  !> This variant uses both the approximations of u and p provided by the TR approximation for computing the waves speed
  !> (WSup algorithm).
  elemental subroutine TR_compute_waves14_up(state1,state4,state23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),  intent(IN)::    state4  !< State 4 (right).
  type(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state 2,3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  call TR_up23(state1=state1,state4=state4,u23=state23%u,p23=state23%p)
  ! computing waves speed 1 and 4
  call state23%compute_w14_from_up(state1=state1,state4=state4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TR_compute_waves14_up
  !> @}
endmodule Lib_Riemann_TR_Solver
