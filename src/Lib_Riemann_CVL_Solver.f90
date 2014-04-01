!> @ingroup Library
!> @{
!> @defgroup Lib_Riemann_CVL_SolverLibrary Lib_Riemann_CVL_Solver
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Riemann_CVL_SolverPublicProcedure Lib_Riemann_CVL_Solver
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Riemann_CVL_SolverPrivateProcedure Lib_Riemann_CVL_Solver
!> @}

!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws my means of
!> CVL (Charactheristic Variables Linearization) approximation.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Riemann_CVL_SolverLibrary
module Lib_Riemann_CVL_Solver
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                    ! precision of integers and reals
USE Data_Type_Riemann_Primitive1D,  only: Type_Riemann_Primitive1D  ! Definition of Type_Riemann_Primitive1D.
USE Data_Type_Riemann_InterState1D, only: Type_Riemann_InterState1D ! Definition of Type_Riemann_InterState1D.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: CVL_u23
public:: CVL_up23
public:: CVL_compute_waves14_u
public:: CVL_compute_waves14_up
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Riemann_CVL_SolverPublicProcedure
  !> @{
  !> @brief Procedure for evaluating speed of intermediate states using CVL (Charactheristic Variables Linearization) approximation.
  elemental function CVL_u23(state1,state4) result(u23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D), intent(IN):: state1 !< State 1 (left).
  type(Type_Riemann_Primitive1D), intent(IN):: state4 !< State 4 (right).
  real(R8P)::                                  u23    !< Velocity of intermediate states.
  real(R8P)::                                  c1     !< Caracteristic of state 1.
  real(R8P)::                                  c4     !< Caracteristic of state 4.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  c1  = state1%r*state1%a
  c4  = state4%r*state4%a
  u23 = (c1*state1%u + c4*state4%u + (state1%p - state4%p))/(c1 + c4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction CVL_u23

  !> @brief Procedure for evaluating speed and pressure of intermediate states using CVL (Charactheristic Variables Linearization)
  !> approximation.
  elemental subroutine CVL_up23(state1,state4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D), intent(IN)::  state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D), intent(IN)::  state4  !< State 4 (right).
  real(R8P),                      intent(OUT):: u23     !< Velocity of intermediate states.
  real(R8P),                      intent(OUT):: p23     !< Pressure of intermediate states.
  real(R8P)::                                   c1      !< Caracteristic of state 1.
  real(R8P)::                                   c4      !< Caracteristic of state 4.
  real(R8P)::                                   c14_inv !< 1/(c1+c4).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  c1      = state1%r*state1%a
  c4      = state4%r*state4%a
  c14_inv = 1._R8P/(c1 + c4)
  u23     = c14_inv*(c1*state1%u + c4*state4%u +       (state1%p - state4%p))
  p23     = c14_inv*(c4*state1%p + c1*state4%p + c1*c4*(state1%u - state4%u))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine CVL_up23

  !> @brief Procedure for evaluating the fastest (1,4) waves speed using CVL (Charactheristic Variables Linearization)
  !> This variant uses only the approximation of u provided by the CVL approximation for computing the waves speed (WSup algorithm).
  elemental subroutine CVL_compute_waves14_u(state1,state4,state23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),  intent(IN)::    state4  !< State 4 (right).
  type(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state 2,3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  state23%u = CVL_u23(state1=state1,state4=state4)
  ! computing waves speed 1 and 4
  call state23%compute_w14_from_u(state1=state1,state4=state4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine CVL_compute_waves14_u

  !> @brief Procedure for evaluating the fastest (1,4) waves speed using CVL (Charactheristic Variables Linearization)
  !> This variant uses both the approximations of u and p provided by the CVL approximation for computing the waves speed
  !> (WSup algorithm).
  elemental subroutine CVL_compute_waves14_up(state1,state4,state23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),  intent(IN)::    state4  !< State 4 (right).
  type(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state 2,3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  call CVL_up23(state1=state1,state4=state4,u23=state23%u,p23=state23%u)
  ! computing waves speed 1 and 4
  call state23%compute_w14_from_up(state1=state1,state4=state4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine CVL_compute_waves14_up
  !> @}
endmodule Lib_Riemann_CVL_Solver
