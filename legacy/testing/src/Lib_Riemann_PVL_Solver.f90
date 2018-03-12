!> @ingroup Library
!> @{
!> @defgroup Lib_Riemann_PVL_SolverLibrary Lib_Riemann_PVL_Solver
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Riemann_PVL_SolverPublicProcedure Lib_Riemann_PVL_Solver
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Riemann_PVL_SolverPrivateProcedure Lib_Riemann_PVL_Solver
!> @}

!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws my means of
!> PVL (Primitive Variables Linearization) approximation.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Riemann_PVL_SolverLibrary
module Lib_Riemann_PVL_Solver
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                    ! precision of integers and reals
USE Data_Type_Riemann_Primitive1D,  only: Type_Riemann_Primitive1D  ! Definition of Type_Riemann_Primitive1D.
USE Data_Type_Riemann_InterState1D, only: Type_Riemann_InterState1D ! Definition of Type_Riemann_InterState1D.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: PVL_u23
public:: PVL_p23
public:: PVL_up23
public:: PVL_compute_waves14_u
public:: PVL_compute_waves14_up
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Riemann_PVL_SolverPublicProcedure
  !> @{
  !> @brief Procedure for evaluating speed of intermediate states using PVL (Primitive Variables Linearization) approximation.
  elemental function PVL_u23(state1,state4) result(u23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D), intent(IN):: state1 !< State 1 (left).
  type(Type_Riemann_Primitive1D), intent(IN):: state4 !< State 4 (right).
  real(R8P)::                                  u23    !< Velocity of intermediate states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  u23 = 0.5_R8P*(state1%u + state4%u) - 2.0_R8P*(state4%p - state1%p)/((state1%r + state4%r)*(state1%a + state4%a))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction PVL_u23

  !> @brief Procedure for evaluating pressure of intermediate states using PVL (Primitive Variables Linearization) approximation.
  elemental function PVL_p23(state1,state4) result(p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D), intent(IN):: state1 !< State 1 (left).
  type(Type_Riemann_Primitive1D), intent(IN):: state4 !< State 4 (right).
  real(R8P)::                                  p23    !< Pressure of intermediate states.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  p23 = 0.5_R8P*((state1%p + state4%p) - 0.25_R8P*(state4%u - state1%u)*(state1%r + state4%r)*(state1%a + state4%a))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction PVL_p23

  !> @brief Procedure for evaluating speed and pressure of intermediate states using PVL (Primitive Variables Linearization)
  !> approximation.
  elemental subroutine PVL_up23(state1,state4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D), intent(IN)::  state1 !< State 1 (left).
  type(Type_Riemann_Primitive1D), intent(IN)::  state4 !< State 4 (right).
  real(R_P),                      intent(OUT):: u23    !< Velocity of intermediate states.
  real(R_P),                      intent(OUT):: p23    !< Pressure of intermediate states.
  real(R_P)::                                   ram    !< Mean value of r*a.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ram = 0.25_R8P*(state1%r + state4%r)*(state1%a + state4%a)
  u23 = 0.5_R8P*((state1%u + state4%u)-(state4%p - state1%p)/ram)
  p23 = 0.5_R8P*((state1%p + state4%p)-(state4%u - state1%u)*ram)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_up23

  !> @brief Procedure for evaluating the fastest (1,4) waves speed using PVL (Primitive Variables Linearization) approximation.
  !> This variant uses only the approximation of u provided by the PVL approximation for computing the waves speed (WSup algorithm).
  pure subroutine PVL_compute_waves14_u(state1,state4,state23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),  intent(IN)::    state4  !< State 4 (right).
  type(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state 2,3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states velocity
  state23%u = PVL_u23(state1=state1,state4=state4)
  ! computing waves speed 1 and 4
  call state23%compute_w14_from_u(state1=state1,state4=state4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_compute_waves14_u

  !> @brief Procedure for evaluating the fastest (1,4) waves speed using PVL (Primitive Variables Linearization) approximation.
  !> This variant uses both the approximations of u and p provided by the PVL approximation for computing the waves speed
  !> (WSup algorithm).
  pure subroutine PVL_compute_waves14_up(state1,state4,state23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),  intent(IN)::    state4  !< State 4 (right).
  type(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state 2,3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  call PVL_up23(state1=state1,state4=state4,u23=state23%u,p23=state23%p)
  ! computing waves speed 1 and 4
  call state23%compute_w14_from_up(state1=state1,state4=state4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine PVL_compute_waves14_up
  !> @}
endmodule Lib_Riemann_PVL_Solver
