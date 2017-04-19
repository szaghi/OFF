!> @ingroup Library
!> @{
!> @defgroup Lib_Riemann_TS_SolverLibrary Lib_Riemann_TS_Solver
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Riemann_TS_SolverPublicProcedure Lib_Riemann_TS_Solver
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Riemann_TS_SolverPrivateProcedure Lib_Riemann_TS_Solver
!> @}

!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws my means of
!> TS (Two Shocks) approximation.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Riemann_TS_SolverLibrary
module Lib_Riemann_TS_Solver
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                    ! precision of integers and reals
USE Data_Type_Riemann_Primitive1D,  only: Type_Riemann_Primitive1D  ! Definition of Type_Riemann_Primitive1D.
USE Data_Type_Riemann_InterState1D, only: Type_Riemann_InterState1D ! Definition of Type_Riemann_InterState1D.
USE Lib_Riemann_PVL_Solver,         only: PVL_p23                   ! Library of PVL Riemann solver.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: TS_up23
public:: TS_compute_waves14_up
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Riemann_TS_SolverPublicProcedure
  !> @{
  ! Subroutine for evaluating speed and pressure of intermediate states using TS (Two Shocks) approximation.
  elemental subroutine TS_up23(state1,state4,u23,p23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D), intent(IN)::  state1      !< State 1 (left).
  type(Type_Riemann_Primitive1D), intent(IN)::  state4      !< State 4 (right).
  real(R8P),                      intent(OUT):: u23         !< Velocity of intermediate states.
  real(R8P),                      intent(OUT):: p23         !< Pressure of intermediate states.
  real(R8P)::                                   g1p,g4p     !< Dummy variables for computing TS pressure of intermediate states.
  real(R8P)::                                   gp1,gm1_gp1 !< g+1, (g-1)/(g+1).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  p23 = PVL_p23(state1=state1,state4=state4)
  gp1 = state1%g + 1._R8P ; gm1_gp1 = (state1%g - 1._R8P)/gp1 ; g1p = sqrt((2._R8P/(gp1*state1%r))/(p23 + gm1_gp1*state1%p))
  gp1 = state4%g + 1._R8P ; gm1_gp1 = (state4%g - 1._R8P)/gp1 ; g4p = sqrt((2._R8P/(gp1*state4%r))/(p23 + gm1_gp1*state4%p))
  p23 = (g1p*state1%p + g4p*state4%p + state1%u - state4%u)/(g1p + g4p)
  u23 = 0.5_R8P*(state1%u + state4%u + (p23 - state4%p)*g4p - (p23 - state1%p)*g1p)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TS_up23

  !> @brief Procedure for evaluating the fastest (1,4) waves speed using TS (Two Shocks) approximation.
  pure subroutine TS_compute_waves14_up(state1,state4,state23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),  intent(IN)::    state4  !< State 4 (right).
  type(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state 2,3.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating intermediate states pressure and velocity
  call TS_up23(state1=state1,state4=state4,u23=state23%u,p23=state23%p)
  ! computing waves speed 1 and 4
  call state23%compute_w14_from_up(state1=state1,state4=state4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine TS_compute_waves14_up
  !> @}
endmodule Lib_Riemann_TS_Solver
