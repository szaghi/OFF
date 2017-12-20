!> @ingroup Library
!> @{
!> @defgroup Lib_Riemann_Z_SolverLibrary Lib_Riemann_Z_Solver
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Riemann_Z_SolverPublicProcedure Lib_Riemann_Z_Solver
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Riemann_Z_SolverPrivateProcedure Lib_Riemann_Z_Solver
!> @}

!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws my means of
!> TS (Two Shocks) approximation.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Riemann_Z_SolverLibrary
module Lib_Riemann_Z_Solver
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                    ! precision of integers and reals
USE Data_Type_Riemann_Primitive1D,  only: Type_Riemann_Primitive1D  ! Definition of Type_Riemann_Primitive1D.
USE Data_Type_Riemann_InterState1D, only: Type_Riemann_InterState1D ! Definition of Type_Riemann_InterState1D.
USE Lib_Riemann_PVL_Solver,         only: PVL_up23                  ! Library of PVL Riemann solver.
USE Lib_Riemann_TR_Solver,          only: TR_up23                   ! Library of TR  Riemann solver.
USE Lib_Riemann_Ts_Solver,          only: TS_up23                   ! Library of TS  Riemann solver.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: Z_compute_waves14_up
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Riemann_Z_SolverPublicProcedure
  !> @{
  !> @brief Procedure for evaluating the fastest (1,4) waves speed using Z algorithm.
  pure subroutine Z_compute_waves14_up(state1,state4,state23)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),  intent(IN)::    state4  !< State 4 (right).
  type(Type_Riemann_InterState1D), intent(INOUT):: state23 !< Intermediate state 2,3.
  real(R8P)::                                      pmin    !< min(p1,p4).
  real(R8P)::                                      pmax    !< max(p1,p4).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pmin = min(state1%p,state4%p)
  pmax = max(state1%p,state4%p)
  ! evaluating intermediate states pressure and velocity
  call PVL_up23(state1=state1,state4=state4,u23=state23%u,p23=state23%p)
  if (pmax/pmin>2._R8P.OR.pmin>state23%p.OR.state23%p>pmax) then
    if (state23%p<pmin) then
      call TR_up23(state1=state1,state4=state4,u23=state23%u,p23=state23%p)
    else
      call TS_up23(state1=state1,state4=state4,u23=state23%u,p23=state23%p)
    endif
  endif
  ! computing waves speed 1 and 4
  call state23%compute_w14_from_up(state1=state1,state4=state4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Z_compute_waves14_up

  !> @}
endmodule Lib_Riemann_Z_Solver
