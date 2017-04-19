!> @ingroup Library
!> @{
!> @defgroup Lib_Riemann_LF_SolverLibrary Lib_Riemann_LF_Solver
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Riemann_LF_SolverPublicProcedure Lib_Riemann_LF_Solver
!> @}

!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws my means of
!> LF (Primitive Variables Linearization) approximation.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Riemann_LF_SolverLibrary
module Lib_Riemann_LF_Solver
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                        ! precision of integers and reals
USE Data_Type_Riemann_Conservative1D, only: Type_Riemann_Conservative1D ! Definition of Type_Riemann_Conservative1D.
USE Data_Type_Riemann_InterState1D,   only: Type_Riemann_InterState1D   ! Definition of Type_Riemann_InterState1D.
USE Data_Type_Riemann_Primitive1D,    only: Type_Riemann_Primitive1D    ! Definition of Type_Riemann_Primitive1D.
USE Lib_Riemann_PVL_Solver                                              ! Library of PVL Riemann solver.
USE Lib_Riemann_Z_Solver                                                ! Library of Z Riemann solver.
USE Lib_Thermodynamic_Laws_Ideal                                        ! Library of thermodynamic laws for ideal gas.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: LF_set_compute_waves14
public:: LF_solver
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
procedure(CW_interface), pointer:: compute_waves14 => null()
!> Abstract interfaces for pointer procedures.
abstract interface
  pure subroutine CW_interface(state1,state4,state23)
  import :: Type_Riemann_InterState1D, Type_Riemann_Primitive1D
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1,state4
  type(Type_Riemann_InterState1D), intent(INOUT):: state23
  endsubroutine CW_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Riemann_LF_SolverPublicProcedure
  !> @{
  !> @brief Procedure for setting the algorithm for computing waves speed accordingly to the Riemann solver selected.
  subroutine LF_set_compute_waves14(solver)
  !---------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN):: solver !< Type of solver selected:
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(trim(adjustl(solver)))
  case("LFpWSu")
    compute_waves14 => PVL_compute_waves14_u
  case("LFpWSup")
    compute_waves14 => PVL_compute_waves14_up
  case("LFz")
    compute_waves14 => Z_compute_waves14_up
  case default
    compute_waves14 => PVL_compute_waves14_up
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine LF_set_compute_waves14

  !> Approximate Riemann solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm.
  pure subroutine LF_solver(state1,state4,flux)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),    intent(IN)::  state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),    intent(IN)::  state4  !< State 4 (right).
  type(Type_Riemann_Conservative1D), intent(OUT):: flux    !< Fluxes of conservative variables.
  type(Type_Riemann_Conservative1D)::              flux1   !< Fluxes of state 1.
  type(Type_Riemann_Conservative1D)::              flux4   !< Fluxes of state 4.
  type(Type_Riemann_InterState1D)::                state23 !< Intermediate state 2,3.
  real(R8P)::                                      lmax    !< Maximum wave speed estimation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! evaluating the waves speed S1, S and S4
  call compute_waves14(state1=state1,state4=state4,state23=state23)
  ! evalutaing the maximum waves speed
  lmax = max(abs(state23%S1),abs(state23%S4))
  ! computing the fluxes of state 1 and 4
  call flux1%compute(state=state1)
  call flux4%compute(state=state4)
  ! computing fluxes
  flux%r=0.5_R8P*(flux1%r+flux4%r-lmax*(state4%r                                               -state1%r))
  flux%u=0.5_R8P*(flux1%u+flux4%u-lmax*(state4%r*state4%u                                      -state1%r*state1%u))
  flux%E=0.5_R8P*(flux1%E+flux4%E-lmax*(state4%r*E(p=state4%p,r=state4%r,u=state4%u,g=state4%g)-&
    state1%r*E(p=state1%p,r=state1%r,u=state1%u,g=state1%g)))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine LF_solver
  !> @}
endmodule Lib_Riemann_LF_Solver
