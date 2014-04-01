!> @ingroup Library
!> @{
!> @defgroup Lib_Riemann_HLLC_SolverLibrary Lib_Riemann_HLLC_Solver
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Riemann_HLLC_SolverPublicProcedure Lib_Riemann_HLLC_Solver
!> @}

!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws my means of
!> HLLC approximation.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Riemann_HLLC_SolverLibrary
module Lib_Riemann_HLLC_Solver
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                        ! precision of integers and reals
USE Data_Type_Riemann_Conservative1D, only: Type_Riemann_Conservative1D ! Definition of Type_Riemann_Conservative1D.
USE Data_Type_Riemann_InterState1D,   only: Type_Riemann_InterState1D   ! Definition of Type_Riemann_InterState1D.
USE Data_Type_Riemann_Primitive1D,    only: Type_Riemann_Primitive1D    ! Definition of Type_Riemann_Primitive1D.
USE Lib_Riemann_CVL_Solver                                              ! Library of CVL Riemann solver.
USE Lib_Riemann_PVL_Solver                                              ! Library of PVL Riemann solver.
USE Lib_Riemann_TR_Solver                                               ! Library of TR Riemann solver.
USE Lib_Riemann_TS_Solver                                               ! Library of TR Riemann solver.
USE Lib_Riemann_Z_Solver                                                ! Library of Z Riemann solver.
USE Lib_Thermodynamic_Laws_Ideal                                        ! Library of thermodynamic laws for ideal gas.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: HLLC_set_compute_waves14
public:: HLLC_solver
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
procedure(CW_interface), pointer:: compute_waves14 => null()
!> Abstract interfaces for pointer procedures.
abstract interface
  elemental subroutine CW_interface(state1,state4,state23)
  use Data_Type_Riemann_Primitive1D
  use Data_Type_Riemann_InterState1D
  type(Type_Riemann_Primitive1D),  intent(IN)::    state1,state4
  type(Type_Riemann_InterState1D), intent(INOUT):: state23
  endsubroutine CW_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Riemann_HLLC_SolverPublicProcedure
  !> @{
  !> @brief Procedure for setting the algorithm for computing waves speed accordingly to the Riemann solver selected.
  subroutine HLLC_set_compute_waves14(solver)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: solver !< Type of solver selected:
                                    !< - "HLLCpWSu"  => approximate HLLC using PVL waves speed estimation based on u;
                                    !< - "HLLCpWSup" => approximate HLLC using PVL waves speed estimation based on u,p;
                                    !< - "HLLCcWSu"  => approximate HLLC using CVL waves speed estimation based on u;
                                    !< - "HLLCcWSup" => approximate HLLC using CVL waves speed estimation based on u,p;
                                    !< - "HLLCtWSu"  => approximate HLLC using TR  waves speed estimation based on u;
                                    !< - "HLLCtWSup" => approximate HLLC using TR  waves speed estimation based on u,p;
                                    !< - "HLLCz"     => approximate HLLC using Z   waves speed estimation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(trim(adjustl(solver)))
  case("HLLCpWSu")
    compute_waves14 => PVL_compute_waves14_u
  case("HLLCpWSup")
    compute_waves14 => PVL_compute_waves14_up
  case("HLLCcWSu")
    compute_waves14 => CVL_compute_waves14_u
  case("HLLCcWSup")
    compute_waves14 => CVL_compute_waves14_up
  case("HLLCtWSu")
    compute_waves14 => TR_compute_waves14_u
  case("HLLCtWSup")
    compute_waves14 => TR_compute_waves14_up
  case("HLLCz")
    compute_waves14 => Z_compute_waves14_up
  case default
    compute_waves14 => PVL_compute_waves14_up
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine HLLC_set_compute_waves14

  !> Approximate Riemann solver based on HLLC algorithm.
  !> @note Five variants are provided: HLLCb, using BCLC estimation of waves speed, HLLCc using CVL algorithm, HLLCp using PVL
  !> algorithm, HLLCt using TR algorithm and HLLCz using Z one.
  elemental subroutine HLLC_solver(state1,state4,flux)
  !--------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Riemann_Primitive1D),    intent(IN)::  state1  !< State 1 (left).
  type(Type_Riemann_Primitive1D),    intent(IN)::  state4  !< State 4 (right).
  type(Type_Riemann_Conservative1D), intent(OUT):: flux    !< Fluxes of conservative variables.
  type(Type_Riemann_InterState1D)::                state23 !< Intermediate state 2,3.
  real(R8P)::                                      U1S     !< Mass conservation.
  real(R8P)::                                      U2S     !< Momentum conservation.
  real(R8P)::                                      U3S     !< Energy conservation.
  real(R8P)::                                      E1      !< Entalpy and internal energy of left state.
  real(R8P)::                                      E4      !< Entalpy and internal energy of right state.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  ! evaluating the waves speed S1, S and S4
  call compute_waves14(state1=state1,state4=state4,state23=state23)
  state23%u = (state4%r*state4%u*(state23%S4-state4%u)-state1%r*state1%u*(state23%S1-state1%u)+state1%p-state4%p)/&
      (state4%r*(state23%S4-state4%u)-state1%r*(state23%S1-state1%u))
  ! computing fluxes
  select case(minloc([-state23%S1,state23%S1*state23%u,state23%u*state23%S4,state23%S4],dim=1))
  case(1)
    call flux%compute(state=state1)
  case(2)
    call flux%compute(state=state1)
    E1  = E(p=state1%p,r=state1%r,u=state1%u,g=state1%g)
    U1S = state1%r*(state23%S1-state1%u)/(state23%S1-state23%u)
    U2S = U1S*state23%u
    U3S = U1S*(E1+(state23%u-state1%u)*(state23%u+state1%p/(state1%r*(state23%S1-state1%u))))

    flux%r = flux%r + state23%S1*(U1S - state1%r)
    flux%u = flux%u + state23%S1*(U2S - state1%r*state1%u)
    flux%E = flux%E + state23%S1*(U3S - state1%r*E1)
  case(3)
    call flux%compute(state=state4)
    E4  = E(p=state4%p,r=state4%r,u=state4%u,g=state4%g)
    U1S = state4%r*(state23%S4-state4%u)/(state23%S4-state23%u)
    U2S = U1S*state23%u
    U3S = U1S*(E4+(state23%u-state4%u)*(state23%u+state4%p/(state4%r*(state23%S4-state4%u))))

    flux%r = flux%r + state23%S4*(U1S - state4%r)
    flux%u = flux%u + state23%S4*(U2S - state4%r*state4%u)
    flux%E = flux%E + state23%S4*(U3S - state4%r*E4)
  case(4)
    call flux%compute(state=state4)
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endsubroutine HLLC_solver
  !> @}
endmodule Lib_Riemann_HLLC_Solver
