!> @ingroup Library
!> @{
!> @defgroup Lib_Riemann_SolversLibrary Lib_Riemann_Solvers
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Riemann_SolversPublicProcedure Lib_Riemann_Solvers
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_Riemann_SolversPrivateProcedure Lib_Riemann_Solvers
!> @}

!> This module contains procedures for computing the solution of the Riemann Problem for the Euler's conservation laws.
!> This is a library module.
!> The Riemann Problem solvers provide as solution the (convective) fluxes
!> \f$(\overline{\overline F} - \overline{\overline {{Q_S}}})_{conv} \f$ normal to
!> the interface direction. The solvers contained into this library have a unique API. They take as input the primitive variables
!> in the left (state 1) and right (state 4) cells with respect the interface and provide the convective fluxes as output:
!> @code
!> ...
!> call Riem_Solver(state1,state4,flux)
!> ...
!> @endcode
!> where \f$p_1,\rho_1,u_1,\gamma_1\f$ define the left state, \f$p_4,\rho_4,u_4,\gamma_4\f$ define the right state and \f$F_r\f$,
!> \f$F_u\f$ and \f$F_E\f$ are the convective fluxes of mass, momentum and energy conservation, respectively.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Riemann_SolversLibrary
module Lib_Riemann_Solvers
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                        ! precision of integers and reals
USE Data_Type_Riemann_Conservative1D, only: Type_Riemann_Conservative1D ! Definition of Type_Riemann_Conservative1D.
USE Data_Type_Riemann_Primitive1D,    only: Type_Riemann_Primitive1D    ! Definition of Type_Riemann_Primitive1D.
USE Lib_Riemann_HLLC_Solver                                             ! Library of HLLC Riemann solver.
USE Lib_Riemann_LF_Solver                                               ! Library of LF Riemann solver.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: set_riemann_solver
public:: riemann_solver
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
procedure(RS_interface), pointer:: riemann_solver  => null()
!> Abstract interfaces for pointer procedures.
abstract interface
  pure subroutine RS_interface(state1,state4,flux)
  import :: Type_Riemann_Conservative1D, Type_Riemann_Primitive1D
  type(Type_Riemann_Primitive1D),    intent(IN)::  state1,state4
  type(Type_Riemann_Conservative1D), intent(OUT):: flux
  endsubroutine RS_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Riemann_SolversPublicProcedure
  !> @{
  !> @brief Procedure for setting the Riemann solver accordingly to the options selected.
  subroutine set_riemann_solver(solver)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: solver !< Type of solver selected:
                                    !< - "HLLCpWSu"  => approximate HLLC using PVL waves speed estimation based on u;
                                    !< - "HLLCpWSup" => approximate HLLC using PVL waves speed estimation based on u,p;
                                    !< - "HLLCcWSu"  => approximate HLLC using CVL waves speed estimation based on u;
                                    !< - "HLLCcWSup" => approximate HLLC using CVL waves speed estimation based on u,p;
                                    !< - "HLLCtWSu"  => approximate HLLC using TR  waves speed estimation based on u;
                                    !< - "HLLCtWSup" => approximate HLLC using TR  waves speed estimation based on u,p;
                                    !< - "HLLCz"     => approximate HLLC using Z   waves speed estimation;
                                    !< - "LFpWSu"    => approximate Lax-Friedrichs using PVL waves speed estimation based on u;
                                    !< - "LFpWSup"   => approximate Lax-Friedrichs using PVL waves speed estimation based on u,p;
                                    !< - "LFz"       => approximate Lax-Friedrichs using Z   waves speed estimation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(trim(adjustl(solver)))
  case("HLLCpWSu","HLLCpWSup","HLLCcWSu","HLLCcWSup","HLLCtWSu","HLLCtWSup","HLLCz")
    riemann_solver  => HLLC_solver ; call HLLC_set_compute_waves14(solver=solver)
  case("LFpWSu","LFpWSup","LFz")
    riemann_solver  => LF_solver   ; call LF_set_compute_waves14(  solver=solver)
  case default
    riemann_solver  => LF_solver   ; call LF_set_compute_waves14(  solver=solver)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_riemann_solver
  !> @}
endmodule Lib_Riemann_Solvers
