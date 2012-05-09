!> This module contains the definition of Thermodynamic Laws for ideal calorically perfect gas.
!> @ingroup Library
module Lib_Thermodynamic_Laws_Ideal
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Definition of integers and reals precision.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: p
public:: r
public:: a
public:: E
public:: H
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> Function for computing the pressure for an ideal calorically perfect gas.
  elemental function p(r,a,g) result(pressure)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: r        !< Density (\f$\rho\f$).
  real(R_P), intent(IN):: a        !< Speed of sound.
  real(R_P), intent(IN):: g        !< Specific heats ratio \f$\gamma=\frac{{c_p}}{{c_v}}\f$.
  real(R_P)               pressure !< Pressure.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pressure = r*a*a/g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction p

  !> Function for computing the density for an ideal calorically perfect gas.
  elemental function r(p,a,g) result(density)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: p       !< Pressure.
  real(R_P), intent(IN):: a       !< Speed of sound.
  real(R_P), intent(IN):: g       !< Specific heats ratio \f$\gamma=\frac{{c_p}}{{c_v}}\f$.
  real(R_P)               density !< Density (\f$\rho\f$).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  density = g*p/(a*a)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction r

  !> Function for computing the speed of sound for an ideal calorically perfect gas.
  elemental function a(p,r,g) result(ss)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: p  !< Pressure.
  real(R_P), intent(IN):: r  !< Density (\f$\rho\f$).
  real(R_P), intent(IN):: g  !< Specific heats ratio \f$\gamma=\frac{{c_p}}{{c_v}}\f$.
  real(R_P)               ss !< Speed of sound.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ss = sqrt(g*p/r)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction a

  !> Function for computing the total specific energy (per unit of mass) for an ideal calorically perfect gas.
  !> @note This law is defined as: \n
  !> \f$ E = \frac{p}{{\left( {\gamma  - 1} \right)\rho }} + \frac{{u^2 }}{2} \f$
  elemental function E(p,r,u,g) result(energy)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: p      !< Pressure.
  real(R_P), intent(IN):: r      !< Density (\f$\rho\f$).
  real(R_P), intent(IN):: u      !< Module of velocity vector.
  real(R_P), intent(IN):: g      !< Specific heats ratio \f$\gamma=\frac{{c_p}}{{c_v}}\f$.
  real(R_P)               energy !< Total specific energy (per unit of mass).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  energy = p/((g-1._R_P)*r)+0.5_R_P*u*u
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction E

  !> Function for computing the total specific entalpy (per unit of mass) for an ideal calorically perfect gas.
  !> @note This law is defined as: \n
  !> \f$ H = \frac{{\gamma p}}{{\left( {\gamma  - 1} \right)\rho }} + \frac{{u^2 }}{2} \f$
  elemental function H(p,r,u,g) result(entalpy)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: p       !< Pressure.
  real(R_P), intent(IN):: r       !< Density (\f$\rho\f$).
  real(R_P), intent(IN):: u       !< Module of velocity vector.
  real(R_P), intent(IN):: g       !< Specific heats ratio \f$\gamma=\frac{{c_p}}{{c_v}}\f$.
  real(R_P)               entalpy !< Total specific entalpy (per unit of mass).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  entalpy = g*p/((g-1._R_P)*r)+0.5_R_P*u*u
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction H
endmodule Lib_Thermodynamic_Laws_Ideal
