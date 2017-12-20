 !> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_AdimensionalDerivedType Data_Type_Adimensional
!> Module definition of Type_Adimensional
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_AdimensionalInterface Data_Type_Adimensional
!> Module definition of Type_Adimensional
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_AdimensionalPrivateProcedure Data_Type_Adimensional
!> Module definition of Type_Adimensional
!> @}

!> @brief Module Data_Type_Adimensional contains the definition of Type_Adimensional and useful procedures for its handling.
!> Type_Adimensional contains all the data defining the main dimensions of block.
module Data_Type_Adimensional
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing global-level non-dimensional numbers and reference values.
!> @ingroup Data_Type_AdimensionalDerivedType
type, public:: Type_Adimensional
  ! non dimensional numbers loaded from input file
  real(R8P):: Re = 1._R8P !< \f$\rm{Re}=\frac{\rho_0 v_0 L_0}{\mu_0}\f$ Reynolds number.
  real(R8P):: Fr = 1._R8P !< \f$\rm{Fr}=\sqrt{\frac{v_0^2}{f_0 L_0}}\f$ Froude number.
  real(R8P):: Pr = 1._R8P !< \f$\rm{Pr}=\frac{\mu_0 c_p}{k_0}\f$ Prandtl number.
  ! reference values loaded from input file
  real(R8P):: L0 = 1._R8P !< Reference length.
  real(R8P):: r0 = 1._R8P !< Reference density.
  real(R8P):: v0 = 1._R8P !< Reference velocity.
  real(R8P):: c0 = 1._R8P !< Reference specific heats (\f$cp_0 = cv_0 = R_0 = c_0\f$).
  ! reference values computed by means of previous values
  real(R8P):: mu0  = 1._R8P !< \f$\mu_0= \frac{\rho_0 v_0 L_0}{\rm{Re}}\f$ Reference dynamic viscosity.
  real(R8P):: f0   = 1._R8P !< \f$f_0= \frac{v_0^2}{L_0 \rm{Fr}^2}\f$ Reference specific force.
  real(R8P):: k0   = 1._R8P !< \f$k_0= \frac{\mu_0 c_0}{\rm{Pr}}\f$ Reference thermal conductivity coefficient.
  real(R8P):: Dt0  = 1._R8P !< \f$Dt_0=\frac{L_0}{v_0}\f$ Reference time interval.
  real(R8P):: p0   = 1._R8P !< \f$p_0=\rho_0 v_0^2\f$ Reference pressure.
  real(R8P):: a0   = 1._R8P !< \f$a_0=v_0\f$ Reference speed of sound.
  real(R8P):: T0   = 1._R8P !< \f$T_0=\frac{v_0^2}{c_0}\f$ Reference temperature.
  real(R8P):: E0   = 1._R8P !< \f$E_0=v_0^2\f$ Reference specific energy.
  real(R8P):: q0   = 1._R8P !< \f$q_0=\frac{v_0^3}{L_0}\f$ Reference specific heat.
  ! equations coefficients computed by means of previous values
  real(R8P):: Re_inv   = 1._R8P !< \f$\frac{1}{\rm{Re}}\f$ Inverse of Reynolds number (coefficient of viscous terms).
  real(R8P):: Fr2_inv  = 1._R8P !< \f$\frac{1}{\rm{Fr}^2}\f$ Inverse of square of Froude number (coefficient of volume forces).
  real(R8P):: PrRe_inv = 1._R8P !< \f$\frac{1}{\rm{Pr Re}}\f$ Inverse of Prandtl and Reynolds numbers (coef. of condution terms).
  contains
    procedure:: compute_values0     ! Procedure for computing non dimensional quantities.
    procedure:: print => print_adim ! Procedure for printing non dimensional quantities with a pretty format.
    ! operators overloading
    generic:: assignment(=) => assign_ad
    ! private procedures
    procedure, pass(ad1), private:: assign_ad
endtype Type_Adimensional
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_AdimensionalPrivateProcedure
  !> @{
  !> @brief Procedure for computing the reference values for non-dimensional quantities.
  subroutine compute_values0(adim)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Adimensional), intent(INOUT):: adim !< Non dimensional data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! reference values
  adim%mu0  = (adim%r0*adim%v0*adim%L0)/adim%Re
  adim%f0   = (adim%v0*adim%v0)/(adim%L0*adim%Fr*adim%Fr)
  adim%k0   = (adim%mu0*adim%c0)/adim%Pr
  adim%Dt0  = adim%L0/adim%v0
  adim%p0   = adim%r0*adim%v0*adim%v0
  adim%a0   = adim%v0
  adim%T0   = (adim%v0*adim%v0)/adim%c0
  adim%E0   = adim%v0*adim%v0
  adim%q0   = (adim%v0*adim%v0*adim%v0)/adim%L0
  ! equations coefficients
  adim%Re_inv   = 1._R_P/adim%Re
  adim%Fr2_inv  = 1._R_P/(adim%Fr*adim%Fr)
  adim%PrRe_inv = 1._R_P/(adim%Pr*adim%Re)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_values0

  !> @brief Procedure for printing non dimensional quantities with a pretty format.
  subroutine print_adim(adim,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Adimensional), intent(IN)::  adim    !< Non dimensional data.
  character(*), optional,   intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional,   intent(OUT):: iostat  !< IO error.
  character(*), optional,   intent(OUT):: iomsg   !< IO error message.
  integer(I4P),             intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::         prefd   !< Prefixing string.
  integer(I4P)::                          iostatd !< IO error.
  character(500)::                        iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Re      : '//trim(str(n=adim%Re      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Fr      : '//trim(str(n=adim%Fr      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Pr      : '//trim(str(n=adim%Pr      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' L0      : '//trim(str(n=adim%L0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' r0      : '//trim(str(n=adim%r0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' v0      : '//trim(str(n=adim%v0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' c0      : '//trim(str(n=adim%c0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' mu0     : '//trim(str(n=adim%mu0     ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' f0      : '//trim(str(n=adim%f0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' k0      : '//trim(str(n=adim%k0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Dt0     : '//trim(str(n=adim%Dt0     ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' p0      : '//trim(str(n=adim%p0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' a0      : '//trim(str(n=adim%a0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' T0      : '//trim(str(n=adim%T0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' E0      : '//trim(str(n=adim%E0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' q0      : '//trim(str(n=adim%q0      ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Re_inv  : '//trim(str(n=adim%Re_inv  ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Fr2_inv : '//trim(str(n=adim%Fr2_inv ))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' PrRe_inv: '//trim(str(n=adim%PrRe_inv))
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_adim

  ! Assignment (=)
  !> @brief Procedure for assignment between two non-dimensional variables.
  elemental subroutine assign_ad(ad1,ad2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Adimensional), intent(INOUT):: ad1
  type(Type_Adimensional),  intent(IN)::    ad2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ad1%Re       = ad2%Re
  ad1%Fr       = ad2%Fr
  ad1%Pr       = ad2%Pr
  ad1%L0       = ad2%L0
  ad1%r0       = ad2%r0
  ad1%v0       = ad2%v0
  ad1%c0       = ad2%c0
  ad1%mu0      = ad2%mu0
  ad1%f0       = ad2%f0
  ad1%k0       = ad2%k0
  ad1%Dt0      = ad2%Dt0
  ad1%p0       = ad2%p0
  ad1%a0       = ad2%a0
  ad1%T0       = ad2%T0
  ad1%E0       = ad2%E0
  ad1%q0       = ad2%q0
  ad1%Re_inv   = ad2%Re_inv
  ad1%Fr2_inv  = ad2%Fr2_inv
  ad1%PrRe_inv = ad2%PrRe_inv
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ad
  !> @}
endmodule Data_Type_Adimensional
