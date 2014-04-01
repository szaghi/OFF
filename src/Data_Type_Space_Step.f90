!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Space_StepDerivedType Data_Type_Space_Step
!> Module definition of Type_Space_Step
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_Space_StepInterface Data_Type_Space_Step
!> Module definition of Type_Space_Step
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_Space_StepPublicProcedure Data_Type_Space_Step
!> Module definition of Type_Space_Step
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Space_StepPrivateProcedure Data_Type_Space_Step
!> Module definition of Type_Space_Step
!> @}

!> @brief Module Data_Type_Space_Step contains the definition of Type_Space_Step, that defines the mains variables concerning with
!> space stepping, namely space integration accuracy.
module Data_Type_Space_Step
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_Space_Step.
!> @ingroup Data_Type_Space_StepDerivedType
type, public:: Type_Space_Step
  integer(I1P):: sp_ord        = 1_I1P    !< Order of space convergence (number of ghost cells).
  logical::      inviscid      = .true.   !< Type of simulation: inviscid (Euler's eq.) or viscous (Navier-Stokes eq.).
  integer(I1P):: gco           = 1_I1P    !< Number of ghost cells necessary to achieve the space reconstruction order.
  real(R8P)::    residual_toll = 0.01_R8P !< Tolerance for residuals vanishing evaluation.
  logical::      residual_stop = .false.  !< Sentinel for stopping steady simulation when residuals vanish.
  contains
    procedure:: set_gc                    ! Procedure for setting the ghost cells number.
    procedure:: print => print_space_step ! Procedure for printing space stepping data with a pretty format.
endtype Type_Space_Step
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Space_StepPrivateProcedure
  !> @{
  !> @brief Procedure for setting the ghost cells number necessary to achieve the space reconstruction order selected.
  elemental subroutine set_gc(space_step)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Space_Step), intent(INOUT):: space_step !< Space stepping data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(space_step%sp_ord)
  case(1_I1P) ! 1st order piecewise constant reconstruction
    space_step%gco = 1_I1P
  case(3_I1P) ! 3rd order WENO reconstruction
    space_step%gco = 2_I1P
  case(5_I1P) ! 5th order WENO reconstruction
    space_step%gco = 3_I1P
  case(7_I1P) ! 7th order WENO reconstruction
    space_step%gco = 4_I1P
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_gc

  !> @brief Procedure for printing space stepping data with a pretty format.
  subroutine print_space_step(space_step,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Space_Step),  intent(IN)::  space_step !< Space stepping data.
  character(*), optional,  intent(IN)::  pref       !< Prefixing string.
  integer(I4P), optional,  intent(OUT):: iostat     !< IO error.
  character(*), optional,  intent(OUT):: iomsg      !< IO error message.
  integer(I4P),            intent(IN)::  unit       !< Logic unit.
  character(len=:), allocatable::        prefd      !< Prefixing string.
  integer(I4P)::                         iostatd    !< IO error.
  character(500)::                       iomsgd     !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  if (space_step%inviscid) then
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' viscosity: inviscid'
  else
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' viscosity: viscous'
  endif
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' sp_ord: '//trim(str(n=space_step%sp_ord))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' gco: '//trim(str(n=space_step%gco))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' residual tollerance: '//trim(str(n=space_step%residual_toll))
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_space_step
  !> @}
endmodule Data_Type_Space_Step
