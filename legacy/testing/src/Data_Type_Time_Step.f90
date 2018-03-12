!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_Time_StepDerivedType Data_Type_Time_Step
!> Module definition of Type_Time_Step
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_Time_StepInterface Data_Type_Time_Step
!> Module definition of Type_Time_Step
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_Time_StepPublicProcedure Data_Type_Time_Step
!> Module definition of Type_Time_Step
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_Time_StepPrivateProcedure Data_Type_Time_Step
!> Module definition of Type_Time_Step
!> @}

!> @brief Module Data_Type_Time_Step contains the definition of Type_Time_Step, that defines the mains variables concerning with
!> time stepping, namely time integration accuracy.
module Data_Type_Time_Step
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_Time_Step.
!> @ingroup Data_Type_Time_StepDerivedType
type, public:: Type_Time_Step
  integer(I8P):: n        = 0_I8P   !< Time steps counter.
  real(R8P)::    t        = 0._R8P  !< Time.
  integer(I8P):: Nmax     = 0_I8P   !< Max number of iterations.
  real(R8P)::    Tmax     = 0._R8P  !< Max time, ignored if Nmax>0.
  integer(I1P):: rk_ord   = 1_I1P   !< Order of time convergence (number of Runge-Kutta stages).
  logical::      unsteady = .true.  !< Type of simulation: unsteady or not.
  real(R8P)::    CFL      = 0.3_R8P !< Value of stability coefficient.
  contains
    procedure:: set_limits               ! Procedure for setting limits for simulation stop condition.
    procedure:: is_the_end               ! Procedure for checking if the end of simulation is reached.
    procedure:: is_to_save               ! Procedure for testing if is time to save final solution accordingly to Nmax or Tmax.
    procedure:: progress                 ! Procedure for evaluating the progress of current simulation.
    procedure:: update                   ! Procedure for updating time stepping data.
    procedure:: print => print_time_step ! Procedure for printing time stepping data with a pretty format.
endtype Type_Time_Step
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_Time_StepPrivateProcedure
  !> @{
  !> @brief Procedure for setting limits for simulation stop condition.
  elemental subroutine set_limits(time_step)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time_Step), intent(INOUT):: time_step !< Time stepping data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (time_step%unsteady) then
    if (time_step%Nmax>0_I8P) then
      time_step%Tmax=-1._R8P ! the value of Tmax is ignored because Nmax>0
    else
      time_step%Nmax=-1_I8P
    endif
  else
    time_step%Tmax=-1._R8P ! the value of Tmax is ignored because steady simulation
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_limits

  !> @brief Procedure for checking if the end of simulation is reached.
  elemental function  is_the_end(time_step) result(finish)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time_Step), intent(IN):: time_step !< Time stepping data.
  logical::                           finish    !< Flag returning the end reaching.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  finish = ((time_step%t==time_step%Tmax).OR.(time_step%n==time_step%Nmax))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_the_end

  !> @brief Procedure for updating time stepping data.
  elemental subroutine update(time_step,gDtmin)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time_Step), intent(INOUT):: time_step !< Time stepping data.
  real(R8P),             intent(INOUT):: gDtmin    !< Global (all processes/all blks) minimum time step.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (time_step%unsteady) then ! for an unsteady accurate simulation each cell is updated by means of global minimum time step
    ! control for the last iterate
    if (time_step%Nmax<=0) then
      if ((time_step%t+gDtmin)>time_step%Tmax) then
        ! the global minimum time step is so high that the last iteration will go over Tmax
        ! it is decreased both for in order to achieve exactly Tmax
        gDtmin=abs(time_step%Tmax-time_step%t)
      endif
    endif
    time_step%t = time_step%t + gDtmin
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine update

  !> @brief Procedure for printing time stepping data with a pretty format.
  subroutine print_time_step(time_step,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time_Step),  intent(IN)::  time_step !< Time stepping data.
  character(*), optional, intent(IN)::  pref      !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat    !< IO error.
  character(*), optional, intent(OUT):: iomsg     !< IO error message.
  integer(I4P),           intent(IN)::  unit      !< Logic unit.
  character(len=:), allocatable::       prefd     !< Prefixing string.
  integer(I4P)::                        iostatd   !< IO error.
  character(500)::                      iomsgd    !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  if (time_step%unsteady) then
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' timing: unsteady'
  else
    write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' timing: steady'
  endif
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Nmax: '//trim(str(n=time_step%Nmax))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Tmax: '//trim(str(n=time_step%Tmax))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' rk_ord: '//trim(str(n=time_step%rk_ord))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' CFL: '//trim(str(n=time_step%CFL))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' n: '//trim(str(n=time_step%n))
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' t: '//trim(str(n=time_step%t))
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_time_step

  !> @brief Procedure for testing if is time to save final solution accordingly to Nmax or Tmax.
  elemental function is_to_save(time_step) result(yes)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time_Step),  intent(IN):: time_step !< Time stepping data.
  logical::                            yes       !< Is to save or not.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  yes=((time_step%t==time_step%Tmax).OR.(time_step%n==time_step%Nmax))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction is_to_save

  !> @brief Procedure for evaluating the progress of current simulation.
  elemental function progress(time_step) result(prog)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time_Step), intent(IN):: time_step !< Time stepping data.
  real(R8P)::                         prog      !< Actual progress value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (time_step%Nmax>0) then
    prog = time_step%n*100/(time_step%Nmax*1._R8P)
  elseif (time_step%Tmax>0._R8P) then
    prog = 100*time_step%t/time_step%Tmax
  else
    prog = 0._R8P
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction progress
  !> @}
endmodule Data_Type_Time_Step
