!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_TimeDerivedType Data_Type_Time
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_TimeInterface Data_Type_Time
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup Data_Type_TimePrivateVarPar Data_Type_Time
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_TimePublicProcedure Data_Type_Time
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_TimePrivateProcedure Data_Type_Time
!> @}

!> This module contains the definition of Type_Time and its procedures.
!> This derived type is useful for handling time and date
!> @todo \b DocComplete: Complete the documentation of internal procedures
module Data_Type_Time
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
#ifdef OPENMP
USE OMP_LIB      ! OpenMP runtime library.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: operator (/=)
public:: operator (<)
public:: operator (<=)
public:: operator (==)
public:: operator (>=)
public:: operator (>)
public:: Get_Date_String
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing time variables.
!> @ingroup Data_Type_TimeDerivedType
type, public:: Type_Time
  integer(I4P):: Days    = 0_I4P  !< Number of days.
  integer(I4P):: Hours   = 0_I4P  !< Number of hours.
  integer(I4P):: Minutes = 0_I4P  !< Number of minutes.
  real(R8P)::    Seconds = 0._R8P !< Number of seconds.
  real(R8P)::    inst0   = 0._R8P !< The inital instant.
  contains
    procedure:: sec2dhms                 ! Procedure for converting seconds to day/hours/minutes/seconds format.
    procedure:: chronos                  ! Procedure for timing the codes.
    procedure:: print => print_time_self ! Procedure for printing time with a pretty format.
endtype Type_Time
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Not-equal-to boolean operator (/=) overloading.
!> @ingroup Data_Type_TimeInterface
interface operator (/=)
  module procedure not_eq
endinterface
!> @brief Lower-than boolean operator (<) overloading.
!> @ingroup Data_Type_TimeInterface
interface operator (<)
  module procedure low
endinterface
!> @brief Lower-equal-than boolean operator (<=) overloading.
!> @ingroup Data_Type_TimeInterface
interface operator (<=)
  module procedure low_eq
endinterface
!> @brief Equal-to boolean operator (==) overloading.
!> @ingroup Data_Type_TimeInterface
interface operator (==)
  module procedure eq
endinterface
!> @brief Higher-equal-than boolean operator (>=) overloading.
!> @ingroup Data_Type_TimeInterface
interface operator (>=)
  module procedure great_eq
endinterface
!> @brief Higher-than boolean operator (>) overloading.
!> @ingroup Data_Type_TimeInterface
interface operator (>)
  module procedure great
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_TimePrivateProcedure
  !> @{
  !> Subroutine for converting seconds to days/hours/minutes/seconds format.
  elemental subroutine sec2dhms(Time,seconds)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time), intent(INOUT):: Time    !< Time data.
  real(R8P),        intent(IN)::    seconds !< Number of seconds.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Time%Days    = int(Seconds/86400)
  Time%Hours   = int((Seconds-Time%Days*86400)/3600)
  Time%Minutes = int((Seconds-Time%Days*86400-Time%Hours*3600)/60)
  Time%Seconds =      Seconds-Time%Days*86400-Time%Hours*3600-Time%Minutes*60
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine sec2dhms

  !> @brief Procedure  for timing the codes.
  subroutine chronos(time,start,instant1)
  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef OPENMP
  USE omp_lib ! OpenMP runtime library.
#elif defined _MPI
  USE MPI ! MPI runtime library.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time),    intent(INOUT):: time     !< Time date.
  logical,   optional, intent(IN)::    start    !< Flag for starting time measurament.
  real(R8P), optional, intent(IN)::    instant1 !< Seconds from instant1 (external supplied, different from inst0).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! getting actual seconds
#ifdef OPENMP
  time%seconds = omp_get_wtime()
#elif defined _MPI
  time%seconds = MPI_Wtime()
#else
  call CPU_TIME(time%seconds)
#endif
  ! updating time
  if (present(start)) then
    time%inst0 = time%seconds                ! Initialize instant0.
  elseif (present(instant1)) then
    time%seconds = time%seconds - instant1   ! Shifting time from instant1 (external supplied).
  else
    time%seconds = time%seconds - time%inst0 ! Shifting time from instant0 (internal supplied).
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine chronos

  !> @brief Procedure for printing time with a pretty format.
  subroutine print_time_self(time,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Time),       intent(IN)::  time    !< Time data.
  character(*), optional, intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::       prefd   !< Prefixing string.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' (D/H/M/S):'//trim(str(.true.,time%Days   ))//'/'//&
                                                                              trim(str(.true.,time%Hours  ))//'/'//&
                                                                              trim(str(.true.,time%Minutes))//'/'//&
                                                                              trim(str('(F6.2)',time%Seconds))
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_time_self

  elemental function not_eq(time1,time2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the the time time1 is /= with respect the time time2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Time), intent(IN):: time1   ! First time.
  type(Type_Time), intent(IN):: time2   ! Second time.
  logical::                     compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = ((time1%Days   /=time2%Days)   .OR. &
             (time1%Hours  /=time2%Hours)  .OR. &
             (time1%Minutes/=time2%Minutes).OR. &
             (time1%Seconds/=time2%Seconds))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction not_eq

  elemental function low(time1,time2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the the time time1 is < with respect the time time2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Time), intent(IN):: time1   ! First time.
  type(Type_Time), intent(IN):: time2   ! Second time.
  logical::                     compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = .false.
  if (time1%Days<time2%Days) then
    compare = .true.
    return
  elseif (time1%Days==time2%Days) then
    if (time1%Hours<time2%Hours) then
      compare = .true.
      return
    elseif (time1%Hours==time2%Hours) then
      if (time1%Minutes<time2%Minutes) then
        compare = .true.
        return
      elseif (time1%Minutes==time2%Minutes) then
        if (time1%Seconds<time2%Seconds) then
          compare = .true.
          return
        else
          compare = .false.
          return
        endif
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction low

  elemental function low_eq(time1,time2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the the time time1 is <= with respect the time time2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Time), intent(IN):: time1   ! First time.
  type(Type_Time), intent(IN):: time2   ! Second time.
  logical::                     compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = .false.
  if (time1%Days<=time2%Days) then
    if (time1%Hours<=time2%Hours) then
      if (time1%Minutes<=time2%Minutes) then
        if (time1%Seconds<=time2%Seconds) then
          compare = .true.
          return
        endif
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction low_eq

  elemental function eq(time1,time2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the the time time1 is = with respect the time time2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Time), intent(IN):: time1   ! First time.
  type(Type_Time), intent(IN):: time2   ! Second time.
  logical::                     compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = ((time1%Days   ==time2%Days)   .AND. &
             (time1%Hours  ==time2%Hours)  .AND. &
             (time1%Minutes==time2%Minutes).AND. &
             (time1%Seconds==time2%Seconds))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction eq

  elemental function great_eq(time1,time2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the the time time1 is >= with respect the time time2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Time), intent(IN):: time1   ! First time.
  type(Type_Time), intent(IN):: time2   ! Second time.
  logical::                     compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = .false.
  if (time1%Days>time2%Days) then
    compare = .true.
    return
  elseif (time1%Days==time2%Days) then
    if (time1%Hours>time2%Hours) then
      compare = .true.
      return
    elseif (time1%Hours==time2%Hours) then
      if (time1%Minutes>time2%Minutes) then
        compare = .true.
        return
      elseif (time1%Minutes==time2%Minutes) then
        if (time1%Seconds>=time2%Seconds) then
          compare = .true.
          return
        else
          compare = .false.
          return
        endif
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction great_eq

  elemental function great(time1,time2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the the time time1 is > with respect the time time2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Time), intent(IN):: time1   ! First time.
  type(Type_Time), intent(IN):: time2   ! Second time.
  logical::                     compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = .false.
  if (time1%Days>time2%Days) then
    compare = .true.
    return
  elseif (time1%Days==time2%Days) then
    if (time1%Hours>time2%Hours) then
      compare = .true.
      return
    elseif (time1%Hours==time2%Hours) then
      if (time1%Minutes>time2%Minutes) then
        compare = .true.
        return
      elseif (time1%Minutes==time2%Minutes) then
        if (time1%Seconds>time2%Seconds) then
          compare = .true.
          return
        else
          compare = .false.
          return
        endif
      endif
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction great
  !> @}

  !> @ingroup Data_Type_TimePublicProcedure
  !> @{
  !> Function for getting actual date and returning it into a string.
  !> @return \b Date_String character(20) variable.
  function Get_Date_String() result(Date_String)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(20):: Date_String                  !< String with actual date.
  integer(I4P)::  Date_Integer(1:8)            !< Integer array for handling the date.
  character(4)::  Year                         !< Dummy variable for year string.
  character(2)::  Day,Month,Hour,Minute,Second !< Dummies variables for day, month, hour, minute and second strings.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! getting the date
  call date_and_time(values=Date_Integer)
  ! integer to string conversion
  Day    = adjustl(trim(strz(2,Date_Integer(3))))
  Month  = adjustl(trim(strz(2,Date_Integer(2))))
  Year   = adjustl(trim(strz(4,Date_Integer(1))))
  Hour   = adjustl(trim(strz(2,Date_Integer(5))))
  Minute = adjustl(trim(strz(2,Date_Integer(6))))
  Second = adjustl(trim(strz(2,Date_Integer(7))))
  ! assembling the final string
  Date_String = Day//'-'//Month//'-'//Year//'  '//Hour//':'//Minute//':'//Second
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Get_Date_String
  !> @}
endmodule Data_Type_Time
