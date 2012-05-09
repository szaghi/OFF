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
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: operator (/=)
public:: operator (<)
public:: operator (<=)
public:: operator (==)
public:: operator (>=)
public:: operator (>)
public:: Get_Date_String
public:: Seconds_To_Time
public:: Crono
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing time variables.
!> @ingroup DerivedType
type, public:: Type_Time
  integer(I_P) Days    !< Number of days.
  integer(I_P) Hours   !< Number of hours.
  integer(I_P) Minutes !< Number of minutes.
  real(R_P)    Seconds !< Number of seconds.
endtype Type_Time
real:: instant0 = 0.0 !< The Crono starting instant (used for timing the code).
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Not-equal-to boolean operator (/=) overloading.
!> @ingroup Interface
interface operator (/=)
  module procedure not_eq
endinterface
!> @brief Lower-than boolean operator (<) overloading.
!> @ingroup Interface
interface operator (<)
  module procedure low
endinterface
!> @brief Lower-equal-than boolean operator (<=) overloading.
!> @ingroup Interface
interface operator (<=)
  module procedure low_eq
endinterface
!> @brief Equal-to boolean operator (==) overloading.
!> @ingroup Interface
interface operator (==)
  module procedure eq
endinterface
!> @brief Higher-equal-than boolean operator (>=) overloading.
!> @ingroup Interface
interface operator (>=)
  module procedure great_eq
endinterface
!> @brief Higher-than boolean operator (>) overloading.
!> @ingroup Interface
interface operator (>)
  module procedure great
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_TimePrivateProcedure
  !> @{
  function not_eq(time1,time2) result(compare)
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

  function low(time1,time2) result(compare)
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

  function low_eq(time1,time2) result(compare)
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

  function eq(time1,time2) result(compare)
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

  function great_eq(time1,time2) result(compare)
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

  function great(time1,time2) result(compare)
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
  integer(I_P)::  Date_Integer(1:8)            !< Integer array for handling the date.
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

  !> Function for converting seconds to Type_Time (days,hours,minutes,seconds) derived type.
  !> @return \b Time type(Type_Time) variable.
  function Seconds_To_Time(seconds) result(Time)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: seconds !< Number of seconds.
  type(Type_Time)::       Time    !< Time in days,hours,minutes,seconds format.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Time%Days    = int(Seconds/86400)
  Time%Hours   = int((Seconds-Time%Days*86400)/3600)
  Time%Minutes = int((Seconds-Time%Days*86400-Time%Hours*3600)/60)
  Time%Seconds =      Seconds-Time%Days*86400-Time%Hours*3600-Time%Minutes*60
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Seconds_To_Time

  !> A simple stop/watch function.
  !> @return \b seconds real variable.
  function Crono(start,instant1) result(seconds)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  logical, intent(IN), optional:: start    !< Flag for starting time measurament.
  real,    intent(IN), optional:: instant1 !< Starting instant1 (external supplied) different from instant0.
  real::                          seconds  !< Seconds from instant0 or instant1.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call CPU_TIME(seconds)         ! Getting seconds.
  if (present(start)) then
    instant0 = seconds           ! Initialize instant0.
  elseif (present(instant1)) then
    seconds = seconds - instant1 ! Shifting time from instant1 (external supplied).
  else
    seconds = seconds - instant0 ! Shifting time from instant0 (internal supplied).
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Crono
  !> @}
endmodule Data_Type_Time
