!> @ingroup Library
!> @{
!> @defgroup Lib_MathLibrary Lib_Math
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Lib_MathInterface Lib_Math
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_MathPublicProcedure Lib_Math
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_MathPrivateProcedure Lib_Math
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Lib_MathGlobalVarPar Lib_Math
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup Lib_MathPrivateVarPar Lib_Math
!> @}

!> This module contains mathematical procedures.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_MathLibrary
module Lib_Math
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                        ! Integers and reals precision definition.
USE Data_Type_SL_List                   ! Definition of Type_SL_List.
USE Data_Type_Vector, only: Type_Vector ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: pi,pi_R16,pi_R8,pi_R4
public:: prime
public:: isort
public:: qsort
public:: spline3
public:: stretchinglside,stretchingrside,stretching2side
public:: average
public:: digit
public:: div2
public:: interpolate1,interpolate2,interpolate3
public:: degree,radiant
public:: delta1_2o,delta2_2o
public:: abs_grad
public:: laplacian
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Pi greek definitions with parametric kind precision.
!> @ingroup Lib_MathGlobalVarPar
!> @{
real(R16P), parameter:: pi_R16 = 2._R16P*asin(1._R16P) !< \f$\pi \f$ with R16P precision.
real(R8P),  parameter:: pi_R8  = 2._R8P* asin(1._R8P)  !< \f$\pi \f$ with R8P  precision.
real(R4P),  parameter:: pi_R4  = 2._R4P* asin(1._R4P)  !< \f$\pi \f$ with R4P  precision.
real(R_P),  parameter:: pi     = 2._R_P* asin(1._R_P)  !< \f$\pi \f$ with R_P  precision.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Subroutine for computing prime factors of input integer number. The results is an array of integers containing the
!> prime factors list.
!> Example of usage:
!> @code
!> ...
!> integer(I4P):: n=246
!> integer(I4P), allocatable:: pf(:)
!> ...
!> call prime(n=n,p=pf)
!> @endcode
!> The input integer n can have all the kinds defined in IR_Precision module.
!> @note The output array containing the factors must allocatable of the same kind of input integer.
!> @ingroup Lib_MathInterface
interface prime
  module procedure prime_I8,prime_I4,prime_I2,prime_I1
endinterface
!> @brief Subroutine for performing InsertionSort with ascending order.
!> Example of usage:
!> @code
!> ...
!> integer(I4P):: a(1:100)
!> ...
!> call isort(a)
!> @endcode
!> The input array can be integer or real with all the kinds defined in IR_Precision module.
!> @ingroup Lib_MathInterface
interface isort
#ifdef r16p
  module procedure isort_R16
#endif
  module procedure isort_R8,isort_R4,isort_I8,isort_I4,isort_I2,isort_I1
endinterface
!> @brief Subroutine for performing QuickSort with ascending order.
!> Example of usage:
!> @code
!> ...
!> integer(I4P):: a(1:100)
!> ...
!> call qsort(a)
!> @endcode
!> The input array can be integer or real with all the kinds defined in IR_Precision module.
!> @ingroup Lib_MathInterface
interface qsort
#ifdef r16p
  module procedure qsort_R16
#endif
  module procedure qsort_R8,qsort_R4,qsort_I8,qsort_I4,qsort_I2,qsort_I1
endinterface
!> @brief Function for computing cubic spline interpolation of an array of values.
!> This function is an interface to three different functions: there is one function for each of the three real kinds supported
!> (R16P, R8P and R4P).
!> @note the input arrays "x", "v" and "vi" must be ordered with increasing values. Details of implementation can be found in
!> "Numerical recipes: The art of scientific computing", third edition, 2007, ISBN 0521880688.
!> @ingroup Lib_MathInterface
interface spline3
#ifdef r16p
  module procedure spline3_R16
#endif
  module procedure spline3_R8,spline3_R4
endinterface
!> Example of usage:
!> @code
!> ...
!> real(R4P):: x(1:100),v(1:100),xi(1:200),vi(1:200)
!> ...
!> vi = spline3(x,v,xi)
!> @endcode
!> The real input arrays can have all the kinds defined in IR_Precision module.

!>Average overloading
!> @ingroup Lib_MathInterface
interface average
#ifdef r16p
  module procedure average_Vectorial1D_R16,average_Vectorial2D_R16,average_Vectorial3D_R16
#endif
  module procedure average_Vectorial1D_R8,average_Vectorial2D_R8,average_Vectorial3D_R8
  module procedure average_Vectorial1D_R4,average_Vectorial2D_R4,average_Vectorial3D_R4
endinterface
!>Digit overloading
!> @ingroup Lib_MathInterface
interface digit
  module procedure digit_I8,digit_I4,digit_I2,digit_I1
endinterface
!>Div2 overloading
!> @ingroup Lib_MathInterface
interface div2
  module procedure div2_I8,div2_I4,div2_I2,div2_I1
endinterface
!>Linear interpolate overloading
!> @ingroup Lib_MathInterface
interface interpolate1
#ifdef r16p
  module procedure interpolate1_R16
#endif
  module procedure interpolate1_R8,interpolate1_R4
endinterface
!>Bi-linear interpolate overloading
!> @ingroup Lib_MathInterface
interface interpolate2
#ifdef r16p
  module procedure interpolate2_R16
#endif
  module procedure interpolate2_R8,interpolate2_R4
endinterface
!>Tri-linear interpolate overloading
!> @ingroup Lib_MathInterface
interface interpolate3
#ifdef r16p
  module procedure interpolate3_R16
#endif
  module procedure interpolate3_R8,interpolate3_R4
endinterface
!>Degree overloading
!> @ingroup Lib_MathInterface
interface degree
#ifdef r16p
  module procedure degree_Scalar_R16,degree_Vectorial1D_R16,degree_Vectorial2D_R16,degree_Vectorial3D_R16
#endif
  module procedure degree_Scalar_R8,degree_Vectorial1D_R8,degree_Vectorial2D_R8,degree_Vectorial3D_R8
  module procedure degree_Scalar_R4,degree_Vectorial1D_R4,degree_Vectorial2D_R4,degree_Vectorial3D_R4
endinterface
!>Radiant overloading
!> @ingroup Lib_MathInterface
interface radiant
#ifdef r16p
  module procedure radiant_Scalar_R16,radiant_Vectorial1D_R16,radiant_Vectorial2D_R16,radiant_Vectorial3D_R16
#endif
  module procedure radiant_Scalar_R8,radiant_Vectorial1D_R8,radiant_Vectorial2D_R8,radiant_Vectorial3D_R8
  module procedure radiant_Scalar_R4,radiant_Vectorial1D_R4,radiant_Vectorial2D_R4,radiant_Vectorial3D_R4
endinterface
!>delta1_2o overloading
!> @ingroup Lib_MathInterface
interface delta1_2o
#ifdef r16p
  module procedure delta1_2o_R16
#endif
  module procedure delta1_2o_R8,delta1_2o_R4
endinterface
!>delta2_2o overloading
!> @ingroup Lib_MathInterface
interface delta2_2o
#ifdef r16p
  module procedure delta2_2o_R16
#endif
  module procedure delta2_2o_R8,delta2_2o_R4
endinterface
!>abs_grad overloading
!> @ingroup Lib_MathInterface
interface abs_grad
#ifdef r16p
  module procedure abs_grad_R16
#endif
  module procedure abs_grad_R8,abs_grad_R4
endinterface
!>laplacian overloading
!> @ingroup Lib_MathInterface
interface laplacian
#ifdef r16p
  module procedure laplacian_R16
#endif
  module procedure laplacian_R8,laplacian_R4
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_MathPrivateProcedure
  !> @{
  !> @brief Subroutine for computing prime factors of input integer number (I8P). The results is an array of integers containing the
  !> prime factors list.
  subroutine prime_I8(n,p)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),              intent(IN)::  n    !< Input number of which prime factors must be computed.
  integer(I8P), allocatable, intent(OUT):: p(:) !< Prime factors of n.
  type(Type_SL_List)::                     pl   !< Prime factors list.
  integer(I8P)::                           nn   !< Copy of Input number.
  integer(I8P)::                           d    !< Divisor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n==1) then
   allocate(p(1:1)) ; p = 1
   return
  endif
  nn=n
  do ! removing all factors of 2
   if (mod(nn,2_I8P)/=0.OR.nn==1) exit
   nn = nn/2_I8P
   call pl%putt(d=2_I8P)
  enddo
  d=3
  do ! removing factor 3, 5, 7, ...
   if (d>nn) exit ! if a factor is too large, exit and done
   do
     if (mod(nn,d)/=0.OR.nn==1) exit
     nn = nn/d ! remove this factor from n
     call pl%putt(d=d)
   enddo
   d = d + 2_I8P ! move to next odd number
  enddo
  call pl%array(p)
  call pl%free()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prime_I8

  !> @brief Subroutine for computing prime factors of input integer number (I4P). The results is an array of integers containing the
  !> prime factors list.
  subroutine prime_I4(n,p)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),              intent(IN)::  n    !< Input number of which prime factors must be computed.
  integer(I4P), allocatable, intent(OUT):: p(:) !< Prime factors of n.
  type(Type_SL_List)::                     pl   !< Prime factors list.
  integer(I4P)::                           nn   !< Copy of Input number.
  integer(I4P)::                           d    !< Divisor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n==1) then
    allocate(p(1:1)) ; p = 1
    return
  endif
  nn=n
  do ! removing all factors of 2
    if (mod(nn,2_I4P)/=0.OR.nn==1) exit
    nn = nn/2_I4P
    call pl%putt(d=2_I4P)
  enddo
  d=3
  do ! removing factor 3, 5, 7, ...
    if (d>nn) exit ! if a factor is too large, exit and done
    do
      if (mod(nn,d)/=0.OR.nn==1) exit
      nn = nn/d ! remove this factor from n
      call pl%putt(d=d)
    enddo
    d = d + 2_I4P ! move to next odd number
  enddo
  call pl%array(p)
  call pl%free()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prime_I4

  !> @brief Subroutine for computing prime factors of input integer number (I2P). The results is an array of integers containing the
  !> prime factors list.
  subroutine prime_I2(n,p)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),              intent(IN)::  n    !< Input number of which prime factors must be computed.
  integer(I2P), allocatable, intent(OUT):: p(:) !< Prime factors of n.
  type(Type_SL_List)::                     pl   !< Prime factors list.
  integer(I2P)::                           nn   !< Copy of Input number.
  integer(I2P)::                           d    !< Divisor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n==1) then
    allocate(p(1:1)) ; p = 1
    return
  endif
  nn=n
  do ! removing all factors of 2
    if (mod(nn,2_I2P)/=0.OR.nn==1) exit
    nn = nn/2_I2P
    call pl%putt(d=2_I2P)
  enddo
  d=3
  do ! removing factor 3, 5, 7, ...
    if (d>nn) exit ! if a factor is too large, exit and done
    do
      if (mod(nn,d)/=0.OR.nn==1) exit
      nn = nn/d ! remove this factor from n
      call pl%putt(d=d)
    enddo
    d = d + 2_I2P ! move to next odd number
  enddo
  call pl%array(p)
  call pl%free()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prime_I2

  !> @brief Subroutine for computing prime factors of input integer number (I1P). The results is an array of integers containing the
  !> prime factors list.
  subroutine prime_I1(n,p)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),              intent(IN)::  n    !< Input number of which prime factors must be computed.
  integer(I1P), allocatable, intent(OUT):: p(:) !< Prime factors of n.
  type(Type_SL_List)::                     pl   !< Prime factors list.
  integer(I1P)::                           nn   !< Copy of Input number.
  integer(I1P)::                           d    !< Divisor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n==1) then
    allocate(p(1:1)) ; p = 1
    return
  endif
  nn=n
  do ! removing all factors of 2
    if (mod(nn,2_I1P)/=0.OR.nn==1) exit
    nn = nn/2_I1P
    call pl%putt(d=2_I1P)
  enddo
  d=3
  do ! removing factor 3, 5, 7, ...
    if (d>nn) exit ! if a factor is too large, exit and done
    do
      if (mod(nn,d)/=0.OR.nn==1) exit
      nn = nn/d ! remove this factor from n
      call pl%putt(d=d)
    enddo
    d = d + 2_I1P ! move to next odd number
  enddo
  call pl%array(p)
  call pl%free()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prime_I1

  ! InsertionSort
  !> @brief Subroutine for performing InsertionSort with ascending order (R16P).
  pure subroutine isort_R16(array)
  !-------------------------------------------------------------------------------------------------------------------------------
  real(R16P), intent(INOUT):: array(1:) !< Array to be sorted.
  real(R16P)::                tmp       !< Temporary array value.
  integer(I4P)::              n,nn      !< Counters.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  do n=2,size(array)
     tmp = array(n)
     if (tmp>=array(n-1)) cycle
     array(n) = array(n-1)
     do nn=n - 2,1,-1
        if (tmp>=array(nn)) exit
        array(nn+1) = array(nn)
     enddo
     array(nn+1) = tmp
  enddo
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endsubroutine isort_R16

  !> @brief Subroutine for performing InsertionSort with ascending order (R8P).
  pure subroutine isort_R8(array)
  !-------------------------------------------------------------------------------------------------------------------------------
  real(R8P), intent(INOUT):: array(1:) !< Array to be sorted.
  real(R8P)::                tmp       !< Temporary array value.
  integer(I4P)::             n,nn      !< Counters.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  do n=2,size(array)
     tmp = array(n)
     if (tmp>=array(n-1)) cycle
     array(n) = array(n-1)
     do nn=n - 2,1,-1
        if (tmp>=array(nn)) exit
        array(nn+1) = array(nn)
     enddo
     array(nn+1) = tmp
  enddo
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endsubroutine isort_R8

  !> @brief Subroutine for performing InsertionSort with ascending order (R4P).
  pure subroutine isort_R4(array)
  !-------------------------------------------------------------------------------------------------------------------------------
  real(R4P), intent(INOUT):: array(1:) !< Array to be sorted.
  real(R4P)::                tmp       !< Temporary array value.
  integer(I4P)::             n,nn      !< Counters.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  do n=2,size(array)
     tmp = array(n)
     if (tmp>=array(n-1)) cycle
     array(n) = array(n-1)
     do nn=n - 2,1,-1
        if (tmp>=array(nn)) exit
        array(nn+1) = array(nn)
     enddo
     array(nn+1) = tmp
  enddo
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endsubroutine isort_R4

  !> @brief Subroutine for performing InsertionSort with ascending order (I8P).
  pure subroutine isort_I8(array)
  !-------------------------------------------------------------------------------------------------------------------------------
  integer(I8P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I8P)::                tmp       !< Temporary array value.
  integer(I4P)::                n,nn      !< Counters.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  do n=2,size(array)
     tmp = array(n)
     if (tmp>=array(n-1)) cycle
     array(n) = array(n-1)
     do nn=n - 2,1,-1
        if (tmp>=array(nn)) exit
        array(nn+1) = array(nn)
     enddo
     array(nn+1) = tmp
  enddo
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endsubroutine isort_I8

  !> @brief Subroutine for performing InsertionSort with ascending order (I4P).
  pure subroutine isort_I4(array)
  !-------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I4P)::                tmp       !< Temporary array value.
  integer(I4P)::                n,nn      !< Counters.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  do n=2,size(array)
     tmp = array(n)
     if (tmp>=array(n-1)) cycle
     array(n) = array(n-1)
     do nn=n - 2,1,-1
        if (tmp>=array(nn)) exit
        array(nn+1) = array(nn)
     enddo
     array(nn+1) = tmp
  enddo
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endsubroutine isort_I4

  !> @brief Subroutine for performing InsertionSort with ascending order (I2P).
  pure subroutine isort_I2(array)
  !-------------------------------------------------------------------------------------------------------------------------------
  integer(I2P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I2P)::                tmp       !< Temporary array value.
  integer(I4P)::                n,nn      !< Counters.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  do n=2,size(array)
     tmp = array(n)
     if (tmp>=array(n-1)) cycle
     array(n) = array(n-1)
     do nn=n - 2,1,-1
        if (tmp>=array(nn)) exit
        array(nn+1) = array(nn)
     enddo
     array(nn+1) = tmp
  enddo
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endsubroutine isort_I2

  !> @brief Subroutine for performing InsertionSort with ascending order (I1P).
  pure subroutine isort_I1(array)
  !-------------------------------------------------------------------------------------------------------------------------------
  integer(I1P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I1P)::                tmp       !< Temporary array value.
  integer(I4P)::                n,nn      !< Counters.
  !-------------------------------------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------------------------------------
  do n=2,size(array)
     tmp = array(n)
     if (tmp>=array(n-1)) cycle
     array(n) = array(n-1)
     do nn=n - 2,1,-1
        if (tmp>=array(nn)) exit
        array(nn+1) = array(nn)
     enddo
     array(nn+1) = tmp
  enddo
  return
  !-------------------------------------------------------------------------------------------------------------------------------
  endsubroutine isort_I1

  ! QuickSort
  !> @brief Subroutine for performing QuickSort with ascending order (R16P).
  !> @note If the number of array elements is "small" the InsertionSort is directly used.
  subroutine qsort_R16(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I4P), parameter::   Nin = 16  !< Maximum dimension for insertion sort.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call subsor(array=array,n1=1,nn=size(array))
  call isort_R16(array)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Sorts array from n1 to nn (QuickSort algorithm).
    recursive subroutine subsor(array,n1,nn)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(R16P),   intent (INOUT):: array(:) !< Array to be sorted.
    integer(I4P), intent (IN)::    n1       !< Initial index of sorting.
    integer(I4P), intent (IN)::    nn       !< Final index of sorting.
    real(R16P)::                   tmp      !< Temporary array value.
    real(R16P)::                   vnM      !< Value of mean (pivot) element.
    integer(I4P)::                 nM       !< Mean (pivot) index.
    integer(I4P)::                 ni1,ni2  !< Index of the two subintervals of array.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if ((nn - n1)>Nin) then ! if there are not enough values for qsort it is left unstorted, using only the final inssor
      nM = (n1+nn)/2
      if (array(nM)<array(n1)) then
        tmp = array(n1)
        array(n1) = array(nM)
        array(nM) = tmp
      endif
      if (array(nM)>array(nn)) then
        tmp = array(nn)
        array(nn) = array(nM)
        array(nM) = tmp
        if (array(nM)<array(n1)) then
          tmp = array(n1)
          array(n1) = array(nM)
          array(nM) = tmp
        endif
      endif
      vnM = array(nM)
      ni1 = n1
      ni2 = nn
      ECH2: do
        do
          ni1 = ni1 + 1
          if (ni1>=ni2) exit ECH2 ! the first > pivot is ni2, the last <= pivot is ni1-1
          if (array(ni1)>vnM) exit
        enddo
        do
          if (array(ni2)<=vnM) exit
          ni2 = ni2 - 1
          if (ni1>=ni2) exit ECH2 ! the last < pivot is always ni1-1
        enddo
        tmp = array(ni2)
        array(ni2) = array(ni1)
        array(ni1) = tmp
      enddo ECH2
      ! sort the two subintervals
      call subsor(array=array,n1=n1, nn=ni1-1)
      call subsor(array=array,n1=ni2,nn=nn   )
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine subsor
  endsubroutine qsort_R16

  !> @brief Subroutine for performing QuickSort with ascending order (R8P).
  !> @note If the number of array elements is "small" the InsertionSort is directly used.
  subroutine qsort_R8(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I4P), parameter::  Nin = 16  !< Maximum dimension for insertion sort.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call subsor(array=array,n1=1,nn=size(array))
  call isort_R8(array)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Sorts array from n1 to nn (QuickSort algorithm).
    recursive subroutine subsor(array,n1,nn)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(R8P),    intent (INOUT):: array(:) !< Array to be sorted.
    integer(I4P), intent (IN)::    n1       !< Initial index of sorting.
    integer(I4P), intent (IN)::    nn       !< Final index of sorting.
    real(R8P)::                    tmp      !< Temporary array value.
    real(R8P)::                    vnM      !< Value of mean (pivot) element.
    integer(I4P)::                 nM       !< Mean (pivot) index.
    integer(I4P)::                 ni1,ni2  !< Index of the two subintervals of array.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if ((nn - n1)>Nin) then ! if there are not enough values for qsort it is left unstorted, using only the final inssor
      nM = (n1+nn)/2
      if (array(nM)<array(n1)) then
        tmp = array(n1)
        array(n1) = array(nM)
        array(nM) = tmp
      endif
      if (array(nM)>array(nn)) then
        tmp = array(nn)
        array(nn) = array(nM)
        array(nM) = tmp
        if (array(nM)<array(n1)) then
          tmp = array(n1)
          array(n1) = array(nM)
          array(nM) = tmp
        endif
      endif
      vnM = array(nM)
      ni1 = n1
      ni2 = nn
      ECH2: do
        do
          ni1 = ni1 + 1
          if (ni1>=ni2) exit ECH2 ! the first > pivot is ni2, the last <= pivot is ni1-1
          if (array(ni1)>vnM) exit
        enddo
        do
          if (array(ni2)<=vnM) exit
          ni2 = ni2 - 1
          if (ni1>=ni2) exit ECH2 ! the last < pivot is always ni1-1
        enddo
        tmp = array(ni2)
        array(ni2) = array(ni1)
        array(ni1) = tmp
      enddo ECH2
      ! sort the two subintervals
      call subsor(array=array,n1=n1, nn=ni1-1)
      call subsor(array=array,n1=ni2,nn=nn   )
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine subsor
  endsubroutine qsort_R8

  !> @brief Subroutine for performing QuickSort with ascending order (R4P).
  !> @note If the number of array elements is "small" the InsertionSort is directly used.
  subroutine qsort_R4(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I4P), parameter::  Nin = 16  !< Maximum dimension for insertion sort.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call subsor(array=array,n1=1,nn=size(array))
  call isort_R4(array)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Sorts array from n1 to nn (QuickSort algorithm).
    recursive subroutine subsor(array,n1,nn)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(R4P),    intent (INOUT):: array(:) !< Array to be sorted.
    integer(I4P), intent (IN)::    n1       !< Initial index of sorting.
    integer(I4P), intent (IN)::    nn       !< Final index of sorting.
    real(R4P)::                    tmp      !< Temporary array value.
    real(R4P)::                    vnM      !< Value of mean (pivot) element.
    integer(I4P)::                 nM       !< Mean (pivot) index.
    integer(I4P)::                 ni1,ni2  !< Index of the two subintervals of array.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if ((nn - n1)>Nin) then ! if there are not enough values for qsort it is left unstorted, using only the final inssor
      nM = (n1+nn)/2
      if (array(nM)<array(n1)) then
        tmp = array(n1)
        array(n1) = array(nM)
        array(nM) = tmp
      endif
      if (array(nM)>array(nn)) then
        tmp = array(nn)
        array(nn) = array(nM)
        array(nM) = tmp
        if (array(nM)<array(n1)) then
          tmp = array(n1)
          array(n1) = array(nM)
          array(nM) = tmp
        endif
      endif
      vnM = array(nM)
      ni1 = n1
      ni2 = nn
      ECH2: do
        do
          ni1 = ni1 + 1
          if (ni1>=ni2) exit ECH2 ! the first > pivot is ni2, the last <= pivot is ni1-1
          if (array(ni1)>vnM) exit
        enddo
        do
          if (array(ni2)<=vnM) exit
          ni2 = ni2 - 1
          if (ni1>=ni2) exit ECH2 ! the last < pivot is always ni1-1
        enddo
        tmp = array(ni2)
        array(ni2) = array(ni1)
        array(ni1) = tmp
      enddo ECH2
      ! sort the two subintervals
      call subsor(array=array,n1=n1, nn=ni1-1)
      call subsor(array=array,n1=ni2,nn=nn   )
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine subsor
  endsubroutine qsort_R4

  !> @brief Subroutine for performing QuickSort with ascending order (I8P).
  !> @note If the number of array elements is "small" the InsertionSort is directly used.
  subroutine qsort_I8(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I4P), parameter::     Nin = 16  !< Maximum dimension for insertion sort.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call subsor(array=array,n1=1,nn=size(array))
  call isort_I8(array)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Sorts array from n1 to nn (QuickSort algorithm).
    recursive subroutine subsor(array,n1,nn)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I8P), intent (INOUT):: array(:) !< Array to be sorted.
    integer(I4P), intent (IN)::    n1       !< Initial index of sorting.
    integer(I4P), intent (IN)::    nn       !< Final index of sorting.
    integer(I8P)::                 tmp      !< Temporary array value.
    integer(I8P)::                 vnM      !< Value of mean (pivot) element.
    integer(I4P)::                 nM       !< Mean (pivot) index.
    integer(I4P)::                 ni1,ni2  !< Index of the two subintervals of array.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if ((nn - n1)>Nin) then ! if there are not enough values for qsort it is left unstorted, using only the final inssor
      nM = (n1+nn)/2
      if (array(nM)<array(n1)) then
        tmp = array(n1)
        array(n1) = array(nM)
        array(nM) = tmp
      endif
      if (array(nM)>array(nn)) then
        tmp = array(nn)
        array(nn) = array(nM)
        array(nM) = tmp
        if (array(nM)<array(n1)) then
          tmp = array(n1)
          array(n1) = array(nM)
          array(nM) = tmp
        endif
      endif
      vnM = array(nM)
      ni1 = n1
      ni2 = nn
      ECH2: do
        do
          ni1 = ni1 + 1
          if (ni1>=ni2) exit ECH2 ! the first > pivot is ni2, the last <= pivot is ni1-1
          if (array(ni1)>vnM) exit
        enddo
        do
          if (array(ni2)<=vnM) exit
          ni2 = ni2 - 1
          if (ni1>=ni2) exit ECH2 ! the last < pivot is always ni1-1
        enddo
        tmp = array(ni2)
        array(ni2) = array(ni1)
        array(ni1) = tmp
      enddo ECH2
      ! sort the two subintervals
      call subsor(array=array,n1=n1, nn=ni1-1)
      call subsor(array=array,n1=ni2,nn=nn   )
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine subsor
  endsubroutine qsort_I8

  !> @brief Subroutine for performing QuickSort with ascending order (I4P).
  !> @note If the number of array elements is "small" the InsertionSort is directly used.
  subroutine qsort_I4(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I4P), parameter::     Nin = 16  !< Maximum dimension for insertion sort.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call subsor(array=array,n1=1,nn=size(array))
  call isort_I4(array)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Sorts array from n1 to nn (QuickSort algorithm).
    recursive subroutine subsor(array,n1,nn)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I4P), intent (INOUT):: array(:) !< Array to be sorted.
    integer(I4P), intent (IN)::    n1       !< Initial index of sorting.
    integer(I4P), intent (IN)::    nn       !< Final index of sorting.
    integer(I4P)::                 tmp      !< Temporary array value.
    integer(I4P)::                 vnM      !< Value of mean (pivot) element.
    integer(I4P)::                 nM       !< Mean (pivot) index.
    integer(I4P)::                 ni1,ni2  !< Index of the two subintervals of array.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if ((nn - n1)>Nin) then ! if there are not enough values for qsort it is left unstorted, using only the final inssor
      nM = (n1+nn)/2
      if (array(nM)<array(n1)) then
        tmp = array(n1)
        array(n1) = array(nM)
        array(nM) = tmp
      endif
      if (array(nM)>array(nn)) then
        tmp = array(nn)
        array(nn) = array(nM)
        array(nM) = tmp
        if (array(nM)<array(n1)) then
          tmp = array(n1)
          array(n1) = array(nM)
          array(nM) = tmp
        endif
      endif
      vnM = array(nM)
      ni1 = n1
      ni2 = nn
      ECH2: do
        do
          ni1 = ni1 + 1
          if (ni1>=ni2) exit ECH2 ! the first > pivot is ni2, the last <= pivot is ni1-1
          if (array(ni1)>vnM) exit
        enddo
        do
          if (array(ni2)<=vnM) exit
          ni2 = ni2 - 1
          if (ni1>=ni2) exit ECH2 ! the last < pivot is always ni1-1
        enddo
        tmp = array(ni2)
        array(ni2) = array(ni1)
        array(ni1) = tmp
      enddo ECH2
      ! sort the two subintervals
      call subsor(array=array,n1=n1, nn=ni1-1)
      call subsor(array=array,n1=ni2,nn=nn   )
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine subsor
  endsubroutine qsort_I4

  !> @brief Subroutine for performing QuickSort with ascending order (I2P).
  !> @note If the number of array elements is "small" the InsertionSort is directly used.
  subroutine qsort_I2(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I4P), parameter::     Nin = 16  !< Maximum dimension for insertion sort.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call subsor(array=array,n1=1,nn=size(array))
  call isort_I2(array)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Sorts array from n1 to nn (QuickSort algorithm).
    recursive subroutine subsor(array,n1,nn)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I2P), intent (INOUT):: array(:) !< Array to be sorted.
    integer(I4P), intent (IN)::    n1       !< Initial index of sorting.
    integer(I4P), intent (IN)::    nn       !< Final index of sorting.
    integer(I2P)::                 tmp      !< Temporary array value.
    integer(I2P)::                 vnM      !< Value of mean (pivot) element.
    integer(I4P)::                 nM       !< Mean (pivot) index.
    integer(I4P)::                 ni1,ni2  !< Index of the two subintervals of array.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if ((nn - n1)>Nin) then ! if there are not enough values for qsort it is left unstorted, using only the final inssor
      nM = (n1+nn)/2
      if (array(nM)<array(n1)) then
        tmp = array(n1)
        array(n1) = array(nM)
        array(nM) = tmp
      endif
      if (array(nM)>array(nn)) then
        tmp = array(nn)
        array(nn) = array(nM)
        array(nM) = tmp
        if (array(nM)<array(n1)) then
          tmp = array(n1)
          array(n1) = array(nM)
          array(nM) = tmp
        endif
      endif
      vnM = array(nM)
      ni1 = n1
      ni2 = nn
      ECH2: do
        do
          ni1 = ni1 + 1
          if (ni1>=ni2) exit ECH2 ! the first > pivot is ni2, the last <= pivot is ni1-1
          if (array(ni1)>vnM) exit
        enddo
        do
          if (array(ni2)<=vnM) exit
          ni2 = ni2 - 1
          if (ni1>=ni2) exit ECH2 ! the last < pivot is always ni1-1
        enddo
        tmp = array(ni2)
        array(ni2) = array(ni1)
        array(ni1) = tmp
      enddo ECH2
      ! sort the two subintervals
      call subsor(array=array,n1=n1, nn=ni1-1)
      call subsor(array=array,n1=ni2,nn=nn   )
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine subsor
  endsubroutine qsort_I2

  !> @brief Subroutine for performing QuickSort with ascending order (I1P).
  !> @note If the number of array elements is "small" the InsertionSort is directly used.
  subroutine qsort_I1(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(INOUT):: array(1:) !< Array to be sorted.
  integer(I4P), parameter::     Nin = 16  !< Maximum dimension for insertion sort.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call subsor(array=array,n1=1,nn=size(array))
  call isort_I1(array)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    !> @brief Sorts array from n1 to nn (QuickSort algorithm).
    recursive subroutine subsor(array,n1,nn)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I1P), intent (INOUT):: array(:) !< Array to be sorted.
    integer(I4P), intent (IN)::    n1       !< Initial index of sorting.
    integer(I4P), intent (IN)::    nn       !< Final index of sorting.
    integer(I1P)::                 tmp      !< Temporary array value.
    integer(I1P)::                 vnM      !< Value of mean (pivot) element.
    integer(I4P)::                 nM       !< Mean (pivot) index.
    integer(I4P)::                 ni1,ni2  !< Index of the two subintervals of array.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    if ((nn - n1)>Nin) then ! if there are not enough values for qsort it is left unstorted, using only the final inssor
      nM = (n1+nn)/2
      if (array(nM)<array(n1)) then
        tmp = array(n1)
        array(n1) = array(nM)
        array(nM) = tmp
      endif
      if (array(nM)>array(nn)) then
        tmp = array(nn)
        array(nn) = array(nM)
        array(nM) = tmp
        if (array(nM)<array(n1)) then
          tmp = array(n1)
          array(n1) = array(nM)
          array(nM) = tmp
        endif
      endif
      vnM = array(nM)
      ni1 = n1
      ni2 = nn
      ECH2: do
        do
          ni1 = ni1 + 1
          if (ni1>=ni2) exit ECH2 ! the first > pivot is ni2, the last <= pivot is ni1-1
          if (array(ni1)>vnM) exit
        enddo
        do
          if (array(ni2)<=vnM) exit
          ni2 = ni2 - 1
          if (ni1>=ni2) exit ECH2 ! the last < pivot is always ni1-1
        enddo
        tmp = array(ni2)
        array(ni2) = array(ni1)
        array(ni1) = tmp
      enddo ECH2
      ! sort the two subintervals
      call subsor(array=array,n1=n1, nn=ni1-1)
      call subsor(array=array,n1=ni2,nn=nn   )
    endif
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endsubroutine subsor
  endsubroutine qsort_I1

  ! cubic spline
  !> @brief Function for computing cubic spline interpolation of an array of values (R16P).
  pure function spline3_R16(bc1,bcn,x,v,xi) result (vi)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN), optional:: bc1            !< Boundary conditions of second derivative at first node.
  real(R16P), intent(IN), optional:: bcn            !< Boundary conditions of second derivative at last node.
  real(R16P), intent(IN)::           x(1:)          !< Independent variable abscissa of original data.
  real(R16P), intent(IN)::           v(1:)          !< Dependent variable (to be interpolated) of original data.
  real(R16P), intent(IN)::           xi(1:)         !< Independent variable abscissa of interpolating points.
  real(R16P)::                       vi(1:size(xi)) !< Dependent variable interpolated at xi points.
  real(R16P)::                       cd2(1:size(x)) !< Coefficients of second derivative.
  integer(I4P)::                     N              !< Number of points of original data.
  integer(I4P)::                     i              !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  N = size(x)
  cd2 = compute_cd2(N=N)
  do i=1,size(xi)
    vi(i) = cspline(xi(i))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure function compute_cd2(N) result(cd2)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I4P), intent(IN):: N               !< Number of points.
    real(R16P)::               cd2(1:N)        !< Coefficients of second derivative.
    real(R16P)::               sig,p,qn,u(1:N) !< Temporary variables.
    integer(I4P)::             i               !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ! first node boundary condition
    if (present(bc1)) then
      cd2(1) = -0.5_R16P
      u(  1) = (3._R16P/(x(2)-x(1)))*((v(2)-v(1))/(x(2)-x(1))-bc1)
    else
      cd2(1) = 0._R16P
      u(  1) = 0._R16P
    endif

    do i=2,N-1
      sig    = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p      = sig*cd2(i-1)+2._R16P
      cd2(i) = (sig-1._R16P)/p
      u(  i) = (6._R16P*((v(i+1)-v(i))/(x(i+1)-x(i))-(v(i)-v(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    ! last node boundary condition
    if (present(bcn)) then
      qn   = 0.5_R16P
      u(N) = (3._R16P/(x(N)-x(N-1)))*(bcn-(v(N)-v(N-1))/(x(N)-x(N-1)))
    else
      qn   = 0._R16P
      u(N) = 0._R16P
    endif
    ! coefficients of second derivative
    cd2(N) = (u(N)-qn*u(N-1))/(qn*cd2(N-1)+1._R16P)
    do i=N-1,1,-1
      cd2(i) = cd2(i)*cd2(i+1)+u(i)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction compute_cd2

    pure function cspline(xint) result(vint)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(R16P), intent(IN):: xint !< Independent variable abscissa of interpolating point.
    real(R16P)::             vint !< Dependent variable interpolated at x point.
    integer(I4P)::           iLo  !< Index of nearest-Low point to xint.
    integer(I4P)::           iHi  !< Index of nearest-High point to xint.
    real(R16P)::             h    !< Distance between the two nearest points to xint.
    real(R16P)::             a    !< Relative distance of nearest-High point to xint.
    real(R16P)::             b    !< Relative distance of nearest-Low point to xint.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    iLo = maxloc(x,mask=(x<xint),dim=1)
    iHi = minloc(x,mask=(x>xint),dim=1)
    h = x(iHi)-x(iLo)
    a = (x(iHi) - xint  )/h
    b = (xint   - x(iLo))/h
    vint = a*v(iLo) + b*v(iHi) + ((a*a*a - a)*cd2(iLo) + (b*b*b - b)*cd2(iHi))*(h*h)/6._R16P
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction cspline
  endfunction spline3_R16

  !> @brief Function for computing cubic spline interpolation of an array of values (R8P).
  pure function spline3_R8(bc1,bcn,x,v,xi) result (vi)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN), optional:: bc1            !< Boundary conditions of second derivative at first node.
  real(R8P), intent(IN), optional:: bcn            !< Boundary conditions of second derivative at last node.
  real(R8P), intent(IN)::           x(1:)          !< Independent variable abscissa of original data.
  real(R8P), intent(IN)::           v(1:)          !< Dependent variable (to be interpolated) of original data.
  real(R8P), intent(IN)::           xi(1:)         !< Independent variable abscissa of interpolating points.
  real(R8P)::                       vi(1:size(xi)) !< Dependent variable interpolated at xi points.
  real(R8P)::                       cd2(1:size(x)) !< Coefficients of second derivative.
  integer(I4P)::                    N              !< Number of points of original data.
  integer(I4P)::                    i              !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  N = size(x)
  cd2 = compute_cd2(N=N)
  do i=1,size(xi)
    vi(i) = cspline(xi(i))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure function compute_cd2(N) result(cd2)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I4P), intent(IN):: N               !< Number of points.
    real(R8P)::                cd2(1:N)        !< Coefficients of second derivative.
    real(R8P)::                sig,p,qn,u(1:N) !< Temporary variables.
    integer(I4P)::             i               !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ! first node boundary condition
    if (present(bc1)) then
      cd2(1) = -0.5_R8P
      u(  1) = (3._R8P/(x(2)-x(1)))*((v(2)-v(1))/(x(2)-x(1))-bc1)
    else
      cd2(1) = 0._R8P
      u(  1) = 0._R8P
    endif

    do i=2,N-1
      sig    = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p      = sig*cd2(i-1)+2._R8P
      cd2(i) = (sig-1._R8P)/p
      u(  i) = (6._R8P*((v(i+1)-v(i))/(x(i+1)-x(i))-(v(i)-v(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    ! last node boundary condition
    if (present(bcn)) then
      qn   = 0.5_R8P
      u(N) = (3._R8P/(x(N)-x(N-1)))*(bcn-(v(N)-v(N-1))/(x(N)-x(N-1)))
    else
      qn   = 0._R8P
      u(N) = 0._R8P
    endif
    ! coefficients of second derivative
    cd2(N) = (u(N)-qn*u(N-1))/(qn*cd2(N-1)+1._R8P)
    do i=N-1,1,-1
      cd2(i) = cd2(i)*cd2(i+1)+u(i)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction compute_cd2

    pure function cspline(xint) result(vint)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(R8P), intent(IN):: xint !< Independent variable abscissa of interpolating point.
    real(R8P)::             vint !< Dependent variable interpolated at x point.
    integer(I4P)::          iLo  !< Index of nearest-Low point to xint.
    integer(I4P)::          iHi  !< Index of nearest-High point to xint.
    real(R8P)::             h    !< Distance between the two nearest points to xint.
    real(R8P)::             a    !< Relative distance of nearest-High point to xint.
    real(R8P)::             b    !< Relative distance of nearest-Low point to xint.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    iLo = maxloc(x,mask=(x<xint),dim=1)
    iHi = minloc(x,mask=(x>xint),dim=1)
    h = x(iHi)-x(iLo)
    a = (x(iHi) - xint  )/h
    b = (xint   - x(iLo))/h
    vint = a*v(iLo) + b*v(iHi) + ((a*a*a - a)*cd2(iLo) + (b*b*b - b)*cd2(iHi))*(h*h)/6._R8P
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction cspline
  endfunction spline3_R8

  !> @brief Function for computing cubic spline interpolation of an array of values (R4P).
  pure function spline3_R4(bc1,bcn,x,v,xi) result (vi)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN), optional:: bc1            !< Boundary conditions of second derivative at first node.
  real(R4P), intent(IN), optional:: bcn            !< Boundary conditions of second derivative at last node.
  real(R4P), intent(IN)::           x(1:)          !< Independent variable abscissa of original data.
  real(R4P), intent(IN)::           v(1:)          !< Dependent variable (to be interpolated) of original data.
  real(R4P), intent(IN)::           xi(1:)         !< Independent variable abscissa of interpolating points.
  real(R4P)::                       vi(1:size(xi)) !< Dependent variable interpolated at xi points.
  real(R4P)::                       cd2(1:size(x)) !< Coefficients of second derivative.
  integer(I4P)::                    N              !< Number of points of original data.
  integer(I4P)::                    i              !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  N = size(x)
  cd2 = compute_cd2(N=N)
  do i=1,size(xi)
    vi(i) = cspline(xi(i))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  contains
    pure function compute_cd2(N) result(cd2)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    integer(I4P), intent(IN):: N               !< Number of points.
    real(R4P)::                cd2(1:N)        !< Coefficients of second derivative.
    real(R4P)::                sig,p,qn,u(1:N) !< Temporary variables.
    integer(I4P)::             i               !< Counter.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    ! first node boundary condition
    if (present(bc1)) then
      cd2(1) = -0.5_R4P
      u(  1) = (3._R4P/(x(2)-x(1)))*((v(2)-v(1))/(x(2)-x(1))-bc1)
    else
      cd2(1) = 0._R4P
      u(  1) = 0._R4P
    endif

    do i=2,N-1
      sig    = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p      = sig*cd2(i-1)+2._R4P
      cd2(i) = (sig-1._R4P)/p
      u(  i) = (6._R4P*((v(i+1)-v(i))/(x(i+1)-x(i))-(v(i)-v(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    ! last node boundary condition
    if (present(bcn)) then
      qn   = 0.5_R4P
      u(N) = (3._R4P/(x(N)-x(N-1)))*(bcn-(v(N)-v(N-1))/(x(N)-x(N-1)))
    else
      qn   = 0._R4P
      u(N) = 0._R4P
    endif
    ! coefficients of second derivative
    cd2(N) = (u(N)-qn*u(N-1))/(qn*cd2(N-1)+1._R4P)
    do i=N-1,1,-1
      cd2(i) = cd2(i)*cd2(i+1)+u(i)
    enddo
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction compute_cd2

    pure function cspline(xint) result(vint)
    !-------------------------------------------------------------------------------------------------------------------------------
    implicit none
    real(R4P), intent(IN):: xint !< Independent variable abscissa of interpolating point.
    real(R4P)::             vint !< Dependent variable interpolated at x point.
    integer(I4P)::          iLo  !< Index of nearest-Low point to xint.
    integer(I4P)::          iHi  !< Index of nearest-High point to xint.
    real(R4P)::             h    !< Distance between the two nearest points to xint.
    real(R4P)::             a    !< Relative distance of nearest-High point to xint.
    real(R4P)::             b    !< Relative distance of nearest-Low point to xint.
    !-------------------------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------------------------
    iLo = maxloc(x,mask=(x<xint),dim=1)
    iHi = minloc(x,mask=(x>xint),dim=1)
    h = x(iHi)-x(iLo)
    a = (x(iHi) - xint  )/h
    b = (xint   - x(iLo))/h
    vint = a*v(iLo) + b*v(iHi) + ((a*a*a - a)*cd2(iLo) + (b*b*b - b)*cd2(iHi))*(h*h)/6._R4P
    return
    !-------------------------------------------------------------------------------------------------------------------------------
    endfunction cspline
  endfunction spline3_R4
  !> @}

  ! stretching functions
  pure function stretchinglside(Ns,s_uniform,stretch) result(s_stretch)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function stretches a uniform grid increasing the density of node to the left side of the grid. It is based on the
  !!hyperbolic tangent function. The stretched abscissa is computed from the uniform by the equation:
  !!\begin{equation}
  !!  s_n^{stretch}  = s_L  + \left( {s_R  - s_L } \right)\left[ {1 + \frac{{\tanh \left[ {strecth\left( {\frac{{n - 1}}{{Ns -
  !!  1}} - 1} \right)} \right]}}{{\tanh \left( {strecth} \right)}}} \right]\quad for\;n = 1,Ns
  !!\end{equation}
  !!\noindent where $stretch$ is the intensity of the stretching, $Ns$ is the number of the nodes grid and $s_L,s_R$ are the
  !!abscissa value at the left and right bounds.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: Ns              ! Number of nodes.
  real(R_P),    intent(IN):: s_uniform(1:Ns) ! Uniform grid.
  real(R_P),    intent(IN):: stretch         ! Intensity of stretch.
  real(R_P)::                s_stretch(1:Ns) ! Stretched grid.
  real(R_P)::                s_L             ! Left value of abscissa.
  real(R_P)::                s_R             ! Right value of abscissa.
  integer(I4P)::             n               ! Nodes counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  s_L     = s_uniform(1)  ! Calculating the left value of abscissa.
  s_R     = s_uniform(Ns) ! Calculating the right value of abscissa.
  do n=1,Ns
    s_stretch(n)=s_L+(s_R-s_L)*(1._R_P+tanh(stretch*((n-1._R_P)/(Ns-1._R_P)-1._R_P))/tanh(stretch)) ! Stretching the unirform grid.
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction stretchinglside

  pure function stretchingrside(Ns,s_uniform,stretch) result(s_stretch)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function stretches a uniform grid increasing the density of node to the right side of the grid. It is based on the
  !!hyperbolic tangent function. The stretched abscissa is computed from the uniform by the equation:
  !!\begin{equation}
  !!  s_n^{stretch}  = s_L  + \left( {s_R  - s_L } \right)\frac{{\tanh \left( {strecth\frac{{n - 1}}{{Ns - 1}}}
  !!  \right)}}{{\tanh \left( {strecth} \right)}}\quad for\;n = 1,Ns
  !!\end{equation}
  !!\noindent where $stretch$ is the intensity of the stretching, $Ns$ is the number of the nodes grid and $s_L,s_R$ are the
  !!abscissa value at the left and right bounds.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: Ns              ! Number of nodes.
  real(R_P),    intent(IN):: s_uniform(1:Ns) ! Uniform grid.
  real(R_P),    intent(IN):: stretch         ! Intensity of stretch.
  real(R_P)::                s_stretch(1:Ns) ! Stretched grid.
  real(R_P)::                s_L             ! Left value of abscissa.
  real(R_P)::                s_R             ! Right value of abscissa.
  integer(I4P)::             n               ! Nodes counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  s_L     = s_uniform(1)  ! Calculating the left value of abscissa.
  s_R     = s_uniform(Ns) ! Calculating the right value of abscissa.
  do n=1,Ns
    s_stretch(n)=s_L+(s_R-s_L)*tanh(stretch*(n-1._R_P)/(Ns-1._R_P))/tanh(stretch) ! Stretching the unirform grid.
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction stretchingrside

  pure function stretching2side(s_sym_in,Ns,s_uniform,Ds_L,Ds_R) result(s_stretch)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function stretching a uniform grid and allowing the imposition of the space steps at the first and at the
  !!last cells (2 sides). It is based on the Vinokur (1980) algorithm:
  !!\begin{equation}
  !!  s_n^{stretch}  = s_L  + \left( {s_R  - s_L } \right)\frac{{u_n }}{{\alpha  + \left( {1 - \alpha } \right)u_n }}\quad
  !!  for\;n = 1,Ns
  !!\end{equation}
  !!\noindent where $u_n$ and $\alpha$ are:
  !!\begin{equation}
  !!  \begin{array}{l}
  !!    u_n  = \frac{{\tanh \left[ {strecth\left( {\frac{{n - 1}}{{Ns - 1}} - s_{sym} } \right)} \right] +
  !!           \tanh \left( {strecth \cdot s_{sym} } \right)}}{{\tanh \left[ {strecth\left( {1 - s_{sym} } \right)} \right] +
  !!           \tanh \left( {strecth \cdot s_{sym} } \right)}} \\
  !!    \alpha  = \sqrt {\frac{{\Delta _R }}{{\Delta _L }}}  \\
  !!  \end{array}
  !!\end{equation}
  !!\noindent where $stretch$ is the intensity of the stretching, $Ns$ is the number of the nodes grid, $s_L,s_R$ are the
  !!abscissa value at the left and right bounds and $s_sym$ is the abscissa where the symmetry point is placed; $s_sym$ falls in
  !!the range $[0,1]$ where $0$ corrisponds to the left bound and $1$ to right one. Note that the
  !!stretch intensity must be calculated according to the space steps imposed at the bounds; the equation for the stretch
  !!intensity is:
  !!\begin{equation}
  !!  \sinh \left( {stretch} \right) - \frac{{stretch}}{{\sqrt {\Delta _L^{adm} \Delta _R^{adm} } }} = 0
  !!\end{equation}
  !!\noindent where:
  !!\begin{equation}
  !! \Delta _L^{adm}  = \frac{{\Delta _L }}{{\left( {s_R  - s_L } \right)}}Ns \quad
  !! \Delta _R^{adm}  = \frac{{\Delta _R }}{{\left( {s_R  - s_L } \right)}}Ns
  !!\end{equation}
  !!\noindent are the adimensional left and right space steps. This equation is non linear and must be solverd with an iterative
  !!method like Newton-Rapson.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),    intent(IN), optional:: s_sym_in              ! Optional input value of abscissa of symmetry.
  integer(I4P), intent(IN)::           Ns                    ! Number of nodes.
  real(R_P),    intent(IN)::           s_uniform(1:Ns)       ! Uniform grid.
  real(R_P),    intent(IN)::           Ds_L                  ! First space step (Left).
  real(R_P),    intent(IN)::           Ds_R                  ! Last  space step (Right).
  real(R_P)::                          s_sym                 ! Abscissa of symmetry.
  real(R_P)::                          s_stretch(1:Ns)       ! Stretched grid.
  real(R_P)::                          s_L                   ! Left value of abscissa.
  real(R_P)::                          s_R                   ! Right value of abscissa.
  real(R_P)::                          c_L                   ! Relative (adimensional) first space step (Left).
  real(R_P)::                          c_R                   ! Relative (adimensional) last  space step (Right).
  real(R_P)::                          stretch               ! Intensity of stretch.
  real(R_P)::                          stretch_n             ! Iterative value of stretch intensity.
  real(R_P), parameter::               toll = 10._R_P**(-10) ! Tollerance for computing the stretch intensity.
  real(R_P)::                          u_n                   ! Node stretch factor.
  integer(I4P)::                       n                     ! Nodes counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(s_sym_in)) then
    s_sym = s_sym_in ! Imposing passed value of symmetry abscissa.
  else
    s_sym = 0.5_R_P  ! Symmetry abscissa is at 50% as default.
  endif
  s_L       = s_uniform(1)                                             ! Compute left value of abscissa.
  s_R       = s_uniform(Ns)                                            ! Compute right value of abscissa.
  c_L       = Ds_L/(s_R-s_L)*Ns                                        ! Compute left adm space step.
  c_R       = Ds_R/(s_R-s_L)*Ns                                        ! Compute right adm space step.
  stretch_n = 1.1_R_P                                                  ! Initialize iterative stretch intensity value.
  stretch   = stretch_n-(sinh(stretch_n)-stretch_n/(sqrt(c_L*c_R)))/ &
                        (cosh(stretch_n)-1._R_P/(sqrt(c_L*c_R)))       ! Update stretch intensity.
  do while(abs(stretch-stretch_n)>toll*abs(stretch))                     ! Evalute iterative convergence.
    stretch_n = stretch                                                  ! Update iterative stretch intensity value.
    stretch   = stretch_n-(sinh(stretch_n)-stretch_n/(sqrt(c_L*c_R)))/ &
                          (cosh(stretch_n)-1._R_P/(sqrt(c_L*c_R)))       ! Update stretch intensity.
  enddo
  do n=1,Ns
    u_n          = (tanh(stretch*((n-1._R_P)/(Ns-1._R_P)-s_sym)) + tanh(stretch*s_sym))/ &
                   (tanh(stretch*(1._R_P-s_sym))                 + tanh(stretch*s_sym))    ! Computing u_n.
    s_stretch(n) = s_L+(s_R-s_L)*u_n/(sqrt(c_R/c_L) + (1._R_P-sqrt(c_R/c_L))*u_n)          ! Stretching the unirform grid.
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction stretching2side

  ! average
  pure function average_Vectorial1D_R16(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R16P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x(:) ! Array.
  real(R16P)::             mean ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R16P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial1D_R16

  pure function average_Vectorial1D_R8(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R8P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x(:) ! Array.
  real(R8P)::             mean ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial1D_R8

  pure function average_Vectorial1D_R4(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R4P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x(:) ! Array.
  real(R4P)::             mean ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial1D_R4

  pure function average_Vectorial2D_R16(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R16P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x(:,:) ! Array.
  real(R16P)::             mean   ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R16P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial2D_R16

  pure function average_Vectorial2D_R8(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R8P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x(:,:) ! Array.
  real(R8P)::             mean   ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial2D_R8

  pure function average_Vectorial2D_R4(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R4P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x(:,:) ! Array.
  real(R4P)::             mean   ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial2D_R4

  pure function average_Vectorial3D_R16(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R16P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x(:,:,:) ! Array.
  real(R16P)::             mean     ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R16p)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial3D_R16

  pure function average_Vectorial3D_R8(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R8P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x(:,:,:) ! Array.
  real(R8P)::             mean     ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial3D_R8

  pure function average_Vectorial3D_R4(x) result(mean)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes arithmetic mean of array elements with parametric precision R4P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x(:,:,:) ! Array.
  real(R4P)::             mean     ! Arithmetic mean.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mean=sum(x)/real(size(x),R4P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial3D_R4

  ! digit
  elemental function digit_I8(n) result(digit)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the number of digits in decimal base of the input integer with parametric precision I8P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: n     ! Input integer.
  character(DI8P)::          str   ! Returned string containing input number plus padding zeros.
  integer(I4P)::             digit ! Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI8P) abs(n)         ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I8

  elemental function digit_I4(n) result(digit)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the number of digits in decimal base of the input integer with parametric precision I4P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: n     ! Input integer.
  character(DI4P)::          str   ! Returned string containing input number plus padding zeros.
  integer(I4P)::             digit ! Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI4P) abs(n)         ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I4

  elemental function digit_I2(n) result(digit)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the number of digits in decimal base of the input integer with parametric precision I2P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P), intent(IN):: n     ! Input integer.
  character(DI2P)::          str   ! Returned string containing input number plus padding zeros.
  integer(I4P)::             digit ! Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI2P) abs(n)         ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I2

  elemental function digit_I1(n) result(digit)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the number of digits in decimal base of the input integer with parametric precision I1P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: n     ! Input integer.
  character(DI1P)::          str   ! Returned string containing input number plus padding zeros.
  integer(I4P)::             digit ! Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI1P) abs(n)         ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I1

  ! div2
  elemental function div2_I1(n) result(d2)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for computing how many times an integer is divisible for 2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: n  ! Input integer.
  integer(I4P)::             d2 ! Number of times n is divisible for 2; is -1 for n=0.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n==0) then
    d2 = -1
    return
  endif
  do d2 = 0_I4P,int(bit_size(n),I4P)-1_I4P
    if (btest(n,d2)) return
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
   endfunction div2_I1

  elemental function div2_I2(n) result(d2)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for computing how many times an integer is divisible for 2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P), intent(IN):: n  ! Input integer.
  integer(I4P)::             d2 ! Number of times n is divisible for 2; is -1 for n=0.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n==0) then
    d2 = -1
    return
  endif
  do d2 = 0_I4P,int(bit_size(n),I4P)-1_I4P
    if (btest(n,d2)) return
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
   endfunction div2_I2

  elemental function div2_I4(n) result(d2)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for computing how many times an integer is divisible for 2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: n  ! Input integer.
  integer(I4P)::             d2 ! Number of times n is divisible for 2; is -1 for n=0.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n==0) then
    d2 = -1
    return
  endif
  do d2 = 0_I4P,int(bit_size(n),I4P)-1_I4P
    if (btest(n,d2)) return
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
   endfunction div2_I4

  elemental function div2_I8(n) result(d2)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for computing how many times an integer is divisible for 2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P), intent(IN):: n  ! Input integer.
  integer(I4P)::             d2 ! Number of times n is divisible for 2; is -1 for n=0.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (n==0) then
    d2 = -1
    return
  endif
  do d2 = 0_I4P,int(bit_size(n),I4P)-1_I4P
    if (btest(n,d2)) return
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
   endfunction div2_I8

  ! interpolate
  elemental function interpolate1_R16(x1,x2,V1,V2,x) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a linear interpolation between 2 values at 2 different points with parametric precision R16P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x1,x2 ! X coordinates of the 2 points.
  real(R16P), intent(IN):: V1,V2 ! Interpolating values.
  real(R16P), intent(IN):: x     ! X coordinate of interpolated value.
  real(R16P)::             V     ! Interpolated value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  V = V1+(V2-V1)*(x-x1)/(x2-x1)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate1_R16

  elemental function interpolate1_R8(x1,x2,V1,V2,x) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a linear interpolation between 2 values at 2 different points with parametric precision R8P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x1,x2 ! X coordinates of the 2 points.
  real(R8P), intent(IN):: V1,V2 ! Interpolating values.
  real(R8P), intent(IN):: x     ! X coordinate of interpolated value.
  real(R8P)::             v     ! Interpolated value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  V = V1+(V2-V1)*(x-x1)/(x2-x1)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate1_R8

  elemental function interpolate1_R4(x1,x2,V1,V2,x) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a linear interpolation between 2 values at 2 different points with parametric precision R4P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x1,x2 ! X coordinates of the 2 points.
  real(R4P), intent(IN):: V1,V2 ! Interpolating values.
  real(R4P), intent(IN):: x     ! X coordinate of interpolated value.
  real(R4P)::             V     ! Interpolated value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  V = V1+(V2-V1)*(x-x1)/(x2-x1)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate1_R4

  elemental function interpolate2_R16(x1,y1,x2,y2,V11,V12,V21,V22,x,y) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a bi-linear interpolation between 4 values evaluated on a 2D plane at 4 different points with parametric
  !!precision R16P.
  !!                            V12(x1,y2).----------.V22(x2,y2)
  !!                                      |          |
  !!                                      |          |
  !!                                      |          |
  !!                                      |          |
  !!                            V11(x1,y1).----------.V21(x2,y1)
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x1,x2 ! X coordinates of the 4 points.
  real(R16P), intent(IN):: y1,y2 ! Y coordinates of the 4 points.
  real(R16P), intent(IN):: V11   ! Interpolating value at the point x1,y1.
  real(R16P), intent(IN):: V12   ! Interpolating value at the point x1,y2.
  real(R16P), intent(IN):: V21   ! Interpolating value at the point x2,y1.
  real(R16P), intent(IN):: V22   ! Interpolating value at the point x2,y2.
  real(R16P), intent(IN):: x     ! X coordinates of interpolated value.
  real(R16P), intent(IN):: y     ! Y coordinates of interpolated value.
  real(R16P)::             V     ! Interpolated value.
  real(R16P)::             Vx1   ! Linear interpolated value along x at y1.
  real(R16P)::             Vx2   ! Linear interpolated value along x at y2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Vx1 = interpolate1_R16(x1=x1,x2=x2,V1=V11,V2=V21,x=x)
  Vx2 = interpolate1_R16(x1=x1,x2=x2,V1=V12,V2=V22,x=x)
  V   = interpolate1_R16(x1=y1,x2=y2,V1=Vx1,V2=Vx2,x=y)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate2_R16

  elemental function interpolate2_R8(x1,y1,x2,y2,V11,V12,V21,V22,x,y) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a bi-linear interpolation between 4 values evaluated on a 2D plane at 4 different points with parametric
  !!precision R8P.
  !!                            V12(x1,y2).----------.V22(x2,y2)
  !!                                      |          |
  !!                                      |          |
  !!                                      |          |
  !!                                      |          |
  !!                            V11(x1,y1).----------.V21(x2,y1)
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x1,x2 ! X coordinates of the 4 points.
  real(R8P), intent(IN):: y1,y2 ! Y coordinates of the 4 points.
  real(R8P), intent(IN):: V11   ! Interpolating value at the point x1,y1.
  real(R8P), intent(IN):: V12   ! Interpolating value at the point x1,y2.
  real(R8P), intent(IN):: V21   ! Interpolating value at the point x2,y1.
  real(R8P), intent(IN):: V22   ! Interpolating value at the point x2,y2.
  real(R8P), intent(IN):: x     ! X coordinates of interpolated value.
  real(R8P), intent(IN):: y     ! Y coordinates of interpolated value.
  real(R8P)::             V     ! Interpolated value.
  real(R8P)::             Vx1   ! Linear interpolated value along x at y1.
  real(R8P)::             Vx2   ! Linear interpolated value along x at y2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Vx1 = interpolate1_R8(x1=x1,x2=x2,V1=V11,V2=V21,x=x)
  Vx2 = interpolate1_R8(x1=x1,x2=x2,V1=V12,V2=V22,x=x)
  V   = interpolate1_R8(x1=y1,x2=y2,V1=Vx1,V2=Vx2,x=y)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate2_R8

  elemental function interpolate2_R4(x1,y1,x2,y2,V11,V12,V21,V22,x,y) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a bi-linear interpolation between 4 values evaluated on a 2D plane at 4 different points with parametric
  !!precision R4P.
  !!                            V12(x1,y2).----------.V22(x2,y2)
  !!                                      |          |
  !!                                      |          |
  !!                                      |          |
  !!                                      |          |
  !!                            V11(x1,y1).----------.V21(x2,y1)
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x1,x2 ! X coordinates of the 4 points.
  real(R4P), intent(IN):: y1,y2 ! Y coordinates of the 4 points.
  real(R4P), intent(IN):: V11   ! Interpolating value at the point x1,y1.
  real(R4P), intent(IN):: V12   ! Interpolating value at the point x1,y2.
  real(R4P), intent(IN):: V21   ! Interpolating value at the point x2,y1.
  real(R4P), intent(IN):: V22   ! Interpolating value at the point x2,y2.
  real(R4P), intent(IN):: x     ! X coordinates of interpolated value.
  real(R4P), intent(IN):: y     ! Y coordinates of interpolated value.
  real(R4P)::             V     ! Interpolated value.
  real(R4P)::             Vx1   ! Linear interpolated value along x at y1.
  real(R4P)::             Vx2   ! Linear interpolated value along x at y2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Vx1 = interpolate1_R4(x1=x1,x2=x2,V1=V11,V2=V21,x=x)
  Vx2 = interpolate1_R4(x1=x1,x2=x2,V1=V12,V2=V22,x=x)
  V   = interpolate1_R4(x1=y1,x2=y2,V1=Vx1,V2=Vx2,x=y)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate2_R4

  elemental function interpolate3_R16(x1,y1,z1,x2,y2,z2,V111,V121,V211,V221,V112,V122,V212,V222,x,y,z) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a tri-linear interpolation between 8 values evaluated on a 3D prism at 8 different points with parametric
  !!precision R16P.
  !!                           V122(x1,y2,z2).----------.V222(x2,y2,z2)
  !!                                        /|         /|
  !!                                       / |        / |
  !!                        V121(x1,y2,z1).----------.V221(x2,y2,z1)
  !!                                      |  |       |  |
  !!                           V112(x1,y1,z2).-------|--.V212(x2,y1,z2)
  !!                                      | /        | /
  !!                                      |/         |/
  !!                        V111(x1,y1,z1).----------.V211(x2,y1,z1)
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x1,x2 ! X coordinates of the 8 points.
  real(R16P), intent(IN):: y1,y2 ! Y coordinates of the 8 points.
  real(R16P), intent(IN):: z1,z2 ! Z coordinates of the 8 points.
  real(R16P), intent(IN):: V111  ! Interpolating value at the point x1,y1,z1.
  real(R16P), intent(IN):: V121  ! Interpolating value at the point x1,y2,z1.
  real(R16P), intent(IN):: V211  ! Interpolating value at the point x2,y1,z1.
  real(R16P), intent(IN):: V221  ! Interpolating value at the point x2,y2,z1.
  real(R16P), intent(IN):: V112  ! Interpolating value at the point x1,y1,z2.
  real(R16P), intent(IN):: V122  ! Interpolating value at the point x1,y2,z2.
  real(R16P), intent(IN):: V212  ! Interpolating value at the point x2,y1,z2.
  real(R16P), intent(IN):: V222  ! Interpolating value at the point x2,y2,z2.
  real(R16P), intent(IN):: x     ! X coordinates of interpolated value.
  real(R16P), intent(IN):: y     ! Y coordinates of interpolated value.
  real(R16P), intent(IN):: z     ! Y coordinates of interpolated value.
  real(R16P)::             V     ! Interpolated value.
  real(R16P)::             Vz1   ! Bi-linear interpolated value along xy at z1.
  real(R16P)::             Vz2   ! Bi-linear interpolated value along xy at z2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Vz1 = interpolate2_R16(x1=x1,y1=y1,x2=x2,y2=y2,V11=V111,V12=V121,V21=V211,V22=V221,x=x,y=y)
  Vz2 = interpolate2_R16(x1=x1,y1=y1,x2=x2,y2=y2,V11=V112,V12=V122,V21=V212,V22=V222,x=x,y=y)
  V   = interpolate1_R16(x1=z1,x2=z2,V1=Vz1,V2=Vz2,x=z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate3_R16

  elemental function interpolate3_R8(x1,y1,z1,x2,y2,z2,V111,V121,V211,V221,V112,V122,V212,V222,x,y,z) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a tri-linear interpolation between 8 values evaluated on a 3D prism at 8 different points with parametric
  !!precision R8P.
  !!                           V122(x1,y2,z2).----------.V222(x2,y2,z2)
  !!                                        /|         /|
  !!                                       / |        / |
  !!                        V121(x1,y2,z1).----------.V221(x2,y2,z1)
  !!                                      |  |       |  |
  !!                           V112(x1,y1,z2).-------|--.V212(x2,y1,z2)
  !!                                      | /        | /
  !!                                      |/         |/
  !!                        V111(x1,y1,z1).----------.V211(x2,y1,z1)
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x1,x2 ! X coordinates of the 8 points.
  real(R8P), intent(IN):: y1,y2 ! Y coordinates of the 8 points.
  real(R8P), intent(IN):: z1,z2 ! Z coordinates of the 8 points.
  real(R8P), intent(IN):: V111  ! Interpolating value at the point x1,y1,z1.
  real(R8P), intent(IN):: V121  ! Interpolating value at the point x1,y2,z1.
  real(R8P), intent(IN):: V211  ! Interpolating value at the point x2,y1,z1.
  real(R8P), intent(IN):: V221  ! Interpolating value at the point x2,y2,z1.
  real(R8P), intent(IN):: V112  ! Interpolating value at the point x1,y1,z2.
  real(R8P), intent(IN):: V122  ! Interpolating value at the point x1,y2,z2.
  real(R8P), intent(IN):: V212  ! Interpolating value at the point x2,y1,z2.
  real(R8P), intent(IN):: V222  ! Interpolating value at the point x2,y2,z2.
  real(R8P), intent(IN):: x     ! X coordinates of interpolated value.
  real(R8P), intent(IN):: y     ! Y coordinates of interpolated value.
  real(R8P), intent(IN):: z     ! Y coordinates of interpolated value.
  real(R8P)::             V     ! Interpolated value.
  real(R8P)::             Vz1   ! Bi-linear interpolated value along xy at z1.
  real(R8P)::             Vz2   ! Bi-linear interpolated value along xy at z2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Vz1 = interpolate2_R8(x1=x1,y1=y1,x2=x2,y2=y2,V11=V111,V12=V121,V21=V211,V22=V221,x=x,y=y)
  Vz2 = interpolate2_R8(x1=x1,y1=y1,x2=x2,y2=y2,V11=V112,V12=V122,V21=V212,V22=V222,x=x,y=y)
  V   = interpolate1_R8(x1=z1,x2=z2,V1=Vz1,V2=Vz2,x=z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate3_R8

  elemental function interpolate3_R4(x1,y1,z1,x2,y2,z2,V111,V121,V211,V221,V112,V122,V212,V222,x,y,z) result(V)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function makes a tri-linear interpolation between 8 values evaluated on a 3D prism at 8 different points with parametric
  !!precision R4P.
  !!                           V122(x1,y2,z2).----------.V222(x2,y2,z2)
  !!                                        /|         /|
  !!                                       / |        / |
  !!                        V121(x1,y2,z1).----------.V221(x2,y2,z1)
  !!                                      |  |       |  |
  !!                           V112(x1,y1,z2).-------|--.V212(x2,y1,z2)
  !!                                      | /        | /
  !!                                      |/         |/
  !!                        V111(x1,y1,z1).----------.V211(x2,y1,z1)
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x1,x2 ! X coordinates of the 8 points.
  real(R4P), intent(IN):: y1,y2 ! Y coordinates of the 8 points.
  real(R4P), intent(IN):: z1,z2 ! Z coordinates of the 8 points.
  real(R4P), intent(IN):: V111  ! Interpolating value at the point x1,y1,z1.
  real(R4P), intent(IN):: V121  ! Interpolating value at the point x1,y2,z1.
  real(R4P), intent(IN):: V211  ! Interpolating value at the point x2,y1,z1.
  real(R4P), intent(IN):: V221  ! Interpolating value at the point x2,y2,z1.
  real(R4P), intent(IN):: V112  ! Interpolating value at the point x1,y1,z2.
  real(R4P), intent(IN):: V122  ! Interpolating value at the point x1,y2,z2.
  real(R4P), intent(IN):: V212  ! Interpolating value at the point x2,y1,z2.
  real(R4P), intent(IN):: V222  ! Interpolating value at the point x2,y2,z2.
  real(R4P), intent(IN):: x     ! X coordinates of interpolated value.
  real(R4P), intent(IN):: y     ! Y coordinates of interpolated value.
  real(R4P), intent(IN):: z     ! Y coordinates of interpolated value.
  real(R4P)::             V     ! Interpolated value.
  real(R4P)::             Vz1   ! Bi-linear interpolated value along xy at z1.
  real(R4P)::             Vz2   ! Bi-linear interpolated value along xy at z2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Vz1 = interpolate2_R4(x1=x1,y1=y1,x2=x2,y2=y2,V11=V111,V12=V121,V21=V211,V22=V221,x=x,y=y)
  Vz2 = interpolate2_R4(x1=x1,y1=y1,x2=x2,y2=y2,V11=V112,V12=V122,V21=V212,V22=V222,x=x,y=y)
  V   = interpolate1_R4(x1=z1,x2=z2,V1=Vz1,V2=Vz2,x=z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction interpolate3_R4

  ! degree
  elemental function degree_Scalar_R16(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R16P (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x_rad ! Radiant value.
  real(R16P)::             x_deg ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R16P/pi_R16
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Scalar_R16

  elemental function degree_Scalar_R8(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R8P (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x_rad ! Radiant value.
  real(R8P)::             x_deg ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R8P/pi_R8
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Scalar_R8

  elemental function degree_Scalar_R4(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R4P (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x_rad ! Radiant value.
  real(R4P)::             x_deg ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R4P/pi_R4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Scalar_R4

  pure function degree_Vectorial1D_R16(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R16P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x_rad(:)                    ! Radiant value.
  real(R16P)::             x_deg(1:size(x_rad,dim=1))  ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R16P/pi_R16
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial1D_R16

  pure function degree_Vectorial1D_R8(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R8P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x_rad(:)                    ! Radiant value.
  real(R8P)::             x_deg(1:size(x_rad,dim=1))  ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R8P/pi_R8
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial1D_R8

  pure function degree_Vectorial1D_R4(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R4P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x_rad(:)                   ! Radiant value.
  real(R4P)::             x_deg(1:size(x_rad,dim=1)) ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R4P/pi_R4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial1D_R4

  pure function degree_Vectorial2D_R16(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R16P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x_rad(:,:)                                      ! Radiant value.
  real(R16P)::             x_deg(1:size(x_rad,dim=1),1:size(x_rad,dim=2))  ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R16P/pi_R16
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial2D_R16

  pure function degree_Vectorial2D_R8(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R8P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x_rad(:,:)                                      ! Radiant value.
  real(R8P)::             x_deg(1:size(x_rad,dim=1),1:size(x_rad,dim=2))  ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R8P/pi_R8
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial2D_R8

  pure function degree_Vectorial2D_R4(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R4P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x_rad(:,:)                                     ! Radiant value.
  real(R4P)::             x_deg(1:size(x_rad,dim=1),1:size(x_rad,dim=2)) ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R4P/pi_R4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial2D_R4

  pure function degree_Vectorial3D_R16(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R16P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x_rad(:,:,:)                                                        ! Radiant value.
  real(R16P)::             x_deg(1:size(x_rad,dim=1),1:size(x_rad,dim=2),1:size(x_rad,dim=3))  ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R16P/pi_R16
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial3D_R16

  pure function degree_Vectorial3D_R8(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R8P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x_rad(:,:,:)                                                        ! Radiant value.
  real(R8P)::             x_deg(1:size(x_rad,dim=1),1:size(x_rad,dim=2),1:size(x_rad,dim=3))  ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R8P/pi_R8
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial3D_R8

  pure function degree_Vectorial3D_R4(x_rad) result(x_deg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting radiants to degrees with parametric precision R4P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x_rad(:,:,:)                                                       ! Radiant value.
  real(R4P)::             x_deg(1:size(x_rad,dim=1),1:size(x_rad,dim=2),1:size(x_rad,dim=3)) ! Degree value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_deg=x_rad*180._R4P/pi_R4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction degree_Vectorial3D_R4

  ! radiant
  elemental function radiant_Scalar_R16(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R16P (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x_deg ! Degree value.
  real(R16P)::             x_rad ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R16/180._R16P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Scalar_R16

  elemental function radiant_Scalar_R8(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R8P (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x_deg ! Degree value.
  real(R8P)::             x_rad ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R8/180._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Scalar_R8

  elemental function radiant_Scalar_R4(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R4P (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x_deg ! Degree value.
  real(R4P)::             x_rad ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R4/180._R4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Scalar_R4

  pure function radiant_Vectorial1D_R16(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R16P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x_deg(:)                   ! Degree value.
  real(R16P)::             x_rad(1:size(x_deg,dim=1)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R16/180._R16P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial1D_R16

  pure function radiant_Vectorial1D_R8(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R8P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x_deg(:)                   ! Degree value.
  real(R8P)::             x_rad(1:size(x_deg,dim=1)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R8/180._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial1D_R8

  pure function radiant_Vectorial1D_R4(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R4P (vectorial 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x_deg(:)                   ! Degree value.
  real(R4P)::             x_rad(1:size(x_deg,dim=1)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R4/180._R4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial1D_R4

  pure function radiant_Vectorial2D_R16(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R16P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x_deg(:,:)                                     ! Degree value.
  real(R16P)::             x_rad(1:size(x_deg,dim=1),1:size(x_deg,dim=2)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R16/180._R16P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial2D_R16

  pure function radiant_Vectorial2D_R8(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R8P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x_deg(:,:)                                     ! Degree value.
  real(R8P)::             x_rad(1:size(x_deg,dim=1),1:size(x_deg,dim=2)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R8/180._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial2D_R8

  pure function radiant_Vectorial2D_R4(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R4P (vectorial 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x_deg(:,:)                                     ! Degree value.
  real(R4P)::             x_rad(1:size(x_deg,dim=1),1:size(x_deg,dim=2)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R4/180._R4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial2D_R4

  pure function radiant_Vectorial3D_R16(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R16P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: x_deg(:,:,:)                                                       ! Degree value.
  real(R16P)::             x_rad(1:size(x_deg,dim=1),1:size(x_deg,dim=2),1:size(x_deg,dim=3)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R16/180._R16P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial3D_R16

  pure function radiant_Vectorial3D_R8(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R8P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: x_deg(:,:,:)                                                       ! Degree value.
  real(R8P)::             x_rad(1:size(x_deg,dim=1),1:size(x_deg,dim=2),1:size(x_deg,dim=3)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R8/180._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial3D_R8

  pure function radiant_Vectorial3D_R4(x_deg) result(x_rad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting degrees to radiants with the parametric precision R4P (vectorial 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: x_deg(:,:,:)                                                       ! Degree value.
  real(R4P)::             x_rad(1:size(x_deg,dim=1),1:size(x_deg,dim=2),1:size(x_deg,dim=3)) ! Radiant value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  x_rad=x_deg*pi_R4/180._R4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction radiant_Vectorial3D_R4

  ! first derivative
  elemental function delta1_2o_R16(Sp1,Sm1,Ds) result(D1)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the first derivative by 2 order central difference with parametric precision R16P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: Sp1 ! value in cell s+1
  real(R16P), intent(IN):: Sm1 ! value in cell s-1
  real(R16P), intent(IN):: Ds  ! step in cell s
  real(R16P)::             D1  ! value of numerical derivative
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  D1 = (Sp1-Sm1)/(2._R16P*Ds)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction delta1_2o_R16

  elemental function delta1_2o_R8(Sp1,Sm1,Ds) result(D1)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the first derivative by 2 order central difference with parametric precision R8P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: Sp1 ! value in cell s+1
  real(R8P), intent(IN):: Sm1 ! value in cell s-1
  real(R8P), intent(IN):: Ds  ! step in cell s
  real(R8P)::             D1  ! value of numerical derivative
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  D1 = (Sp1-Sm1)/(2._R8P*Ds)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction delta1_2o_R8

  elemental function delta1_2o_R4(Sp1,Sm1,Ds) result(D1)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the first derivative by 2 order central difference with parametric precision R4P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: Sp1 ! value in cell s+1
  real(R4P), intent(IN):: Sm1 ! value in cell s-1
  real(R4P), intent(IN):: Ds  ! step in cell s
  real(R4P)::             D1  ! value of numerical derivative
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  D1 = (Sp1-Sm1)/(2._R4P*Ds)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction delta1_2o_R4

  ! second derivative
  elemental function delta2_2o_R16(Sp1,S,Sm1,Ds) result(D2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the second derivative by 2 order central difference with parametric precision R16P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: Sp1 ! Value of node s+1.
  real(R16P), intent(IN):: S   ! Value of node s.
  real(R16P), intent(IN):: Sm1 ! Value of node s-1.
  real(R16P), intent(IN):: Ds  ! Step value.
  real(R16P)::             D2  ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  D2 = (Sp1-2._R16P*S+Sm1)/(Ds*Ds)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction delta2_2o_R16

  elemental function delta2_2o_R8(Sp1,S,Sm1,Ds) result(D2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the second derivative by 2 order central difference with parametric precision R8P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: Sp1 ! Value of node s+1.
  real(R8P), intent(IN):: S   ! Value of node s.
  real(R8P), intent(IN):: Sm1 ! Value of node s-1.
  real(R8P), intent(IN):: Ds  ! Step value.
  real(R8P)::             D2  ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  D2 = (Sp1-2._R8P*S+Sm1)/(Ds*Ds)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction delta2_2o_R8

  elemental function delta2_2o_R4(Sp1,S,Sm1,Ds) result(D2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the second derivative by 2 order central difference with parametric precision R4P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: Sp1 ! Value of node s+1.
  real(R4P), intent(IN):: S   ! Value of node s.
  real(R4P), intent(IN):: Sm1 ! Value of node s-1.
  real(R4P), intent(IN):: Ds  ! Step value.
  real(R4P)::             D2  ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  D2 = (Sp1-2._R4P*S+Sm1)/(Ds*Ds)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction delta2_2o_R4

  ! absolute value of gradient
  elemental function abs_grad_R16(Vip1,Vim1,Di,Vjp1,Vjm1,Dj,Vkp1,Vkm1,Dk) result(abs_grad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the absolute value of the gradient of a scalar by 2 order central difference with
  !!parametric precision R16P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: Vip1     ! Value of node i+1.
  real(R16P), intent(IN):: Vim1     ! Value of node i-1.
  real(R16P), intent(IN):: Di       ! I step value.
  real(R16P), intent(IN):: Vjp1     ! Value of node j+1.
  real(R16P), intent(IN):: Vjm1     ! Value of node j-1.
  real(R16P), intent(IN):: Dj       ! J step value.
  real(R16P), intent(IN):: Vkp1     ! Value of node k+1.
  real(R16P), intent(IN):: Vkm1     ! Value of node k-1.
  real(R16P), intent(IN):: Dk       ! K step value.
  real(R16P)::             abs_grad ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  abs_grad = sqrt(delta1_2o_R16(Sp1 = Vip1,Sm1 = Vim1,Ds = Di)*delta1_2o_R16(Sp1 = Vip1,Sm1 = Vim1,Ds = Di) + &
                  delta1_2o_R16(Sp1 = Vjp1,Sm1 = Vjm1,Ds = Dj)*delta1_2o_R16(Sp1 = Vjp1,Sm1 = Vjm1,Ds = Dj) + &
                  delta1_2o_R16(Sp1 = Vkp1,Sm1 = Vkm1,Ds = Dk)*delta1_2o_R16(Sp1 = Vkp1,Sm1 = Vkm1,Ds = Dk))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction abs_grad_R16

  elemental function abs_grad_R8(Vip1,Vim1,Di,Vjp1,Vjm1,Dj,Vkp1,Vkm1,Dk) result(abs_grad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the absolute value of the gradient of a scalar by 2 order central difference with
  !!parametric precision R8P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: Vip1     ! Value of node i+1.
  real(R8P), intent(IN):: Vim1     ! Value of node i-1.
  real(R8P), intent(IN):: Di       ! I step value.
  real(R8P), intent(IN):: Vjp1     ! Value of node j+1.
  real(R8P), intent(IN):: Vjm1     ! Value of node j-1.
  real(R8P), intent(IN):: Dj       ! J step value.
  real(R8P), intent(IN):: Vkp1     ! Value of node k+1.
  real(R8P), intent(IN):: Vkm1     ! Value of node k-1.
  real(R8P), intent(IN):: Dk       ! K step value.
  real(R8P)::             abs_grad ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  abs_grad = sqrt(delta1_2o_R8(Sp1 = Vip1,Sm1 = Vim1,Ds = Di)*delta1_2o_R8(Sp1 = Vip1,Sm1 = Vim1,Ds = Di) + &
                  delta1_2o_R8(Sp1 = Vjp1,Sm1 = Vjm1,Ds = Dj)*delta1_2o_R8(Sp1 = Vjp1,Sm1 = Vjm1,Ds = Dj) + &
                  delta1_2o_R8(Sp1 = Vkp1,Sm1 = Vkm1,Ds = Dk)*delta1_2o_R8(Sp1 = Vkp1,Sm1 = Vkm1,Ds = Dk))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction abs_grad_R8

  elemental function abs_grad_R4(Vip1,Vim1,Di,Vjp1,Vjm1,Dj,Vkp1,Vkm1,Dk) result(abs_grad)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the absolute value of the gradient of a scalar by 2 order central difference with
  !!parametric precision R4P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: Vip1     ! Value of node i+1.
  real(R4P), intent(IN):: Vim1     ! Value of node i-1.
  real(R4P), intent(IN):: Di       ! I step value.
  real(R4P), intent(IN):: Vjp1     ! Value of node j+1.
  real(R4P), intent(IN):: Vjm1     ! Value of node j-1.
  real(R4P), intent(IN):: Dj       ! J step value.
  real(R4P), intent(IN):: Vkp1     ! Value of node k+1.
  real(R4P), intent(IN):: Vkm1     ! Value of node k-1.
  real(R4P), intent(IN):: Dk       ! K step value.
  real(R4P)::             abs_grad ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  abs_grad = sqrt(delta1_2o_R4(Sp1 = Vip1,Sm1 = Vim1,Ds = Di)*delta1_2o_R4(Sp1 = Vip1,Sm1 = Vim1,Ds = Di) + &
                  delta1_2o_R4(Sp1 = Vjp1,Sm1 = Vjm1,Ds = Dj)*delta1_2o_R4(Sp1 = Vjp1,Sm1 = Vjm1,Ds = Dj) + &
                  delta1_2o_R4(Sp1 = Vkp1,Sm1 = Vkm1,Ds = Dk)*delta1_2o_R4(Sp1 = Vkp1,Sm1 = Vkm1,Ds = Dk))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction abs_grad_R4

  ! laplacian operator
  elemental function laplacian_R16(Vip1,Vi,Vim1,Di,Vjp1,Vj,Vjm1,Dj,Vkp1,Vk,Vkm1,Dk) result(laplace)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the laplacian operator of a scalar by 2 order central difference with parametric
  !!precision R16P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P), intent(IN):: Vip1    ! Value of node i+1.
  real(R16P), intent(IN):: Vi      ! Value of node i.
  real(R16P), intent(IN):: Vim1    ! Value of node i-1.
  real(R16P), intent(IN):: Di      ! I step value.
  real(R16P), intent(IN):: Vjp1    ! Value of node j+1.
  real(R16P), intent(IN):: Vj      ! Value of node j.
  real(R16P), intent(IN):: Vjm1    ! Value of node j-1.
  real(R16P), intent(IN):: Dj      ! J step value.
  real(R16P), intent(IN):: Vkp1    ! Value of node k+1.
  real(R16P), intent(IN):: Vk      ! Value of node k.
  real(R16P), intent(IN):: Vkm1    ! Value of node k-1.
  real(R16P), intent(IN):: Dk      ! K step value.
  real(R16P)::             laplace ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  laplace = delta2_2o_R16(Sp1 = Vip1,S = Vi,Sm1 = Vim1,Ds = Di) + &
            delta2_2o_R16(Sp1 = Vjp1,S = Vj,Sm1 = Vjm1,Ds = Dj) + &
            delta2_2o_R16(Sp1 = Vkp1,S = Vk,Sm1 = Vkm1,Ds = Dk)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction laplacian_R16

  elemental function laplacian_R8(Vip1,Vi,Vim1,Di,Vjp1,Vj,Vjm1,Dj,Vkp1,Vk,Vkm1,Dk) result(laplace)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the laplacian operator of a scalar by 2 order central difference with parametric
  !!precision R8P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P), intent(IN):: Vip1    ! Value of node i+1.
  real(R8P), intent(IN):: Vi      ! Value of node i.
  real(R8P), intent(IN):: Vim1    ! Value of node i-1.
  real(R8P), intent(IN):: Di      ! I step value.
  real(R8P), intent(IN):: Vjp1    ! Value of node j+1.
  real(R8P), intent(IN):: Vj      ! Value of node j.
  real(R8P), intent(IN):: Vjm1    ! Value of node j-1.
  real(R8P), intent(IN):: Dj      ! J step value.
  real(R8P), intent(IN):: Vkp1    ! Value of node k+1.
  real(R8P), intent(IN):: Vk      ! Value of node k.
  real(R8P), intent(IN):: Vkm1    ! Value of node k-1.
  real(R8P), intent(IN):: Dk      ! K step value.
  real(R8P)::             laplace ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  laplace = delta2_2o_R8(Sp1 = Vip1,S = Vi,Sm1 = Vim1,Ds = Di) + &
            delta2_2o_R8(Sp1 = Vjp1,S = Vj,Sm1 = Vjm1,Ds = Dj) + &
            delta2_2o_R8(Sp1 = Vkp1,S = Vk,Sm1 = Vkm1,Ds = Dk)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction laplacian_R8

  elemental function laplacian_R4(Vip1,Vi,Vim1,Di,Vjp1,Vj,Vjm1,Dj,Vkp1,Vk,Vkm1,Dk) result(laplace)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes numerically the laplacian operator of a scalar by 2 order central difference with parametric
  !!precision R4P.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P), intent(IN):: Vip1    ! Value of node i+1.
  real(R4P), intent(IN):: Vi      ! Value of node i.
  real(R4P), intent(IN):: Vim1    ! Value of node i-1.
  real(R4P), intent(IN):: Di      ! I step value.
  real(R4P), intent(IN):: Vjp1    ! Value of node j+1.
  real(R4P), intent(IN):: Vj      ! Value of node j.
  real(R4P), intent(IN):: Vjm1    ! Value of node j-1.
  real(R4P), intent(IN):: Dj      ! J step value.
  real(R4P), intent(IN):: Vkp1    ! Value of node k+1.
  real(R4P), intent(IN):: Vk      ! Value of node k.
  real(R4P), intent(IN):: Vkm1    ! Value of node k-1.
  real(R4P), intent(IN):: Dk      ! K step value.
  real(R4P)::             laplace ! Value of numerical derivative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  laplace = delta2_2o_R4(Sp1 = Vip1,S = Vi,Sm1 = Vim1,Ds = Di) + &
            delta2_2o_R4(Sp1 = Vjp1,S = Vj,Sm1 = Vjm1,Ds = Dj) + &
            delta2_2o_R4(Sp1 = Vkp1,S = Vk,Sm1 = Vkm1,Ds = Dk)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction laplacian_R4
endmodule Lib_Math
