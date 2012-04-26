module Data_Type_Conservative
!-----------------------------------------------------------------------------------------------------------------------------------
!!This module contains Data_Type_Conservative, a derived type that handles conservative fluidynamic variables.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                              ! Integers and reals precision definition.
USE Data_Type_Vector,                       & ! Definition of Type_Vector.
                      init_vector => init,  & ! Function for initializing Type_Vector.
                      set_vector  => set,   & ! Function for setting Type_Vector.
                      get_vector  => get      ! Function for getting Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: Type_Conservative
public:: init,set,free
!public:: write,read
public:: assignment (=)
public:: operator (*)
public:: operator (/)
public:: operator (+)
public:: operator (-)
public:: cons2array,array2cons
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Type_Conservative definition:
type:: Type_Conservative
  sequence
  real(R_P), allocatable:: rs(:)       ! Density of single species [1:Ns].
  type(Type_Vector)::      rv          ! Momentum vector.
  real(R_P)::              re = 0._R_P ! Product of density for total internal energy (sum(r)*E).
endtype Type_Conservative
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Init overloading.
interface init
  module procedure Init_Scalar,  & ! scalar
                   Init_Array1D, & ! array1D
                   Init_Array2D, & ! array2D
                   Init_Array3D, & ! array3D
                   Init_Array4D    ! array4D
endinterface
!!Set overloading.
interface set
  module procedure Set_Scalar,  & ! scalar
                   Set_Array1D, & ! array1D
                   Set_Array2D, & ! array2D
                   Set_Array3D, & ! array3D
                   Set_Array4D    ! array4D
endinterface
!!Free overloading.
interface free
  module procedure Free_Scalar,  & ! scalar
                   Free_Array1D, & ! array1D
                   Free_Array2D, & ! array2D
                   Free_Array3D, & ! array3D
                   Free_Array4D    ! array4D
endinterface
!!Assigment (=) overloading.
interface assignment (=)
  module procedure assign_cons,     &
#ifdef r16p
                   assign_ScalR16P, &
#endif
                   assign_ScalR8P,  &
                   assign_ScalR4P,  &
                   assign_ScalI8P,  &
                   assign_ScalI4P,  &
                   assign_ScalI2P,  &
                   assign_ScalI1P
end interface
!!Multiplication (*) overloading.
interface operator (*)
  module procedure cons_mul_cons,     &
#ifdef r16p
                   ScalR16P_mul_cons, &
#endif
                   ScalR8P_mul_cons,  &
                   ScalR4P_mul_cons,  &
                   ScalI8P_mul_cons,  &
                   ScalI4P_mul_cons,  &
                   ScalI2P_mul_cons,  &
                   ScalI1P_mul_cons,  &
#ifdef r16p
                   cons_mul_ScalR16P, &
#endif
                   cons_mul_ScalR8P,  &
                   cons_mul_ScalR4P,  &
                   cons_mul_ScalI8P,  &
                   cons_mul_ScalI4P,  &
                   cons_mul_ScalI2P,  &
                   cons_mul_ScalI1P
endinterface
!!Division (/) overloading.
interface operator (/)
  module procedure cons_div_cons,     &
#ifdef r16p
                   cons_div_ScalR16P, &
#endif
                   cons_div_ScalR8P,  &
                   cons_div_ScalR4P,  &
                   cons_div_ScalI8P,  &
                   cons_div_ScalI4P,  &
                   cons_div_ScalI2P,  &
                   cons_div_ScalI1P
endinterface
!!Sum (+) overloading.
interface operator (+)
  module procedure positive_cons,     &
                   cons_sum_cons,     &
#ifdef r16p
                   ScalR16P_sum_cons, &
#endif
                   ScalR8P_sum_cons,  &
                   ScalR4P_sum_cons,  &
                   ScalI8P_sum_cons,  &
                   ScalI4P_sum_cons,  &
                   ScalI2P_sum_cons,  &
                   ScalI1P_sum_cons,  &
#ifdef r16p
                   cons_sum_ScalR16P, &
#endif
                   cons_sum_ScalR8P,  &
                   cons_sum_ScalR4P,  &
                   cons_sum_ScalI8P,  &
                   cons_sum_ScalI4P,  &
                   cons_sum_ScalI2P,  &
                   cons_sum_ScalI1P
endinterface
!!Subtraction (-) overloading.
interface operator (-)
  module procedure negative_cons,     &
                   cons_sub_cons,     &
#ifdef r16p
                   ScalR16P_sub_cons, &
#endif
                   ScalR8P_sub_cons,  &
                   ScalR4P_sub_cons,  &
                   ScalI8P_sub_cons,  &
                   ScalI4P_sub_cons,  &
                   ScalI2P_sub_cons,  &
                   ScalI1P_sub_cons,  &
#ifdef r16p
                   cons_sub_ScalR16P, &
#endif
                   cons_sub_ScalR8P,  &
                   cons_sub_ScalR4P,  &
                   cons_sub_ScalI8P,  &
                   cons_sub_ScalI4P,  &
                   cons_sub_ScalI2P,  &
                   cons_sub_ScalI1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! init
  pure subroutine Init_Scalar(rs,rv,re,Ns,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Conservative (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  integer(I_P),            intent(IN)::           Ns
  type(Type_Conservative), intent(INOUT)::        cons
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) deallocate(cons%rs) ; allocate(cons%rs(1:Ns)) ; cons%rs = 0._R_P
  if (present(rs)) cons%rs = rs
  if (present(rv)) cons%rv = rv
  if (present(re)) cons%re = re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Scalar

  pure subroutine Init_Array1D(rs,rv,re,Ns,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Conservative (array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  integer(I_P),            intent(IN)::           Ns
  type(Type_Conservative), intent(INOUT)::        cons(:)
  integer(I4P)::                                  i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=lbound(cons,dim=1),ubound(cons,dim=1)
    if (allocated(cons(i)%rs)) deallocate(cons(i)%rs) ; allocate(cons(i)%rs(1:Ns)) ; cons(i)%rs=0._R_P
    if (present(rs)) cons(i)%rs = rs
    if (present(rv)) cons(i)%rv = rv
    if (present(re)) cons(i)%re = re
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Array1D

  pure subroutine Init_Array2D(rs,rv,re,Ns,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Conservative (array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  integer(I_P),            intent(IN)::           Ns
  type(Type_Conservative), intent(INOUT)::        cons(:,:)
  integer(I4P)::                                  i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do j=lbound(cons,dim=2),ubound(cons,dim=2)
    do i=lbound(cons,dim=1),ubound(cons,dim=1)
      if (allocated(cons(i,j)%rs)) deallocate(cons(i,j)%rs) ; allocate(cons(i,j)%rs(1:Ns)) ; cons(i,j)%rs=0._R_P
      if (present(rs)) cons(i,j)%rs = rs
      if (present(rv)) cons(i,j)%rv = rv
      if (present(re)) cons(i,j)%re = re
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Array2D

  pure subroutine Init_Array3D(rs,rv,re,Ns,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Conservative (array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  integer(I_P),            intent(IN)::           Ns
  type(Type_Conservative), intent(INOUT)::        cons(:,:,:)
  integer(I4P)::                                  i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do k=lbound(cons,dim=3),ubound(cons,dim=3)
    do j=lbound(cons,dim=2),ubound(cons,dim=2)
      do i=lbound(cons,dim=1),ubound(cons,dim=1)
        if (allocated(cons(i,j,k)%rs)) deallocate(cons(i,j,k)%rs) ; allocate(cons(i,j,k)%rs(1:Ns)) ; cons(i,j,k)%rs=0._R_P
        if (present(rs)) cons(i,j,k)%rs = rs
        if (present(rv)) cons(i,j,k)%rv = rv
        if (present(re)) cons(i,j,k)%re = re
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Array3D

  pure subroutine Init_Array4D(rs,rv,re,Ns,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Conservative (array 4D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  integer(I_P),            intent(IN)::           Ns
  type(Type_Conservative), intent(INOUT)::        cons(:,:,:,:)
  integer(I4P)::                                  i,j,k,p
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do p=lbound(cons,dim=4),ubound(cons,dim=4)
    do k=lbound(cons,dim=3),ubound(cons,dim=3)
      do j=lbound(cons,dim=2),ubound(cons,dim=2)
        do i=lbound(cons,dim=1),ubound(cons,dim=1)
          if (allocated(cons(i,j,k,p)%rs)) deallocate(cons(i,j,k,p)%rs)
          allocate(cons(i,j,k,p)%rs(1:Ns)) ; cons(i,j,k,p)%rs=0._R_P
          if (present(rs)) cons(i,j,k,p)%rs = rs
          if (present(rv)) cons(i,j,k,p)%rv = rv
          if (present(re)) cons(i,j,k,p)%re = re
        enddo
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Array4D

  ! set
  pure subroutine Set_Scalar(rs,rv,re,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Conservative (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  type(Type_Conservative), intent(INOUT)::        cons
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(rs)) then
    if (allocated(cons%rs)) deallocate(cons%rs) ; allocate(cons%rs(1:size(rs,dim=1))) ; cons%rs = rs
  endif
  if (present(rv)) cons%rv = rv
  if (present(re)) cons%re = re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Scalar

  pure subroutine Set_Array1D(rs,rv,re,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Conservative (array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  type(Type_Conservative), intent(INOUT)::        cons(:)
  integer(I4P)::                                  i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=lbound(cons,dim=1),ubound(cons,dim=1)
    if (present(rs)) then
      if (allocated(cons(i)%rs)) deallocate(cons(i)%rs) ; allocate(cons(i)%rs(1:size(rs,dim=1))) ; cons(i)%rs = rs
    endif
    if (present(rv)) cons(i)%rv = rv
    if (present(re)) cons(i)%re = re
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Array1D

  pure subroutine Set_Array2D(rs,rv,re,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Conservative (array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  type(Type_Conservative), intent(INOUT)::        cons(:,:)
  integer(I4P)::                                  i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do j=lbound(cons,dim=2),ubound(cons,dim=2)
    do i=lbound(cons,dim=1),ubound(cons,dim=1)
      if (present(rs)) then
        if (allocated(cons(i,j)%rs)) deallocate(cons(i,j)%rs) ; allocate(cons(i,j)%rs(1:size(rs,dim=1))) ; cons(i,j)%rs = rs
      endif
      if (present(rv)) cons(i,j)%rv = rv
      if (present(re)) cons(i,j)%re = re
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Array2D

  pure subroutine Set_Array3D(rs,rv,re,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Conservative (array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  type(Type_Conservative), intent(INOUT)::        cons(:,:,:)
  integer(I4P)::                                  i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do k=lbound(cons,dim=3),ubound(cons,dim=3)
    do j=lbound(cons,dim=2),ubound(cons,dim=2)
      do i=lbound(cons,dim=1),ubound(cons,dim=1)
        if (present(rs)) then
          if (allocated(cons(i,j,k)%rs)) deallocate(cons(i,j,k)%rs) ; allocate(cons(i,j,k)%rs(1:size(rs,dim=1))) ; cons(i,j,k)%rs=rs
        endif
        if (present(rv)) cons(i,j,k)%rv = rv
        if (present(re)) cons(i,j,k)%re = re
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Array3D

  pure subroutine Set_Array4D(rs,rv,re,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Conservative (array 4D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN), optional:: rs(:)
  type(Type_Vector),       intent(IN), optional:: rv
  real(R_P),               intent(IN), optional:: re
  type(Type_Conservative), intent(INOUT)::        cons(:,:,:,:)
  integer(I4P)::                                  i,j,k,p
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do p=lbound(cons,dim=4),ubound(cons,dim=4)
    do k=lbound(cons,dim=3),ubound(cons,dim=3)
      do j=lbound(cons,dim=2),ubound(cons,dim=2)
        do i=lbound(cons,dim=1),ubound(cons,dim=1)
          if (present(rs)) then
            if (allocated(cons(i,j,k,p)%rs)) deallocate(cons(i,j,k,p)%rs)
            allocate(cons(i,j,k,p)%rs(1:size(rs,dim=1))) ; cons(i,j,k,p)%rs=rs
          endif
          if (present(rv)) cons(i,j,k,p)%rv = rv
          if (present(re)) cons(i,j,k,p)%re = re
        enddo
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Array4D

  ! free dynamic memory
  function Free_Scalar(cons) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Conservative (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I4P)::                           err  ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  if (allocated(cons%rs)) deallocate(cons%rs,stat=err)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Scalar

  function Free_Array1D(cons) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Conservative (array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons(:)
  integer(I4P)::                           err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                           i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do i=lbound(cons,dim=1),ubound(cons,dim=1)
    if (allocated(cons(i)%rs)) deallocate(cons(i)%rs,stat=err)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array1D

  function Free_Array2D(cons) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Conservative (array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons(:,:)
  integer(I4P)::                           err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                           i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do j=lbound(cons,dim=2),ubound(cons,dim=2)
    do i=lbound(cons,dim=1),ubound(cons,dim=1)
      if (allocated(cons(i,J)%rs)) deallocate(cons(i,J)%rs,stat=err)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array2D

  function Free_Array3D(cons) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Conservative (array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons(:,:,:)
  integer(I4P)::                           err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                           i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do k=lbound(cons,dim=3),ubound(cons,dim=3)
    do j=lbound(cons,dim=2),ubound(cons,dim=2)
      do i=lbound(cons,dim=1),ubound(cons,dim=1)
        if (allocated(cons(i,j,k)%rs)) deallocate(cons(i,j,k)%rs,stat=err)
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array3D

  function Free_Array4D(cons) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Conservative (array 4D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons(:,:,:,:)
  integer(I4P)::                           err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                           i,j,k,p
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do p=lbound(cons,dim=4),ubound(cons,dim=4)
    do k=lbound(cons,dim=3),ubound(cons,dim=3)
      do j=lbound(cons,dim=2),ubound(cons,dim=2)
        do i=lbound(cons,dim=1),ubound(cons,dim=1)
          if (allocated(cons(i,j,k,p)%rs)) deallocate(cons(i,j,k,p)%rs,stat=err)
        enddo
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array4D

  ! Assignment (=)
  elemental subroutine assign_cons(cons1,cons2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between two cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons1
  type(Type_Conservative), intent(IN)::    cons2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons1%rs).and.allocated(cons2%rs)) cons1%rs = cons2%rs
                                                   cons1%rv = cons2%rv
                                                   cons1%re = cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_cons

#ifdef r16p
  elemental subroutine assign_ScalR16P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R16P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  real(R16P),              intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) cons%rs = real(scal,R_P)
                          cons%rv = real(scal,R_P)
                          cons%re = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR16P
#endif

  elemental subroutine assign_ScalR8P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R8P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  real(R8P),               intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) cons%rs = real(scal,R_P)
                          cons%rv = real(scal,R_P)
                          cons%re = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR8P

  elemental subroutine assign_ScalR4P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R4P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  real(R4P),               intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) cons%rs = real(scal,R_P)
                          cons%rv = real(scal,R_P)
                          cons%re = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR4P

  elemental subroutine assign_ScalI8P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I8P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I8P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) cons%rs = real(scal,R_P)
                          cons%rv = real(scal,R_P)
                          cons%re = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI8P

  elemental subroutine assign_ScalI4P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I4P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I4P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) cons%rs = real(scal,R_P)
                          cons%rv = real(scal,R_P)
                          cons%re = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI4P

  elemental subroutine assign_ScalI2P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I2P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I2P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) cons%rs = real(scal,R_P)
                          cons%rv = real(scal,R_P)
                          cons%re = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI2P

  elemental subroutine assign_ScalI1P(cons,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I1P) and cons.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(INOUT):: cons
  integer(I1P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(cons%rs)) cons%rs = real(scal,R_P)
                          cons%rv = real(scal,R_P)
                          cons%re = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI1P

  ! Multiplication (*)
  elemental function cons_mul_cons(cons1,cons2) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply (by components) conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::  cons1 ! First cons obj.
  type(Type_Conservative), intent(IN)::  cons2 ! Second cons obj.
  type(Type_Conservative)::              mul   ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons1%rs)))
  mul%rs = cons1%rs * cons2%rs
  mul%rv = cons1%rv * cons2%rv
  mul%re = cons1%re * cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_cons

#ifdef r16p
  elemental function ScalR16P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R16P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),              intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_mul_cons

  elemental function cons_mul_ScalR16P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R16P),              intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalR16P
#endif

  elemental function ScalR8P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_mul_cons

  elemental function cons_mul_ScalR8P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R8P),               intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalR8P

  elemental function ScalR4P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_mul_cons

  elemental function cons_mul_ScalR4P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R4P),               intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalR4P

  elemental function ScalI8P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_mul_cons

  elemental function cons_mul_ScalI8P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalI8P

  elemental function ScalI4P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_mul_cons

  elemental function cons_mul_ScalI4P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalI4P

  elemental function ScalI2P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I2P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_mul_cons

  elemental function cons_mul_ScalI2P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalI2P

  elemental function ScalI1P_mul_cons(scal,cons) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I1P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_mul_cons

  elemental function cons_mul_ScalI1P(cons,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply conservative object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%rs)) allocate(mul%rs(1:size(cons%rs)))
  mul%rs = real(scal,R_P) * cons%rs
  mul%rv = real(scal,R_P) * cons%rv
  mul%re = real(scal,R_P) * cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_mul_ScalI1P

  ! Division (/)
  elemental function cons_div_cons(cons1,cons2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide (by components) conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::  cons1 ! First cons obj.
  type(Type_Conservative), intent(IN)::  cons2 ! Second cons obj.
  type(Type_Conservative)::              div   ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%rs)) allocate(div%rs(1:size(cons1%rs)))
  div%rs = cons1%rs / cons2%rs
  div%rv = cons1%rv / cons2%rv
  div%re = cons1%re / cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_cons

#ifdef r16p
  elemental function cons_div_ScalR16P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R16P),              intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%rs)) allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R_P)
  div%rv = cons%rv / real(scal,R_P)
  div%re = cons%re / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalR16P
#endif

  elemental function cons_div_ScalR8P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R8P),               intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%rs)) allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R_P)
  div%rv = cons%rv / real(scal,R_P)
  div%re = cons%re / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalR8P

  elemental function cons_div_ScalR4P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R4P),               intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%rs)) allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R_P)
  div%rv = cons%rv / real(scal,R_P)
  div%re = cons%re / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalR4P

  elemental function cons_div_ScalI8P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%rs)) allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R_P)
  div%rv = cons%rv / real(scal,R_P)
  div%re = cons%re / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalI8P

  elemental function cons_div_ScalI4P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%rs)) allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R_P)
  div%rv = cons%rv / real(scal,R_P)
  div%re = cons%re / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalI4P

  elemental function cons_div_ScalI2P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%rs)) allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R_P)
  div%rv = cons%rv / real(scal,R_P)
  div%re = cons%re / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalI2P

  elemental function cons_div_ScalI1P(cons,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide conservative object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%rs)) allocate(div%rs(1:size(cons%rs)))
  div%rs = cons%rs / real(scal,R_P)
  div%rv = cons%rv / real(scal,R_P)
  div%re = cons%re / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_div_ScalI1P

  ! Sum (+)
  elemental function positive_cons(cons) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary + to a conservative objecy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(pos%rs)) allocate(pos%rs(1:size(cons%rs)))
  pos%rs =  + cons%rs
  pos%rv =  + cons%rv
  pos%re =  + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_cons

  elemental function cons_sum_cons(cons1,cons2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum (by components) conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::  cons1 ! First cons obj.
  type(Type_Conservative), intent(IN)::  cons2 ! Second cons obj.
  type(Type_Conservative)::              summ  ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons1%rs)))
  summ%rs = cons1%rs + cons2%rs
  summ%rv = cons1%rv + cons2%rv
  summ%re = cons1%re + cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_cons

#ifdef r16p
  elemental function ScalR16P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R16P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),              intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sum_cons

  elemental function cons_sum_ScalR16P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R16P),              intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalR16P
#endif

  elemental function ScalR8P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sum_cons

  elemental function cons_sum_ScalR8P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R8P),               intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalR8P

  elemental function ScalR4P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sum_cons

  elemental function cons_sum_ScalR4P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R4P),               intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalR4P

  elemental function ScalI8P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sum_cons

  elemental function cons_sum_ScalI8P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalI8P

  elemental function ScalI4P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sum_cons

  elemental function cons_sum_ScalI4P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalI4P

  elemental function ScalI2P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I2P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sum_cons

  elemental function cons_sum_ScalI2P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalI2P

  elemental function ScalI1P_sum_cons(scal,cons) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I1P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sum_cons

  elemental function cons_sum_ScalI1P(cons,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum conservative object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%rs)) allocate(summ%rs(1:size(cons%rs)))
  summ%rs = real(scal,R_P) + cons%rs
  summ%rv = real(scal,R_P) + cons%rv
  summ%re = real(scal,R_P) + cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sum_ScalI1P

  ! Subtraction (-)
    elemental function negative_cons(cons) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary - to a conservative objecy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(neg%rs)) allocate(neg%rs(1:size(cons%rs)))
  neg%rs =  - cons%rs
  neg%rv =  - cons%rv
  neg%re =  - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_cons

  elemental function cons_sub_cons(cons1,cons2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract (by components) conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::  cons1 ! First cons obj.
  type(Type_Conservative), intent(IN)::  cons2 ! Second cons obj.
  type(Type_Conservative)::              sub  ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons1%rs)))
  sub%rs = cons1%rs - cons2%rs
  sub%rv = cons1%rv - cons2%rv
  sub%re = cons1%re - cons2%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_cons

#ifdef r16p
  elemental function ScalR16P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R16P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),              intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R_P) - cons%rs
  sub%rv = real(scal,R_P) - cons%rv
  sub%re = real(scal,R_P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sub_cons

  elemental function cons_sub_ScalR16P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R16P),              intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R_P)
  sub%rv = cons%rv - real(scal,R_P)
  sub%re = cons%re - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalR16P
#endif

  elemental function ScalR8P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R_P) - cons%rs
  sub%rv = real(scal,R_P) - cons%rv
  sub%re = real(scal,R_P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sub_cons

  elemental function cons_sub_ScalR8P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R8P),               intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R_P)
  sub%rv = cons%rv - real(scal,R_P)
  sub%re = cons%re - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalR8P

  elemental function ScalR4P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),               intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R_P) - cons%rs
  sub%rv = real(scal,R_P) - cons%rv
  sub%re = real(scal,R_P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sub_cons

  elemental function cons_sub_ScalR4P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R4P),               intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R_P)
  sub%rv = cons%rv - real(scal,R_P)
  sub%re = cons%re - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalR4P

  elemental function ScalI8P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I8P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R_P) - cons%rs
  sub%rv = real(scal,R_P) - cons%rv
  sub%re = real(scal,R_P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sub_cons

  elemental function cons_sub_ScalI8P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I8P),            intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R_P)
  sub%rv = cons%rv - real(scal,R_P)
  sub%re = cons%re - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalI8P

  elemental function ScalI4P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I4P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R_P) - cons%rs
  sub%rv = real(scal,R_P) - cons%rv
  sub%re = real(scal,R_P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sub_cons

  elemental function cons_sub_ScalI4P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I4P),            intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R_P)
  sub%rv = cons%rv - real(scal,R_P)
  sub%re = cons%re - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalI4P

  elemental function ScalI2P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I2P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R_P) - cons%rs
  sub%rv = real(scal,R_P) - cons%rv
  sub%re = real(scal,R_P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sub_cons

  elemental function cons_sub_ScalI2P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I2P),            intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R_P)
  sub%rv = cons%rv - real(scal,R_P)
  sub%re = cons%re - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalI2P

  elemental function ScalI1P_sub_cons(scal,cons) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I1P) for conservative object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative), intent(IN):: cons
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = real(scal,R_P) - cons%rs
  sub%rv = real(scal,R_P) - cons%rv
  sub%re = real(scal,R_P) - cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sub_cons

  elemental function cons_sub_ScalI1P(cons,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract conservative object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  integer(I1P),            intent(IN):: scal
  type(Type_Conservative)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%rs)) allocate(sub%rs(1:size(cons%rs)))
  sub%rs = cons%rs - real(scal,R_P)
  sub%rv = cons%rv - real(scal,R_P)
  sub%re = cons%re - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons_sub_ScalI1P

  pure function cons2array(cons) result(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting derived type Type_Conservative to array.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN):: cons
  real(R_P)::                           array(1:size(cons%rs)+4)
  integer(I_P)::                        Ns
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(cons%rs)
  array(1:Ns) = cons%rs
  array(Ns+1) = cons%rv%x
  array(Ns+2) = cons%rv%y
  array(Ns+3) = cons%rv%z
  array(Ns+4) = cons%re
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction cons2array

  pure function array2cons(array) result(cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting array to derived type Type_Conservative.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN)::   array(:)
  type(Type_Conservative):: cons
  integer(I_P)::            Ns
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = ubound(array,dim=1)-4
  if (allocated(cons%rs)) deallocate(cons%rs) ; allocate(cons%rs(1:Ns))
  cons%rs   = array(1:Ns)
  cons%rv%x = array(Ns+1)
  cons%rv%y = array(Ns+2)
  cons%rv%z = array(Ns+3)
  cons%re   = array(Ns+4)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction array2cons
endmodule Data_Type_Conservative
