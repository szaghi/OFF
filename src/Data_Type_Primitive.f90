module Data_Type_Primitive
!-----------------------------------------------------------------------------------------------------------------------------------
!!This module contains Data_Type_Primitive, a derived type that handles primitive fluidynamic variables.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                ! Integers and reals precision definition.
USE Data_Type_Vector,                         & ! Definition of Type_Vector.
                      init_vector   => init,  & ! Function for initializing Type_Vector.
                      set_vector    => set,   & ! Function for setting Type_Vector.
                      get_vector    => get,   & ! Function for getting Type_Vector.
                      write_vector  => write, & ! Function for writing Type_Vector.
                      read_vector   => read     ! Function for reading Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: Type_Primitive
public:: init,set,free
public:: write,read
public:: assignment (=)
public:: operator (*)
public:: operator (/)
public:: operator (+)
public:: operator (-)
public:: prim2array,array2prim
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Type_Primivite definition:
type:: Type_Primitive
  sequence
  real(R_P), allocatable:: r(:)       ! Density of single species [1:Ns].
  type(Type_Vector)::      v          ! Velocity vector.
  real(R_P)::              p = 0._R_P ! Pressure.
  real(R_P)::              d = 0._R_P ! Density = sum(r(1:Ns)).
  real(R_P)::              g = 0._R_P ! Specific heats ratio cp/cv.
endtype Type_Primitive
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Init overloading.
interface init
  module procedure Init_Scalar,  & ! scalar
                   Init_Array1D, & ! array1D
                   Init_Array2D, & ! array2D
                   Init_Array3D    ! array3D
endinterface
!!Set overloading.
interface set
  module procedure Set_Scalar,  & ! scalar
                   Set_Array1D, & ! array1D
                   Set_Array2D, & ! array2D
                   Set_Array3D    ! array3D
endinterface
!!Free overloading.
interface free
  module procedure Free_Scalar,  & ! scalar
                   Free_Array1D, & ! array1D
                   Free_Array2D, & ! array2D
                   Free_Array3D    ! array3D
endinterface
!!Write overloading.
interface write
  module procedure Write_Bin_Scalar,    & ! binary scalar
                   Write_Ascii_Scalar,  & ! ascii scalar
                   Write_Bin_Array1D,   & ! binary array1D
                   Write_Ascii_Array1D, & ! ascii array1D
                   Write_Bin_Array2D,   & ! binary array2D
                   Write_Ascii_Array2D, & ! ascii array2D
                   Write_Bin_Array3D,   & ! binary array3D
                   Write_Ascii_Array3D    ! ascii array3D
endinterface
!!Read overloading.
interface read
  module procedure Read_Bin_Scalar,    & ! binary scalar
                   Read_Ascii_Scalar,  & ! ascii scalar
                   Read_Bin_Array1D,   & ! binary array1D
                   Read_Ascii_Array1D, & ! ascii array1D
                   Read_Bin_Array2D,   & ! binary array2D
                   Read_Ascii_Array2D, & ! ascii array2D
                   Read_Bin_Array3D,   & ! binary array3D
                   Read_Ascii_Array3D    ! ascii array3D
endinterface
!!Assigment (=) overloading.
interface assignment (=)
  module procedure assign_prim,     &
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
  module procedure prim_mul_prim,     &
#ifdef r16p
                   ScalR16P_mul_prim, &
#endif
                   ScalR8P_mul_prim,  &
                   ScalR4P_mul_prim,  &
                   ScalI8P_mul_prim,  &
                   ScalI4P_mul_prim,  &
                   ScalI2P_mul_prim,  &
                   ScalI1P_mul_prim,  &
#ifdef r16p
                   prim_mul_ScalR16P, &
#endif
                   prim_mul_ScalR8P,  &
                   prim_mul_ScalR4P,  &
                   prim_mul_ScalI8P,  &
                   prim_mul_ScalI4P,  &
                   prim_mul_ScalI2P,  &
                   prim_mul_ScalI1P
endinterface
!!Division (/) overloading.
interface operator (/)
  module procedure prim_div_prim,     &
#ifdef r16p
                   prim_div_ScalR16P, &
#endif
                   prim_div_ScalR8P,  &
                   prim_div_ScalR4P,  &
                   prim_div_ScalI8P,  &
                   prim_div_ScalI4P,  &
                   prim_div_ScalI2P,  &
                   prim_div_ScalI1P
endinterface
!!Sum (+) overloading.
interface operator (+)
  module procedure positive_prim,     &
                   prim_sum_prim,     &
#ifdef r16p
                   ScalR16P_sum_prim, &
#endif
                   ScalR8P_sum_prim,  &
                   ScalR4P_sum_prim,  &
                   ScalI8P_sum_prim,  &
                   ScalI4P_sum_prim,  &
                   ScalI2P_sum_prim,  &
                   ScalI1P_sum_prim,  &
#ifdef r16p
                   prim_sum_ScalR16P, &
#endif
                   prim_sum_ScalR8P,  &
                   prim_sum_ScalR4P,  &
                   prim_sum_ScalI8P,  &
                   prim_sum_ScalI4P,  &
                   prim_sum_ScalI2P,  &
                   prim_sum_ScalI1P
endinterface
!!Subtraction (-) overloading.
interface operator (-)
  module procedure negative_prim,     &
                   prim_sub_prim,     &
#ifdef r16p
                   ScalR16P_sub_prim, &
#endif
                   ScalR8P_sub_prim,  &
                   ScalR4P_sub_prim,  &
                   ScalI8P_sub_prim,  &
                   ScalI4P_sub_prim,  &
                   ScalI2P_sub_prim,  &
                   ScalI1P_sub_prim,  &
#ifdef r16p
                   prim_sub_ScalR16P, &
#endif
                   prim_sub_ScalR8P,  &
                   prim_sub_ScalR4P,  &
                   prim_sub_ScalI8P,  &
                   prim_sub_ScalI4P,  &
                   prim_sub_ScalI2P,  &
                   prim_sub_ScalI1P
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! init
  pure subroutine Init_Scalar(r,v,p,d,g,Ns,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Primitive (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),            intent(IN), optional:: r(:)
  type(Type_Vector),    intent(IN), optional:: v
  real(R_P),            intent(IN), optional:: p
  real(R_P),            intent(IN), optional:: d
  real(R_P),            intent(IN), optional:: g
  integer(I_P),         intent(IN)::           Ns
  type(Type_Primitive), intent(INOUT)::        prim
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r)) deallocate(prim%r) ; allocate(prim%r(1:Ns)) ; prim%r = 0._R_P
  if (present(r)) prim%r = r
  if (present(v)) prim%v = v
  if (present(p)) prim%p = p
  if (present(d)) prim%d = d
  if (present(g)) prim%g = g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Scalar

  pure subroutine Init_Array1D(r,v,p,d,g,Ns,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Primitive (array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),            intent(IN), optional:: r(:)
  type(Type_Vector),    intent(IN), optional:: v
  real(R_P),            intent(IN), optional:: p
  real(R_P),            intent(IN), optional:: d
  real(R_P),            intent(IN), optional:: g
  integer(I_P),         intent(IN)::           Ns
  type(Type_Primitive), intent(INOUT)::        prim(:)
  integer(I4P)::                               i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=lbound(prim,dim=1),ubound(prim,dim=1)
    if (allocated(prim(i)%r)) deallocate(prim(i)%r) ; allocate(prim(i)%r(1:Ns)) ; prim(i)%r = 0._R_P
    if (present(r)) prim(i)%r = r
    if (present(v)) prim(i)%v = v
    if (present(p)) prim(i)%p = p
    if (present(d)) prim(i)%d = d
    if (present(g)) prim(i)%g = g
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Array1D

  pure subroutine Init_Array2D(r,v,p,d,g,Ns,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Primitive (array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),            intent(IN), optional:: r(:)
  type(Type_Vector),    intent(IN), optional:: v
  real(R_P),            intent(IN), optional:: p
  real(R_P),            intent(IN), optional:: d
  real(R_P),            intent(IN), optional:: g
  integer(I_P),         intent(IN)::           Ns
  type(Type_Primitive), intent(INOUT)::        prim(:,:)
  integer(I4P)::                               i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do j=lbound(prim,dim=2),ubound(prim,dim=2)
    do i=lbound(prim,dim=1),ubound(prim,dim=1)
      if (allocated(prim(i,j)%r)) deallocate(prim(i,j)%r) ; allocate(prim(i,j)%r(1:Ns)) ; prim(i,j)%r = 0._R_P
      if (present(r)) prim(i,j)%r = r
      if (present(v)) prim(i,j)%v = v
      if (present(p)) prim(i,j)%p = p
      if (present(d)) prim(i,j)%d = d
      if (present(g)) prim(i,j)%g = g
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Array2D

  pure subroutine Init_Array3D(r,v,p,d,g,Ns,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for initializing Type_Primitive (array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),            intent(IN), optional:: r(:)
  type(Type_Vector),    intent(IN), optional:: v
  real(R_P),            intent(IN), optional:: p
  real(R_P),            intent(IN), optional:: d
  real(R_P),            intent(IN), optional:: g
  integer(I_P),         intent(IN)::           Ns
  type(Type_Primitive), intent(INOUT)::        prim(:,:,:)
  integer(I4P)::                               i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do k=lbound(prim,dim=3),ubound(prim,dim=3)
    do j=lbound(prim,dim=2),ubound(prim,dim=2)
      do i=lbound(prim,dim=1),ubound(prim,dim=1)
        if (allocated(prim(i,j,k)%r)) deallocate(prim(i,j,k)%r) ; allocate(prim(i,j,k)%r(1:Ns)) ; prim(i,j,k)%r = 0._R_P
        if (present(r)) prim(i,j,k)%r = r
        if (present(v)) prim(i,j,k)%v = v
        if (present(p)) prim(i,j,k)%p = p
        if (present(d)) prim(i,j,k)%d = d
        if (present(g)) prim(i,j,k)%g = g
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Init_Array3D

  ! set
  pure subroutine Set_Scalar(r,v,p,d,g,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Primitive (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),            intent(IN), optional:: r(:)
  type(Type_Vector),    intent(IN), optional:: v
  real(R_P),            intent(IN), optional:: p
  real(R_P),            intent(IN), optional:: d
  real(R_P),            intent(IN), optional:: g
  type(Type_Primitive), intent(INOUT)::        prim
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(r)) then
    if (allocated(prim%r)) deallocate(prim%r) ; allocate(prim%r(1:size(r,dim=1))) ; prim%r = r
  endif
  if (present(v)) prim%v = v
  if (present(p)) prim%p = p
  if (present(d)) prim%d = d
  if (present(g)) prim%g = g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Scalar

  pure subroutine Set_Array1D(r,v,p,d,g,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Primitive (array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),            intent(IN), optional:: r(:)
  type(Type_Vector),    intent(IN), optional:: v
  real(R_P),            intent(IN), optional:: p
  real(R_P),            intent(IN), optional:: d
  real(R_P),            intent(IN), optional:: g
  type(Type_Primitive), intent(INOUT)::        prim(:)
  integer(I4P)::                               i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do i=lbound(prim,dim=1),ubound(prim,dim=1)
    if (present(r)) then
      if (allocated(prim(i)%r)) deallocate(prim(i)%r) ; allocate(prim(i)%r(1:size(r,dim=1))) ; prim(i)%r = r
    endif
    if (present(v)) prim(i)%v = v
    if (present(p)) prim(i)%p = p
    if (present(d)) prim(i)%d = d
    if (present(g)) prim(i)%g = g
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Array1D

  pure subroutine Set_Array2D(r,v,p,d,g,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Primitive (array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),            intent(IN), optional:: r(:)
  type(Type_Vector),    intent(IN), optional:: v
  real(R_P),            intent(IN), optional:: p
  real(R_P),            intent(IN), optional:: d
  real(R_P),            intent(IN), optional:: g
  type(Type_Primitive), intent(INOUT)::        prim(:,:)
  integer(I4P)::                               i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do j=lbound(prim,dim=2),ubound(prim,dim=2)
    do i=lbound(prim,dim=1),ubound(prim,dim=1)
      if (present(r)) then
        if (allocated(prim(i,j)%r)) deallocate(prim(i,j)%r) ; allocate(prim(i,j)%r(1:size(r,dim=1))) ; prim(i,j)%r = r
      endif
      if (present(v)) prim(i,j)%v = v
      if (present(p)) prim(i,j)%p = p
      if (present(d)) prim(i,j)%d = d
      if (present(g)) prim(i,j)%g = g
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Array2D

  pure subroutine Set_Array3D(r,v,p,d,g,prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for setting Type_Primitive (array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),            intent(IN), optional:: r(:)
  type(Type_Vector),    intent(IN), optional:: v
  real(R_P),            intent(IN), optional:: p
  real(R_P),            intent(IN), optional:: d
  real(R_P),            intent(IN), optional:: g
  type(Type_Primitive), intent(INOUT)::        prim(:,:,:)
  integer(I4P)::                               i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do k=lbound(prim,dim=3),ubound(prim,dim=3)
    do j=lbound(prim,dim=2),ubound(prim,dim=2)
      do i=lbound(prim,dim=1),ubound(prim,dim=1)
        if (present(r)) then
          if (allocated(prim(i,j,k)%r)) deallocate(prim(i,j,k)%r) ; allocate(prim(i,j,k)%r(1:size(r,dim=1))) ; prim(i,j,k)%r = r
        endif
        if (present(v)) prim(i,j,k)%v = v
        if (present(p)) prim(i,j,k)%p = p
        if (present(d)) prim(i,j,k)%d = d
        if (present(g)) prim(i,j,k)%g = g
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine Set_Array3D

  ! free dynamic memory
  function Free_Scalar(prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Primitive (scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I4P)::                        err  ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  if (allocated(prim%r)) deallocate(prim%r,stat=err)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Scalar

  function Free_Array1D(prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Primitive (array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim(:)
  integer(I4P)::                        err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do i=lbound(prim,dim=1),ubound(prim,dim=1)
    if (allocated(prim(i)%r)) deallocate(prim(i)%r,stat=err)
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array1D

  function Free_Array2D(prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Primitive (array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim(:,:)
  integer(I4P)::                        err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do j=lbound(prim,dim=2),ubound(prim,dim=2)
    do i=lbound(prim,dim=1),ubound(prim,dim=1)
      if (allocated(prim(i,J)%r)) deallocate(prim(i,J)%r,stat=err)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array2D

  function Free_Array3D(prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Primitive (array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim(:,:,:)
  integer(I4P)::                        err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do k=lbound(prim,dim=3),ubound(prim,dim=3)
    do j=lbound(prim,dim=2),ubound(prim,dim=2)
      do i=lbound(prim,dim=1),ubound(prim,dim=1)
        if (allocated(prim(i,j,k)%r)) deallocate(prim(i,j,k)%r,stat=err)
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array3D

  function Free_Array4D(prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for freeing the memory of Type_Primitive (array 4D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim(:,:,:,:)
  integer(I4P)::                        err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        i,j,k,p
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 0_I4P
  do p=lbound(prim,dim=4),ubound(prim,dim=4)
    do k=lbound(prim,dim=3),ubound(prim,dim=3)
      do j=lbound(prim,dim=2),ubound(prim,dim=2)
        do i=lbound(prim,dim=1),ubound(prim,dim=1)
          if (allocated(prim(i,j,k,p)%r)) deallocate(prim(i,j,k,p)%r,stat=err)
        enddo
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Free_Array4D

  ! write
  function Write_Bin_Scalar(unit,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Primitive (binary, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: unit ! logic unit
  type(Type_Primitive), intent(IN):: prim
  integer(I4P)::                     err  ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)prim%r(:)
  err = write_vector(unit,prim%v)
  write(unit,iostat=err)prim%p
  write(unit,iostat=err)prim%d
  write(unit,iostat=err)prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Scalar

  function Write_Ascii_Scalar(unit,format,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Primitive (ascii, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: unit   ! logic unit
  character(*),         intent(IN):: format ! format specifier
  type(Type_Primitive), intent(IN):: prim
  integer(I4P)::                     err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)prim%r(:)
    err = write_vector(unit,format,prim%v)
    write(unit,*,iostat=err)prim%p
    write(unit,*,iostat=err)prim%d
    write(unit,*,iostat=err)prim%g
  case default
    write(unit,adjustl(trim(format)),iostat=err)prim%r(:)
    err = write_vector(unit,format,prim%v)
    write(unit,adjustl(trim(format)),iostat=err)prim%p
    write(unit,adjustl(trim(format)),iostat=err)prim%d
    write(unit,adjustl(trim(format)),iostat=err)prim%g
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Scalar

  function Write_Bin_Array1D(unit,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Primitive (binary, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: unit    ! logic unit
  type(Type_Primitive), intent(IN):: prim(:)
  integer(I4P)::                     err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                     Ni,i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  write(unit,iostat=err)(prim(i)%r(:),i=1,Ni)
  err = write_vector(unit,prim%v)
  write(unit,iostat=err)prim%p
  write(unit,iostat=err)prim%d
  write(unit,iostat=err)prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array1D

  function Write_Ascii_Array1D(unit,format,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Primitive (ascii, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: unit    ! logic unit
  character(*),         intent(IN):: format
  type(Type_Primitive), intent(IN):: prim(:)
  integer(I4P)::                     err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                     Ni,i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)(prim(i)%r(:),i=1,Ni)
    err = write_vector(unit,format,prim%v)
    write(unit,*,iostat=err)prim%p
    write(unit,*,iostat=err)prim%d
    write(unit,*,iostat=err)prim%g
  case default
    write(unit,adjustl(trim(format)),iostat=err)(prim(i)%r(:),i=1,Ni)
    err = write_vector(unit,format,prim%v)
    write(unit,adjustl(trim(format)),iostat=err)prim%p
    write(unit,adjustl(trim(format)),iostat=err)prim%d
    write(unit,adjustl(trim(format)),iostat=err)prim%g
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array1D

  function Write_Bin_Array2D(unit,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Primitive (binary, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: unit     ! logic unit
  type(Type_Primitive), intent(IN):: prim(:,:)
  integer(I4P)::                     err      ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                     Ni,Nj,i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  Nj = size(prim,dim=2)
  write(unit,iostat=err)((prim(i,j)%r(:),i=1,Ni),j=1,Nj)
  err = write_vector(unit,prim%v)
  write(unit,iostat=err)prim%p
  write(unit,iostat=err)prim%d
  write(unit,iostat=err)prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array2D

  function Write_Ascii_Array2D(unit,format,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Primitive (ascii, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: unit      ! logic unit
  character(*),         intent(IN):: format
  type(Type_Primitive), intent(IN):: prim(:,:)
  integer(I4P)::                     err       ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                     Ni,Nj,i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  Nj = size(prim,dim=2)
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)((prim(i,j)%r(:),i=1,Ni),j=1,Nj)
    err = write_vector(unit,format,prim%v)
    write(unit,*,iostat=err)prim%p
    write(unit,*,iostat=err)prim%d
    write(unit,*,iostat=err)prim%g
  case default
    write(unit,adjustl(trim(format)),iostat=err)((prim(i,j)%r(:),i=1,Ni),j=1,Nj)
    err = write_vector(unit,format,prim%v)
    write(unit,adjustl(trim(format)),iostat=err)prim%p
    write(unit,adjustl(trim(format)),iostat=err)prim%d
    write(unit,adjustl(trim(format)),iostat=err)prim%g
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array2D

  function Write_Bin_Array3D(unit,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Primitive (binary, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: unit        ! logic unit
  type(Type_Primitive), intent(IN):: prim(:,:,:)
  integer(I4P)::                     err         ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                     Ni,Nj,Nk,i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  Nj = size(prim,dim=2)
  Nk = size(prim,dim=3)
  write(unit,iostat=err)(((prim(i,j,k)%r(:),i=1,Ni),j=1,Nj),k=1,Nk)
  err = write_vector(unit,prim%v)
  write(unit,iostat=err)prim%p
  write(unit,iostat=err)prim%d
  write(unit,iostat=err)prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array3D

  function Write_Ascii_Array3D(unit,format,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Primitive (ascii, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: unit        ! logic unit
  character(*),         intent(IN):: format
  type(Type_Primitive), intent(IN):: prim(:,:,:)
  integer(I4P)::                     err         ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                     Ni,Nj,Nk,i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  Nj = size(prim,dim=2)
  Nk = size(prim,dim=3)
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)(((prim(i,j,k)%r(:),i=1,Ni),j=1,Nj),k=1,Nk)
    err = write_vector(unit,format,prim%v)
    write(unit,*,iostat=err)prim%p
    write(unit,*,iostat=err)prim%d
    write(unit,*,iostat=err)prim%g
  case default
    write(unit,adjustl(trim(format)),iostat=err)(((prim(i,j,k)%r(:),i=1,Ni),j=1,Nj),k=1,Nk)
    err = write_vector(unit,format,prim%v)
    write(unit,adjustl(trim(format)),iostat=err)prim%p
    write(unit,adjustl(trim(format)),iostat=err)prim%d
    write(unit,adjustl(trim(format)),iostat=err)prim%g
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array3D

  ! read
  function Read_Bin_Scalar(unit,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Primitive (binary, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN)::    unit ! logic unit
  type(Type_Primitive), intent(INOUT):: prim
  integer(I4P)::                        err  ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)prim%r(:)
  err = read_vector(unit,prim%v)
  read(unit,iostat=err)prim%p
  read(unit,iostat=err)prim%d
  read(unit,iostat=err)prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Scalar

  function Read_Ascii_Scalar(unit,format,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Primitive (ascii, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN)::    unit   ! logic unit
  character(*),         intent(IN)::    format ! format specifier
  type(Type_Primitive), intent(INOUT):: prim
  integer(I4P)::                        err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)prim%r(:)
    err = read_vector(unit,format,prim%v)
    read(unit,*,iostat=err)prim%p
    read(unit,*,iostat=err)prim%d
    read(unit,*,iostat=err)prim%g
  case default
    read(unit,adjustl(trim(format)),iostat=err)prim%r(:)
    err = read_vector(unit,format,prim%v)
    read(unit,adjustl(trim(format)),iostat=err)prim%p
    read(unit,adjustl(trim(format)),iostat=err)prim%d
    read(unit,adjustl(trim(format)),iostat=err)prim%g
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Scalar

  function Read_Bin_Array1D(unit,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Primitive (binary, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN)::    unit    ! logic unit
  type(Type_Primitive), intent(INOUT):: prim(:)
  integer(I4P)::                        err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        Ni,i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  read(unit,iostat=err)(prim(i)%r(:),i=1,Ni)
  err = read_vector(unit,prim%v)
  read(unit,iostat=err)prim%p
  read(unit,iostat=err)prim%d
  read(unit,iostat=err)prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array1D

  function Read_Ascii_Array1D(unit,format,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Primitive (ascii, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN)::    unit    ! logic unit
  character(*),         intent(IN)::    format  ! format specifier
  type(Type_Primitive), intent(INOUT):: prim(:)
  integer(I4P)::                        err     ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        Ni,i
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)(prim(i)%r(:),i=1,Ni)
    err = read_vector(unit,format,prim%v)
    read(unit,*,iostat=err)prim%p
    read(unit,*,iostat=err)prim%d
    read(unit,*,iostat=err)prim%g
  case default
    read(unit,adjustl(trim(format)),iostat=err)(prim(i)%r(:),i=1,Ni)
    err = read_vector(unit,format,prim%v)
    read(unit,adjustl(trim(format)),iostat=err)prim%p
    read(unit,adjustl(trim(format)),iostat=err)prim%d
    read(unit,adjustl(trim(format)),iostat=err)prim%g
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array1D

  function Read_Bin_Array2D(unit,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Primitive (binary, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN)::    unit      ! logic unit
  type(Type_Primitive), intent(INOUT):: prim(:,:)
  integer(I4P)::                        err       ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        Ni,Nj,i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  Nj = size(prim,dim=2)
  read(unit,iostat=err)((prim(i,j)%r(:),i=1,Ni),j=1,Nj)
  err = read_vector(unit,prim%v)
  read(unit,iostat=err)prim%p
  read(unit,iostat=err)prim%d
  read(unit,iostat=err)prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array2D

  function Read_Ascii_Array2D(unit,format,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Primitive (ascii, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN)::    unit     ! logic unit
  character(*),         intent(IN)::    format   ! format specifier
  type(Type_Primitive), intent(INOUT):: prim(:,:)
  integer(I4P)::                        err      ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        Ni,Nj,i,j
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  Nj = size(prim,dim=2)
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)((prim(i,j)%r(:),i=1,Ni),j=1,Nj)
    err = read_vector(unit,format,prim%v)
    read(unit,*,iostat=err)prim%p
    read(unit,*,iostat=err)prim%d
    read(unit,*,iostat=err)prim%g
  case default
    read(unit,adjustl(trim(format)),iostat=err)((prim(i,j)%r(:),i=1,Ni),j=1,Nj)
    err = read_vector(unit,format,prim%v)
    read(unit,adjustl(trim(format)),iostat=err)prim%p
    read(unit,adjustl(trim(format)),iostat=err)prim%d
    read(unit,adjustl(trim(format)),iostat=err)prim%g
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array2D

  function Read_Bin_Array3D(unit,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Primitive (binary, Array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN)::    unit        ! logic unit
  type(Type_Primitive), intent(INOUT):: prim(:,:,:)
  integer(I4P)::                        err        ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        Ni,Nj,Nk,Ns,i,j,k,s
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  Nj = size(prim,dim=2)
  Nk = size(prim,dim=3)
  Ns = size(prim(1,1,1)%r(:),dim=1)
  read(unit,iostat=err)((((prim(i,j,k)%r(s),s=1,Ns),i=1,Ni),j=1,Nj),k=1,Nk)
  err = read_vector(unit,prim%v)
  read(unit,iostat=err)prim%p
  read(unit,iostat=err)prim%d
  read(unit,iostat=err)prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array3D

  function Read_Ascii_Array3D(unit,format,prim) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for reading Type_Primitive (ascii, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN)::    unit        ! logic unit
  character(*),         intent(IN)::    format      ! format specifier
  type(Type_Primitive), intent(INOUT):: prim(:,:,:)
  integer(I4P)::                        err         ! Error traping flag: 0 no errors, >0 error occours.
  integer(I4P)::                        Ni,Nj,Nk,i,j,k
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ni = size(prim,dim=1)
  Nj = size(prim,dim=2)
  Nk = size(prim,dim=3)
  !Ns = size(prim(1,1,1)%r(:))
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)(((prim(i,j,k)%r(:),i=1,Ni),j=1,Nj),k=1,Nk)
    err = read_vector(unit,format,prim%v)
    read(unit,*,iostat=err)prim%p
    read(unit,*,iostat=err)prim%d
    read(unit,*,iostat=err)prim%g
  case default
    read(unit,adjustl(trim(format)),iostat=err)(((prim(i,j,k)%r(:),i=1,Ni),j=1,Nj),k=1,Nk)
    err = read_vector(unit,format,prim%v)
    read(unit,adjustl(trim(format)),iostat=err)prim%p
    read(unit,adjustl(trim(format)),iostat=err)prim%d
    read(unit,adjustl(trim(format)),iostat=err)prim%g
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array3D

  ! Assignment (=)
  elemental subroutine assign_prim(prim1,prim2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between two prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim1
  type(Type_Primitive), intent(IN)::    prim2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim1%r).and.allocated(prim2%r)) prim1%r = prim2%r
                                                 prim1%v = prim2%v
                                                 prim1%p = prim2%p
                                                 prim1%d = prim2%d
                                                 prim1%g = prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_prim

#ifdef r16p
  elemental subroutine assign_ScalR16P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R16P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  real(R16P),           intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r))  prim%r = real(scal,R_P)
                          prim%v = real(scal,R_P)
                          prim%p = real(scal,R_P)
                          prim%d = real(scal,R_P)
                          prim%g = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR16P
#endif

  elemental subroutine assign_ScalR8P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R8P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  real(R8P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r))  prim%r = real(scal,R_P)
                          prim%v = real(scal,R_P)
                          prim%p = real(scal,R_P)
                          prim%d = real(scal,R_P)
                          prim%g = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR8P

  elemental subroutine assign_ScalR4P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R4P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  real(R4P),            intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r))  prim%r = real(scal,R_P)
                          prim%v = real(scal,R_P)
                          prim%p = real(scal,R_P)
                          prim%d = real(scal,R_P)
                          prim%g = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR4P

  elemental subroutine assign_ScalI8P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I8P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I8P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r))  prim%r = real(scal,R_P)
                          prim%v = real(scal,R_P)
                          prim%p = real(scal,R_P)
                          prim%d = real(scal,R_P)
                          prim%g = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI8P

  elemental subroutine assign_ScalI4P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I4P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I4P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r))  prim%r = real(scal,R_P)
                          prim%v = real(scal,R_P)
                          prim%p = real(scal,R_P)
                          prim%d = real(scal,R_P)
                          prim%g = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI4P

  elemental subroutine assign_ScalI2P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I2P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I2P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r))  prim%r = real(scal,R_P)
                          prim%v = real(scal,R_P)
                          prim%p = real(scal,R_P)
                          prim%d = real(scal,R_P)
                          prim%g = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI2P

  elemental subroutine assign_ScalI1P(prim,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I1P) and prim.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(INOUT):: prim
  integer(I1P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(prim%r))  prim%r = real(scal,R_P)
                          prim%v = real(scal,R_P)
                          prim%p = real(scal,R_P)
                          prim%d = real(scal,R_P)
                          prim%g = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI1P

  ! Multiplication (*)
  elemental function prim_mul_prim(prim1,prim2) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply (by components) Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN)::  prim1 ! First prim obj.
  type(Type_Primitive), intent(IN)::  prim2 ! Second prim obj.
  type(Type_Primitive)::              mul   ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim1%r)))
  mul%r = prim1%r * prim2%r
  mul%v = prim1%v * prim2%v
  mul%p = prim1%p * prim2%p
  mul%d = prim1%d * prim2%d
  mul%g = prim1%g * prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_prim

#ifdef r16p
  elemental function ScalR16P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R16P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),           intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_mul_prim

  elemental function prim_mul_ScalR16P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R16P),           intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalR16P
#endif

  elemental function ScalR8P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_mul_prim

  elemental function prim_mul_ScalR8P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R8P),            intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalR8P

  elemental function ScalR4P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_mul_prim

  elemental function prim_mul_ScalR4P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R4P),            intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalR4P

  elemental function ScalI8P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_mul_prim

  elemental function prim_mul_ScalI8P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalI8P

  elemental function ScalI4P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_mul_prim

  elemental function prim_mul_ScalI4P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalI4P

  elemental function ScalI2P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I2P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_mul_prim

  elemental function prim_mul_ScalI2P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalI2P

  elemental function ScalI1P_mul_prim(scal,prim) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I1P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_mul_prim

  elemental function prim_mul_ScalI1P(prim,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply Primitive object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(mul%r)) allocate(mul%r(1:size(prim%r)))
  mul%r = real(scal,R_P) * prim%r
  mul%v = real(scal,R_P) * prim%v
  mul%p = real(scal,R_P) * prim%p
  mul%d = real(scal,R_P) * prim%d
  mul%g = real(scal,R_P) * prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_mul_ScalI1P

  ! Division (/)
  elemental function prim_div_prim(prim1,prim2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide (by components) Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN)::  prim1 ! First prim obj.
  type(Type_Primitive), intent(IN)::  prim2 ! Second prim obj.
  type(Type_Primitive)::              div   ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%r)) allocate(div%r(1:size(prim1%r)))
  div%r = prim1%r / prim2%r
  div%v = prim1%v / prim2%v
  div%p = prim1%p / prim2%p
  div%d = prim1%d / prim2%d
  div%g = prim1%g / prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_prim

#ifdef r16p
  elemental function prim_div_ScalR16P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R16P),           intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%r)) allocate(div%r(1:size(prim1%r)))
  div%r = prim%r / real(scal,R_P)
  div%v = prim%v / real(scal,R_P)
  div%p = prim%p / real(scal,R_P)
  div%d = prim%d / real(scal,R_P)
  div%g = prim%g / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalR16P
#endif

  elemental function prim_div_ScalR8P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R8P),            intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%r)) allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R_P)
  div%v = prim%v / real(scal,R_P)
  div%p = prim%p / real(scal,R_P)
  div%d = prim%d / real(scal,R_P)
  div%g = prim%g / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalR8P

  elemental function prim_div_ScalR4P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R4P),            intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%r)) allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R_P)
  div%v = prim%v / real(scal,R_P)
  div%p = prim%p / real(scal,R_P)
  div%d = prim%d / real(scal,R_P)
  div%g = prim%g / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalR4P

  elemental function prim_div_ScalI8P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%r)) allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R_P)
  div%v = prim%v / real(scal,R_P)
  div%p = prim%p / real(scal,R_P)
  div%d = prim%d / real(scal,R_P)
  div%g = prim%g / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalI8P

  elemental function prim_div_ScalI4P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%r)) allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R_P)
  div%v = prim%v / real(scal,R_P)
  div%p = prim%p / real(scal,R_P)
  div%d = prim%d / real(scal,R_P)
  div%g = prim%g / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalI4P

  elemental function prim_div_ScalI2P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%r)) allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R_P)
  div%v = prim%v / real(scal,R_P)
  div%p = prim%p / real(scal,R_P)
  div%d = prim%d / real(scal,R_P)
  div%g = prim%g / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalI2P

  elemental function prim_div_ScalI1P(prim,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide Primitive object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(div%r)) allocate(div%r(1:size(prim%r)))
  div%r = prim%r / real(scal,R_P)
  div%v = prim%v / real(scal,R_P)
  div%p = prim%p / real(scal,R_P)
  div%d = prim%d / real(scal,R_P)
  div%g = prim%g / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_div_ScalI1P

  ! Sum (+)
  elemental function positive_prim(prim) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary + to a Primitive objecy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(pos%r)) allocate(pos%r(1:size(prim%r)))
  pos%r =  + prim%r
  pos%v =  + prim%v
  pos%p =  + prim%p
  pos%d =  + prim%d
  pos%g =  + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_prim

  elemental function prim_sum_prim(prim1,prim2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum (by components) Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN)::  prim1 ! First prim obj.
  type(Type_Primitive), intent(IN)::  prim2 ! Second prim obj.
  type(Type_Primitive)::              summ  ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim1%r)))
  summ%r = prim1%r + prim2%r
  summ%v = prim1%v + prim2%v
  summ%p = prim1%p + prim2%p
  summ%d = prim1%d + prim2%d
  summ%g = prim1%g + prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_prim

#ifdef r16p
  elemental function ScalR16P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R16P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),           intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sum_prim

  elemental function prim_sum_ScalR16P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R16P),           intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalR16P
#endif

  elemental function ScalR8P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sum_prim

  elemental function prim_sum_ScalR8P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R8P),            intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalR8P

  elemental function ScalR4P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sum_prim

  elemental function prim_sum_ScalR4P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R4P),            intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalR4P

  elemental function ScalI8P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sum_prim

  elemental function prim_sum_ScalI8P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalI8P

  elemental function ScalI4P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sum_prim

  elemental function prim_sum_ScalI4P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalI4P

  elemental function ScalI2P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I2P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sum_prim

  elemental function prim_sum_ScalI2P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalI2P

  elemental function ScalI1P_sum_prim(scal,prim) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I1P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sum_prim

  elemental function prim_sum_ScalI1P(prim,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum Primitive object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(summ%r)) allocate(summ%r(1:size(prim%r)))
  summ%r = real(scal,R_P) + prim%r
  summ%v = real(scal,R_P) + prim%v
  summ%p = real(scal,R_P) + prim%p
  summ%d = real(scal,R_P) + prim%d
  summ%g = real(scal,R_P) + prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sum_ScalI1P

  ! Subtraction (-)
  elemental function negative_prim(prim) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary - to a Primitive objecy.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(neg%r)) allocate(neg%r(1:size(prim%r)))
  neg%r =  - prim%r
  neg%v =  - prim%v
  neg%p =  - prim%p
  neg%d =  - prim%d
  neg%g =  - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_prim

  elemental function prim_sub_prim(prim1,prim2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract (by components) Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN)::  prim1 ! First prim obj.
  type(Type_Primitive), intent(IN)::  prim2 ! Second prim obj.
  type(Type_Primitive)::              sub  ! Resulting obj.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim1%r)))
  sub%r = prim1%r - prim2%r
  sub%v = prim1%v - prim2%v
  sub%p = prim1%p - prim2%p
  sub%d = prim1%d - prim2%d
  sub%g = prim1%g - prim2%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_prim

#ifdef r16p
  elemental function ScalR16P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R16P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),           intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R_P) - prim%r
  sub%v = real(scal,R_P) - prim%v
  sub%p = real(scal,R_P) - prim%p
  sub%d = real(scal,R_P) - prim%d
  sub%g = real(scal,R_P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sub_prim

  elemental function prim_sub_ScalR16P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R16P),           intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R_P)
  sub%v = prim%v - real(scal,R_P)
  sub%p = prim%p - real(scal,R_P)
  sub%d = prim%d - real(scal,R_P)
  sub%g = prim%g - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalR16P
#endif

  elemental function ScalR8P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R_P) - prim%r
  sub%v = real(scal,R_P) - prim%v
  sub%p = real(scal,R_P) - prim%p
  sub%d = real(scal,R_P) - prim%d
  sub%g = real(scal,R_P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sub_prim

  elemental function prim_sub_ScalR8P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R8P),            intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R_P)
  sub%v = prim%v - real(scal,R_P)
  sub%p = prim%p - real(scal,R_P)
  sub%d = prim%d - real(scal,R_P)
  sub%g = prim%g - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalR8P

  elemental function ScalR4P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),            intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R_P) - prim%r
  sub%v = real(scal,R_P) - prim%v
  sub%p = real(scal,R_P) - prim%p
  sub%d = real(scal,R_P) - prim%d
  sub%g = real(scal,R_P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sub_prim

  elemental function prim_sub_ScalR4P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R4P),            intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R_P)
  sub%v = prim%v - real(scal,R_P)
  sub%p = prim%p - real(scal,R_P)
  sub%d = prim%d - real(scal,R_P)
  sub%g = prim%g - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalR4P

  elemental function ScalI8P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I8P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R_P) - prim%r
  sub%v = real(scal,R_P) - prim%v
  sub%p = real(scal,R_P) - prim%p
  sub%d = real(scal,R_P) - prim%d
  sub%g = real(scal,R_P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sub_prim

  elemental function prim_sub_ScalI8P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I8P),         intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R_P)
  sub%v = prim%v - real(scal,R_P)
  sub%p = prim%p - real(scal,R_P)
  sub%d = prim%d - real(scal,R_P)
  sub%g = prim%g - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalI8P

  elemental function ScalI4P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I4P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R_P) - prim%r
  sub%v = real(scal,R_P) - prim%v
  sub%p = real(scal,R_P) - prim%p
  sub%d = real(scal,R_P) - prim%d
  sub%g = real(scal,R_P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sub_prim

  elemental function prim_sub_ScalI4P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I4P),         intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R_P)
  sub%v = prim%v - real(scal,R_P)
  sub%p = prim%p - real(scal,R_P)
  sub%d = prim%d - real(scal,R_P)
  sub%g = prim%g - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalI4P

  elemental function ScalI2P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I2P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R_P) - prim%r
  sub%v = real(scal,R_P) - prim%v
  sub%p = real(scal,R_P) - prim%p
  sub%d = real(scal,R_P) - prim%d
  sub%g = real(scal,R_P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sub_prim

  elemental function prim_sub_ScalI2P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I2P),         intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R_P)
  sub%v = prim%v - real(scal,R_P)
  sub%p = prim%p - real(scal,R_P)
  sub%d = prim%d - real(scal,R_P)
  sub%g = prim%g - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalI2P

  elemental function ScalI1P_sub_prim(scal,prim) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I1P) for Primitive object.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive), intent(IN):: prim
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = real(scal,R_P) - prim%r
  sub%v = real(scal,R_P) - prim%v
  sub%p = real(scal,R_P) - prim%p
  sub%d = real(scal,R_P) - prim%d
  sub%g = real(scal,R_P) - prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sub_prim

  elemental function prim_sub_ScalI1P(prim,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract Primitive object for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  integer(I1P),         intent(IN):: scal
  type(Type_Primitive)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.allocated(sub%r)) allocate(sub%r(1:size(prim%r)))
  sub%r = prim%r - real(scal,R_P)
  sub%v = prim%v - real(scal,R_P)
  sub%p = prim%p - real(scal,R_P)
  sub%d = prim%d - real(scal,R_P)
  sub%g = prim%g - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim_sub_ScalI1P

  pure function prim2array(prim) result(array)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting derived type Type_Primitive to array.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive), intent(IN):: prim
  real(R_P)::                        array(1:size(prim%r)+6)
  integer(I_P)::                     Ns
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(prim%r)
  array(1:Ns) = prim%r
  array(Ns+1) = prim%v%x
  array(Ns+2) = prim%v%y
  array(Ns+3) = prim%v%z
  array(Ns+4) = prim%p
  array(Ns+5) = prim%d
  array(Ns+6) = prim%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction prim2array

  pure function array2prim(array) result(prim)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for converting array to derived type Type_Primitive.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN):: array(:)
  type(Type_Primitive)::  prim
  integer(I_P)::          Ns
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = ubound(array,dim=1)-6
  if (allocated(prim%r)) deallocate(prim%r)  ; allocate(prim%r(1:Ns))
  prim%r   = array(1:Ns)
  prim%v%x = array(Ns+1)
  prim%v%y = array(Ns+2)
  prim%v%z = array(Ns+3)
  prim%p   = array(Ns+4)
  prim%d   = array(Ns+5)
  prim%g   = array(Ns+6)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction array2prim
endmodule Data_Type_Primitive
