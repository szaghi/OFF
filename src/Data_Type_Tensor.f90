module Data_Type_Tensor
!-----------------------------------------------------------------------------------------------------------------------------------
!!This module contains Data_Type_Tensor, a "class" for manipulating second order tensors in 3D space. The components of the tensors
!!are derived type Type_Vector. The components are defined in a three-dimensional cartesian frame of reference.
!!All the tensorial math procedures (cross, dot products, normL2...) assume a three-dimensional cartesian frame of reference.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                       ! Integers and reals precision definition.
USE Data_Type_Vector, init_vec => init, set_vec => set, get_vec => get ! Definition of type Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: Type_Tensor
public:: unity
public:: init,set,get
public:: write,read
public:: assignment (=)
public:: operator (*)
public:: operator (/)
public:: operator (+)
public:: operator (-)
public:: operator (.ddot.)
public:: operator (.dot.)
public:: operator (.diad.)
public:: sq_norm
public:: normL2
public:: normalize
public:: transpose
public:: determinant
public:: invert,invertible
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Type_Tensor definition:
type:: Type_Tensor
  sequence
  type(Type_Vector):: x ! Cartesian vector component in x direction    |       |T11 T12 T13|   |x%x x%y x%z|.
  type(Type_Vector):: y ! Cartesian vector component in y direction => | Tij = |T21 T22 T23| = |y%x y%y y%z|.
  type(Type_Vector):: z ! Cartesian vector component in z direction    |       |T31 T32 T33|   |z%x z%y z%z|.
endtype Type_Tensor
type(Type_Tensor), parameter:: unity = Type_Tensor(ex,ey,ez)
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Write overloading.
interface write
  module procedure Write_Bin_Scalar,    & ! binary scalar
                   Write_Ascii_Scalar,  & ! ascii scalar
                   Write_Bin_Array1D,   & ! binary array 1D
                   Write_Ascii_Array1D, & ! ascii array 1D
                   Write_Bin_Array2D,   & ! binary array 2D
                   Write_Ascii_Array2D, & ! ascii array 2D
                   Write_Bin_Array3D,   & ! binary array 3D
                   Write_Ascii_Array3D    ! ascii array 3D
endinterface
!!Read overloading.
interface read
  module procedure Read_Bin_Scalar,    & ! binary scalar
                   Read_Ascii_Scalar,  & ! ascii scalar
                   Read_Bin_Array1D,   & ! binary array 1D
                   Read_Ascii_Array1D, & ! ascii array 1D
                   Read_Bin_Array2D,   & ! binary array 2D
                   Read_Ascii_Array2D, & ! ascii array 2D
                   Read_Bin_Array3D,   & ! binary array 3D
                   Read_Ascii_Array3D    ! ascii array 3D
endinterface
!!Assigment (=) overloading.
interface assignment (=)
 module procedure                  &
                  assign_Vec,      &
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
  module procedure ten_mul_ten,      &
#ifdef r16p
                   ScalR16P_mul_ten, &
#endif
                   ScalR8P_mul_ten,  &
                   ScalR4P_mul_ten,  &
                   ScalI8P_mul_ten,  &
                   ScalI4P_mul_ten,  &
                   ScalI2P_mul_ten,  &
                   ScalI1P_mul_ten,  &
#ifdef r16p
                   ten_mul_ScalR16P, &
#endif
                   ten_mul_ScalR8P,  &
                   ten_mul_ScalR4P,  &
                   ten_mul_ScalI8P,  &
                   ten_mul_ScalI4P,  &
                   ten_mul_ScalI2P,  &
                   ten_mul_ScalI1P
endinterface
!!Division (/) overloading.
interface operator (/)
  module procedure ten_div_ten,      &
#ifdef r16p
                   ten_div_ScalR16P, &
#endif
                   ten_div_ScalR8P,  &
                   ten_div_ScalR4P,  &
                   ten_div_ScalI8P,  &
                   ten_div_ScalI4P,  &
                   ten_div_ScalI2P,  &
                   ten_div_ScalI1P
endinterface
!!Sum (+) overloading.
interface operator (+)
  module procedure positive_ten,     &
                   ten_sum_ten,      &
#ifdef r16p
                   ScalR16P_sum_ten, &
#endif
                   ScalR8P_sum_ten,  &
                   ScalR4P_sum_ten,  &
                   ScalI8P_sum_ten,  &
                   ScalI4P_sum_ten,  &
                   ScalI2P_sum_ten,  &
                   ScalI1P_sum_ten,  &
#ifdef r16p
                   ten_sum_ScalR16P, &
#endif
                   ten_sum_ScalR8P,  &
                   ten_sum_ScalR4P,  &
                   ten_sum_ScalI8P,  &
                   ten_sum_ScalI4P,  &
                   ten_sum_ScalI2P,  &
                   ten_sum_ScalI1P
endinterface
!!Subtraction (-) overloading.
interface operator (-)
  module procedure negative_ten,     &
                   ten_sub_ten,      &
#ifdef r16p
                   ScalR16P_sub_ten, &
#endif
                   ScalR8P_sub_ten,  &
                   ScalR4P_sub_ten,  &
                   ScalI8P_sub_ten,  &
                   ScalI4P_sub_ten,  &
                   ScalI2P_sub_ten,  &
                   ScalI1P_sub_ten,  &
#ifdef r16p
                   ten_sub_ScalR16P, &
#endif
                   ten_sub_ScalR8P,  &
                   ten_sub_ScalR4P,  &
                   ten_sub_ScalI8P,  &
                   ten_sub_ScalI4P,  &
                   ten_sub_ScalI2P,  &
                   ten_sub_ScalI1P
endinterface
!!Double dot product operator definition.
interface operator (.ddot.)
  module procedure ddotproduct
endinterface
!!Dot product operator definition.
interface operator (.dot.)
  module procedure ten_dot_vec,vec_dot_ten
endinterface
!!Diadic product operator definition.
interface operator (.diad.)
  module procedure diadicproduct
endinterface
!!sq_norm overloading.
interface sq_norm
  module procedure sq_norm,    & ! vector sq_norm
                   sq_norm_ten   ! tensor sq_norm
endinterface
!!normL2 overloading.
interface normL2
  module procedure normL2,    & ! vector normL2
                   normL2_ten   ! tensor normL2
endinterface
!!normalize overloading.
interface normalize
  module procedure normalize,    & ! vector normalize
                   normalize_ten   ! tensor normalize
endinterface
!!transpose overloading.
interface transpose
  module procedure transpose_ten
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function init(x,y,z) result(ten)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for initializing Type_Tensor: all components are initialized to zero and if there is a dummy argument the corresponding
  !!component is set to dummy value passed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN), optional:: x,y,z ! Vector's components.
  type(Type_Tensor)::                       ten   ! Tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(x)) ten%x = x
  if (present(y)) ten%y = y
  if (present(z)) ten%z = z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction init

  elemental subroutine set(x,y,z,ten)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment Type_Tensor: if there is a dummy argument the corresponding component is set to dummy value passed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN), optional:: x,y,z ! Vector's components.
  type(Type_Tensor), intent(INOUT)::        ten   ! Tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(x)) ten%x = x
  if (present(y)) ten%y = y
  if (present(z)) ten%z = z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  elemental subroutine get(x,y,z,ten)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for extraction Type_Tensor components: if there is a dummy argument it is set to the corresponding component value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(OUT), optional:: x,y,z ! Vector's components.
  type(Type_Tensor), intent(IN)::            ten   ! Tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(x)) x = ten%x
  if (present(y)) y = ten%y
  if (present(z)) z = ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get

  ! write
  function Write_Bin_Scalar(unit,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Tensor (binary, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit ! Logic unit.
  type(Type_Tensor), intent(IN):: ten  ! Tensor.
  integer(I_P)::                  err  ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)ten
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Scalar

  function Write_Ascii_Scalar(unit,format,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Tensor (ascii, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit   ! Logic unit.
  character(*),      intent(IN):: format ! Format specifier.
  type(Type_Tensor), intent(IN):: ten    ! Tensor.
  integer(I_P)::                  err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)ten
  case default
    write(unit,adjustl(trim(format)),iostat=err)ten
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Scalar

  function Write_Bin_Array1D(unit,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Tensor (binary, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit   ! Logic unit.
  type(Type_Tensor), intent(IN):: ten(:) ! Tensor.
  integer(I_P)::                  err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)ten
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array1D

  function Write_Ascii_Array1D(unit,format,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Tensor (ascii, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit   ! Logic unit.
  character(*),      intent(IN):: format ! Format specifier.
  type(Type_Tensor), intent(IN):: ten(:) ! Tensor.
  integer(I_P)::                  err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)ten
  case default
    write(unit,adjustl(trim(format)),iostat=err)ten
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array1D

  function Write_Bin_Array2D(unit,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Tensor (binary, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit     ! Logic unit.
  type(Type_Tensor), intent(IN):: ten(:,:) ! Tensor.
  integer(I_P)::                  err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)ten
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array2D

  function Write_Ascii_Array2D(unit,format,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Tensor (ascii, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit     ! Logic unit.
  character(*),      intent(IN):: format   ! Format specifier.
  type(Type_Tensor), intent(IN):: ten(:,:) ! Tensor.
  integer(I_P)::                  err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)ten
  case default
    write(unit,adjustl(trim(format)),iostat=err)ten
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array2D

  function Write_Bin_Array3D(unit,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Tensor (binary, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit       ! Logic unit.
  type(Type_Tensor), intent(IN):: ten(:,:,:) ! Tensor.
  integer(I_P)::                  err        ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)ten
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array3D

  function Write_Ascii_Array3D(unit,format,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Tensor (ascii, Array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit       ! Logic unit
  character(*),      intent(IN):: format     ! Format specifier.
  type(Type_Tensor), intent(IN):: ten(:,:,:) ! Tensor.
  integer(I_P)::                  err        ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)ten
  case default
    write(unit,adjustl(trim(format)),iostat=err)ten
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array3D

  ! read
  function Read_Bin_Scalar(unit,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Tensor (binary, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit  ! Logic unit.
  type(Type_Tensor), intent(INOUT):: ten   ! Tensor.
  integer(I_P)::                     err   ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)ten
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Scalar

  function Read_Ascii_Scalar(unit,format,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Tensor (ascii, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit   ! Logic unit.
  character(*),      intent(IN)::    format ! Format specifier.
  type(Type_Tensor), intent(INOUT):: ten    ! Tensor.
  integer(I_P)::                     err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)ten
  case default
    read(unit,adjustl(trim(format)),iostat=err)ten
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Scalar

  function Read_Bin_Array1D(unit,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Tensor (binary, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit   ! Logic unit.
  type(Type_Tensor), intent(INOUT):: ten(:) ! Tensor.
  integer(I_P)::                     err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)ten
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array1D

  function Read_Ascii_Array1D(unit,format,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Tensor (ascii, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit   ! logic unit
  character(*),      intent(IN)::    format ! format specifier
  type(Type_Tensor), intent(INOUT):: ten(:) ! Tensor.
  integer(I_P)::                     err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)ten
  case default
    read(unit,adjustl(trim(format)),iostat=err)ten
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array1D

  function Read_Bin_Array2D(unit,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Tensor (binary, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit     ! Logic unit.
  type(Type_Tensor), intent(INOUT):: ten(:,:) ! Tensor.
  integer(I_P)::                     err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)ten
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array2D

  function Read_Ascii_Array2D(unit,format,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Tensor (ascii, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit     ! Logic unit.
  character(*),      intent(IN)::    format   ! Format specifier.
  type(Type_Tensor), intent(INOUT):: ten(:,:) ! Tensor.
  integer(I_P)::                     err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)ten
  case default
    read(unit,adjustl(trim(format)),iostat=err)ten
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array2D

  function Read_Bin_Array3D(unit,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Tensor (binary, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit       ! Logic unit.
  type(Type_Tensor), intent(INOUT):: ten(:,:,:) ! Tensor.
  integer(I_P)::                     err        ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)ten
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array3D

  function Read_Ascii_Array3D(unit,format,ten) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for reading Type_Tensor (ascii, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit       ! Logic unit.
  character(*),      intent(IN)::    format     ! Format specifier.
  type(Type_Tensor), intent(INOUT):: ten(:,:,:) ! Tensor.
  integer(I_P)::                     err        ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)ten
  case default
    read(unit,adjustl(trim(format)),iostat=err)ten
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array3D

  ! Assignment (=)
  elemental subroutine assign_Vec(ten,vec)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a vector and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT):: ten
  type(Type_Vector), intent(IN)::    vec
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = vec
  ten%y = vec
  ten%z = vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_Vec

#ifdef r16p
  elemental subroutine assign_ScalR16P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R16P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT):: ten
  real(R16P),        intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R_P)
  ten%y = real(scal,R_P)
  ten%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR16P
#endif

  elemental subroutine assign_ScalR8P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R8P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT):: ten
  real(R8P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R_P)
  ten%y = real(scal,R_P)
  ten%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR8P

  elemental subroutine assign_ScalR4P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R4P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT):: ten
  real(R4P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R_P)
  ten%y = real(scal,R_P)
  ten%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR4P

  elemental subroutine assign_ScalI8P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I8P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT):: ten
  integer(I8P),      intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R_P)
  ten%y = real(scal,R_P)
  ten%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI8P

  elemental subroutine assign_ScalI4P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I4P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT):: ten
  integer(I4P),      intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R_P)
  ten%y = real(scal,R_P)
  ten%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI4P

  elemental subroutine assign_ScalI2P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I2P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT):: ten
  integer(I2P),      intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R_P)
  ten%y = real(scal,R_P)
  ten%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI2P

  elemental subroutine assign_ScalI1P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I1P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT):: ten
  integer(I1P),      intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R_P)
  ten%y = real(scal,R_P)
  ten%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI1P

  ! Multiplication (*)
  elemental function ten_mul_ten(ten1,ten2) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply (by components) tensors.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten1 ! First tensor.
  type(Type_Tensor), intent(IN):: ten2 ! Second tensor.
  type(Type_Tensor)::             mul  ! Resulting tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = ten1%x * ten2%x
  mul%y = ten1%y * ten2%y
  mul%z = ten1%z * ten2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ten

#ifdef r16p
  elemental function ScalR16P_mul_ten(scal,ten) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R16P) for ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_mul_ten

  elemental function ten_mul_ScalR16P(ten,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply ten for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R16P),        intent(IN):: scal
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ScalR16P
#endif

  elemental function ScalR8P_mul_ten(scal,ten) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R8P) for ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_mul_ten

  elemental function ten_mul_ScalR8P(ten,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply ten for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R8P),         intent(IN):: scal
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ScalR8P

  elemental function ScalR4P_mul_ten(scal,ten) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R4P) for ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_mul_ten

  elemental function ten_mul_ScalR4P(ten,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply ten for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R4P),         intent(IN):: scal
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ScalR4P

  elemental function ScalI8P_mul_ten(scal,ten) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I8P) for ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_mul_ten

  elemental function ten_mul_ScalI8P(ten,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply ten for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I8P),      intent(IN):: scal
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ScalI8P

  elemental function ScalI4P_mul_ten(scal,ten) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I4P) for ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_mul_ten

  elemental function ten_mul_ScalI4P(ten,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply ten for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I4P),      intent(IN):: scal
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ScalI4P

  elemental function ScalI2P_mul_ten(scal,ten) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I2P) for ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_mul_ten

  elemental function ten_mul_ScalI2P(ten,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply ten for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I2P),      intent(IN):: scal
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ScalI2P

  elemental function ScalI1P_mul_ten(scal,ten) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I1P) for ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_mul_ten

  elemental function ten_mul_ScalI1P(ten,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply ten for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I1P),      intent(IN):: scal
  type(Type_Tensor)::             mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * ten%x
  mul%y = real(scal,R_P) * ten%y
  mul%z = real(scal,R_P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ScalI1P

  ! Division (/)
  elemental function ten_div_ten(ten1,ten2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide ten for ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten1
  type(Type_Tensor), intent(IN):: ten2
  type(Type_Tensor)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = ten1%x / ten2%x
  div%y = ten1%y / ten2%y
  div%z = ten1%z / ten2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ten

#ifdef r16p
  elemental function ten_div_ScalR16P(ten,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide ten for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R16P),        intent(IN):: scal
  type(Type_Tensor)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = ten%x / real(scal,R_P)
  div%y = ten%y / real(scal,R_P)
  div%z = ten%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ScalR16P
#endif

  elemental function ten_div_ScalR8P(ten,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide ten for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R8P),         intent(IN):: scal
  type(Type_Tensor)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = ten%x / real(scal,R_P)
  div%y = ten%y / real(scal,R_P)
  div%z = ten%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ScalR8P

  elemental function ten_div_ScalR4P(ten,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide ten for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R4P),         intent(IN):: scal
  type(Type_Tensor)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = ten%x / real(scal,R_P)
  div%y = ten%y / real(scal,R_P)
  div%z = ten%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ScalR4P

  elemental function ten_div_ScalI8P(ten,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide ten for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I8P),      intent(IN):: scal
  type(Type_Tensor)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = ten%x / real(scal,R_P)
  div%y = ten%y / real(scal,R_P)
  div%z = ten%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ScalI8P

  elemental function ten_div_ScalI4P(ten,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide ten for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I4P),      intent(IN):: scal
  type(Type_Tensor)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = ten%x / real(scal,R_P)
  div%y = ten%y / real(scal,R_P)
  div%z = ten%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ScalI4P

  elemental function ten_div_ScalI2P(ten,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide ten for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I2P),      intent(IN):: scal
  type(Type_Tensor)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = ten%x / real(scal,R_P)
  div%y = ten%y / real(scal,R_P)
  div%z = ten%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ScalI2P

  elemental function ten_div_ScalI1P(ten,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide ten for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I1P),      intent(IN):: scal
  type(Type_Tensor)::             div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = ten%x / real(scal,R_P)
  div%y = ten%y / real(scal,R_P)
  div%z = ten%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ScalI1P

  ! Sum (+)
  elemental function positive_ten(ten) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary + to an ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pos%x =  + ten%x
  pos%y =  + ten%y
  pos%z =  + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_ten

  elemental function ten_sum_ten(ten1,ten2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten1
  type(Type_Tensor), intent(IN):: ten2
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = ten1%x + ten2%x
  summ%y = ten1%y + ten2%y
  summ%z = ten1%z + ten2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ten

#ifdef r16p
  elemental function ScalR16P_sum_ten(scal,ten) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R16P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sum_ten

  elemental function ten_sum_ScalR16P(ten,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum ten and scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R16P),        intent(IN):: scal
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ScalR16P
#endif

  elemental function ScalR8P_sum_ten(scal,ten) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R8P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sum_ten

  elemental function ten_sum_ScalR8P(ten,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum ten and scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R8P),         intent(IN):: scal
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ScalR8P

  elemental function ScalR4P_sum_ten(scal,ten) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R4P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sum_ten

  elemental function ten_sum_ScalR4P(ten,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum ten and scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R4P),         intent(IN):: scal
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ScalR4P

  elemental function ScalI8P_sum_ten(scal,ten) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I8P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sum_ten

  elemental function ten_sum_ScalI8P(ten,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum ten and scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I8P),      intent(IN):: scal
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ScalI8P

  elemental function ScalI4P_sum_ten(scal,ten) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I4P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sum_ten

  elemental function ten_sum_ScalI4P(ten,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum ten and scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I4P),      intent(IN):: scal
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ScalI4P

  elemental function ScalI2P_sum_ten(scal,ten) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I2P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sum_ten

  elemental function ten_sum_ScalI2P(ten,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum ten and scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I2P),      intent(IN):: scal
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ScalI2P

  elemental function ScalI1P_sum_ten(scal,ten) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I1P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sum_ten

  elemental function ten_sum_ScalI1P(ten,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum ten and scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I1P),      intent(IN):: scal
  type(Type_Tensor)::             summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + ten%x
  summ%y = real(scal,R_P) + ten%y
  summ%z = real(scal,R_P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ScalI1P

  ! Subtraction (-)
  elemental function negative_ten(ten) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary - to an ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neg%x =  - ten%x
  neg%y =  - ten%y
  neg%z =  - ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_ten

  elemental function ten_sub_ten(ten1,ten2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten1
  type(Type_Tensor), intent(IN):: ten2
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = ten1%x - ten2%x
  sub%y = ten1%y - ten2%y
  sub%z = ten1%z - ten2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ten

#ifdef r16p
  elemental function ScalR16P_sub_ten(scal,ten) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R16P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - ten%x
  sub%y = real(scal,R_P) - ten%y
  sub%z = real(scal,R_P) - ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sub_ten

  elemental function ten_sub_ScalR16P(ten,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract ten and scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R16P),        intent(IN):: scal
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = ten%x - real(scal,R_P)
  sub%y = ten%y - real(scal,R_P)
  sub%z = ten%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ScalR16P
#endif

  elemental function ScalR8P_sub_ten(scal,ten) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R8P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - ten%x
  sub%y = real(scal,R_P) - ten%y
  sub%z = real(scal,R_P) - ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sub_ten

  elemental function ten_sub_ScalR8P(ten,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract ten and scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R8P),         intent(IN):: scal
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = ten%x - real(scal,R_P)
  sub%y = ten%y - real(scal,R_P)
  sub%z = ten%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ScalR8P

  elemental function ScalR4P_sub_ten(scal,ten) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R4P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - ten%x
  sub%y = real(scal,R_P) - ten%y
  sub%z = real(scal,R_P) - ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sub_ten

  elemental function ten_sub_ScalR4P(ten,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract ten and scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  real(R4P),         intent(IN):: scal
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = ten%x - real(scal,R_P)
  sub%y = ten%y - real(scal,R_P)
  sub%z = ten%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ScalR4P

  elemental function ScalI8P_sub_ten(scal,ten) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I8P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - ten%x
  sub%y = real(scal,R_P) - ten%y
  sub%z = real(scal,R_P) - ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sub_ten

  elemental function ten_sub_ScalI8P(ten,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract ten and scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I8P),      intent(IN):: scal
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = ten%x - real(scal,R_P)
  sub%y = ten%y - real(scal,R_P)
  sub%z = ten%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ScalI8P

  elemental function ScalI4P_sub_ten(scal,ten) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I4P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - ten%x
  sub%y = real(scal,R_P) - ten%y
  sub%z = real(scal,R_P) - ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sub_ten

  elemental function ten_sub_ScalI4P(ten,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract ten and scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I4P),      intent(IN):: scal
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = ten%x - real(scal,R_P)
  sub%y = ten%y - real(scal,R_P)
  sub%z = ten%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ScalI4P

  elemental function ScalI2P_sub_ten(scal,ten) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I2P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - ten%x
  sub%y = real(scal,R_P) - ten%y
  sub%z = real(scal,R_P) - ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sub_ten

  elemental function ten_sub_ScalI2P(ten,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract ten and scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I2P),      intent(IN):: scal
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = ten%x - real(scal,R_P)
  sub%y = ten%y - real(scal,R_P)
  sub%z = ten%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ScalI2P

  elemental function ScalI1P_sub_ten(scal,ten) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I1P) and ten.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal
  type(Type_Tensor), intent(IN):: ten
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - ten%x
  sub%y = real(scal,R_P) - ten%y
  sub%z = real(scal,R_P) - ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sub_ten

  elemental function ten_sub_ScalI1P(ten,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract ten and scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten
  integer(I1P),      intent(IN):: scal
  type(Type_Tensor)::             sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = ten%x - real(scal,R_P)
  sub%y = ten%y - real(scal,R_P)
  sub%z = ten%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ScalI1P

  ! dot product (.dot.)
  elemental function ten_dot_vec(ten,vec) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the vector (dot) product of a tensor and a vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten ! Tensor.
  type(Type_Vector), intent(IN):: vec ! Vector.
  type(Type_Vector)::             dot ! Dot product (vector).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call set_vec(x=(ten%x.dot.vec),y=(ten%y.dot.vec),z=(ten%z.dot.vec),vec=dot)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_dot_vec

  elemental function vec_dot_ten(vec,ten) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the vector (dot) product of a vector and a tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec ! Vector.
  type(Type_Tensor), intent(IN):: ten ! Tensor.
  type(Type_Vector)::             dot ! Dot product (vector).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call set_vec(x=((ten.dot.ex).dot.vec),y=((ten.dot.ey).dot.vec),z=((ten.dot.ez).dot.vec),vec=dot)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_dot_ten

  ! double dot product (.ddot.)
  elemental function ddotproduct(ten1,ten2) result(ddot)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the scalar (double dot) product of 2 tensors.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten1 ! First tensor.
  type(Type_Tensor), intent(IN):: ten2 ! Second tensor.
  real(R_P)::                     ddot ! Double dot product.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ddot = (ten1%x.dot.(ten2.dot.ex)) + (ten1%y.dot.(ten2.dot.ey)) + (ten1%z.dot.(ten2.dot.ez))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ddotproduct

  ! diadic product (.diad.)
  elemental function diadicproduct(vec1,vec2) result(ten)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the diadic product of 2 vectors producing a second order tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1 ! First vector.
  type(Type_Vector), intent(IN):: vec2 ! Second vector.
  type(Type_Tensor)::             ten  ! Tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = vec1%x*vec2
  ten%y = vec1%y*vec2
  ten%z = vec1%z*vec2
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction diadicproduct

  ! sq_norm
  elemental function sq_norm_ten(ten) result(sq)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the square of the norm of a tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten ! Tensor.
  real(R_P)::                     sq  ! Square of the Norm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sq = sq_norm(ten%x) + sq_norm(ten%y) + sq_norm(ten%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sq_norm_ten

  ! normL2
  elemental function normL2_ten(ten) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the norm L2 of a tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten  ! Tensor.
  real(R_P)::                     norm ! Norm L2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  norm = sqrt(sq_norm(ten%x) + sq_norm(ten%y) + sq_norm(ten%z))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normL2_ten

  ! normalize
  elemental function normalize_ten(ten) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function normalize a tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten  ! Tensor to be normalized.
  type(Type_Tensor)::             norm ! Tensor normalized.
  real(R_P)::                     nm   ! Norm L2 of tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  nm = normL2_ten(ten)
  if (nm < smallR_P) then
    nm = nm + smallR_P
  endif
  norm%x = ten%x/nm
  norm%y = ten%y/nm
  norm%z = ten%z/nm
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normalize_ten

  ! transpose
  elemental function transpose_ten(ten) result(tran)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function transpose a tensor.
  !       |x%x x%y x%z|           |x%x y%x z%x|
  ! ten = |y%x y%y y%z| => tran = |x%y y%y z%y|
  !       |z%x z%y z%z|           |x%z y%z z%z|
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten  ! Tensor to be transposed.
  type(Type_Tensor)::             tran ! Tensor transposed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  tran%x = ten.dot.ex
  tran%y = ten.dot.ey
  tran%z = ten.dot.ez
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction transpose_ten

  ! determinant
  elemental function determinant(ten) result(det)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function compute the determinant of a tensor.
  !!      |x%x x%y x%z|
  !!ten = |y%x y%y y%z| => det = x%x*(z%z*y%y-z%y*y%z)-y%x*(z%z*x%y-z%y*x%z)+z%x*(y%z*x%y-y%y*x%z)
  !!      |z%x z%y z%z|
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten ! Tensor.
  real(R_P)::                     det ! Determinant.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  det=ten%x%x*(ten%z%z*ten%y%y-ten%z%y*ten%y%z)-ten%y%x*(ten%z%z*ten%x%y-ten%z%y*ten%x%z)+ten%z%x*(ten%y%z*ten%x%y-ten%y%y*ten%x%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction determinant

  ! invert
  elemental function invert(ten) result(inv)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function invert a tensor.
  !!   -1  |x%x x%y x%z|-1                 |   z%z*y%y-z%y*y%z  -(z%z*x%y-z%y*x%z)   y%z*x%y-y%y*x%z  |
  !!ten  = |y%x y%y y%z|  => inv = 1/DET * | -(z%z*y%x-z%x*y%z)   z%z*x%x-z%x*x%z  -(y%z*x%x-y%x*x%z) |
  !!       |z%x z%y z%z|                   |   z%y*y%x-z%x*y%y  -(z%y*x%x-z%x*x%y)   y%y*x%x-y%x*x%y  |
  !!with DET = x%x*(z%z*y%y-z%y*y%z)-y%x*(z%z*x%y-z%y*x%z)+z%x*(y%z*x%y-y%y*x%z)
  !!If the tensor is not invertible a null tensor is returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten    ! Tensor to be inverted.
  type(Type_Tensor)::             inv    ! Tensor inverted.
  real(R_P)::                     det,di ! Determinant and 1/Determinant.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  det = determinant(ten)
  if (det/=0._R_P) then
    di = 1._R_P/det
   inv%x=di*init_vec(x= (ten%z%z*ten%y%y-ten%z%y*ten%y%z),y=-(ten%z%z*ten%x%y-ten%z%y*ten%x%z),z= (ten%y%z*ten%x%y-ten%y%y*ten%x%z))
   inv%y=di*init_vec(x=-(ten%z%z*ten%y%x-ten%z%x*ten%y%z),y= (ten%z%z*ten%x%x-ten%z%x*ten%x%z),z=-(ten%y%z*ten%x%x-ten%y%x*ten%x%z))
   inv%z=di*init_vec(x= (ten%z%y*ten%y%x-ten%z%x*ten%y%y),y=-(ten%z%y*ten%x%x-ten%z%x*ten%x%y),z= (ten%y%y*ten%x%x-ten%y%x*ten%x%y))
  else
    inv%x=0._R_P
    inv%y=0._R_P
    inv%z=0._R_P
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction invert

  elemental function invertible(ten) result(inv)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function check if a tensor is invertible (determinant /=0, not singular tensor).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten ! Tensor to be inverted.
  logical::                       inv ! True if the tensor is not singular.
  real(R_P)::                     det ! Determinant.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  det = determinant(ten) ; inv = (det/=0._R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction invertible
endmodule Data_Type_Tensor
