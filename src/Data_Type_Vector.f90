module Data_Type_Vector
!-----------------------------------------------------------------------------------------------------------------------------------
!!This module contains Data_Type_Vector, a "class" for manipulating vectors in 3D space. The components of the vectors are real with
!!R_P kind as defined by the IR_Precision module. The components are defined in a three-dimensional cartesian frame of reference.
!!All the vectorial math procedures (cross, dot products, parallel...) assume a three-dimensional cartesian frame of reference.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: Type_Vector
public:: ex,ey,ez
public:: init,set,get
public:: pprint
public:: write,read
public:: assignment (=)
public:: operator (*)
public:: operator (/)
public:: operator (+)
public:: operator (-)
public:: operator (/=)
public:: operator (<)
public:: operator (<=)
public:: operator (==)
public:: operator (>=)
public:: operator (>)
public:: operator (.cross.)
public:: operator (.dot.)
public:: operator (.paral.)
public:: operator (.ortho.)
public:: sq_norm
public:: normL2
public:: normalize
public:: face_normal3,face_normal4
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!!Type_Vector definition:
type:: Type_Vector
  sequence
  real(R_P):: x = 0._R_P ! Cartesian component in x direction.
  real(R_P):: y = 0._R_P ! Cartesian component in y direction.
  real(R_P):: z = 0._R_P ! Cartesian component in z direction.
endtype Type_Vector
type(Type_Vector), parameter:: ex = Type_Vector(1._R_P,0._R_P,0._R_P) ! x versor.
type(Type_Vector), parameter:: ey = Type_Vector(0._R_P,1._R_P,0._R_P) ! x versor.
type(Type_Vector), parameter:: ez = Type_Vector(0._R_P,0._R_P,1._R_P) ! x versor.
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
  module procedure vec_mul_vec,      &
#ifdef r16p
                   ScalR16P_mul_vec, &
#endif
                   ScalR8P_mul_vec,  &
                   ScalR4P_mul_vec,  &
                   ScalI8P_mul_vec,  &
                   ScalI4P_mul_vec,  &
                   ScalI2P_mul_vec,  &
                   ScalI1P_mul_vec,  &
#ifdef r16p
                   vec_mul_ScalR16P, &
#endif
                   vec_mul_ScalR8P,  &
                   vec_mul_ScalR4P,  &
                   vec_mul_ScalI8P,  &
                   vec_mul_ScalI4P,  &
                   vec_mul_ScalI2P,  &
                   vec_mul_ScalI1P
endinterface
!!Division (/) overloading.
interface operator (/)
  module procedure vec_div_vec,      &
#ifdef r16p
                   vec_div_ScalR16P, &
#endif
                   vec_div_ScalR8P,  &
                   vec_div_ScalR4P,  &
                   vec_div_ScalI8P,  &
                   vec_div_ScalI4P,  &
                   vec_div_ScalI2P,  &
                   vec_div_ScalI1P
endinterface
!!Sum (+) overloading.
interface operator (+)
  module procedure positive_vec,     &
                   vec_sum_vec,      &
#ifdef r16p
                   ScalR16P_sum_vec, &
#endif
                   ScalR8P_sum_vec,  &
                   ScalR4P_sum_vec,  &
                   ScalI8P_sum_vec,  &
                   ScalI4P_sum_vec,  &
                   ScalI2P_sum_vec,  &
                   ScalI1P_sum_vec,  &
#ifdef r16p
                   vec_sum_ScalR16P, &
#endif
                   vec_sum_ScalR8P,  &
                   vec_sum_ScalR4P,  &
                   vec_sum_ScalI8P,  &
                   vec_sum_ScalI4P,  &
                   vec_sum_ScalI2P,  &
                   vec_sum_ScalI1P
endinterface
!!Subtraction (-) overloading.
interface operator (-)
  module procedure negative_vec,     &
                   vec_sub_vec,      &
#ifdef r16p
                   ScalR16P_sub_vec, &
#endif
                   ScalR8P_sub_vec,  &
                   ScalR4P_sub_vec,  &
                   ScalI8P_sub_vec,  &
                   ScalI4P_sub_vec,  &
                   ScalI2P_sub_vec,  &
                   ScalI1P_sub_vec,  &
#ifdef r16p
                   vec_sub_ScalR16P, &
#endif
                   vec_sub_ScalR8P,  &
                   vec_sub_ScalR4P,  &
                   vec_sub_ScalI8P,  &
                   vec_sub_ScalI4P,  &
                   vec_sub_ScalI2P,  &
                   vec_sub_ScalI1P
endinterface
!!Conditional operators overloading.
interface operator (/=)
  module procedure vec_not_eq_vec,  &
#ifdef r16p
                   R16P_not_eq_vec, &
                   vec_not_eq_R16P, &
#endif
                   R8P_not_eq_vec,  &
                   vec_not_eq_R8P,  &
                   R4P_not_eq_vec,  &
                   vec_not_eq_R4P,  &
                   I8P_not_eq_vec,  &
                   vec_not_eq_I8P,  &
                   I4P_not_eq_vec,  &
                   vec_not_eq_I4P,  &
                   I2P_not_eq_vec,  &
                   vec_not_eq_I2P,  &
                   I1P_not_eq_vec,  &
                   vec_not_eq_I1P
endinterface
interface operator (<)
  module procedure vec_low_vec,  &
#ifdef r16p
                   R16P_low_vec, &
                   vec_low_R16P, &
#endif
                   R8P_low_vec,  &
                   vec_low_R8P,  &
                   R4P_low_vec,  &
                   vec_low_R4P,  &
                   I8P_low_vec,  &
                   vec_low_I8P,  &
                   I4P_low_vec,  &
                   vec_low_I4P,  &
                   I2P_low_vec,  &
                   vec_low_I2P,  &
                   I1P_low_vec,  &
                   vec_low_I1P
endinterface
interface operator (<=)
  module procedure vec_low_eq_vec,  &
#ifdef r16p
                   R16P_low_eq_vec, &
                   vec_low_eq_R16P, &
#endif
                   R8P_low_eq_vec,  &
                   vec_low_eq_R8P,  &
                   R4P_low_eq_vec,  &
                   vec_low_eq_R4P,  &
                   I8P_low_eq_vec,  &
                   vec_low_eq_I8P,  &
                   I4P_low_eq_vec,  &
                   vec_low_eq_I4P,  &
                   I2P_low_eq_vec,  &
                   vec_low_eq_I2P,  &
                   I1P_low_eq_vec,  &
                   vec_low_eq_I1P
endinterface
interface operator (==)
  module procedure vec_eq_vec,  &
#ifdef r16p
                   R16P_eq_vec, &
                   vec_eq_R16P, &
#endif
                   R8P_eq_vec,  &
                   vec_eq_R8P,  &
                   R4P_eq_vec,  &
                   vec_eq_R4P,  &
                   I8P_eq_vec,  &
                   vec_eq_I8P,  &
                   I4P_eq_vec,  &
                   vec_eq_I4P,  &
                   I2P_eq_vec,  &
                   vec_eq_I2P,  &
                   I1P_eq_vec,  &
                   vec_eq_I1P
endinterface
interface operator (>=)
  module procedure vec_great_eq_vec,  &
#ifdef r16p
                   R16P_great_eq_vec, &
                   vec_great_eq_R16P, &
#endif
                   R8P_great_eq_vec,  &
                   vec_great_eq_R8P,  &
                   R4P_great_eq_vec,  &
                   vec_great_eq_R4P,  &
                   I8P_great_eq_vec,  &
                   vec_great_eq_I8P,  &
                   I4P_great_eq_vec,  &
                   vec_great_eq_I4P,  &
                   I2P_great_eq_vec,  &
                   vec_great_eq_I2P,  &
                   I1P_great_eq_vec,  &
                   vec_great_eq_I1P
endinterface
interface operator (>)
  module procedure vec_great_vec,  &
#ifdef r16p
                   R16P_great_vec, &
                   vec_great_R16P, &
#endif
                   R8P_great_vec,  &
                   vec_great_R8P,  &
                   R4P_great_vec,  &
                   vec_great_R4P,  &
                   I8P_great_vec,  &
                   vec_great_I8P,  &
                   I4P_great_vec,  &
                   vec_great_I4P,  &
                   I2P_great_vec,  &
                   vec_great_I2P,  &
                   I1P_great_vec,  &
                   vec_great_I1P
endinterface
!!Cross product operator definition.
interface operator (.cross.)
  module procedure crossproduct
endinterface
!!Dot product operator definition.
interface operator (.dot.)
  module procedure dotproduct
endinterface
!!Parallel operator definition.
interface operator (.paral.)
  module procedure parallel
endinterface
!!Orthogonal operator definition.
interface operator (.ortho.)
  module procedure orthogonal
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  elemental function init(x,y,z) result(vec)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for initializing Type_Vector: all components are initialized to zero and if there is a dummy argument the corresponding
  !!component is set to dummy value passed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P), intent(IN), optional:: x,y,z ! Vector's components.
  type(Type_Vector)::               vec   ! Vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(x)) vec%x = x
  if (present(y)) vec%y = y
  if (present(z)) vec%z = z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction init

  elemental subroutine set(x,y,z,vec)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!subroutine for assignment Type_Vector: if there is a dummy argument the corresponding component is set to dummy value passed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),         intent(IN), optional:: x,y,z ! Vector's components.
  type(Type_Vector), intent(INOUT)::        vec   ! Vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(x)) vec%x = x
  if (present(y)) vec%y = y
  if (present(z)) vec%z = z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  elemental subroutine get(x,y,z,vec)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for extraction Type_Vector components: if there is a dummy argument it is set to the corresponding component value.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),         intent(OUT), optional:: x,y,z ! Vector's components.
  type(Type_Vector), intent(IN)::            vec   ! Vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(x)) x = vec%x
  if (present(y)) y = vec%y
  if (present(z)) z = vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine get

  ! write
  function pprint(vec,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function that write to unit the components of vec with a "pretty" format.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec  ! Vector.
  integer(I_P),      intent(IN):: unit ! Logic unit.
  integer(I_P)::                  err  ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,'(A)',iostat=err)' Component x '//str(n=vec%x)
  write(unit,'(A)',iostat=err)' Component y '//str(n=vec%y)
  write(unit,'(A)',iostat=err)' Component z '//str(n=vec%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction pprint

  function Write_Bin_Scalar(unit,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Vector (binary, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit ! Logic unit.
  type(Type_Vector), intent(IN):: vec  ! Vector.
  integer(I_P)::                  err  ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Scalar

  function Write_Ascii_Scalar(unit,format,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Vector (ascii, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit   ! Logic unit.
  character(*),      intent(IN):: format ! Format specifier.
  type(Type_Vector), intent(IN):: vec    ! Vector.
  integer(I_P)::                  err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)vec
  case default
    write(unit,adjustl(trim(format)),iostat=err)vec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Scalar

  function Write_Bin_Array1D(unit,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Vector (binary, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit   ! Logic unit.
  type(Type_Vector), intent(IN):: vec(:) ! Vector.
  integer(I_P)::                  err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array1D

  function Write_Ascii_Array1D(unit,format,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Vector (ascii, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit   ! Logic unit.
  character(*),      intent(IN):: format ! Format specifier.
  type(Type_Vector), intent(IN):: vec(:) ! Vector.
  integer(I_P)::                  err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)vec
  case default
    write(unit,adjustl(trim(format)),iostat=err)vec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array1D

  function Write_Bin_Array2D(unit,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Vector (binary, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit     ! Logic unit.
  type(Type_Vector), intent(IN):: vec(:,:) ! Vector.
  integer(I_P)::                  err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array2D

  function Write_Ascii_Array2D(unit,format,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Vector (ascii, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit     ! Logic unit.
  character(*),      intent(IN):: format   ! Format specifier.
  type(Type_Vector), intent(IN):: vec(:,:) ! Vector.
  integer(I_P)::                  err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)vec
  case default
    write(unit,adjustl(trim(format)),iostat=err)vec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array2D

  function Write_Bin_Array3D(unit,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Vector (binary, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit       ! Logic unit.
  type(Type_Vector), intent(IN):: vec(:,:,:) ! Vector.
  integer(I_P)::                  err        ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(unit,iostat=err)vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Bin_Array3D

  function Write_Ascii_Array3D(unit,format,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for writing Type_Vector (ascii, Array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: unit       ! Logic unit
  character(*),      intent(IN):: format     ! Format specifier.
  type(Type_Vector), intent(IN):: vec(:,:,:) ! Vector.
  integer(I_P)::                  err        ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    write(unit,*,iostat=err)vec
  case default
    write(unit,adjustl(trim(format)),iostat=err)vec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Write_Ascii_Array3D

  ! read
  function Read_Bin_Scalar(unit,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Vector (binary, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit  ! Logic unit.
  type(Type_Vector), intent(INOUT):: vec   ! Vector.
  integer(I_P)::                     err   ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Scalar

  function Read_Ascii_Scalar(unit,format,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Vector (ascii, scalar).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit   ! Logic unit.
  character(*),      intent(IN)::    format ! Format specifier.
  type(Type_Vector), intent(INOUT):: vec    ! Vector.
  integer(I_P)::                     err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)vec
  case default
    read(unit,adjustl(trim(format)),iostat=err)vec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Scalar

  function Read_Bin_Array1D(unit,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Vector (binary, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit   ! Logic unit.
  type(Type_Vector), intent(INOUT):: vec(:) ! Vector.
  integer(I_P)::                     err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array1D

  function Read_Ascii_Array1D(unit,format,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Vector (ascii, array 1D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit   ! logic unit
  character(*),      intent(IN)::    format ! format specifier
  type(Type_Vector), intent(INOUT):: vec(:) ! Vector.
  integer(I_P)::                     err    ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)vec
  case default
    read(unit,adjustl(trim(format)),iostat=err)vec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array1D

  function Read_Bin_Array2D(unit,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Vector (binary, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit     ! Logic unit.
  type(Type_Vector), intent(INOUT):: vec(:,:) ! Vector.
  integer(I_P)::                     err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array2D

  function Read_Ascii_Array2D(unit,format,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Vector (ascii, array 2D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit     ! Logic unit.
  character(*),      intent(IN)::    format   ! Format specifier.
  type(Type_Vector), intent(INOUT):: vec(:,:) ! Vector.
  integer(I_P)::                     err      ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)vec
  case default
    read(unit,adjustl(trim(format)),iostat=err)vec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array2D

  function Read_Bin_Array3D(unit,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for reading Type_Vector (binary, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit       ! Logic unit.
  type(Type_Vector), intent(INOUT):: vec(:,:,:) ! Vector.
  integer(I_P)::                     err        ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  read(unit,iostat=err)vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Bin_Array3D

  function Read_Ascii_Array3D(unit,format,vec) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for reading Type_Vector (ascii, array 3D).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::    unit       ! Logic unit.
  character(*),      intent(IN)::    format     ! Format specifier.
  type(Type_Vector), intent(INOUT):: vec(:,:,:) ! Vector.
  integer(I_P)::                     err        ! Error traping flag: 0 no errors, >0 error occours.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(adjustl(trim(format)))
  case('*')
    read(unit,*,iostat=err)vec
  case default
    read(unit,adjustl(trim(format)),iostat=err)vec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Read_Ascii_Array3D

  ! Assignment (=)
#ifdef r16p
  elemental subroutine assign_ScalR16P(vec,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R16P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(INOUT):: vec
  real(R16P),        intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vec%x = real(scal,R_P)
  vec%y = real(scal,R_P)
  vec%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR16P
#endif

  elemental subroutine assign_ScalR8P(vec,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R8P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(INOUT):: vec
  real(R8P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vec%x = real(scal,R_P)
  vec%y = real(scal,R_P)
  vec%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR8P

  elemental subroutine assign_ScalR4P(vec,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (real R4P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(INOUT):: vec
  real(R4P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vec%x = real(scal,R_P)
  vec%y = real(scal,R_P)
  vec%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR4P

  elemental subroutine assign_ScalI8P(vec,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I8P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(INOUT):: vec
  integer(I8P),      intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vec%x = real(scal,R_P)
  vec%y = real(scal,R_P)
  vec%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI8P

  elemental subroutine assign_ScalI4P(vec,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I4P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(INOUT):: vec
  integer(I4P),      intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vec%x = real(scal,R_P)
  vec%y = real(scal,R_P)
  vec%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI4P

  elemental subroutine assign_ScalI2P(vec,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I2P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(INOUT):: vec
  integer(I2P),      intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vec%x = real(scal,R_P)
  vec%y = real(scal,R_P)
  vec%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI2P

  elemental subroutine assign_ScalI1P(vec,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Subroutine for assignment between a scalar (integer I1P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(INOUT):: vec
  integer(I1P),      intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vec%x = real(scal,R_P)
  vec%y = real(scal,R_P)
  vec%z = real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI1P

  ! Multiplication (*)
  elemental function vec_mul_vec(vec1,vec2) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply (by components) vectors.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec1 ! First vector.
  type(Type_Vector), intent(IN)::  vec2 ! Second vector.
  type(Type_Vector)::              mul  ! Resulting vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = vec1%x * vec2%x
  mul%y = vec1%y * vec2%y
  mul%z = vec1%z * vec2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_mul_vec

#ifdef r16p
  elemental function ScalR16P_mul_vec(scal,vec) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R16P) for vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_mul_vec

  elemental function vec_mul_ScalR16P(vec,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply vec for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R16P),        intent(IN)::  scal
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_mul_ScalR16P
#endif

  elemental function ScalR8P_mul_vec(scal,vec) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R8P) for vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_mul_vec

  elemental function vec_mul_ScalR8P(vec,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply vec for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R8P),         intent(IN)::  scal
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_mul_ScalR8P

  elemental function ScalR4P_mul_vec(scal,vec) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (real R4P) for vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_mul_vec

  elemental function vec_mul_ScalR4P(vec,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply vec for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R4P),         intent(IN)::  scal
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_mul_ScalR4P

  elemental function ScalI8P_mul_vec(scal,vec) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I8P) for vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_mul_vec

  elemental function vec_mul_ScalI8P(vec,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply vec for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I8P),      intent(IN)::  scal
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_mul_ScalI8P

  elemental function ScalI4P_mul_vec(scal,vec) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I4P) for vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_mul_vec

  elemental function vec_mul_ScalI4P(vec,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply vec for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I4P),      intent(IN)::  scal
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_mul_ScalI4P

  elemental function ScalI2P_mul_vec(scal,vec) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I2P) for vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_mul_vec

  elemental function vec_mul_ScalI2P(vec,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply vec for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I2P),      intent(IN)::  scal
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_mul_ScalI2P

  elemental function ScalI1P_mul_vec(scal,vec) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply scalar (integer I1P) for vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_mul_vec

  elemental function vec_mul_ScalI1P(vec,scal) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for multiply vec for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I1P),      intent(IN)::  scal
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = real(scal,R_P) * vec%x
  mul%y = real(scal,R_P) * vec%y
  mul%z = real(scal,R_P) * vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_mul_ScalI1P

  ! Division (/)
  elemental function vec_div_vec(vec1,vec2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide vec for vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec1
  type(Type_Vector), intent(IN)::  vec2
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = vec1%x / vec2%x
  div%y = vec1%y / vec2%y
  div%z = vec1%z / vec2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_div_vec

#ifdef r16p
  elemental function vec_div_ScalR16P(vec,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide vec for scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R16P),        intent(IN)::  scal
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = vec%x / real(scal,R_P)
  div%y = vec%y / real(scal,R_P)
  div%z = vec%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_div_ScalR16P
#endif

  elemental function vec_div_ScalR8P(vec,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide vec for scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R8P),         intent(IN)::  scal
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = vec%x / real(scal,R_P)
  div%y = vec%y / real(scal,R_P)
  div%z = vec%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_div_ScalR8P

  elemental function vec_div_ScalR4P(vec,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide vec for scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R4P),         intent(IN)::  scal
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = vec%x / real(scal,R_P)
  div%y = vec%y / real(scal,R_P)
  div%z = vec%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_div_ScalR4P

  elemental function vec_div_ScalI8P(vec,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide vec for scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I8P),      intent(IN)::  scal
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = vec%x / real(scal,R_P)
  div%y = vec%y / real(scal,R_P)
  div%z = vec%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_div_ScalI8P

  elemental function vec_div_ScalI4P(vec,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide vec for scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I4P),      intent(IN)::  scal
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = vec%x / real(scal,R_P)
  div%y = vec%y / real(scal,R_P)
  div%z = vec%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_div_ScalI4P

  elemental function vec_div_ScalI2P(vec,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide vec for scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I2P),      intent(IN)::  scal
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = vec%x / real(scal,R_P)
  div%y = vec%y / real(scal,R_P)
  div%z = vec%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_div_ScalI2P

  elemental function vec_div_ScalI1P(vec,scal) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for divide vec for scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I1P),      intent(IN)::  scal
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = vec%x / real(scal,R_P)
  div%y = vec%y / real(scal,R_P)
  div%z = vec%z / real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_div_ScalI1P

  ! Sum (+)
  elemental function positive_vec(vec) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary + to an vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pos%x =  + vec%x
  pos%y =  + vec%y
  pos%z =  + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_vec

  elemental function vec_sum_vec(vec1,vec2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec1
  type(Type_Vector), intent(IN)::  vec2
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = vec1%x + vec2%x
  summ%y = vec1%y + vec2%y
  summ%z = vec1%z + vec2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sum_vec

#ifdef r16p
  elemental function ScalR16P_sum_vec(scal,vec) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R16P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sum_vec

  elemental function vec_sum_ScalR16P(vec,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R16P),        intent(IN)::  scal
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sum_ScalR16P
#endif

  elemental function ScalR8P_sum_vec(scal,vec) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R8P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sum_vec

  elemental function vec_sum_ScalR8P(vec,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R8P),         intent(IN)::  scal
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sum_ScalR8P

  elemental function ScalR4P_sum_vec(scal,vec) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (real R4P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sum_vec

  elemental function vec_sum_ScalR4P(vec,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R4P),         intent(IN)::  scal
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sum_ScalR4P

  elemental function ScalI8P_sum_vec(scal,vec) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I8P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sum_vec

  elemental function vec_sum_ScalI8P(vec,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I8P),      intent(IN)::  scal
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sum_ScalI8P

  elemental function ScalI4P_sum_vec(scal,vec) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I4P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sum_vec

  elemental function vec_sum_ScalI4P(vec,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I4P),      intent(IN)::  scal
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sum_ScalI4P

  elemental function ScalI2P_sum_vec(scal,vec) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I2P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sum_vec

  elemental function vec_sum_ScalI2P(vec,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I2P),      intent(IN)::  scal
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sum_ScalI2P

  elemental function ScalI1P_sum_vec(scal,vec) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum scalar (integer I1P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sum_vec

  elemental function vec_sum_ScalI1P(vec,scal) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for sum vec and scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I1P),      intent(IN)::  scal
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = real(scal,R_P) + vec%x
  summ%y = real(scal,R_P) + vec%y
  summ%z = real(scal,R_P) + vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sum_ScalI1P

  ! Subtraction (-)
  elemental function negative_vec(vec) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for applay unary - to an vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neg%x =  - vec%x
  neg%y =  - vec%y
  neg%z =  - vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_vec

  elemental function vec_sub_vec(vec1,vec2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec1
  type(Type_Vector), intent(IN)::  vec2
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = vec1%x - vec2%x
  sub%y = vec1%y - vec2%y
  sub%z = vec1%z - vec2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sub_vec

#ifdef r16p
  elemental function ScalR16P_sub_vec(scal,vec) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R16P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - vec%x
  sub%y = real(scal,R_P) - vec%y
  sub%z = real(scal,R_P) - vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR16P_sub_vec

  elemental function vec_sub_ScalR16P(vec,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and scalar (real R16P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R16P),        intent(IN)::  scal
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = vec%x - real(scal,R_P)
  sub%y = vec%y - real(scal,R_P)
  sub%z = vec%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sub_ScalR16P
#endif

  elemental function ScalR8P_sub_vec(scal,vec) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R8P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - vec%x
  sub%y = real(scal,R_P) - vec%y
  sub%z = real(scal,R_P) - vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR8P_sub_vec

  elemental function vec_sub_ScalR8P(vec,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and scalar (real R8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R8P),         intent(IN)::  scal
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = vec%x - real(scal,R_P)
  sub%y = vec%y - real(scal,R_P)
  sub%z = vec%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sub_ScalR8P

  elemental function ScalR4P_sub_vec(scal,vec) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (real R4P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - vec%x
  sub%y = real(scal,R_P) - vec%y
  sub%z = real(scal,R_P) - vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalR4P_sub_vec

  elemental function vec_sub_ScalR4P(vec,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and scalar (real R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  real(R4P),         intent(IN)::  scal
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = vec%x - real(scal,R_P)
  sub%y = vec%y - real(scal,R_P)
  sub%z = vec%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sub_ScalR4P

  elemental function ScalI8P_sub_vec(scal,vec) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I8P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - vec%x
  sub%y = real(scal,R_P) - vec%y
  sub%z = real(scal,R_P) - vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI8P_sub_vec

  elemental function vec_sub_ScalI8P(vec,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and scalar (integer I8P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I8P),      intent(IN)::  scal
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = vec%x - real(scal,R_P)
  sub%y = vec%y - real(scal,R_P)
  sub%z = vec%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sub_ScalI8P

  elemental function ScalI4P_sub_vec(scal,vec) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I4P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - vec%x
  sub%y = real(scal,R_P) - vec%y
  sub%z = real(scal,R_P) - vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI4P_sub_vec

  elemental function vec_sub_ScalI4P(vec,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and scalar (integer I4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I4P),      intent(IN)::  scal
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = vec%x - real(scal,R_P)
  sub%y = vec%y - real(scal,R_P)
  sub%z = vec%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sub_ScalI4P

  elemental function ScalI2P_sub_vec(scal,vec) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I2P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - vec%x
  sub%y = real(scal,R_P) - vec%y
  sub%z = real(scal,R_P) - vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI2P_sub_vec

  elemental function vec_sub_ScalI2P(vec,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and scalar (integer I2P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I2P),      intent(IN)::  scal
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = vec%x - real(scal,R_P)
  sub%y = vec%y - real(scal,R_P)
  sub%z = vec%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sub_ScalI2P

  elemental function ScalI1P_sub_vec(scal,vec) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract scalar (integer I1P) and vec.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN)::  scal
  type(Type_Vector), intent(IN)::  vec
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = real(scal,R_P) - vec%x
  sub%y = real(scal,R_P) - vec%y
  sub%z = real(scal,R_P) - vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ScalI1P_sub_vec

  elemental function vec_sub_ScalI1P(vec,scal) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!Function for subtract vec and scalar (integer I1P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN)::  vec
  integer(I1P),      intent(IN)::  scal
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = vec%x - real(scal,R_P)
  sub%y = vec%y - real(scal,R_P)
  sub%z = vec%z - real(scal,R_P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_sub_ScalI1P

  ! Conditional operators
  ! /=
  elemental function vec_not_eq_vec(vec1,vec2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is /= with respect the normL2 of vec2 or if the directions of vec1 and
  !!vec2 are different, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1    ! First vector.
  type(Type_Vector), intent(IN):: vec2    ! Second vector.
  logical::                       compare ! The result of the comparison.
  type(Type_Vector)::             n1,n2   ! Normalizations of vec1 and vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec1)/=normL2(vec2))
  if (.not.compare) then ! the normL2 are the same, checking the directions
    n1 = normalize(vec1)
    n2 = normalize(vec2)
    compare = ((n1%x/=n2%x).OR.(n1%y/=n2%y).OR.(n1%z/=n2%z))
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_not_eq_vec

#ifdef r16p
  elemental function R16P_not_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_not_eq_vec

  elemental function vec_not_eq_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_not_eq_R16P
#endif

  elemental function R8P_not_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_not_eq_vec

  elemental function vec_not_eq_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_not_eq_R8P

  elemental function R4P_not_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_not_eq_vec

  elemental function vec_not_eq_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_not_eq_R4P

  elemental function I8P_not_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_not_eq_vec

  elemental function vec_not_eq_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_not_eq_I8P

  elemental function I4P_not_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_not_eq_vec

  elemental function vec_not_eq_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_not_eq_I4P

  elemental function I2P_not_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_not_eq_vec

  elemental function vec_not_eq_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_not_eq_I2P

  elemental function I1P_not_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_not_eq_vec

  elemental function vec_not_eq_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_not_eq_I1P

  ! <
  elemental function vec_low_vec(vec1,vec2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the normL2 of vec2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1    ! First vector.
  type(Type_Vector), intent(IN):: vec2    ! Second vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec1)<normL2(vec2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_vec

#ifdef r16p
  elemental function R16P_low_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_low_vec

  elemental function vec_low_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_R16P
#endif

  elemental function R8P_low_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_low_vec

  elemental function vec_low_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_R8P

  elemental function R4P_low_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_low_vec

  elemental function vec_low_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_R4P

  elemental function I8P_low_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_low_vec

  elemental function vec_low_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_I8P

  elemental function I4P_low_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_low_vec

  elemental function vec_low_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_I4P

  elemental function I2P_low_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_low_vec

  elemental function vec_low_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_I2P

  elemental function I1P_low_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_low_vec

  elemental function vec_low_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_I1P

  ! <=
  elemental function vec_low_eq_vec(vec1,vec2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the normL2 of vec2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1    ! First vector.
  type(Type_Vector), intent(IN):: vec2    ! Second vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec1)<=normL2(vec2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_eq_vec

#ifdef r16p
  elemental function R16P_low_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_low_eq_vec

  elemental function vec_low_eq_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_eq_R16P
#endif

  elemental function R8P_low_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_low_eq_vec

  elemental function vec_low_eq_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_eq_R8P

  elemental function R4P_low_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_low_eq_vec

  elemental function vec_low_eq_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_eq_R4P

  elemental function I8P_low_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_low_eq_vec

  elemental function vec_low_eq_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_eq_I8P

  elemental function I4P_low_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_low_eq_vec

  elemental function vec_low_eq_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_eq_I4P

  elemental function I2P_low_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_low_eq_vec

  elemental function vec_low_eq_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_eq_I2P

  elemental function I1P_low_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_low_eq_vec

  elemental function vec_low_eq_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_low_eq_I1P

  ! ==
  elemental function vec_eq_vec(vec1,vec2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is = with respect the normL2 of vec2 and the directions of vec1 and vec2
  !!are the same, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1    ! First vector.
  type(Type_Vector), intent(IN):: vec2    ! Second vector.
  logical::                       compare ! The result of the comparison.
  type(Type_Vector)::             n1,n2   ! Normalizations of vec1 and vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec1)==normL2(vec2))
  if (compare) then ! the normL2 are the same, checking the directions
    n1 = normalize(vec1)
    n2 = normalize(vec2)
    compare = ((n1%x==n2%x).AND.(n1%y==n2%y).AND.(n1%z==n2%z))
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_eq_vec

#ifdef r16p
  elemental function R16P_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_eq_vec

  elemental function vec_eq_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_eq_R16P
#endif

  elemental function R8P_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_eq_vec

  elemental function vec_eq_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_eq_R8P

  elemental function R4P_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_eq_vec

  elemental function vec_eq_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_eq_R4P

  elemental function I8P_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_eq_vec

  elemental function vec_eq_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_eq_I8P

  elemental function I4P_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_eq_vec

  elemental function vec_eq_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_eq_I4P

  elemental function I2P_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_eq_vec

  elemental function vec_eq_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_eq_I2P

  elemental function I1P_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_eq_vec

  elemental function vec_eq_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_eq_I1P

  ! >=
  elemental function vec_great_eq_vec(vec1,vec2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the normL2 of vec2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1    ! First vector.
  type(Type_Vector), intent(IN):: vec2    ! Second vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec1)>=normL2(vec2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_eq_vec

#ifdef r16p
  elemental function R16P_great_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_great_eq_vec

  elemental function vec_great_eq_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_eq_R16P
#endif

  elemental function R8P_great_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_great_eq_vec

  elemental function vec_great_eq_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_eq_R8P

  elemental function R4P_great_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_great_eq_vec

  elemental function vec_great_eq_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_eq_R4P

  elemental function I8P_great_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_great_eq_vec

  elemental function vec_great_eq_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_eq_I8P

  elemental function I4P_great_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_great_eq_vec

  elemental function vec_great_eq_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_eq_I4P

  elemental function I2P_great_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_great_eq_vec

  elemental function vec_great_eq_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_eq_I2P

  elemental function I1P_great_eq_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_great_eq_vec

  elemental function vec_great_eq_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_eq_I1P

  ! >
  elemental function vec_great_vec(vec1,vec2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the normL2 of vec2, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1    ! First vector.
  type(Type_Vector), intent(IN):: vec2    ! Second vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec1)>normL2(vec2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_vec

#ifdef r16p
  elemental function R16P_great_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),        intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_great_vec

  elemental function vec_great_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_R16P
#endif

  elemental function R8P_great_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_great_vec

  elemental function vec_great_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_R8P

  elemental function R4P_great_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),         intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_great_vec

  elemental function vec_great_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_R4P

  elemental function I8P_great_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_great_vec

  elemental function vec_great_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_I8P

  elemental function I4P_great_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_great_vec

  elemental function vec_great_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_I4P

  elemental function I2P_great_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_great_vec

  elemental function vec_great_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_I2P

  elemental function I1P_great_vec(scal,vec) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),      intent(IN):: scal    ! Scalar.
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R_P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_great_vec

  elemental function vec_great_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! 3D vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R_P))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction vec_great_I1P

  ! cross product (.cross.)
  elemental function crossproduct(vec1,vec2) result(cross)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the cross product of 2 vectors (rectangular coordinates).
  !!\begin{equation}
  !!\vec V=\left({y_1 z_2 - z_1 y_2}\right)\vec i + \left({z_1 x_2 - x_1 z_2}\right)\vec j + \left({x_1 y_2 - y_1 x_2}\right)\vec k
  !!\end{equation}
  !!\noindent where $x_i$, $y_i$ and $z_i$ $i=1,2$ are the components of the vectors.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1  ! First vector.
  type(Type_Vector), intent(IN):: vec2  ! Second vector.
  type(Type_Vector)::             cross ! Cross product vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cross%x = (vec1%y*vec2%z) - (vec1%z*vec2%y)
  cross%y = (vec1%z*vec2%x) - (vec1%x*vec2%z)
  cross%z = (vec1%x*vec2%y) - (vec1%y*vec2%x)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction crossproduct

  ! dot product (.dot.)
  elemental function dotproduct(vec1,vec2) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the scalar (dot) product of 2 vectors.
  !!\begin{equation}
  !!{\rm D}= x_1 \cdot x_2 + y_1 \cdot y_2 + z_1 \cdot z_2
  !!\end{equation}
  !!\noindent where $x_i$, $y_i$ and $z_i$ $i=1,2$ are the components of the vectors.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1 ! First vector.
  type(Type_Vector), intent(IN):: vec2 ! Second vector.
  real(R_P)::                     dot  ! Dot product.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  dot = (vec1%x*vec2%x) + (vec1%y*vec2%y) + (vec1%z*vec2%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dotproduct

  ! sq_norm
  elemental function sq_norm(vec) result(sq)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the square of the norm of a 3D vector.
  !!\begin{equation}
  !!N = x^2  + y^2  + z^2
  !!\end{equation}
  !!\noindent where $x$, $y$ and $z$ are the components of the vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec ! Vector.
  real(R_P)::                     sq  ! Square of the Norm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sq = (vec%x*vec%x) + (vec%y*vec%y) + (vec%z*vec%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sq_norm

  ! normL2
  elemental function normL2(vec) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function computes the norm L2 of a 3D vector.
  !!\begin{equation}
  !!N = \sqrt {x^2  + y^2  + z^2 }
  !!\end{equation}
  !!\noindent where $x$, $y$ and $z$ are the components of the vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec  ! Vector.
  real(R_P)::                     norm ! Norm L2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  norm = sqrt((vec%x*vec%x) + (vec%y*vec%y) + (vec%z*vec%z))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normL2

  ! normalize
  elemental function normalize(vec) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function normalize a 3D vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec  ! Vector to be normalized.
  type(Type_Vector)::             norm ! Vector normalized.
  real(R_P)::                     nm   ! Norm L2 of 3D vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  nm = normL2(vec)
  if (nm < smallR_P) then
    nm = nm + smallR_P
  endif
  norm%x = vec%x/nm
  norm%y = vec%y/nm
  norm%z = vec%z/nm
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normalize

  ! othogonal & parallel
  elemental function parallel(vec1,vec2) result(paral)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function calculates the component of vec1 parallel to vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1  ! First vector.
  type(Type_Vector), intent(IN):: vec2  ! Second vector.
  type(Type_Vector)::             paral ! vector parallel to vec2 with the module of vec1.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  paral = (vec1.dot.normalize(vec2))*normalize(vec2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction parallel

  elemental function orthogonal(vec1,vec2) result(ortho)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function calculates the vector orthogonal to vec2 with the module of vec1.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1  ! First vector.
  type(Type_Vector), intent(IN):: vec2  ! Second vector.
  type(Type_Vector)::             ortho ! vector orthogonal to vec2 with the module of vec1.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ortho = vec1 - (vec1.paral.vec2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction orthogonal

  ! face normal
  elemental function face_normal4(norm,pt1,pt2,pt3,pt4) result(fnormal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function calculates the normal of the face defined by the 4 points pt1, pt2, pt3 and pt4. The convention for the points
  !!numeration is the following:
  !!                                     1.----------.2  The normal is calculated by the cross product of the diagonal d13 for
  !!                                      |          |   the diagonal d24: d13 x d24.
  !!                                      |          |   The normal is normalized if the variable 'norm' is passed.
  !!                                      |          |
  !!                                      |          |
  !!                                     4.----------.3
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(1),      intent(IN), optional:: norm            ! If 'norm' is passed as argument the normal is normalized.
  type(Type_Vector), intent(IN)::           pt1,pt2,pt3,pt4 ! Face points.
  type(Type_Vector)::                       fnormal         ! Face normal.
  type(Type_Vector)::                       d13,d24         ! Face diagonals.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d13 = pt3 - pt1
  d24 = pt4 - pt2
  if (present(norm)) then
    fnormal = normalize(d13.cross.d24)
  else
    fnormal = 0.5_R_P*(d13.cross.d24)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction face_normal4

  elemental function face_normal3(norm,pt1,pt2,pt3) result(fnormal)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function calculates the normal of the face defined by the 3 points pt1, pt2 and pt3. The convention for the points
  !!numeration is the following:
  !!                                     1.----.2  The normal is calculated by the cross product of the side s12 for the side s13:
  !!                                       \   |   s12 x s13.
  !!                                        \  |   The normal is normalized if the variable 'norm' is passed.
  !!                                         \ |
  !!                                          \|
  !!                                           .3
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(1),      intent(IN), optional:: norm        ! If 'norm' is passed as argument the normal is normalized.
  type(Type_Vector), intent(IN)::           pt1,pt2,pt3 ! Face points.
  type(Type_Vector)::                       fnormal     ! Face normal.
  type(Type_Vector)::                       s12,s13     ! Face diagonals.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  s12 = pt2 - pt1
  s13 = pt3 - pt1
  if (present(norm)) then
    fnormal = normalize(s12.cross.s13)
  else
    fnormal = 0.5_R_P*(s12.cross.s13)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction face_normal3
endmodule Data_Type_Vector
