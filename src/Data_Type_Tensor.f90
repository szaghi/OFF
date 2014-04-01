!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_TensorDerivedType Data_Type_Tensor
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_TensorInterface Data_Type_Tensor
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_TensorGlobalVarPar Data_Type_Tensor
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_TensorPublicProcedure Data_Type_Tensor
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_TensorPrivateProcedure Data_Type_Tensor
!> @}

!>This module contains the definition of Type_Tensor and its procedures.
!>This derived type is useful for manipulating second order tensors in 3D space. The components of the tensors
!>are derived type of Type_Vector. The components are defined in a three-dimensional cartesian frame of reference.
!>All the tensorial math procedures (cross, dot products, normL2...) assume a three-dimensional cartesian frame of reference.
!> @note The operators of assignment (=), multiplication (*), division (/), sum (+) and subtraction (-) have been overloaded.
!> Furthermore the \em dot, \em double \em dot and \em diadic products have been defined.
!> Therefore this module provides a far-complete algebra based on Type_Tensor derived type. This algebra simplifies the
!> tensorial operations of Partial Differential Equations (PDE) systems.
!> @todo \b DocComplete: Complete the documentation of internal procedures
module Data_Type_Tensor
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                             ! Integers and reals precision definition.
USE Data_Type_Vector, sq_norm_vec => sq_norm ! Definition of type Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: unity
public:: sq_norm
public:: normL2
public:: normalize
public:: transpose
public:: determinant
public:: invert,invertible
public:: write_tensor,read_tensor
public:: operator (*)
public:: operator (/)
public:: operator (+)
public:: operator (-)
public:: operator (.ddot.)
public:: operator (.dot.)
public:: operator (.diad.)
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type defining tensors.
!> @note The components of the tensor are 3 vectors arranged as following: \n
!> \f$ T_{ij}=\left[{\begin{array}{*{20}{c}}t_{11}&t_{12}&t_{13}\\ t_{21}&t_{22}&t_{23}\\ t_{31}&t_{32}&t_{33}\end{array}}\right]
!> =\left[{\begin{array}{*{20}{c}} {x\% x}&{x\% y}&{x\% z}\\ {y\% x}&{y\% y}&{y\% z}\\ {z\% x}&{z\% y}&{z\% z}\end{array}}\right]\f$
!> @ingroup Data_Type_TensorDerivedType
type, public:: Type_Tensor
  type(Type_Vector):: x !< Cartesian vector component in x direction.
  type(Type_Vector):: y !< Cartesian vector component in y direction.
  type(Type_Vector):: z !< Cartesian vector component in z direction.
  contains
    procedure:: set                         ! Procedure for setting tensor components.
    procedure:: sq_norm => sq_norm_ten      ! Procedure for computing the square of the norm of a tensor.
    procedure:: normL2 => normL2_ten        ! Procedure for computing the norm L2 of a tensor.
    procedure:: normalize => normalize_self ! Procedure for normalizing a tensor.
    procedure:: transpose => transpose_self ! Procedure for transposing a tensor.
    procedure:: determinant                 ! Procedure for computing the determinant of a tensor.
    procedure:: invert => invert_self       ! Procedure for inverting a tensor.
    procedure:: invertible                  ! Procedure for checking if a tensor is invertible.
    procedure:: rotox                       ! Procedure for computing the rotation tensor along x axis.
    procedure:: rotoy                       ! Procedure for computing the rotation tensor along y axis.
    procedure:: rotoz                       ! Procedure for computing the rotation tensor along z axis.
    procedure:: rotou                       ! Procedure for computing the rotation tensor along a vector axis.
    ! operators overloading
    generic:: assignment(=) => assign_Vec,                   &
#ifdef r16p
                               assign_ScalR16P,              &
#endif
                               assign_ScalR8P,assign_ScalR4P,&
                               assign_ScalI8P,assign_ScalI4P,assign_ScalI2P,assign_ScalI1P
    ! private procedures
    procedure, pass(ten ), private:: assign_Vec
#ifdef r16p
    procedure, pass(ten ), private:: assign_ScalR16P
#endif
    procedure, pass(ten ), private:: assign_ScalR8P
    procedure, pass(ten ), private:: assign_ScalR4P
    procedure, pass(ten ), private:: assign_ScalI8P
    procedure, pass(ten ), private:: assign_ScalI4P
    procedure, pass(ten ), private:: assign_ScalI2P
    procedure, pass(ten ), private:: assign_ScalI1P
endtype Type_Tensor
!> @ingroup Data_Type_TensorGlobalVarPar
!> @{
type(Type_Tensor), parameter:: unity = Type_Tensor(ex,ey,ez) !< Unity (identity) tensor
                                                             !< (see \ref data_type_tensor::type_tensor "definition").
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Square norm function \em sq_norm overloading.
!> The function \em sq_norm defined for Type_Vector is overloaded for handling also Type_Tensor.
!> @ingroup Data_Type_TensorInterface
interface sq_norm
  module procedure sq_norm_ten
endinterface
!> @brief L2 norm function \em normL2 overloading.
!> The function \em normL2 defined for Type_Vector is overloaded for handling also Type_Tensor.
!> @ingroup Data_Type_TensorInterface
interface normL2
  module procedure normL2,normL2_ten
endinterface
!> @brief Normalize function \em normalize overloading.
!> The function \em normalize defined for Type_Vector is overloaded for handling also Type_Tensor.
!> @ingroup Data_Type_TensorInterface
interface normalize
  module procedure normalize,normalize_ten
endinterface
!> @brief Transpose function \em transpose overloading.
!> The built in function \em transpose defined for rank 2 arrays is overloaded for handling also Type_Tensor.
!> @ingroup Data_Type_TensorInterface
interface transpose
  module procedure transpose_ten
endinterface
!> @brief Multiplication operator (*) overloading.
!> @note The admissible multiplications are:
!>       - Type_Tensor * Type_Tensor: each component of first tensor variable (ten1) is multiplied for the
!>         corresponding component of the second one (ten2), i.e. \n
!>         \f$ {\rm result\%x = ten1\%x*ten2\%x} \f$ \n
!>         \f$ {\rm result\%y = ten1\%y*ten2\%y} \f$ \n
!>         \f$ {\rm result\%z = ten1\%z*ten2\%z} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) * Type_Tensor: each component of
!>         Type_Tensor is multiplied for the scalar, i.e. \n
!>         \f$ {\rm result\%x = ten\%x*scalar} \f$ \n
!>         \f$ {\rm result\%y = ten\%y*scalar} \f$ \n
!>         \f$ {\rm result\%z = ten\%z*scalar} \f$ \n
!>       - Type_Tensor * scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Tensor is multiplied for the scalar, i.e. \n
!>         \f$ {\rm result\%x = ten\%x*scalar} \f$ \n
!>         \f$ {\rm result\%y = ten\%y*scalar} \f$ \n
!>         \f$ {\rm result\%z = ten\%z*scalar} \f$ \n
!> @ingroup Data_Type_TensorInterface
interface operator (*)
  module procedure ten_mul_ten
#ifdef r16p
  module procedure ScalR16P_mul_ten
#endif
  module procedure ScalR8P_mul_ten
  module procedure ScalR4P_mul_ten
  module procedure ScalI8P_mul_ten
  module procedure ScalI4P_mul_ten
  module procedure ScalI2P_mul_ten
  module procedure ScalI1P_mul_ten
#ifdef r16p
  module procedure ten_mul_ScalR16P
#endif
  module procedure ten_mul_ScalR8P
  module procedure ten_mul_ScalR4P
  module procedure ten_mul_ScalI8P
  module procedure ten_mul_ScalI4P
  module procedure ten_mul_ScalI2P
  module procedure ten_mul_ScalI1P
endinterface
!> @brief Division operator (/) overloading.
!> @note The admissible divisions are:
!>       - Type_Tensor / Type_Tensor: each component of first tensor variable (ten1) is divided for the
!>         corresponding component of the second one (ten2), i.e. \n
!>         \f$ {\rm result\%x = \frac{ten1\%x}{ten2\%x}} \f$ \n
!>         \f$ {\rm result\%y = \frac{ten1\%y}{ten2\%y}} \f$ \n
!>         \f$ {\rm result\%z = \frac{ten1\%z}{ten2\%z}} \f$ \n
!>       - Type_Tensor / scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Tensor is divided for the scalar, i.e. \n
!>         \f$ {\rm result\%x = \frac{ten\%x}{scalar}} \f$ \n
!>         \f$ {\rm result\%y = \frac{ten\%y}{scalar}} \f$ \n
!>         \f$ {\rm result\%z = \frac{ten\%z}{scalar}} \f$ \n
!> @ingroup Data_Type_TensorInterface
interface operator (/)
  module procedure ten_div_ten
#ifdef r16p
  module procedure ten_div_ScalR16P
#endif
  module procedure ten_div_ScalR8P
  module procedure ten_div_ScalR4P
  module procedure ten_div_ScalI8P
  module procedure ten_div_ScalI4P
  module procedure ten_div_ScalI2P
  module procedure ten_div_ScalI1P
endinterface
!> @brief Sum operator (+) overloading.
!> @note The admissible summations are:
!>       - Type_Tensor + Type_Tensor: each component of first tensor variable (ten1) is summed with the
!>         corresponding component of the second one (ten2), i.e. \n
!>         \f$ {\rm result\%x = ten1\%x+ten2\%x} \f$ \n
!>         \f$ {\rm result\%y = ten1\%y+ten2\%y} \f$ \n
!>         \f$ {\rm result\%z = ten1\%z+ten2\%z} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) + Type_Tensor: each component of
!>         Type_Tensor is summed with the scalar, i.e. \n
!>         \f$ {\rm result\%x = ten\%x+scalar} \f$ \n
!>         \f$ {\rm result\%y = ten\%y+scalar} \f$ \n
!>         \f$ {\rm result\%z = ten\%z+scalar} \f$ \n
!>       - Type_Tensor + scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Tensor is summed with the scalar, i.e. \n
!>         \f$ {\rm result\%x = ten\%x+scalar} \f$ \n
!>         \f$ {\rm result\%y = ten\%y+scalar} \f$ \n
!>         \f$ {\rm result\%z = ten\%z+scalar} \f$ \n
!> @ingroup Data_Type_TensorInterface
interface operator (+)
  module procedure positive_ten
  module procedure ten_sum_ten
#ifdef r16p
  module procedure ScalR16P_sum_ten
#endif
  module procedure ScalR8P_sum_ten
  module procedure ScalR4P_sum_ten
  module procedure ScalI8P_sum_ten
  module procedure ScalI4P_sum_ten
  module procedure ScalI2P_sum_ten
  module procedure ScalI1P_sum_ten
#ifdef r16p
  module procedure ten_sum_ScalR16P
#endif
  module procedure ten_sum_ScalR8P
  module procedure ten_sum_ScalR4P
  module procedure ten_sum_ScalI8P
  module procedure ten_sum_ScalI4P
  module procedure ten_sum_ScalI2P
  module procedure ten_sum_ScalI1P
endinterface
!> @brief Subtraction operator (-) overloading.
!> @note The admissible subtractions are:
!>       - Type_Tensor - Type_Tensor: each component of first tensor variable (ten1) is subtracted with the
!>         corresponding component of the second one (ten2), i.e. \n
!>         \f$ {\rm result\%x = ten1\%x-ten2\%x} \f$ \n
!>         \f$ {\rm result\%y = ten1\%y-ten2\%y} \f$ \n
!>         \f$ {\rm result\%z = ten1\%z-ten2\%z} \f$ \n
!>       - scalar number (real or integer of any kinds defined in IR_Precision module) - Type_Tensor: each component of
!>         Type_Tensor is subtracted with the scalar, i.e. \n
!>         \f$ {\rm result\%x = scalar-ten\%x} \f$ \n
!>         \f$ {\rm result\%y = scalar-ten\%y} \f$ \n
!>         \f$ {\rm result\%z = scalar-ten\%z} \f$ \n
!>       - Type_Tensor - scalar number (real or integer of any kinds defined in IR_Precision module): each component of
!>         Type_Tensor is subtracted with the scalar, i.e. \n
!>         \f$ {\rm result\%x = ten\%x-scalar} \f$ \n
!>         \f$ {\rm result\%y = ten\%y-scalar} \f$ \n
!>         \f$ {\rm result\%z = ten\%z-scalar} \f$ \n
!> @ingroup Data_Type_TensorInterface
interface operator (-)
  module procedure negative_ten
  module procedure ten_sub_ten
#ifdef r16p
  module procedure ScalR16P_sub_ten
#endif
  module procedure ScalR8P_sub_ten
  module procedure ScalR4P_sub_ten
  module procedure ScalI8P_sub_ten
  module procedure ScalI4P_sub_ten
  module procedure ScalI2P_sub_ten
  module procedure ScalI1P_sub_ten
#ifdef r16p
  module procedure ten_sub_ScalR16P
#endif
  module procedure ten_sub_ScalR8P
  module procedure ten_sub_ScalR4P
  module procedure ten_sub_ScalI8P
  module procedure ten_sub_ScalI4P
  module procedure ten_sub_ScalI2P
  module procedure ten_sub_ScalI1P
endinterface
!> @brief Double dot product operator (.ddot.) definition.
!> @ingroup Data_Type_TensorInterface
interface operator (.ddot.)
  module procedure ddotproduct
endinterface
!> @brief Dot product operator (.dot.) definition.
!> @ingroup Data_Type_TensorInterface
interface operator (.dot.)
  module procedure ten_dot_vec,vec_dot_ten
endinterface
!> @brief Diadic product operator (.diad.) definition.
!> @ingroup Data_Type_TensorInterface
interface operator (.diad.)
  module procedure diadicproduct
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_TensorPublicProcedure
  !> @{
  !> @brief Function for computing the determinant of a tensor.
  !> @return \b det real(R8P) variable
  !> @note The determinant is computed according the following equation: \n
  !> \f$ \det  = \left| {\begin{array}{*{20}{c}} {x\% x}&{x\% y}&{x\% z}\\ {y\% x}&{y\% y}&{y\% z}\\ {z\% x}&{z\% y}&{z\% z}
  !> \end{array}} \right| = \f$ \n
  !> \f$=x\%x(z\%z\cdot y\%y-z\%y\cdot y\%z)-\f$
  !> \f$ y\%x(z\%z\cdot x\%y-z\%y\cdot x\%z)+\f$
  !> \f$ z\%x(y\%z\cdot x\%y-y\%y\cdot x\%z) \f$
  elemental function determinant(ten) result(det)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(IN):: ten !< Tensor.
  real(R8P)::                      det !< Determinant.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  det=ten%x%x*(ten%z%z*ten%y%y-ten%z%y*ten%y%z)-ten%y%x*(ten%z%z*ten%x%y-ten%z%y*ten%x%z)+ten%z%x*(ten%y%z*ten%x%y-ten%y%y*ten%x%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction determinant

  !> @brief Function for computing the inverse of a tensor.
  !> If the tensor is not invertible a null tensor is returned.
  !> @return \b inv type(Type_Tensor) variable
  !> @note The inverse tensor is computed according the following equation: \n
  !> \f$ T_{ij}^{ - 1} = \frac{1}{{\det (T)}}\left[ {\begin{array}{*{20}{c}}
  !> {z\%z\cdot y\%y-z\%y\cdot y\%z}&{-(z\%z\cdot x\%y-z\%y\cdot x\%z)}&{y\%z\cdot x\%y-y\%y\cdot x\%z}\\{-(z\%z\cdot y\%x-
  !> z\%x\cdot y\%z)}&{z\%z\cdot x\%x-z\%x\cdot x\%z}&{-(y\%z\cdot x\%x-y\%x\cdot x\%z)}\\{z\%y\cdot y\%x-z\%x\cdot y\%y}&{-(z\%y
  !> \cdot x\%x-z\%x\cdot x\%y)}&{y\%y\cdot x\%x-y\%x\cdot x\%y} \end{array}} \right]\f$ \n
  !> where det(T) is the determinant computed by means of the function determinant.
  elemental function invert(ten) result(inv)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten !< Tensor to be inverted.
  type(Type_Tensor)::             inv !< Tensor inverted.
  real(R8P)::                     det !< Determinant and 1/Determinant.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  det = determinant(ten)
  if (det/=0._R8P) then
    det = 1._R8P/det
   inv%x=det*( (ten%z%z*ten%y%y-ten%z%y*ten%y%z)*ex - (ten%z%z*ten%x%y-ten%z%y*ten%x%z)*ey + (ten%y%z*ten%x%y-ten%y%y*ten%x%z)*ez)
   inv%y=det*(-(ten%z%z*ten%y%x-ten%z%x*ten%y%z)*ex + (ten%z%z*ten%x%x-ten%z%x*ten%x%z)*ey - (ten%y%z*ten%x%x-ten%y%x*ten%x%z)*ez)
   inv%z=det*( (ten%z%y*ten%y%x-ten%z%x*ten%y%y)*ex - (ten%z%y*ten%x%x-ten%z%x*ten%x%y)*ey + (ten%y%y*ten%x%x-ten%y%x*ten%x%y)*ez)
  else
    inv%x=0._R8P
    inv%y=0._R8P
    inv%z=0._R8P
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction invert

  !> @brief Function for checking if a tensor is invertible (determinant /=0, not singular tensor).
  !> @return \b inv logical variable
  elemental function invertible(ten) result(inv)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(IN):: ten !< Tensor to be inverted.
  logical::                        inv !< True if the tensor is not singular.
  real(R8P)::                      det !< Determinant.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  det = determinant(ten) ; inv = (det/=0._R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction invertible

  !> @brief Function for writing Type_Tensor data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function write_tensor(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN), optional:: scalar         !< Scalar tensor data.
  type(Type_Tensor), intent(IN), optional:: array1D(:)     !< One dimensional array tensor data.
  type(Type_Tensor), intent(IN), optional:: array2D(:,:)   !< Two dimensional array tensor data.
  type(Type_Tensor), intent(IN), optional:: array3D(:,:,:) !< Three dimensional array tensor data.
  character(*),      intent(IN), optional:: format         !< Format specifier.
  integer(I4P),      intent(IN)::           unit           !< Logic unit.
  integer(I_P)::                            err            !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        write(unit,*,iostat=err)scalar
      elseif (present(array1D)) then
        write(unit,*,iostat=err)array1D
      elseif (present(array2D)) then
        write(unit,*,iostat=err)array2D
      elseif (present(array3D)) then
        write(unit,*,iostat=err)array3D
      endif
    case default
      if (present(scalar)) then
        write(unit,adjustl(trim(format)),iostat=err)scalar
      elseif (present(array1D)) then
        write(unit,adjustl(trim(format)),iostat=err)array1D
      elseif (present(array2D)) then
        write(unit,adjustl(trim(format)),iostat=err)array2D
      elseif (present(array3D)) then
        write(unit,adjustl(trim(format)),iostat=err)array3D
      endif
    endselect
  else
    if (present(scalar)) then
      write(unit,iostat=err)scalar
    elseif (present(array1D)) then
      write(unit,iostat=err)array1D
    elseif (present(array2D)) then
      write(unit,iostat=err)array2D
    elseif (present(array3D)) then
      write(unit,iostat=err)array3D
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction write_tensor

  !> @brief Function for reading Type_Tensor data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function read_tensor(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(INOUT), optional:: scalar         !< Scalar tensor data.
  type(Type_Tensor), intent(INOUT), optional:: array1D(:)     !< One dimensional array tensor data.
  type(Type_Tensor), intent(INOUT), optional:: array2D(:,:)   !< Two dimensional array tensor data.
  type(Type_Tensor), intent(INOUT), optional:: array3D(:,:,:) !< Three dimensional array tensor data.
  character(*),      intent(IN),    optional:: format         !< Format specifier.
  integer(I4P),      intent(IN)::              unit           !< Logic unit.
  integer(I_P)::                               err            !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(format)) then
    select case(adjustl(trim(format)))
    case('*')
      if (present(scalar)) then
        read(unit,*,iostat=err)scalar
      elseif (present(array1D)) then
        read(unit,*,iostat=err)array1D
      elseif (present(array2D)) then
        read(unit,*,iostat=err)array2D
      elseif (present(array3D)) then
        read(unit,*,iostat=err)array3D
      endif
    case default
      if (present(scalar)) then
        read(unit,adjustl(trim(format)),iostat=err)scalar
      elseif (present(array1D)) then
        read(unit,adjustl(trim(format)),iostat=err)array1D
      elseif (present(array2D)) then
        read(unit,adjustl(trim(format)),iostat=err)array2D
      elseif (present(array3D)) then
        read(unit,adjustl(trim(format)),iostat=err)array3D
      endif
    endselect
  else
    if (present(scalar)) then
      read(unit,iostat=err)scalar
    elseif (present(array1D)) then
      read(unit,iostat=err)array1D
    elseif (present(array2D)) then
      read(unit,iostat=err)array2D
    elseif (present(array3D)) then
      read(unit,iostat=err)array3D
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction read_tensor
  !> @}

  !> @ingroup Data_Type_TensorPrivateProcedure
  !> @{
  !> Subroutine for setting components of Type_Tensor variable.
  elemental subroutine set(ten,x,y,z)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT)::        ten !< Tensor.
  type(Type_Vector),  intent(IN), optional:: x   !< Cartesian vector component in x direction.
  type(Type_Vector),  intent(IN), optional:: y   !< Cartesian vector component in y direction.
  type(Type_Vector),  intent(IN), optional:: z   !< Cartesian vector component in z direction.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(x)) ten%x = x
  if (present(y)) ten%y = y
  if (present(z)) ten%z = z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> @brief Function for computing the square of the norm of a tensor.
  !> The square norm if defined as \f$ N = x_x^2  + x_y^2  + x_z^2 + y_x^2 + y_y^2  + y_z^2 +...\f$.
  !> @return \b sq square norm
  elemental function sq_norm_ten(ten) result(sq)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(IN):: ten ! Tensor.
  real(R8P)::                      sq  ! Square of the Norm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sq = sq_norm_vec(ten%x) + sq_norm_vec(ten%y) + sq_norm_vec(ten%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sq_norm_ten

  !> @brief Function for computing the norm L2 of a tensor.
  !> The norm L2 if defined as \f$N = \sqrt {x_x^2  + x_y^2  + x_z^2 + y_y^2  + y_z^2 +...}\f$.
  !> @return \b norm norm L2
  elemental function normL2_ten(ten) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(IN):: ten  ! Tensor.
  real(R8P)::                      norm ! Norm L2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  norm = sqrt(sq_norm_vec(ten%x) + sq_norm_vec(ten%y) + sq_norm_vec(ten%z))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normL2_ten

  !> @brief Subroutine for normalizing a tensor.
  !> The normalization is made by means of norm L2. If the norm L2 of the tensor is less than the parameter smallR8P the
  !> normalization value is set to normL2(ten)+smallR8P.
  elemental subroutine normalize_self(ten)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten ! Tensor to be normalized.
  real(R8P)::                         nm  ! Norm L2 of tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  nm = normL2_ten(ten)
  if (nm < smallR8P) then
    nm = nm + smallR8P
  endif
  ten%x = ten%x/nm
  ten%y = ten%y/nm
  ten%z = ten%z/nm
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine normalize_self

  !> @brief Function for normalizing a tensor.
  !> The normalization is made by means of norm L2. If the norm L2 of the tensor is less than the parameter smallR8P the
  !> normalization value is set to normL2(ten)+smallR8P.
  !> @return \b norm normalized tensor
  elemental function normalize_ten(ten) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten  ! Tensor to be normalized.
  type(Type_Tensor)::             norm ! Tensor normalized.
  real(R8P)::                     nm   ! Norm L2 of tensor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  nm = normL2_ten(ten)
  if (nm < smallR8P) then
    nm = nm + smallR8P
  endif
  norm%x = ten%x/nm
  norm%y = ten%y/nm
  norm%z = ten%z/nm
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normalize_ten

  !> @brief Subroutine for transposing a tensor.
  !> The transpositions is done as following:
  !> @code
  !>       |x%x x%y x%z|           |x%x y%x z%x|
  !> ten = |y%x y%y y%z| => tran = |x%y y%y z%y|
  !>       |z%x z%y z%z|           |x%z y%z z%z|
  !> @endcode
  elemental subroutine transpose_self(ten)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten  !< Tensor to be transposed.
  type(Type_Tensor)::                 tran !< Temporary tensor transposed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! the following do not compile with gnu gfortran...
  !tran%x = ten.dot.ex
  !tran%y = ten.dot.ey
  !tran%z = ten.dot.ez
  ! used for correct compiling with gnu gfortran...
  select type(ten)
  type is(Type_Tensor)
    tran%x = ten_dot_vec(ten,ex)
    tran%y = ten_dot_vec(ten,ey)
    tran%z = ten_dot_vec(ten,ez)
  endselect
  ten%x = tran%x
  ten%y = tran%y
  ten%z = tran%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine transpose_self

  !> @brief Function for transposing a tensor.
  !> The transpositions is done as following:
  !> @code
  !>       |x%x x%y x%z|           |x%x y%x z%x|
  !> ten = |y%x y%y y%z| => tran = |x%y y%y z%y|
  !>       |z%x z%y z%z|           |x%z y%z z%z|
  !> @endcode
  elemental function transpose_ten(ten) result(tran)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Tensor), intent(IN):: ten  !< Tensor to be transposed.
  type(Type_Tensor)::             tran !< Tensor transposed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !tran%x = ten.dot.ex
  !tran%y = ten.dot.ey
  !tran%z = ten.dot.ez
  tran%x = ten_dot_vec(ten,ex)
  tran%y = ten_dot_vec(ten,ey)
  tran%z = ten_dot_vec(ten,ez)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction transpose_ten

  !> @brief Subroutine for computing the inverse of a tensor.
  !> If the tensor is not invertible a null tensor is returned.
  !> @return \b inv type(Type_Tensor) variable
  !> @note The inverse tensor is computed according the following equation: \n
  !> \f$ T_{ij}^{ - 1} = \frac{1}{{\det (T)}}\left[ {\begin{array}{*{20}{c}}
  !> {z\%z\cdot y\%y-z\%y\cdot y\%z}&{-(z\%z\cdot x\%y-z\%y\cdot x\%z)}&{y\%z\cdot x\%y-y\%y\cdot x\%z}\\{-(z\%z\cdot y\%x-
  !> z\%x\cdot y\%z)}&{z\%z\cdot x\%x-z\%x\cdot x\%z}&{-(y\%z\cdot x\%x-y\%x\cdot x\%z)}\\{z\%y\cdot y\%x-z\%x\cdot y\%y}&{-(z\%y
  !> \cdot x\%x-z\%x\cdot x\%y)}&{y\%y\cdot x\%x-y\%x\cdot x\%y} \end{array}} \right]\f$ \n
  !> where det(T) is the determinant computed by means of the function determinant.
  elemental subroutine invert_self(ten)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten !< Tensor to be inverted.
  type(Type_Tensor)::                 inv !< Temporary tensor inverted.
  real(R8P)::                         det !< Determinant and 1/Determinant.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  det = determinant(ten)
  if (det/=0._R8P) then
   det = 1._R8P/det
   inv%x=det*( (ten%z%z*ten%y%y-ten%z%y*ten%y%z)*ex - (ten%z%z*ten%x%y-ten%z%y*ten%x%z)*ey + (ten%y%z*ten%x%y-ten%y%y*ten%x%z)*ez)
   inv%y=det*(-(ten%z%z*ten%y%x-ten%z%x*ten%y%z)*ex + (ten%z%z*ten%x%x-ten%z%x*ten%x%z)*ey - (ten%y%z*ten%x%x-ten%y%x*ten%x%z)*ez)
   inv%z=det*( (ten%z%y*ten%y%x-ten%z%x*ten%y%y)*ex - (ten%z%y*ten%x%x-ten%z%x*ten%x%y)*ey + (ten%y%y*ten%x%x-ten%y%x*ten%x%y)*ez)
  else
    inv%x=0._R8P
    inv%y=0._R8P
    inv%z=0._R8P
  endif
  ten%x = inv%x
  ten%y = inv%y
  ten%z = inv%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine invert_self

  !> @brief Subroutine for computing the rotation tensor along x.
  elemental subroutine rotox(ten,ang)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten !< Rotating tensor.
  real(R8P),          intent(IN)::    ang !< Angle (radians) of rotation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call ten%x%set(x = 1._R8P,y = 0._R8P  ,z =  0._R8P  )
  call ten%y%set(x = 0._R8P,y = cos(ang),z = -sin(ang))
  call ten%z%set(x = 0._R8P,y = sin(ang),z =  cos(ang))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rotox

  !> @brief Subroutine for computing the rotation tensor along y.
  elemental subroutine rotoy(ten,ang)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten !< Rotating tensor.
  real(R8P),          intent(IN)::    ang !< Angle (radians) of rotation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call ten%x%set(x =  cos(ang),y = 0._R8P,z = sin(ang))
  call ten%y%set(x =  0._R8P  ,y = 1._R8P,z = 0._R8P  )
  call ten%z%set(x = -sin(ang),y = 0._R8P,z = cos(ang))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rotoy

  !> @brief Subroutine for computing the rotation tensor along z.
  elemental subroutine rotoz(ten,ang)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten !< Rotating tensor.
  real(R8P),          intent(IN)::    ang !< Angle (radians) of rotation.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call ten%x%set(x = cos(ang),y = -sin(ang),z = 0._R8P)
  call ten%y%set(x = sin(ang),y =  cos(ang),z = 0._R8P)
  call ten%z%set(x = 0._R8P  ,y =  0._R8P  ,z = 1._R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rotoz

  !> @brief Subroutine for computing the rotation tensor along a generic vector axis.
  elemental subroutine rotou(ten,u,ang)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten !< Rotating tensor.
  type(Type_Vector),  intent(IN)::    u   !< Vector axis.
  real(R8P),          intent(IN)::    ang !< Angle (radians) of rotation.
  type(Type_Vector)::                 n   !< Normalized vector axis.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  n = u
  call n%normalize
  call ten%x%set(x=cos(ang)+n%x*n%x*(1._R8P-cos(ang))    ,&
                 y=n%x*n%y*(1._R8P-cos(ang))-u%z*sin(ang),&
                 z=n%x*n%z*(1._R8P-cos(ang))+u%y*sin(ang))
  call ten%y%set(x=n%y*n%x*(1._R8P-cos(ang))+u%z*sin(ang),&
                 y=cos(ang)+n%y*n%y*(1._R8P-cos(ang))    ,&
                 z=n%y*n%z*(1._R8P-cos(ang))-u%x*sin(ang))
  call ten%z%set(x=n%z*n%x*(1._R8P-cos(ang))-u%y*sin(ang),&
                 y=n%z*n%y*(1._R8P-cos(ang))+u%x*sin(ang),&
                 z=cos(ang)+n%z*n%z*(1._R8P-cos(ang)))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rotou

  ! Assignment (=)
  !!Subroutine for assignment between a vector and ten.
  elemental subroutine assign_Vec(ten,vec)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten
  type(Type_Vector),  intent(IN)::    vec
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = vec
  ten%y = vec
  ten%z = vec
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_Vec

  !!Subroutine for assignment between a scalar (real R16P) and ten.
  elemental subroutine assign_ScalR16P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten
  real(R16P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R8P)
  ten%y = real(scal,R8P)
  ten%z = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR16P

  !!Subroutine for assignment between a scalar (real R8P) and ten.
  elemental subroutine assign_ScalR8P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten
  real(R8P),          intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R8P)
  ten%y = real(scal,R8P)
  ten%z = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR8P

  !!Subroutine for assignment between a scalar (real R4P) and ten.
  elemental subroutine assign_ScalR4P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten
  real(R4P),          intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R8P)
  ten%y = real(scal,R8P)
  ten%z = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR4P

  !!Subroutine for assignment between a scalar (integer I8P) and ten.
  elemental subroutine assign_ScalI8P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten
  integer(I8P),       intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R8P)
  ten%y = real(scal,R8P)
  ten%z = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI8P

  !!Subroutine for assignment between a scalar (integer I4P) and ten.
  elemental subroutine assign_ScalI4P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten
  integer(I4P),       intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R8P)
  ten%y = real(scal,R8P)
  ten%z = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI4P

  !!Subroutine for assignment between a scalar (integer I2P) and ten.
  elemental subroutine assign_ScalI2P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten
  integer(I2P),       intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R8P)
  ten%y = real(scal,R8P)
  ten%z = real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI2P

  !!Subroutine for assignment between a scalar (integer I1P) and ten.
  elemental subroutine assign_ScalI1P(ten,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Tensor), intent(INOUT):: ten
  integer(I1P),       intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ten%x = real(scal,R8P)
  ten%y = real(scal,R8P)
  ten%z = real(scal,R8P)
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_mul_ScalR16P

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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  mul%x = real(scal,R8P) * ten%x
  mul%y = real(scal,R8P) * ten%y
  mul%z = real(scal,R8P) * ten%z
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
  div%x = ten%x / real(scal,R8P)
  div%y = ten%y / real(scal,R8P)
  div%z = ten%z / real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_div_ScalR16P

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
  div%x = ten%x / real(scal,R8P)
  div%y = ten%y / real(scal,R8P)
  div%z = ten%z / real(scal,R8P)
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
  div%x = ten%x / real(scal,R8P)
  div%y = ten%y / real(scal,R8P)
  div%z = ten%z / real(scal,R8P)
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
  div%x = ten%x / real(scal,R8P)
  div%y = ten%y / real(scal,R8P)
  div%z = ten%z / real(scal,R8P)
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
  div%x = ten%x / real(scal,R8P)
  div%y = ten%y / real(scal,R8P)
  div%z = ten%z / real(scal,R8P)
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
  div%x = ten%x / real(scal,R8P)
  div%y = ten%y / real(scal,R8P)
  div%z = ten%z / real(scal,R8P)
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
  div%x = ten%x / real(scal,R8P)
  div%y = ten%y / real(scal,R8P)
  div%z = ten%z / real(scal,R8P)
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sum_ScalR16P

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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  summ%x = real(scal,R8P) + ten%x
  summ%y = real(scal,R8P) + ten%y
  summ%z = real(scal,R8P) + ten%z
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
  sub%x = real(scal,R8P) - ten%x
  sub%y = real(scal,R8P) - ten%y
  sub%z = real(scal,R8P) - ten%z
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
  sub%x = ten%x - real(scal,R8P)
  sub%y = ten%y - real(scal,R8P)
  sub%z = ten%z - real(scal,R8P)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction ten_sub_ScalR16P

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
  sub%x = real(scal,R8P) - ten%x
  sub%y = real(scal,R8P) - ten%y
  sub%z = real(scal,R8P) - ten%z
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
  sub%x = ten%x - real(scal,R8P)
  sub%y = ten%y - real(scal,R8P)
  sub%z = ten%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - ten%x
  sub%y = real(scal,R8P) - ten%y
  sub%z = real(scal,R8P) - ten%z
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
  sub%x = ten%x - real(scal,R8P)
  sub%y = ten%y - real(scal,R8P)
  sub%z = ten%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - ten%x
  sub%y = real(scal,R8P) - ten%y
  sub%z = real(scal,R8P) - ten%z
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
  sub%x = ten%x - real(scal,R8P)
  sub%y = ten%y - real(scal,R8P)
  sub%z = ten%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - ten%x
  sub%y = real(scal,R8P) - ten%y
  sub%z = real(scal,R8P) - ten%z
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
  sub%x = ten%x - real(scal,R8P)
  sub%y = ten%y - real(scal,R8P)
  sub%z = ten%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - ten%x
  sub%y = real(scal,R8P) - ten%y
  sub%z = real(scal,R8P) - ten%z
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
  sub%x = ten%x - real(scal,R8P)
  sub%y = ten%y - real(scal,R8P)
  sub%z = ten%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - ten%x
  sub%y = real(scal,R8P) - ten%y
  sub%z = real(scal,R8P) - ten%z
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
  sub%x = ten%x - real(scal,R8P)
  sub%y = ten%y - real(scal,R8P)
  sub%z = ten%z - real(scal,R8P)
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
  call dot%set(x=(ten%x.dot.vec),y=(ten%y.dot.vec),z=(ten%z.dot.vec))
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
  !call dot%set(x=((ten.dot.ex).dot.vec),y=((ten.dot.ey).dot.vec),z=((ten.dot.ez).dot.vec))
  call dot%set(x=((ten_dot_vec(ten,ex)).dot.vec),y=((ten_dot_vec(ten,ey)).dot.vec),z=((ten_dot_vec(ten,ez)).dot.vec))
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
  real(R8P)::                     ddot ! Double dot product.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !ddot = (ten1%x.dot.(ten2.dot.ex)) + (ten1%y.dot.(ten2.dot.ey)) + (ten1%z.dot.(ten2.dot.ez))
  ddot = (ten1%x.dot.(ten_dot_vec(ten2,ex))) + (ten1%y.dot.(ten_dot_vec(ten2,ey))) + (ten1%z.dot.(ten_dot_vec(ten2,ez)))
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
  !> @}
endmodule Data_Type_Tensor
