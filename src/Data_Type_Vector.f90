!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_VectorDerivedType Data_Type_Vector
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_VectorGlobalVarPar Data_Type_Vector
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_VectorInterface Data_Type_Vector
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_VectorPublicProcedure Data_Type_Vector
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_VectorPrivateProcedure Data_Type_Vector
!> @}

!> This module contains the definition of Type_Vector and its procedures.
!> This derived type is useful for manipulating vectors in 3D space. The components of the vectors are reals with
!> R_P kind as defined by the IR_Precision module. The components are defined in a three-dimensional cartesian frame of reference.
!> All the vectorial math procedures (cross, dot products, parallel...) assume a three-dimensional cartesian frame of reference.
!> @note The operators of assignment (=), multiplication (*), division (/), sum (+) and subtraction (-) have been overloaded.
!> Furthermore the \em dot and \em cross products have been defined.
!> Therefore this module provides a far-complete algebra based on Type_Vector derived type. This algebra simplifies the
!> vectorial operations of Partial Differential Equations (PDE) systems.
!> @todo \b DocComplete: Complete the documentation
module Data_Type_Vector
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision ! Integers and reals precision definition.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: ex,ey,ez
public:: sq_norm
public:: normL2
public:: normalize
public:: face_normal3,face_normal4
public:: write_vector,read_vector
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
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_VectorDerivedType
!> @{
!> Derived type defining vectors.
type, public:: Type_Vector
  real(R8P):: x = 0._R8P !< Cartesian component in x direction.
  real(R8P):: y = 0._R8P !< Cartesian component in y direction.
  real(R8P):: z = 0._R8P !< Cartesian component in z direction.
  contains
    procedure, non_overridable:: set                         ! Procedure for setting vector components.
    procedure, non_overridable:: pprint                      ! Procedure for printing vector components with a "pretty" format.
    procedure, non_overridable:: sq_norm                     ! Procedure for computing the square of the norm of a vector.
    procedure, non_overridable:: normL2                      ! Procedure for computing the norm L2 of a vector.
    procedure, non_overridable:: normalize => normalize_self ! Procedure for normalizing a vector.
endtype Type_Vector
!> Pointer of Type_Vector for creating array of pointers of Type_Vector.
type, public:: Type_Vector_Ptr
  type(Type_Vector), pointer:: p => null()
endtype Type_Vector_Ptr
!> @}
!> @ingroup Data_Type_VectorGlobalVarPar
!> @{
type(Type_Vector), parameter:: ex = Type_Vector(1._R8P,0._R8P,0._R8P) !< X direction versor
                                                                      !< (see \ref data_type_vector::type_vector "definition").
type(Type_Vector), parameter:: ey = Type_Vector(0._R8P,1._R8P,0._R8P) !< Y direction versor
                                                                      !< (see \ref data_type_vector::type_vector "definition").
type(Type_Vector), parameter:: ez = Type_Vector(0._R8P,0._R8P,1._R8P) !< Z direction versor
                                                                      !< (see \ref data_type_vector::type_vector "definition").
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Assignment operator (=) overloading.
!> @ingroup Data_Type_VectorInterface
interface assignment (=)
#ifdef r16p
  module procedure assign_ScalR16P
#endif
  module procedure assign_ScalR8P
  module procedure assign_ScalR4P
  module procedure assign_ScalI8P
  module procedure assign_ScalI4P
  module procedure assign_ScalI2P
  module procedure assign_ScalI1P
endinterface
!> @brief Multiplication operator (*) overloading.
!> @note The admissible multiplications are:
!>       - Type_Vector * Type_Vector: each component of first vector variable (vec1) is multiplied for the
!>         corresponding component of the second one (vec2), i.e. \n
!>         \f$ {\rm result\%x = vec1\%x*vec2\%x} \f$ \n
!>         \f$ {\rm result\%y = vec1\%y*vec2\%y} \f$ \n
!>         \f$ {\rm result\%z = vec1\%z*vec2\%z} \f$ \n
!>       - scalar number (real or integer of any kind defined in IR_Precision module) * Type_Vector: each component of
!>         Type_Vector is multiplied for the scalar, i.e. \n
!>         \f$ {\rm result\%x = vec\%x*scalar} \f$ \n
!>         \f$ {\rm result\%y = vec\%y*scalar} \f$ \n
!>         \f$ {\rm result\%z = vec\%z*scalar} \f$ \n
!>       - Type_Vector * scalar number (real or integer of any kind defined in IR_Precision module): each component of
!>         Type_Vector is multiplied for the scalar, i.e. \n
!>         \f$ {\rm result\%x = vec\%x*scalar} \f$ \n
!>         \f$ {\rm result\%y = vec\%y*scalar} \f$ \n
!>         \f$ {\rm result\%z = vec\%z*scalar} \f$ \n
!> @ingroup Data_Type_VectorInterface
interface operator (*)
  module procedure vec_mul_vec
#ifdef r16p
  module procedure ScalR16P_mul_vec
#endif
  module procedure ScalR8P_mul_vec
  module procedure ScalR4P_mul_vec
  module procedure ScalI8P_mul_vec
  module procedure ScalI4P_mul_vec
  module procedure ScalI2P_mul_vec
  module procedure ScalI1P_mul_vec
#ifdef r16p
  module procedure vec_mul_ScalR16P
#endif
  module procedure vec_mul_ScalR8P
  module procedure vec_mul_ScalR4P
  module procedure vec_mul_ScalI8P
  module procedure vec_mul_ScalI4P
  module procedure vec_mul_ScalI2P
  module procedure vec_mul_ScalI1P
endinterface
!> @brief Division operator (/) overloading.
!> @note The admissible divisions are:
!>       - Type_Vector / Type_Vector: each component of first vector variable (vec1) is divided for the
!>         corresponding component of the second one (vec2), i.e. \n
!>         \f$ {\rm result\%x = \frac{vec1\%x}{vec2\%x}} \f$ \n
!>         \f$ {\rm result\%y = \frac{vec1\%y}{vec2\%y}} \f$ \n
!>         \f$ {\rm result\%z = \frac{vec1\%z}{vec2\%z}} \f$ \n
!>       - Type_Vector / scalar number (real or integer of any kind defined in IR_Precision module): each component of
!>         Type_Vector is divided for the scalar, i.e. \n
!>         \f$ {\rm result\%x = \frac{vec\%x}{scalar}} \f$ \n
!>         \f$ {\rm result\%y = \frac{vec\%y}{scalar}} \f$ \n
!>         \f$ {\rm result\%z = \frac{vec\%z}{scalar}} \f$ \n
!> @ingroup Data_Type_VectorInterface
interface operator (/)
  module procedure vec_div_vec
#ifdef r16p
  module procedure vec_div_ScalR16P
#endif
  module procedure vec_div_ScalR8P
  module procedure vec_div_ScalR4P
  module procedure vec_div_ScalI8P
  module procedure vec_div_ScalI4P
  module procedure vec_div_ScalI2P
  module procedure vec_div_ScalI1P
endinterface
!> @brief Sum operator (+) overloading.
!> @note The admissible summations are:
!>       - Type_Vector + Type_Vector: each component of first Vector variable (vec1) is summed with the
!>         corresponding component of the second one (vec2), i.e. \n
!>         \f$ {\rm result\%x = vec1\%x+vec2\%x} \f$ \n
!>         \f$ {\rm result\%y = vec1\%y+vec2\%y} \f$ \n
!>         \f$ {\rm result\%z = vec1\%z+vec2\%z} \f$ \n
!>       - scalar number (real or integer of any kind defined in IR_Precision module) + Type_Vector: each component of
!>         Type_Vector is summed with the scalar, i.e. \n
!>         \f$ {\rm result\%x = vec\%x+scalar} \f$ \n
!>         \f$ {\rm result\%y = vec\%y+scalar} \f$ \n
!>         \f$ {\rm result\%z = vec\%z+scalar} \f$ \n
!>       - Type_Vector + scalar number (real or integer of any kind defined in IR_Precision module): each component of
!>         Type_Vector is summed with the scalar, i.e. \n
!>         \f$ {\rm result\%x = vec\%x+scalar} \f$ \n
!>         \f$ {\rm result\%y = vec\%y+scalar} \f$ \n
!>         \f$ {\rm result\%z = vec\%z+scalar} \f$ \n
!> @ingroup Data_Type_VectorInterface
interface operator (+)
  module procedure positive_vec
  module procedure vec_sum_vec
#ifdef r16p
  module procedure ScalR16P_sum_vec
#endif
  module procedure ScalR8P_sum_vec
  module procedure ScalR4P_sum_vec
  module procedure ScalI8P_sum_vec
  module procedure ScalI4P_sum_vec
  module procedure ScalI2P_sum_vec
  module procedure ScalI1P_sum_vec
#ifdef r16p
  module procedure vec_sum_ScalR16P
#endif
  module procedure vec_sum_ScalR8P
  module procedure vec_sum_ScalR4P
  module procedure vec_sum_ScalI8P
  module procedure vec_sum_ScalI4P
  module procedure vec_sum_ScalI2P
  module procedure vec_sum_ScalI1P
endinterface
!> @brief Subtraction operator (-) overloading.
!> @note The admissible subtractions are:
!>       - Type_Vector - Type_Vector: each component of first vector variable (vec1) is subtracted with the
!>         corresponding component of the second one (vec2), i.e. \n
!>         \f$ {\rm result\%x = vec1\%x-vec2\%x} \f$ \n
!>         \f$ {\rm result\%y = vec1\%y-vec2\%y} \f$ \n
!>         \f$ {\rm result\%z = vec1\%z-vec2\%z} \f$ \n
!>       - scalar number (real or integer of any kind defined in IR_Precision module) - Type_Vector: each component of
!>         Type_Vector is subtracted with the scalar, i.e. \n
!>         \f$ {\rm result\%x = scalar-vec\%x} \f$ \n
!>         \f$ {\rm result\%y = scalar-vec\%y} \f$ \n
!>         \f$ {\rm result\%z = scalar-vec\%z} \f$ \n
!>       - Type_Vector - scalar number (real or integer of any kind defined in IR_Precision module): each component of
!>         Type_Vector is subtracted with the scalar, i.e. \n
!>         \f$ {\rm result\%x = vec\%x-scalar} \f$ \n
!>         \f$ {\rm result\%y = vec\%y-scalar} \f$ \n
!>         \f$ {\rm result\%z = vec\%z-scalar} \f$ \n
!> @ingroup Data_Type_VectorInterface
interface operator (-)
  module procedure negative_vec
  module procedure vec_sub_vec
#ifdef r16p
  module procedure ScalR16P_sub_vec
#endif
  module procedure ScalR8P_sub_vec
  module procedure ScalR4P_sub_vec
  module procedure ScalI8P_sub_vec
  module procedure ScalI4P_sub_vec
  module procedure ScalI2P_sub_vec
  module procedure ScalI1P_sub_vec
#ifdef r16p
  module procedure vec_sub_ScalR16P
#endif
  module procedure vec_sub_ScalR8P
  module procedure vec_sub_ScalR4P
  module procedure vec_sub_ScalI8P
  module procedure vec_sub_ScalI4P
  module procedure vec_sub_ScalI2P
  module procedure vec_sub_ScalI1P
endinterface
!> @brief Not-equal-to boolean operator (/=) overloading.
!> @note The boolean comparison between two vectors is made on normL2 and direction, while the comparison between scalar number
!> (real or integer of any kind as defined in IR_Precision module) is made on only normL2 of the vector.
!> @ingroup Data_Type_VectorInterface
interface operator (/=)
  module procedure vec_not_eq_vec
#ifdef r16p
  module procedure R16P_not_eq_vec,vec_not_eq_R16P
#endif
  module procedure R8P_not_eq_vec,vec_not_eq_R8P
  module procedure R4P_not_eq_vec,vec_not_eq_R4P
  module procedure I8P_not_eq_vec,vec_not_eq_I8P
  module procedure I4P_not_eq_vec,vec_not_eq_I4P
  module procedure I2P_not_eq_vec,vec_not_eq_I2P
  module procedure I1P_not_eq_vec,vec_not_eq_I1P
endinterface
!> @brief Lower-than boolean operator (<) overloading.
!> @note The boolean comparison between two vectors is made on normL2 and direction, while the comparison between scalar number
!> (real or integer of any kind as defined in IR_Precision module) is made on only normL2 of the vector.
!> @ingroup Data_Type_VectorInterface
interface operator (<)
  module procedure vec_low_vec
#ifdef r16p
  module procedure R16P_low_vec,vec_low_R16P
#endif
  module procedure R8P_low_vec,vec_low_R8P
  module procedure R4P_low_vec,vec_low_R4P
  module procedure I8P_low_vec,vec_low_I8P
  module procedure I4P_low_vec,vec_low_I4P
  module procedure I2P_low_vec,vec_low_I2P
  module procedure I1P_low_vec,vec_low_I1P
endinterface
!> @brief Lower-equal-than boolean operator (<=) overloading.
!> @note The boolean comparison between two vectors is made on normL2 and direction, while the comparison between scalar number
!> (real or integer of any kind as defined in IR_Precision module) is made on only normL2 of the vector.
!> @ingroup Data_Type_VectorInterface
interface operator (<=)
  module procedure vec_low_eq_vec
#ifdef r16p
  module procedure R16P_low_eq_vec,vec_low_eq_R16P
#endif
  module procedure R8P_low_eq_vec,vec_low_eq_R8P
  module procedure R4P_low_eq_vec,vec_low_eq_R4P
  module procedure I8P_low_eq_vec,vec_low_eq_I8P
  module procedure I4P_low_eq_vec,vec_low_eq_I4P
  module procedure I2P_low_eq_vec,vec_low_eq_I2P
  module procedure I1P_low_eq_vec,vec_low_eq_I1P
endinterface
!> @brief Equal-to boolean operator (==) overloading.
!> @note The boolean comparison between two vectors is made on normL2 and direction, while the comparison between scalar number
!> (real or integer of any kind as defined in IR_Precision module) is made on only normL2 of the vector.
!> @ingroup Data_Type_VectorInterface
interface operator (==)
  module procedure vec_eq_vec
#ifdef r16p
  module procedure R16P_eq_vec,vec_eq_R16P
#endif
  module procedure R8P_eq_vec,vec_eq_R8P
  module procedure R4P_eq_vec,vec_eq_R4P
  module procedure I8P_eq_vec,vec_eq_I8P
  module procedure I4P_eq_vec,vec_eq_I4P
  module procedure I2P_eq_vec,vec_eq_I2P
  module procedure I1P_eq_vec,vec_eq_I1P
endinterface
!> @brief Higher-equal-than boolean operator (>=) overloading.
!> @note The boolean comparison between two vectors is made on normL2 and direction, while the comparison between scalar number
!> (real or integer of any kind as defined in IR_Precision module) is made on only normL2 of the vector.
!> @ingroup Data_Type_VectorInterface
interface operator (>=)
  module procedure vec_great_eq_vec
#ifdef r16p
  module procedure R16P_great_eq_vec,vec_great_eq_R16P
#endif
  module procedure R8P_great_eq_vec,vec_great_eq_R8P
  module procedure R4P_great_eq_vec,vec_great_eq_R4P
  module procedure I8P_great_eq_vec,vec_great_eq_I8P
  module procedure I4P_great_eq_vec,vec_great_eq_I4P
  module procedure I2P_great_eq_vec,vec_great_eq_I2P
  module procedure I1P_great_eq_vec,vec_great_eq_I1P
endinterface
!> @brief Higher-than boolean operator (>) overloading.
!> @note The boolean comparison between two vectors is made on normL2 and direction, while the comparison between scalar number
!> (real or integer of any kind as defined in IR_Precision module) is made on only normL2 of the vector.
!> @ingroup Data_Type_VectorInterface
interface operator (>)
  module procedure vec_great_vec
#ifdef r16p
  module procedure R16P_great_vec,vec_great_R16P
#endif
  module procedure R8P_great_vec,vec_great_R8P
  module procedure R4P_great_vec,vec_great_R4P
  module procedure I8P_great_vec,vec_great_I8P
  module procedure I4P_great_vec,vec_great_I4P
  module procedure I2P_great_vec,vec_great_I2P
  module procedure I1P_great_vec,vec_great_I1P
endinterface
!> @brief Cross product operator (.cross.) definition.
!> @ingroup Data_Type_VectorInterface
interface operator (.cross.)
  module procedure crossproduct
endinterface
!> @brief Dot product operator (.dot.) definition.
!> @ingroup Data_Type_VectorInterface
interface operator (.dot.)
  module procedure dotproduct
endinterface
!> @brief Parallel product operator (.paral.) definition.
!> @note This operator produces a vector with the normL2 of first vector and parallel to the second vector.
!> @ingroup Data_Type_VectorInterface
interface operator (.paral.)
  module procedure parallel
endinterface
!> @brief Orthogonal product operator (.ortho.) definition.
!> @note This operator produces a vector with the normL2 of first vector and orthogonal to the second vector.
!> @ingroup Data_Type_VectorInterface
interface operator (.ortho.)
  module procedure orthogonal
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_VectorPublicProcedure
  !> @{
  !> @brief Function for computing the square of the norm of a vector.
  !> The square norm if defined as \f$ N = x^2  + y^2  + z^2\f$.
  !> @return \b sq square norm
  elemental function sq_norm(vec) result(sq)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec !< Vector.
  real(R8P)::                      sq  !< Square of the Norm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sq = (vec%x*vec%x) + (vec%y*vec%y) + (vec%z*vec%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sq_norm

  !> @brief Function for computing the norm L2 of a vector.
  !> The norm L2 if defined as \f$N = \sqrt {x^2  + y^2  + z^2 }\f$.
  !> @return \b norm norm L2
  elemental function normL2(vec) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec  !< Vector.
  real(R8P)::                      norm !< Norm L2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  norm = sqrt((vec%x*vec%x) + (vec%y*vec%y) + (vec%z*vec%z))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normL2

  !> @brief Function for normalizing a vector.
  !> The normalization is made by means of norm L2. If the norm L2 of the vector is less than the parameter smallR8P the
  !> normalization value is set to normL2(vec)+smallR8P.
  !> @return \b norm normalized vector
  elemental function normalize(vec) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec  !< Vector to be normalized.
  type(Type_Vector)::             norm !< Vector normalized.
  real(R8P)::                     nm   !< Norm L2 of vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  nm = normL2(vec)
  if (nm < smallR8P) then
    nm = nm + smallR8P
  endif
  norm%x = vec%x/nm
  norm%y = vec%y/nm
  norm%z = vec%z/nm
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normalize

  !> @brief Function for calculating the normal of the face defined by 4 points vector pt1, pt2, pt3 and pt4.
  !> The convention for the points numeration is the following:
  !> @code
  !> 1.----------.2
  !>  |          |
  !>  |          |
  !>  |          |
  !>  |          |
  !> 4.----------.3
  !> @endcode
  !> The normal is calculated by the cross product of the diagonal d13 for the diagonal d24: d13 x d24.
  !> The normal is normalized if the variable 'norm' is passed (with any value).
  !> @return \b fnormal face normal
  elemental function face_normal4(norm,pt1,pt2,pt3,pt4) result(fnormal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(1),      intent(IN), optional:: norm    !< If 'norm' is passed as argument the normal is normalized.
  type(Type_Vector), intent(IN)::           pt1     !< First face point.
  type(Type_Vector), intent(IN)::           pt2     !< Second face point.
  type(Type_Vector), intent(IN)::           pt3     !< Third face point.
  type(Type_Vector), intent(IN)::           pt4     !< Fourth face point.
  type(Type_Vector)::                       fnormal !< Face normal.
  type(Type_Vector)::                       d13,d24 !< Face diagonals.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  d13 = pt3 - pt1
  d24 = pt4 - pt2
  if (present(norm)) then
    fnormal = normalize(d13.cross.d24)
  else
    fnormal = 0.5_R8P*(d13.cross.d24)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction face_normal4

  !> @brief Function for calculating the normal of the face defined by the 3 points vector pt1, pt2 and pt3.
  ! The convention for the points numeration is the following:
  !> @code
  !> 1.----.2
  !>   \   |
  !>    \  |
  !>     \ |
  !>      \|
  !>       .3
  !> @endcode
  !> The normal is calculated by the cross product of the side s12 for the side s13: s12 x s13.
  !> The normal is normalized if the variable 'norm' is passed (with any value).
  !> @return \b fnormal face normal
  elemental function face_normal3(norm,pt1,pt2,pt3) result(fnormal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(1),      intent(IN), optional:: norm    !< If 'norm' is passed as argument the normal is normalized.
  type(Type_Vector), intent(IN)::           pt1     !< First face point.
  type(Type_Vector), intent(IN)::           pt2     !< Second face point.
  type(Type_Vector), intent(IN)::           pt3     !< Third face point.
  type(Type_Vector)::                       fnormal !< Face normal.
  type(Type_Vector)::                       s12,s13 !< Face diagonals.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  s12 = pt2 - pt1
  s13 = pt3 - pt1
  if (present(norm)) then
    fnormal = normalize(s12.cross.s13)
  else
    fnormal = 0.5_R8P*(s12.cross.s13)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction face_normal3

  !> @brief Function for writing Type_Vector data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function write_vector(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN), optional:: scalar         !< Scalar vector data.
  type(Type_Vector), intent(IN), optional:: array1D(:)     !< One dimensional array vector data.
  type(Type_Vector), intent(IN), optional:: array2D(:,:)   !< Two dimensional array vector data.
  type(Type_Vector), intent(IN), optional:: array3D(:,:,:) !< Three dimensional array vector data.
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
  endfunction write_vector

  !> @brief Function for reading Type_Vector data.
  !> The vector data could be scalar, one, two and three dimensional array. The format could be ascii or binary.
  !> @return \b err integer(I_P) variable for error trapping.
  function read_vector(scalar,array1D,array2D,array3D,format,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(INOUT), optional:: scalar         !< Scalar vector data.
  type(Type_Vector), intent(INOUT), optional:: array1D(:)     !< One dimensional array vector data.
  type(Type_Vector), intent(INOUT), optional:: array2D(:,:)   !< Two dimensional array vector data.
  type(Type_Vector), intent(INOUT), optional:: array3D(:,:,:) !< Three dimensional array vector data.
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
  endfunction read_vector
  !> @}

  !> @ingroup Data_Type_VectorPrivateProcedure
  !> @{
  !> @brief Subroutine for setting components of Type_Vector variable.
  elemental subroutine set(vec,x,y,z)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(INOUT)::        vec !< Vector.
  real(R8P),          intent(IN), optional:: x   !< Cartesian component in x direction.
  real(R8P),          intent(IN), optional:: y   !< Cartesian component in y direction.
  real(R8P),          intent(IN), optional:: z   !< Cartesian component in z direction.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(x)) vec%x = x
  if (present(y)) vec%y = y
  if (present(z)) vec%z = z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set

  !> @brief Function for printing in a pretty ascii format the components of type Type_Vector.
  !> @return \b err integer(I_P) variable for error trapping.
  function pprint(vec,myrank,unit) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN)::           vec    !< Vector.
  integer(I_P),       intent(IN), optional:: myrank !< Actual rank process.
  integer(I_P),       intent(IN)::           unit   !< Logic unit.
  integer(I_P)::                             err    !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DI_P)::                          rks    !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,0_I_P)) ; if (present(myrank)) rks = 'rank'//trim(str(.true.,myrank))
  write(unit,'(A)',iostat=err)trim(rks)//' Component x '//str(n=vec%x)
  write(unit,'(A)',iostat=err)trim(rks)//' Component y '//str(n=vec%y)
  write(unit,'(A)',iostat=err)trim(rks)//' Component z '//str(n=vec%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction pprint

  !> @brief Subroutine for normalizing a vector.
  !> The normalization is made by means of norm L2. If the norm L2 of the vector is less than the parameter smallR8P the
  !> normalization value is set to normL2(vec)+smallR8P.
  elemental subroutine normalize_self(vec)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(INOUT):: vec !< Vector to be normalized.
  real(R8P)::                         nm  !< Norm L2 of vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  nm = normL2(vec)
  if (nm < smallR8P) then
    nm = nm + smallR8P
  endif
  vec%x = vec%x/nm
  vec%y = vec%y/nm
  vec%z = vec%z/nm
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine normalize_self

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
  vec%x = real(scal,R8P)
  vec%y = real(scal,R8P)
  vec%z = real(scal,R8P)
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
  vec%x = real(scal,R8P)
  vec%y = real(scal,R8P)
  vec%z = real(scal,R8P)
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
  vec%x = real(scal,R8P)
  vec%y = real(scal,R8P)
  vec%z = real(scal,R8P)
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
  vec%x = real(scal,R8P)
  vec%y = real(scal,R8P)
  vec%z = real(scal,R8P)
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
  vec%x = real(scal,R8P)
  vec%y = real(scal,R8P)
  vec%z = real(scal,R8P)
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
  vec%x = real(scal,R8P)
  vec%y = real(scal,R8P)
  vec%z = real(scal,R8P)
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
  vec%x = real(scal,R8P)
  vec%y = real(scal,R8P)
  vec%z = real(scal,R8P)
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  mul%x = real(scal,R8P) * vec%x
  mul%y = real(scal,R8P) * vec%y
  mul%z = real(scal,R8P) * vec%z
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
  div%x = vec%x / real(scal,R8P)
  div%y = vec%y / real(scal,R8P)
  div%z = vec%z / real(scal,R8P)
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
  div%x = vec%x / real(scal,R8P)
  div%y = vec%y / real(scal,R8P)
  div%z = vec%z / real(scal,R8P)
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
  div%x = vec%x / real(scal,R8P)
  div%y = vec%y / real(scal,R8P)
  div%z = vec%z / real(scal,R8P)
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
  div%x = vec%x / real(scal,R8P)
  div%y = vec%y / real(scal,R8P)
  div%z = vec%z / real(scal,R8P)
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
  div%x = vec%x / real(scal,R8P)
  div%y = vec%y / real(scal,R8P)
  div%z = vec%z / real(scal,R8P)
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
  div%x = vec%x / real(scal,R8P)
  div%y = vec%y / real(scal,R8P)
  div%z = vec%z / real(scal,R8P)
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
  div%x = vec%x / real(scal,R8P)
  div%y = vec%y / real(scal,R8P)
  div%z = vec%z / real(scal,R8P)
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  summ%x = real(scal,R8P) + vec%x
  summ%y = real(scal,R8P) + vec%y
  summ%z = real(scal,R8P) + vec%z
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
  sub%x = real(scal,R8P) - vec%x
  sub%y = real(scal,R8P) - vec%y
  sub%z = real(scal,R8P) - vec%z
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
  sub%x = vec%x - real(scal,R8P)
  sub%y = vec%y - real(scal,R8P)
  sub%z = vec%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - vec%x
  sub%y = real(scal,R8P) - vec%y
  sub%z = real(scal,R8P) - vec%z
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
  sub%x = vec%x - real(scal,R8P)
  sub%y = vec%y - real(scal,R8P)
  sub%z = vec%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - vec%x
  sub%y = real(scal,R8P) - vec%y
  sub%z = real(scal,R8P) - vec%z
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
  sub%x = vec%x - real(scal,R8P)
  sub%y = vec%y - real(scal,R8P)
  sub%z = vec%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - vec%x
  sub%y = real(scal,R8P) - vec%y
  sub%z = real(scal,R8P) - vec%z
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
  sub%x = vec%x - real(scal,R8P)
  sub%y = vec%y - real(scal,R8P)
  sub%z = vec%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - vec%x
  sub%y = real(scal,R8P) - vec%y
  sub%z = real(scal,R8P) - vec%z
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
  sub%x = vec%x - real(scal,R8P)
  sub%y = vec%y - real(scal,R8P)
  sub%z = vec%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - vec%x
  sub%y = real(scal,R8P) - vec%y
  sub%z = real(scal,R8P) - vec%z
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
  sub%x = vec%x - real(scal,R8P)
  sub%y = vec%y - real(scal,R8P)
  sub%z = vec%z - real(scal,R8P)
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
  sub%x = real(scal,R8P) - vec%x
  sub%y = real(scal,R8P) - vec%y
  sub%z = real(scal,R8P) - vec%z
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
  sub%x = vec%x - real(scal,R8P)
  sub%y = vec%y - real(scal,R8P)
  sub%z = vec%z - real(scal,R8P)
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_not_eq_vec

  elemental function vec_not_eq_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_not_eq_vec

  elemental function vec_not_eq_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_not_eq_vec

  elemental function vec_not_eq_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_not_eq_vec

  elemental function vec_not_eq_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_not_eq_vec

  elemental function vec_not_eq_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_not_eq_vec

  elemental function vec_not_eq_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)/=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_not_eq_vec

  elemental function vec_not_eq_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is /= with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)/=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_low_vec

  elemental function vec_low_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_low_vec

  elemental function vec_low_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_low_vec

  elemental function vec_low_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_low_vec

  elemental function vec_low_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_low_vec

  elemental function vec_low_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_low_vec

  elemental function vec_low_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_low_vec

  elemental function vec_low_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is < with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_low_eq_vec

  elemental function vec_low_eq_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_low_eq_vec

  elemental function vec_low_eq_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_low_eq_vec

  elemental function vec_low_eq_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_low_eq_vec

  elemental function vec_low_eq_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_low_eq_vec

  elemental function vec_low_eq_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_low_eq_vec

  elemental function vec_low_eq_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)<=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_low_eq_vec

  elemental function vec_low_eq_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is <= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)<=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_eq_vec

  elemental function vec_eq_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_eq_vec

  elemental function vec_eq_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_eq_vec

  elemental function vec_eq_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_eq_vec

  elemental function vec_eq_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_eq_vec

  elemental function vec_eq_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_eq_vec

  elemental function vec_eq_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)==normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_eq_vec

  elemental function vec_eq_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec is = with respect the value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)==real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_great_eq_vec

  elemental function vec_great_eq_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_great_eq_vec

  elemental function vec_great_eq_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_great_eq_vec

  elemental function vec_great_eq_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_great_eq_vec

  elemental function vec_great_eq_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_great_eq_vec

  elemental function vec_great_eq_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_great_eq_vec

  elemental function vec_great_eq_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>=normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_great_eq_vec

  elemental function vec_great_eq_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is >= with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>=real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_great_vec

  elemental function vec_great_R16P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R16P),        intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_great_vec

  elemental function vec_great_R8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R8P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_great_vec

  elemental function vec_great_R4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  real(R4P),         intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_great_vec

  elemental function vec_great_I8P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I8P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_great_vec

  elemental function vec_great_I4P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I4P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_great_vec

  elemental function vec_great_I2P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I2P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R8P))
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
  type(Type_Vector), intent(IN):: vec     ! Vector.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (real(scal,R8P)>normL2(vec))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_great_vec

  elemental function vec_great_I1P(vec,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function returns .true. if the normL2 of the vec1 is > with respect the  value of scalar scal, .false. otherwise.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec     ! Vector.
  integer(I1P),      intent(IN):: scal    ! Scalar.
  logical::                       compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(vec)>real(scal,R8P))
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
  real(R8P)::                     dot  ! Dot product.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  dot = (vec1%x*vec2%x) + (vec1%y*vec2%y) + (vec1%z*vec2%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dotproduct

  ! othogonal & parallel
  elemental function parallel(vec1,vec2) result(paral)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function calculates the component of vec1 parallel to vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1  ! First vector.
  type(Type_Vector), intent(IN):: vec2  ! Second vector.
  type(Type_Vector)::             paral ! Component of of vec1 parallel to vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  paral = (vec1.dot.normalize(vec2))*normalize(vec2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction parallel

  elemental function orthogonal(vec1,vec2) result(ortho)
  !---------------------------------------------------------------------------------------------------------------------------------
  !!This function calculates the component of vec1 orthogonal to vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec1  ! First vector.
  type(Type_Vector), intent(IN):: vec2  ! Second vector.
  type(Type_Vector)::             ortho ! Component of of vec1 orthogonal to vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ortho = vec1 - (vec1.paral.vec2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction orthogonal
  !> @}
endmodule Data_Type_Vector
