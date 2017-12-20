!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_VectorDerivedType Data_Type_Vector
!> Module definition of Type_Vector
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_VectorGlobalVarPar Data_Type_Vector
!> Module definition of Type_Vector
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_VectorPublicProcedure Data_Type_Vector
!> Module definition of Type_Vector
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_VectorPrivateProcedure Data_Type_Vector
!> Module definition of Type_Vector
!> @}

!> This module contains the definition of Type_Vector and its procedures.
!> This derived type is useful for manipulating selftors in 3D space. The components of the selftors are reals with
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
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type defining vectors.
!> @ingroup Data_Type_VectorDerivedType
type, public:: Type_Vector
  real(R8P):: x = 0._R8P !< Cartesian component in x direction.
  real(R8P):: y = 0._R8P !< Cartesian component in y direction.
  real(R8P):: z = 0._R8P !< Cartesian component in z direction.
  contains
    procedure:: init            => init_vector_self  ! Procedure for initializing vector components.
    procedure:: set             => set_vector_self   ! Procedure for setting vector components.
    procedure:: iolen           => iolen_vector_self ! Procedure for computing IO length.
    procedure:: load            => load_vector_self  ! Procedure for loading Type_Vector data.
    procedure:: save            => save_vector_self  ! Procedure for saving Type_Vector data.
    procedure:: print           => print_vector_self ! Procedure for printing vector components with a "pretty" format.
    procedure:: sq_norm         => sq_norm_self      ! Procedure for computing the square of the norm of a vector.
    procedure:: normL2          => normL2_self       ! Procedure for computing the norm L2 of a vector.
    procedure:: normalize       => normalize_self    ! Procedure for normalizing a vector.
    procedure:: normalized      => normalized_self   ! Procedure for obtaining a normalized copy of a vector.
    procedure:: face_normal4    => face_normal4_self ! Procedure for calculating the normal of the face defined by 4 points vector.
    procedure:: face_normal3    => face_normal3_self ! Procedure for calculating the normal of the face defined by 3 points vector.
    generic:: operator(.cross.) => crossproduct      ! Procedure for computing the cross product of 2 vectors.
    generic:: operator(.dot.)   => dotproduct        ! Procedure for computing the scalar (dot) product of 2 vectors.
    generic:: operator(.paral.) => parallel          ! Procedure for computing the component of vec1 parallel to vec2.
    generic:: operator(.ortho.) => orthogonal        ! Procedure for computign the component of vec1 orthogonal to vec2.
    ! private procedures
    procedure, pass(vec1), private:: crossproduct
    procedure, pass(vec1), private:: dotproduct
    procedure, pass(vec1), private:: parallel
    procedure, pass(vec1), private:: orthogonal
    ! operators overloading
#include "Data_Type_Bounds_Proc_OpOverloading.inc"
    ! conditional operators overloading
#include "Data_Type_Bounds_Proc_CondOpOverloading.inc"
endtype Type_Vector
!> Pointer of Type_Vector for creating array of pointers of Type_Vector.
!> @ingroup Data_Type_VectorDerivedType
type, public:: Type_Vector_Ptr
  type(Type_Vector), pointer:: p => null()
endtype Type_Vector_Ptr
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
contains
  !> @ingroup Data_Type_VectorPublicProcedure
  !> @{
  !> @brief Function for computing the square of the norm of a vector.
  !> The square norm if defined as \f$ N = x^2  + y^2  + z^2\f$.
  elemental function sq_norm(vec) result(sq)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec !< Vector.
  real(R8P)::                     sq  !< Square of the Norm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sq = (vec%x*vec%x) + (vec%y*vec%y) + (vec%z*vec%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sq_norm

  !> @brief Function for computing the norm L2 of a vector.
  !> The norm L2 if defined as \f$N = \sqrt {x^2  + y^2  + z^2 }\f$.
  elemental function normL2(vec) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Vector), intent(IN):: vec  !< Vector.
  real(R8P)::                     norm !< Norm L2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  norm = sqrt((vec%x*vec%x) + (vec%y*vec%y) + (vec%z*vec%z))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normL2

  !> @brief Function for normalizing a vector.
  !> The normalization is made by means of norm L2. If the norm L2 of the vector is less than the parameter smallR8P the
  !> normalization value is set to normL2(vec)+smallR8P.
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
  !> @}

  !> @ingroup Data_Type_VectorPrivateProcedure
  !> @{
  !> @brief Subroutine for initializing components of Type_Vector variable.
  elemental subroutine init_vector_self(vec)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(INOUT):: vec !< Vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vec%x = 0._R8P
  vec%y = 0._R8P
  vec%z = 0._R8P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_vector_self

  !> @brief Subroutine for setting components of Type_Vector variable.
  elemental subroutine set_vector_self(vec,x,y,z)
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
  endsubroutine set_vector_self

  !> @brief Procedure for computing IO length.
  function iolen_vector_self(vec) result(iolen)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec   !< Vector.
  integer(I4P)::                   iolen !< IO length.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(iolength=iolen) vec%x,vec%y,vec%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction iolen_vector_self

  !> @brief Procedure for loading Type_Vector data.
  subroutine load_vector_self(vec,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector),     intent(INOUT):: vec     !< Vector data.
  integer(I8P), optional, intent(IN)::    pos     !< Position specifier.
  integer(I4P), optional, intent(OUT)::   iostat  !< IO error.
  character(*), optional, intent(OUT)::   iomsg   !< IO error message.
  integer(I4P),           intent(IN)::    unit    !< Logic unit.
  integer(I4P)::                          iostatd !< IO error.
  character(500)::                        iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    read(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)vec%x,vec%y,vec%z
  else
    read(unit=unit,        iostat=iostatd,iomsg=iomsgd)vec%x,vec%y,vec%z
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine load_vector_self

  !> @brief Procedure for saving Type_Vector data.
  subroutine save_vector_self(vec,pos,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector),     intent(IN)::  vec     !< Vector data.
  integer(I8P), optional, intent(IN)::  pos     !< Position specifier.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(pos)) then
    write(unit=unit,pos=pos,iostat=iostatd,iomsg=iomsgd)vec%x,vec%y,vec%z
  else
    write(unit=unit,        iostat=iostatd,iomsg=iomsgd)vec%x,vec%y,vec%z
  endif
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine save_vector_self

  !> @brief Function for printing in a pretty ascii format the components of type Type_Vector.
  subroutine print_vector_self(vec,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector),     intent(IN)::  vec     !< Vector.
  character(*), optional, intent(IN)::  pref    !< Prefixing string for outputs.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::       prefd   !< Prefixing string.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)pref//' Component x '//str(n=vec%x)
  write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)pref//' Component y '//str(n=vec%y)
  write(unit,'(A)',iostat=iostatd,iomsg=iomsgd)pref//' Component z '//str(n=vec%z)
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = trim(adjustl(iomsgd))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_vector_self

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

  !> @brief Subroutine for obtaining a normalized copy of a vector.
  !> The normalization is made by means of norm L2. If the norm L2 of the vector is less than the parameter smallR8P the
  !> normalization value is set to normL2(vec)+smallR8P.
  elemental function normalized_self(vec) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec  !< Vector to be normalized.
  type(Type_Vector)::              norm !< Normalized copy.
  real(R8P)::                      nm   !< Norm L2 of vector.
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
  endfunction normalized_self

  !> @brief Function for computing the square of the norm of a vector.
  !> The square norm if defined as \f$ N = x^2  + y^2  + z^2\f$.
  elemental function sq_norm_self(vec) result(sq)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec !< Vector.
  real(R8P)::                      sq  !< Square of the Norm.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sq = (vec%x*vec%x) + (vec%y*vec%y) + (vec%z*vec%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction sq_norm_self

  !> @brief Function for computing the norm L2 of a vector.
  !> The norm L2 if defined as \f$N = \sqrt {x^2  + y^2  + z^2 }\f$.
  elemental function normL2_self(vec) result(norm)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec  !< Vector.
  real(R8P)::                      norm !< Norm L2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  norm = sqrt((vec%x*vec%x) + (vec%y*vec%y) + (vec%z*vec%z))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction normL2_self

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
  elemental subroutine face_normal4_self(fnormal,norm,pt1,pt2,pt3,pt4)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector),     intent(INOUT):: fnormal !< Face normal.
  character(1), optional, intent(IN)::    norm    !< If 'norm' is passed as argument the normal is normalized.
  type(Type_Vector),      intent(IN)::    pt1     !< First face point.
  type(Type_Vector),      intent(IN)::    pt2     !< Second face point.
  type(Type_Vector),      intent(IN)::    pt3     !< Third face point.
  type(Type_Vector),      intent(IN)::    pt4     !< Fourth face point.
  type(Type_Vector)::                     d13,d24 !< Face diagonals.
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
  endsubroutine face_normal4_self

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
  elemental subroutine face_normal3_self(fnormal,norm,pt1,pt2,pt3)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector),     intent(INOUT):: fnormal !< Face normal.
  character(1), optional, intent(IN)::    norm    !< If 'norm' is passed as argument the normal is normalized.
  type(Type_Vector),      intent(IN)::    pt1     !< First face point.
  type(Type_Vector),      intent(IN)::    pt2     !< Second face point.
  type(Type_Vector),      intent(IN)::    pt3     !< Third face point.
  type(Type_Vector)::                     s12,s13 !< Face diagonals.
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
  endsubroutine face_normal3_self

  !> @brief Procedure for computing the cross product of 2 vectors.
  !> \f$\vec V=\left({y_1 z_2 - z_1 y_2}\right)\vec i +
  !>           \left({z_1 x_2 - x_1 z_2}\right)\vec j +
  !>           \left({x_1 y_2 - y_1 x_2}\right)\vec k\f$
  !> where \f$x_i\f$, \f$y_i\f$ and \f$z_i\f$ \f$i=1,2\f$ are the components of the vectors.
  elemental function crossproduct(vec1,vec2) result(cross)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec1  ! First vector.
  type(Type_Vector),  intent(IN):: vec2  ! Second vector.
  type(Type_Vector)::              cross ! Cross product vector.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cross%x = (vec1%y*vec2%z) - (vec1%z*vec2%y)
  cross%y = (vec1%z*vec2%x) - (vec1%x*vec2%z)
  cross%z = (vec1%x*vec2%y) - (vec1%y*vec2%x)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction crossproduct

  !> @brief Procedure for computing the scalar (dot) product of 2 vectors.
  !> \f${\rm D}= x_1 \cdot x_2 + y_1 \cdot y_2 + z_1 \cdot z_2\f$
  !> where \f$x_i\f$, \f$y_i\f$ and \f$z_i\f$ \f$i=1,2\f$ are the components of the vectors.
  elemental function dotproduct(vec1,vec2) result(dot)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec1 ! First vector.
  type(Type_Vector),  intent(IN):: vec2 ! Second vector.
  real(R8P)::                      dot  ! Dot product.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  dot = (vec1%x*vec2%x) + (vec1%y*vec2%y) + (vec1%z*vec2%z)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction dotproduct

  !> @brief Procedure for computing the component of vec1 parallel to vec2.
  elemental function parallel(vec1,vec2) result(paral)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec1  ! First vector.
  type(Type_Vector),  intent(IN):: vec2  ! Second vector.
  type(Type_Vector)::              paral ! Component of of vec1 parallel to vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  paral = (vec1.dot.normalize(vec2))*normalize(vec2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction parallel

  !> @brief Procedure for computign the component of vec1 orthogonal to vec2.
  elemental function orthogonal(vec1,vec2) result(ortho)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: vec1  ! First vector.
  type(Type_Vector),  intent(IN):: vec2  ! Second vector.
  type(Type_Vector)::              ortho ! Component of of vec1 orthogonal to vec2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ortho = vec1 - (vec1.paral.vec2)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction orthogonal

  ! Operators overloading.
  ! Operator (=)
#define self_type_ class(Type_Vector), intent(INOUT):: self
#define ass_scal_ self%x = real(scal,R8P) ; self%y = real(scal,R8P) ; self%z = real(scal,R8P)
  !> @brief Procedure for assignment between two selfs.
  pure subroutine assign_self(self1,self2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(INOUT):: self1
  class(Type_Vector), intent(IN)::    self2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self1%x = self2%x
  self1%y = self2%y
  self1%z = self2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_self
#include "Data_Type_Bounds_Proc_AssDefinitions.inc"
#undef self_type_
#undef ass_scal_

  ! Operator (*)
#define self_type_ class(Type_Vector), intent(IN):: self
#define mul_type_ type(Type_Vector):: mul
#define mul_scal_ mul%x = real(scal,R8P)*self%x ; mul%y = real(scal,R8P)*self%y ; mul%z = real(scal,R8P)*self%z
  !> @brief Procedure for multiply (by components) two selfs.
  elemental function self_mul_self(self1,self2) result(mul)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1
  type(Type_Vector),  intent(IN):: self2
  type(Type_Vector)::              mul
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  mul%x = self1%x * self2%x
  mul%y = self1%y * self2%y
  mul%z = self1%z * self2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_mul_self
#include "Data_Type_Bounds_Proc_MulDefinitions.inc"
#undef self_type_
#undef mul_type_
#undef mul_scal_

  ! Operator (/)
#define self_type_ class(Type_Vector), intent(IN):: self
#define div_type_ type(Type_Vector):: div
#define div_scal_ div%x = self%x/real(scal,R8P) ; div%y = self%y/real(scal,R8P) ; div%z = self%z/real(scal,R8P)
  !> @brief Procedure for divide self for self.
  elemental function self_div_self(self1,self2) result(div)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1
  type(Type_Vector),  intent(IN):: self2
  type(Type_Vector)::              div
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  div%x = self1%x / self2%x
  div%y = self1%y / self2%y
  div%z = self1%z / self2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_div_self
#include "Data_Type_Bounds_Proc_DivDefinitions.inc"
#undef self_type_
#undef div_type_
#undef div_scal_

  ! Operator (+)
#define self_type_ class(Type_Vector), intent(IN):: self
#define summ_type_ type(Type_Vector):: summ
#define sum_scal_ summ%x = real(scal,R8P)+self%x ; summ%y = real(scal,R8P)+self%y ; summ%z = real(scal,R8P)+self%z
  !> @brief Procedure for applay unary + to a self.
  elemental function positive_self(self) result(pos)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  type(Type_Vector)::              pos
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  pos%x = + self%x
  pos%y = + self%y
  pos%z = + self%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction positive_self

  !> @brief Procedure for sum self and self.
  elemental function self_sum_self(self1,self2) result(summ)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1
  type(Type_Vector),  intent(IN):: self2
  type(Type_Vector)::              summ
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  summ%x = self1%x + self2%x
  summ%y = self1%y + self2%y
  summ%z = self1%z + self2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sum_self
#include "Data_Type_Bounds_Proc_SumDefinitions.inc"
#undef self_type_
#undef summ_type_
#undef sum_scal_

  ! Operator (-)
#define self_type_ class(Type_Vector), intent(IN):: self
#define sub_type_ type(Type_Vector):: sub
#define self_sub_scal_ sub%x = self%x-real(scal,R8P) ; sub%y = self%y-real(scal,R8P) ; sub%z = self%z-real(scal,R8P)
#define scal_sub_self_ sub%x = real(scal,R8P)-self%x ; sub%y = real(scal,R8P)-self%y ; sub%z = real(scal,R8P)-self%z
  !> @brief Procedure for applay unary - to a self.
  elemental function negative_self(self) result(neg)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  type(Type_Vector)::              neg
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neg%x = - self%x
  neg%y = - self%y
  neg%z = - self%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction negative_self

  !> @brief Procedure for subtract self and self.
  elemental function self_sub_self(self1,self2) result(sub)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1
  type(Type_Vector),  intent(IN):: self2
  type(Type_Vector)::              sub
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  sub%x = self1%x - self2%x
  sub%y = self1%y - self2%y
  sub%z = self1%z - self2%z
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_sub_self
#include "Data_Type_Bounds_Proc_SubDefinitions.inc"
#undef self_type_
#undef sub_type_
#undef self_sub_scal_
#undef scal_sub_self_

  ! Conditional operators
  ! Oprator /=
#define neq_scal_ compare = (real(scal,R8P)/=normL2(self))
  !> @brief Procedure returns .true. if the normL2 of the self1 is /= with respect the normL2 of self2 or if the directions of self1
  !> and self2 are different, .false. otherwise.
  elemental function self_not_eq_self(self1,self2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1   !< First selftor.
  type(Type_Vector),  intent(IN):: self2   !< Second selftor.
  logical::                        compare !< The result of the comparison.
  type(Type_Vector)::              n1,n2   !< Normalizations of self1 and self2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(self1)/=normL2(self2))
  if (.not.compare) then ! the normL2 are the same, checking the directions
    n1 = normalize(self1)
    n2 = normalize(self2)
    compare = ((n1%x/=n2%x).OR.(n1%y/=n2%y).OR.(n1%z/=n2%z))
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_not_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function R16P_not_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),         intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_not_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function self_not_eq_R16P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  real(R16P),         intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_not_eq_R16P

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function R8P_not_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),          intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_not_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function self_not_eq_R8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  real(R8P),          intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_not_eq_R8P

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function R4P_not_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),          intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_not_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function self_not_eq_R4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  real(R4P),          intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_not_eq_R4P

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function I8P_not_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),       intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_not_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function self_not_eq_I8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  integer(I8P),       intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_not_eq_I8P

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function I4P_not_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),       intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_not_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function self_not_eq_I4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  integer(I4P),       intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_not_eq_I4P

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function I2P_not_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),       intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_not_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function self_not_eq_I2P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  integer(I2P),       intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_not_eq_I2P

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function I1P_not_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),       intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_not_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is /= with respect the value of scalar scal, .false. otherwise.
  elemental function self_not_eq_I1P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  integer(I1P),       intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  neq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_not_eq_I1P
#undef neq_scal_

  ! Oprator <
#define self_lt_scal_ compare = (normL2(self)<real(scal,R8P))
#define scal_lt_self_ compare = (real(scal,R8P)<normL2(self))
  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the normL2 of self2, .false. otherwise.
  elemental function self_low_self(self1,self2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1
  type(Type_Vector),  intent(IN):: self2
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(self1)<normL2(self2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function R16P_low_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),         intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_lt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_low_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_R16P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  real(R16P),         intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_lt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_R16P

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function R8P_low_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),          intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_lt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_low_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_R8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  real(R8P),          intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_lt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_R8P

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function R4P_low_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),          intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_lt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_low_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_R4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  real(R4P),          intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_lt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_R4P

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function I8P_low_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),       intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_lt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_low_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_I8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  integer(I8P),       intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_lt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_I8P

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function I4P_low_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),       intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_lt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_low_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_I4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  integer(I4P),       intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_lt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_I4P

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function I2P_low_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),       intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_lt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_low_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_I2P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  integer(I2P),       intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_lt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_I2P

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function I1P_low_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),       intent(IN):: scal
  class(Type_Vector), intent(IN):: self
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_lt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_low_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is < with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_I1P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self
  integer(I1P),       intent(IN):: scal
  logical::                        compare
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_lt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_I1P
#undef self_lt_scal_
#undef scal_lt_self_

  ! Oprator <=
#define self_le_scal_ compare = (normL2(self)<=real(scal,R8P))
#define scal_le_self_ compare = (real(scal,R8P)<=normL2(self))
  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the normL2 of self2, .false. otherwise.
  elemental function self_low_eq_self(self1,self2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1    ! First selftor.
  type(Type_Vector),  intent(IN):: self2    ! Second selftor.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(self1)<=normL2(self2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function R16P_low_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),         intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_le_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_low_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_eq_R16P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R16P),         intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_le_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_eq_R16P

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function R8P_low_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),          intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_le_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_low_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_eq_R8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R8P),          intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_le_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_eq_R8P

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function R4P_low_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),          intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_le_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_low_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_eq_R4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R4P),          intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_le_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_eq_R4P

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function I8P_low_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_le_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_low_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_eq_I8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I8P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_le_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_eq_I8P

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function I4P_low_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_le_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_low_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_eq_I4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I4P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_le_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_eq_I4P

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function I2P_low_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_le_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_low_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_eq_I2P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I2P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_le_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_eq_I2P

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function I1P_low_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_le_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_low_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is <= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_low_eq_I1P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I1P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_le_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_low_eq_I1P
#undef self_le_scal_
#undef scal_le_self_

  ! Oprator ==
#define eq_scal_ compare = (normL2(self)==real(scal,R8P))
  !> @brief Procedure returns .true. if the normL2 of the self1 is = with respect the normL2 of self2 and the directions of
  !> self1 and self2 are the same, .false. otherwise.
  elemental function self_eq_self(self1,self2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1    ! First selftor.
  type(Type_Vector),  intent(IN):: self2    ! Second selftor.
  logical::                        compare ! The result of the comparison.
  type(Type_Vector)::              n1,n2   ! Normalizations of self1 and self2.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(self1)==normL2(self2))
  if (compare) then ! the normL2 are the same, checking the directions
    n1 = normalize(self1)
    n2 = normalize(self2)
    compare = ((n1%x==n2%x).AND.(n1%y==n2%y).AND.(n1%z==n2%z))
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function R16P_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),         intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function self_eq_R16P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R16P),         intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_eq_R16P

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function R8P_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),          intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function self_eq_R8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R8P),          intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_eq_R8P

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function R4P_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),          intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function self_eq_R4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R4P),          intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_eq_R4P

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function I8P_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function self_eq_I8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I8P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_eq_I8P

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function I4P_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function self_eq_I4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I4P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_eq_I4P

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function I2P_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function self_eq_I2P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I2P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_eq_I2P

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function I1P_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self is = with respect the value of scalar scal, .false. otherwise.
  elemental function self_eq_I1P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I1P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  eq_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_eq_I1P
#undef eq_scal_

  ! Oprator >=
#define self_ge_scal_ compare = (normL2(self)>=real(scal,R8P))
#define scal_ge_self_ compare = (real(scal,R8P)>=normL2(self))
  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the normL2 of self2, .false. otherwise.
  elemental function self_great_eq_self(self1,self2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1    ! First selftor.
  type(Type_Vector),  intent(IN):: self2    ! Second selftor.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(self1)>=normL2(self2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function R16P_great_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),         intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_ge_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_great_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_eq_R16P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R16P),         intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_ge_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_eq_R16P

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function R8P_great_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),          intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_ge_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_great_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_eq_R8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R8P),          intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_ge_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_eq_R8P

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function R4P_great_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),          intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_ge_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_great_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_eq_R4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R4P),          intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_ge_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_eq_R4P

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function I8P_great_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_ge_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_great_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_eq_I8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I8P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_ge_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_eq_I8P

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function I4P_great_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_ge_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_great_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_eq_I4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I4P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_ge_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_eq_I4P

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function I2P_great_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_ge_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_great_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_eq_I2P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I2P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_ge_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_eq_I2P

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function I1P_great_eq_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_ge_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_great_eq_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is >= with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_eq_I1P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I1P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_ge_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_eq_I1P
#undef self_ge_scal_
#undef scal_ge_self_

  ! Oprator >
#define self_gt_scal_ compare = (normL2(self)>real(scal,R8P))
#define scal_gt_self_ compare = (real(scal,R8P)>normL2(self))
  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the normL2 of self2, .false. otherwise.
  elemental function self_great_self(self1,self2) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self1    ! First selftor.
  type(Type_Vector),  intent(IN):: self2    ! Second selftor.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  compare = (normL2(self1)>normL2(self2))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function R16P_great_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R16P),         intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_gt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R16P_great_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_R16P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R16P),         intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_gt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_R16P

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function R8P_great_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),          intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_gt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R8P_great_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_R8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R8P),          intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_gt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_R8P

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function R4P_great_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R4P),          intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_gt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction R4P_great_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_R4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  real(R4P),          intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_gt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_R4P

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function I8P_great_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I8P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_gt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I8P_great_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_I8P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I8P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_gt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_I8P

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function I4P_great_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_gt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I4P_great_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_I4P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I4P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_gt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_I4P

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function I2P_great_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I2P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_gt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I2P_great_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_I2P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I2P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_gt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_I2P

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function I1P_great_self(scal,self) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P),       intent(IN):: scal    ! Scalar.
  class(Type_Vector), intent(IN):: self     ! Vector.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  scal_gt_self_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction I1P_great_self

  !> @brief Procedure returns .true. if the normL2 of the self1 is > with respect the  value of scalar scal, .false. otherwise.
  elemental function self_great_I1P(self,scal) result(compare)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Vector), intent(IN):: self     ! Vector.
  integer(I1P),       intent(IN):: scal    ! Scalar.
  logical::                        compare ! The result of the comparison.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self_gt_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction self_great_I1P
#undef self_gt_scal_
#undef scal_gt_self_
  !> @}
endmodule Data_Type_Vector
