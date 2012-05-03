!> This module contains mathematical procedures.
!> @todo \b DocComplete: Complete the documentation of internal procedures
module Lib_Math
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                        ! Integers and reals precision definition.
USE Data_Type_Vector, only: Type_Vector ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: stretchinglside
public:: stretchingrside
public:: stretching2side
public:: average
public:: digit
public:: interpolate1
public:: interpolate2
public:: interpolate3
public:: pi
#ifdef r16p
public:: pi_R16
#endif
public:: pi_R8
public:: pi_R4
public:: degree
public:: radiant
public:: delta1_2o
public:: delta2_2o
public:: abs_grad
public:: laplacian
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! Pi greek definitions with parametric kind precision.
#ifdef r16p
real(R16P), parameter:: pi_R16 = 2._R16P*asin(1._R16P) !< Pi greek with R16P precision.
#endif
real(R8P),  parameter:: pi_R8  = 2._R8P* asin(1._R8P)  !< Pi greek with R8P  precision.
real(R4P),  parameter:: pi_R4  = 2._R4P* asin(1._R4P)  !< Pi greek with R4P  precision.
real(R_P),  parameter:: pi     = 2._R_P* asin(1._R_P)  !< Pi greek with R_P  precision.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!>Average overloading
interface average
#ifdef r16p
  module procedure average_Vectorial1D_R16,average_Vectorial2D_R16,average_Vectorial3D_R16
#endif
  module procedure average_Vectorial1D_R8,average_Vectorial2D_R8,average_Vectorial3D_R8
  module procedure average_Vectorial1D_R4,average_Vectorial2D_R4,average_Vectorial3D_R4
endinterface
!>Digit overloading
interface digit
  module procedure digit_I8,digit_I4,digit_I2,digit_I1
endinterface
!>Linear interpolate overloading
interface interpolate1
#ifdef r16p
  module procedure interpolate1_R16
#endif
  module procedure interpolate1_R8,interpolate1_R4
endinterface
!>Bi-linear interpolate overloading
interface interpolate2
#ifdef r16p
  module procedure interpolate2_R16
#endif
  module procedure interpolate2_R8,interpolate2_R4
endinterface
!>Tri-linear interpolate overloading
interface interpolate3
#ifdef r16p
  module procedure interpolate3_R16
#endif
  module procedure interpolate3_R8,interpolate3_R4
endinterface
!>Degree overloading
interface degree
#ifdef r16p
  module procedure degree_Scalar_R16,degree_Vectorial1D_R16,degree_Vectorial2D_R16,degree_Vectorial3D_R16
#endif
  module procedure degree_Scalar_R8,degree_Vectorial1D_R8,degree_Vectorial2D_R8,degree_Vectorial3D_R8
  module procedure degree_Scalar_R4,degree_Vectorial1D_R4,degree_Vectorial2D_R4,degree_Vectorial3D_R4
endinterface
!>Radiant overloading
interface radiant
#ifdef r16p
  module procedure radiant_Scalar_R16,radiant_Vectorial1D_R16,radiant_Vectorial2D_R16,radiant_Vectorial3D_R16
#endif
  module procedure radiant_Scalar_R8,radiant_Vectorial1D_R8,radiant_Vectorial2D_R8,radiant_Vectorial3D_R8
  module procedure radiant_Scalar_R4,radiant_Vectorial1D_R4,radiant_Vectorial2D_R4,radiant_Vectorial3D_R4
endinterface
!>delta1_2o overloading
interface delta1_2o
#ifdef r16p
  module procedure delta1_2o_R16
#endif
  module procedure delta1_2o_R8,delta1_2o_R4
endinterface
!>delta2_2o overloading
interface delta2_2o
#ifdef r16p
  module procedure delta2_2o_R16
#endif
  module procedure delta2_2o_R8,delta2_2o_R4
endinterface
!>abs_grad overloading
interface abs_grad
#ifdef r16p
  module procedure abs_grad_R16
#endif
  module procedure abs_grad_R8,abs_grad_R4
endinterface
!>laplacian overloading
interface laplacian
#ifdef r16p
  module procedure laplacian_R16
#endif
  module procedure laplacian_R8,laplacian_R4
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
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
#ifdef r16p
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
  mean=sum(x)/size(x)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial1D_R16
#endif

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
  mean=sum(x)/size(x)
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
  mean=sum(x)/size(x)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial1D_R4

#ifdef r16p
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
  mean=sum(x)/size(x)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial2D_R16
#endif

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
  mean=sum(x)/size(x)
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
  mean=sum(x)/size(x)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial2D_R4

#ifdef r16p
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
  mean=sum(x)/size(x)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction average_Vectorial3D_R16
#endif

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
  mean=sum(x)/size(x)
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
  mean=sum(x)/size(x)
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
  integer::                  digit ! Number of digits.
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
  integer::                  digit ! Number of digits.
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
  integer::                  digit ! Number of digits.
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
  integer::                  digit ! Number of digits.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(str,FI1P) abs(n)         ! Casting of n to string.
  digit = len_trim(adjustl(str)) ! Calculating the digits number of n.
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction digit_I1

  ! interpolate
#ifdef r16p
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
#endif

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

#ifdef r16p
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
#endif

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

#ifdef r16p
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
#endif

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
#ifdef r16p
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
#endif

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

#ifdef r16p
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
#endif

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

#ifdef r16p
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
#endif

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

#ifdef r16p
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
#endif

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
#ifdef r16p
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
#endif

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

#ifdef r16p
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
#endif

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

#ifdef r16p
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
#endif

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

#ifdef r16p
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
#endif

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
#ifdef r16p
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
#endif

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
#ifdef r16p
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
#endif

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
#ifdef r16p
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
#endif

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
#ifdef r16p
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
#endif

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
