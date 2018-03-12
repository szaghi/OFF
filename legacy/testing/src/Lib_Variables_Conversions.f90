!> @ingroup Library
!> @{
!> @defgroup Lib_Variables_ConversionsLibrary Lib_Variables_Conversions
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Variables_ConversionsPublicProcedure Lib_Variables_Conversions
!> @}

!> This module contains the definition of variables set (conservative, primitive, ecc...) conversions.
!> This is a library module.
!> @todo \b DocComplete: Complete the documentation of internal procedures
!> @ingroup Lib_Variables_ConversionsLibrary
module Lib_Variables_Conversions
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                  ! Integers and reals precision definition.
USE Data_Type_Conservative,        only: Type_Conservative        ! Definition of Type_Conservative.
USE Data_Type_Primitive,           only: Type_Primitive           ! Definition of Type_Primitive.
USE Data_Type_Primitive1D,         only: Type_Primitive1D         ! Definition of Type_Primitive1D.
USE Data_Type_Riemann_Primitive1D, only: Type_Riemann_Primitive1D ! Definition of Type_Riemann_Primitive1D.
USE Data_Type_Species,             only: Type_Species             ! Definition of Type_Species.
USE Data_Type_Vector,              only: Type_Vector              ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: prim2cons,cons2prim
public:: prim3D2prim1D,prim1D2prim3D
public:: charac2prim1D,prim1D2charac
public:: prim1D2riemprim1D
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Variables_ConversionsPublicProcedure
  !> @{
  !> @brief Procedure for converting primitive variables to conservative variables.
  !> The conversion laws are: \n
  !> - Partial densities\f$|_{conservative} =\left| \rho_s\right|_{primitive}\f$
  !> - Momentum\f$|_{conservative} = \left| \rho \vec V \right|_{primitive}\f$
  !> - Energy\f$|_{conservative}=\left|\frac{p}{\gamma-1}+\frac{1}{2}\rho
  !>   \left(V_x^2+V_y^2+V_z^2\right)\right|_{primitive}\f$
  elemental subroutine prim2cons(prim,cons)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive),    intent(IN)::    prim !< Primitive variables.
  type(Type_Conservative), intent(INOUT):: cons !< Conservative variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  cons%rs = prim%r
  cons%rv = prim%d*prim%v
  cons%re = prim%p/(prim%g-1._R8P) + 0.5_R8P*prim%d*prim%v%sq_norm()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim2cons

  !> @brief Procedure for converting conservative variables to primitive variables.
  !> The conversion laws are: \n
  !> - Partial densities\f$|_{primitive} =\left| \rho_s\right|_{conservative}\f$
  !> - Density\f$|_{primitive} =\left|\sum_{s=1}^{Ns}\rho_s\right|_{conservative}\f$
  !> - Specific heats ratio\f$|_{primitive}=\left|\frac{\sum_{s=1}^{Ns}{\frac{\rho_s}{\rho}cp_s^0}}
  !>   {\sum_{s=1}^{Ns}{\frac{\rho_s}{\rho}cv_s^0}}\right|_{conservative}\f$
  !> - Velocity\f$|_{primitive}=\left|\frac{\vec momentum}{\rho}\right|_{conservative}\f$
  !> - Pressure\f$|_{primitive}=\left|(\gamma-1)\left[energy-\frac{1}{2}\rho
  !>   \left(V_x^2+V_y^2+V_z^2\right)\right]\right|_{conservative}\f$ \n
  !> where Ns is the number of initial species and \f$cp_s^0,cv_s^0\f$ are the initial specific heats.
  elemental subroutine cons2prim(cons,prim,species0)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Conservative), intent(IN)::    cons     !< Conservative variables.
  type(Type_Primitive),    intent(INOUT):: prim     !< Primitive variables.
  type(Type_Species),      intent(IN)::    species0 !< Initial species.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim%r = cons%rs
  call prim%compute_d
  call prim%compute_g(species0=species0)
  prim%v = cons%rv/prim%d
  prim%p = (prim%g - 1._R8P)*(cons%re - 0.5_R8P*prim%d*prim%v%sq_norm())
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine cons2prim

  !> @brief Subroutine for converting (3D) primitive variables to 1D primitive ones.
  elemental subroutine prim3D2prim1D(prim3D,N,prim1D)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive),   intent(IN)::    prim3D !< 3D primitive data.
  type(Type_Vector),      intent(IN)::    N      !< Normal vector for 1D projection.
  type(Type_Primitive1D), intent(INOUT):: prim1D !< 1D primitive data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim1D%r =  prim3D%r
  prim1D%u = (prim3D%v.dot.N)
  prim1D%p =  prim3D%p
  prim1D%d =  prim3D%d
  prim1D%g =  prim3D%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim3D2prim1D

  !> @brief Subroutine for converting 1D primitive variables to (3D) primitive ones.
  elemental subroutine prim1D2prim3D(prim1D,N,prim3D)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive1D), intent(IN)::    prim1D !< 1D primitive data.
  type(Type_Vector),      intent(IN)::    N      !< Normal vector for 1D projection.
  type(Type_Primitive),   intent(INOUT):: prim3D !< 3D primitive data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prim3D%r =  prim1D%r
  prim3D%v = (prim1D%u*N)
  prim3D%p =  prim1D%p
  prim3D%d =  prim1D%d
  prim3D%g =  prim1D%g
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim1D2prim3D

  !> @brief Subroutine for converting (local) characteristic variables to 1D primitive ones.
  pure subroutine charac2prim1D(prim,eigenvectR,charac)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R8P),              intent(IN)::    charac(1:)                                !< Charactetistic data.
  real(R8P),              intent(IN)::    eigenvectR(1:size(charac),1:size(charac)) !< Eigenvectors matrix.
  type(Type_Primitive1D), intent(INOUT):: prim                                      !< Primitive data.
  real(R8P)::                             prim_array(1:size(charac))                !< Temporary array for the primitive data.
  integer(I4P)::                          s,Ns                                      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(charac,dim=1)-4
  do s=1,Ns+4
    prim_array(s) = dot_product(eigenvectR(s,1:Ns+4),charac(1:Ns+4))
  enddo
  call prim%array2prim(prim_array)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine charac2prim1D

  !> @brief Subroutine for converting 1D primitive variables to (local) characteristic ones.
  pure subroutine prim1D2charac(prim,eigenvectL,charac)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive1D), intent(IN)::    prim                                          !< Primitive data.
  real(R8P),              intent(IN)::    eigenvectL(1:size(prim%r)+4,1:size(prim%r)+4) !< Eigenvectors matrix.
  real(R8P),              intent(INOUT):: charac(1:size(prim%r)+4)                      !< Charactetistic data.
  real(R8P)::                             prim_array(1:size(prim%r)+4)                  !< Temporary array for the primitive data.
  integer(I4P)::                          s,Ns                                          !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Ns = size(prim%r,dim=1)
  prim_array = prim%prim2array()
  do s=1,Ns+4
    charac(s) = dot_product(eigenvectL(s,1:Ns+4),prim_array(1:Ns+4))
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim1D2charac

  !> @brief Subroutine for converting 1D primitive variables to Riemann primitive 1D ones.
  elemental subroutine prim1D2riemprim1D(prim1D,riemprim1D)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Primitive1D),         intent(IN)::    prim1D     !< 1D primitive data.
  type(Type_Riemann_Primitive1D), intent(INOUT):: riemprim1D !< Riemann 1D primitive data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call riemprim1D%set(r=prim1D%d,u=prim1D%u,p=prim1D%p,g=prim1D%g) ; call riemprim1D%compute_dependent
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine prim1D2riemprim1D
  !> @}
endmodule Lib_Variables_Conversions
