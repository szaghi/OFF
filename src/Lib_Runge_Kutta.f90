!> @ingroup Library
!> @{
!> @defgroup Lib_Runge_KuttaLibrary Lib_Runge_Kutta
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup GvarLib_Runge_Kutta Lib_Runge_Kutta
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_Runge_KuttaPublicProcedure Lib_Runge_Kutta
!> @}

!> This module contains the definition of Runge-Kutta procedures for time integration.
!> @note
!> This is a library module. This module provides procedure for performing ODE time integration by means numerical scheme belongs to
!> the explicit, Strong Stability Preserving Runge-Kutta family. The implementation is explicit without low-storage form.
!> The scheme is written using the Butcher table (see Butcher, J.C., "Coefficients for the study of Runge-Kutta integration
!> processes", J. Austral. Math. Soc., Vol. 3, pages: 185--201, 1963). Considering the following ODE (Ordinary Differential
!> Equation): \n
!> \f$ u_t = R(t,u) \f$ \n
!> where R is the residual function. The scheme adopted can written in the form: \n
!> \f$ u^{n+1} = u^n +\Delta t \sum_{s=1}^{Ns}c_1^s K^s \f$ \n
!> where Ns is the number of stages used and \f$K^s\f$ is the \f$s^{th}\f$ stage computed as: \n
!> \f$ K^s = R\left( t^n+c_3^s \Delta t, u^n +\Delta t \sum_{i=1}^{s-1}c_2^{s,i} K^i \right) \f$ \n
!> The scheme is explicit thus the above summation is up to \f$s-1\f$. The coefficients \f$c_1\f$, \f$c_2\f$ and \f$c_3\f$ are given
!> in the Butcher table form: \n
!> \f$\begin{array}{l} \begin{array}{*{20}{c}} {c_3^1}\\ {c_3^2}\\ \vdots \\ {c_3^{Ns}}
!> \end{array}\left| \!{\underline {\, {\begin{array}{*{20}{c}} {c_2^{1,1}}&{c_2^{1,2}}&
!> \cdots &{c_2^{1,Ns}}\\ {c_2^{2,1}}&{c_2^{2,2}}& \cdots &{c_2^{1,Ns}}\\ \vdots & \vdots & \ddots &
!> \vdots \\ {c_2^{Ns,1}}&{c_2^{Ns,2}}& \ldots &{c_2^{Ns,Ns}} \end{array}} \,}} \right. \\ \begin{array}{*{20}{c}}
!> {\;\;\;\;\;\;\;}&{c_1^1}&{\;\;\;\;\;c_1^2}&{\;\;\;\;\; \cdots }&{c_1^{Ns}} \end{array} \end{array}\f$ \n
!> Because only explicit scheme is considered the Butcher table reduces to diagonal matrix: \n
!> \f$\begin{array}{l} \begin{array}{*{20}{c}} {c_3^1}\\ {c_3^1}\\ \vdots \\ {c_3^{Ns}}
!> \end{array}\left| \!{\underline {\, {\begin{array}{*{20}{c}} 0&0& \cdots &0\\ {c_2^{2,1}}&0& \cdots &0\\ \vdots & \vdots &
!> \ddots & \vdots \\ {c_2^{Ns,1}}&{c_2^{Ns,2}}& \ldots &0 \end{array}} \,}} \right. \\ \begin{array}{*{20}{c}}
!> {\;\;\;\;\;\;\;}&{c_1^1}&{\;\;\;\;\;c_1^2}&{\;\;\;\;\; \cdots }&{c_1^{Ns}} \end{array} \end{array}\f$ \n
!> Moreover the following relation always holds:
!> \f$ c_3^s = \sum_{i=1}^{Ns}c_2^{s,i} \f$ \n
!> According to the number of stages used different schemes are available. See the \ref rk_init "procedure rk_init" for more
!> details.
!>
!> @todo \b RK-implicit: Implement (semi)-implicit schemes.
!> @ingroup Lib_Runge_KuttaLibrary
module Lib_Runge_Kutta
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision           ! Integers and reals precision definition.
USE Data_Type_Conservative ! Definition of Type_Conservative.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
save
private
public:: rk_init
public:: rk_stage
public:: rk_time_integ
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup GvarLib_Runge_Kutta
!> @{
real(R_P), allocatable:: c1(:)   !< c1 coefficients of Runge-Kutta integration [1:S]
real(R_P), allocatable:: c2(:,:) !< c2 coefficients of Runge-Kutta integration [1:S,1:S]
real(R_P), allocatable:: c3(:)   !< c3 coefficients of Runge-Kutta integration [    1:S]
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Lib_Runge_KuttaPublicProcedure
  !> @{
  !> Subroutine for initializing Runge-Kutta coefficients.
  !> @note
  !> Currently 5 different schemes are available whose coefficients are (see \ref Lib_Runge_Kutta "definition" in Runge-Kutta
  !> library):
  !> - <b> 1 stage</b>: Forward Euler, one stage, \f$1^{st}\f$ order \n
  !>   - \f$ c_1 = \left[1\right]\f$
  !>   - \f$ c_2 = \left[ {\begin{array}{*{20}{c}} 0&0 \end{array}} \right] \f$
  !>   - \f$ c_3 = \left[0\right]\f$
  !> - <b> 2 stages</b>: Strong Stability Preserving, two stages, \f$2^{nd}\f$ order \n
  !>   - \f$ c_1 = \left[ {\begin{array}{*{20}{c}} \frac{1}{2}&\frac{1}{2} \end{array}} \right] \f$
  !>   - \f$ c_2 = \left[ {\begin{array}{*{20}{c}} 0&0\\ 1&0 \end{array}} \right] \f$
  !>   - \f$ c_3 = \left[ {\begin{array}{*{20}{c}} 0 \\ 1 \end{array}} \right] \f$
  !> - <b> 3 stages</b>: Strong Stability Preserving, three stages, \f$3^{rd}\f$ order \n
  !>   - \f$ c_1 = \left[ {\begin{array}{*{20}{c}} \frac{1}{6}&\frac{1}{6}&\frac{1}{3}  \end{array}} \right] \f$
  !>   - \f$ c_2 = \left[ {\begin{array}{*{20}{c}} 0&0&0\\ 1&0&0\\ \frac{1}{4}&\frac{1}{4}&0 \end{array}} \right] \f$
  !>   - \f$ c_3 = \left[ {\begin{array}{*{20}{c}} 0 \\ 1 \\ \frac{1}{2} \end{array}} \right] \f$
  !> - <b> 4 stages</b>: four stages, \f$4^{rd}\f$ order \n
  !>   - \f$ c_1 = \left[ {\begin{array}{*{20}{c}} \frac{55}{108}&\frac{1}{3}&\frac{1}{3}&\frac{1}{6}  \end{array}} \right] \f$
  !>   - \f$ c_2 = \left[ {\begin{array}{*{20}{c}} 0&0&0&0\\ 1&0&0&0\\ \frac{1}{4}&\frac{1}{2}&0&0\\ \frac{5}{18}&0&1&0
  !>               \end{array}} \right] \f$
  !>   - \f$ c_3 = \left[ {\begin{array}{*{20}{c}} 0 \\ 1 \\ \frac{3}{4} \\ \frac{23}{18} \end{array}} \right] \f$
  !> - <b> 5 stages</b>: Strong Stability Preserving, five stages, \f$4^{rd}\f$ order \n \n
  !>   - \f$ c_1 = \left[ {\begin{array}{*{20}{c}} 0.14681187618661&0.24848290924556&0.10425883036650&0.27443890091960&
  !>               0.22600748319395 \end{array}} \right] \f$
  !>   - \f$ c_2 = \left[ {\begin{array}{*{20}{c}}
  !> 0&0&0&0&0 \\ 0.39175222700392&0&0&0&0 \\ 0.21766909633821&0.36841059262959&0&0&0 \\ 0.08269208670950&
  !> 0.13995850206999&0.25189177424738&0&0 \\ 0.06796628370320&0.11503469844438&0.20703489864929&0.54497475021237&0
  !>               \end{array}} \right] \f$
  !>   - \f$ c_3 = \left[ {\begin{array}{*{20}{c}} 0\\0.39175222700392\\0.58607968896780\\0.47454236302687\\0.93501063100924
  !>               \end{array}} \right] \f$
  subroutine rk_init(S)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I1P), intent(IN):: S !< Number of stages used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! allocating variables
  if (allocated(c1)) deallocate(c1) ; allocate(c1(1:S    )) ; c1 = 0._R_P
  if (allocated(c2)) deallocate(c2) ; allocate(c2(1:S,1:S)) ; c2 = 0._R_P
  if (allocated(c3)) deallocate(c3) ; allocate(c3(    1:S)) ; c3 = 0._R_P
  ! initializing the coefficients
  select case(S)
  case(1_I1P)
    ! RK(1,1) Forward-Euler
    c1(1) = 1._R_P
  case(2_I1P)
    ! SSPRK(2,2)
    c1(1) = 0.5_R_P
    c1(2) = 0.5_R_P

    c2(2,1) = 1._R_P

    c3(2) = 1._R_P
  case(3_I1P)
    ! SSPRK(3,3)
    c1(1) = 1._R_P/6._R_P
    c1(2) = 1._R_P/6._R_P
    c1(3) = 2._R_P/3._R_P

    c2(2,1) = 1._R_P
    c2(3,1) = 0.25_R_P ; c2(3,2) = 0.25_R_P

    c3(2) = 1._R_P
    c3(3) = 0.5_R_P
  case(4_I1P)
    ! NSSPRK(4,4)
    c1(1) = 55._R_P/108._R_P
    c1(2) = 1._R_P/3._R_P
    c1(3) = 1._R_P/3._R_P
    c1(4) = 1._R_P/6._R_P

    c2(2,1) = 1._R_P
    c2(3,1) = 0.25_R_P       ; c2(3,2) = 0.5_R_P
    c2(4,1) = 5._R_P/18._R_P ; c2(4,2) = 0._R_P  ; c2(4,3) = 1._R_P

    c3(2) = 1._R_P
    c3(3) = 3._R_P/4._R_P
    c3(4) = 23._R_P/18._R_P
  case(5_I1P)
    ! SSPRK(5,4)
    c1(1) = 0.14681187618661_R_P
    c1(2) = 0.24848290924556_R_P
    c1(3) = 0.10425883036650_R_P
    c1(4) = 0.27443890091960_R_P
    c1(5) = 0.22600748319395_R_P

    c2(2,1)=0.39175222700392_R_P
    c2(3,1)=0.21766909633821_R_P;c2(3,2)=0.36841059262959_R_P
    c2(4,1)=0.08269208670950_R_P;c2(4,2)=0.13995850206999_R_P;c2(4,3)=0.25189177424738_R_P
    c2(5,1)=0.06796628370320_R_P;c2(5,2)=0.11503469844438_R_P;c2(5,3)=0.20703489864929_R_P;c2(5,4)=0.54497475021237_R_P

    c3(2) = 0.39175222700392_R_P
    c3(3) = 0.58607968896780_R_P
    c3(4) = 0.47454236302687_R_P
    c3(5) = 0.93501063100924_R_P
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rk_init

  !> Subroutine for computing \f$s1^{th}\f$ Runge-Kutta stage.
  !> @note For avoid the creation of temporary arrays (improving the efficiency) the array \b KS is declared as assumed-shape with
  !> only the lower bound defined. Its extention is [1:s1-1].
  pure subroutine rk_stage(s1,Dt,Un,KS,KS1)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),            intent(IN)::    s1     !< Current Runge-Kutta stage number.
  real(R_P),               intent(IN)::    Dt     !< Current time step.
  type(Type_Conservative), intent(IN)::    Un     !< Current integrating variable.
  type(Type_Conservative), intent(IN)::    KS(1:) !< Previous Runge-Kutta stages [1:s1-1].
  type(Type_Conservative), intent(INOUT):: KS1    !< Current Runge-Kutta stage.
  real(R_P)::                              ssum   !< Stages sum of partial densities.
  integer(I_P)::                           s      !< Species counter.
  integer(I_P)::                           ss     !< Stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !Ks1 = Un + Dt*(c2(s1,1:s1-1).dot.KS(1:s1-1)) ! overloaded operators form: non efficient!
  do s=1,size(Un%rs) ! loop over species
    ssum = 0._R_P
    do ss=1,s1-1 ! loop over stages
      ssum = ssum + c2(s1,ss)*KS(ss)%rs(s)
    enddo
    Ks1%rs(s) = Un%rs(s) + Dt*ssum
  enddo
  Ks1%rv%x = Un%rv%x + Dt*(dot_product(c2(s1,1:s1-1),KS(1:s1-1)%rv%x))
  Ks1%rv%y = Un%rv%y + Dt*(dot_product(c2(s1,1:s1-1),KS(1:s1-1)%rv%y))
  Ks1%rv%z = Un%rv%z + Dt*(dot_product(c2(s1,1:s1-1),KS(1:s1-1)%rv%z))
  Ks1%re = Un%re + Dt*(dot_product(c2(s1,1:s1-1),KS(1:s1-1)%re))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rk_stage

  !> Subroutine for computing Runge-Kutta one time step integration \f$u^{n+1}\f$.
  pure subroutine rk_time_integ(Dt,Un,KS,Unp1)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(R_P),               intent(IN)::    Dt    !< Current time step.
  type(Type_Conservative), intent(IN)::    Un    !< Current integrating variable.
  type(Type_Conservative), intent(IN)::    KS(:) !< Runge-Kutta stages.
  type(Type_Conservative), intent(INOUT):: Unp1  !< Time-integrated variable.
  real(R_P)::                              ssum  !< Stages sum of partial densities.
  integer(I_P)::                           s     !< Species counter.
  integer(I_P)::                           ss    !< Stages counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !Unp1 = Un + Dt*(c1.dot.KS) ! overloaded operators form: non efficient!
  do s=1,size(Un%rs) ! loop over species
    ssum = 0._R_P
    do ss=1,size(KS) ! loop over stages
      ssum = ssum + c1(ss)*KS(ss)%rs(s)
    enddo
    Unp1%rs(s) = Un%rs(s) + Dt*ssum
  enddo
  Unp1%rv%x = Un%rv%x + Dt*(dot_product(c1,KS%rv%x))
  Unp1%rv%y = Un%rv%y + Dt*(dot_product(c1,KS%rv%y))
  Unp1%rv%z = Un%rv%z + Dt*(dot_product(c1,KS%rv%z))
  Unp1%re = Un%re + Dt*(dot_product(c1,KS%re))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine rk_time_integ
  !> @}
endmodule Lib_Runge_Kutta
