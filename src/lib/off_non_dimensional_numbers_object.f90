!< OFF non dimensional numbers object definition and implementation.

module off_non_dimensional_numbers_object
!< OFF non dimensional numbers object definition and implementation.

use off_error_object, only : error_object
use finer, only : file_ini
use penf, only : I4P, R8P, str

implicit none
private
public :: non_dimensional_numbers_object

character(len=23), parameter :: INI_SECTION_NAME='non_dimensional_numbers' !< INI (config) file section name containing the values
                                                                           !< of non dimensional numbers.

type :: non_dimensional_numbers_object
   !< Non dimensional numbers object class.
   type(error_object) :: error !< Errors handler.
   ! non dimensional numbers imposed
   real(R8P) :: Re=1._R8P !< \(\rm{Re}=\frac{\rho_0 v_0 L_0}{\mu_0}\) Reynolds number.
   real(R8P) :: Fr=1._R8P !< \(\rm{Fr}=\sqrt{\frac{v_0^2}{f_0 L_0}}\) Froude number.
   real(R8P) :: We=1._R8P !< \(\rm{We}=\sqrt{\frac{\rho_0 v_0^2 L_0}{\sigma_0}}\) Weber number.
   real(R8P) :: Ma=1._R8P !< \(\rm{We}=\sqrt{\frac{v_0}{a_0}}\) Mach number.
   real(R8P) :: Pr=1._R8P !< \(\rm{Pr}=\frac{\mu_0 c_p}{k_0}\) Prandtl number.
   ! reference values imposed
   real(R8P) :: L0=1._R8P !< Reference length.
   real(R8P) :: r0=1._R8P !< Reference density.
   real(R8P) :: v0=1._R8P !< Reference velocity.
   real(R8P) :: c0=1._R8P !< Reference specific heats (\f$cp_0 = cv_0 = R_0 = c_0\f$).
   ! derived reference values
   real(R8P) :: mu0=1._R8P !< \(\mu_0= \frac{\rho_0 v_0 L_0}{\rm{Re}}\) Reference dynamic viscosity.
   real(R8P) :: f0 =1._R8P !< \(f_0= \frac{v_0^2}{L_0 \rm{Fr}^2}\) Reference specific force.
   real(R8P) :: k0 =1._R8P !< \(k_0= \frac{\mu_0 c_0}{\rm{Pr}}\) Reference thermal conductivity coefficient.
   real(R8P) :: Dt0=1._R8P !< \(Dt_0=\frac{L_0}{v_0}\) Reference time interval.
   real(R8P) :: p0 =1._R8P !< \(p_0=\rho_0 v_0^2\) Reference pressure.
   real(R8P) :: a0 =1._R8P !< \(a_0=v_0\) Reference speed of sound.
   real(R8P) :: T0 =1._R8P !< \(T_0=\frac{v_0^2}{c_0}\) Reference temperature.
   real(R8P) :: E0 =1._R8P !< \(E_0=v_0^2\) Reference specific energy.
   real(R8P) :: q0 =1._R8P !< \(q_0=\frac{v_0^3}{L_0}\) Reference specific heat.
   ! equations coefficients
   real(R8P) :: Re_inv=1._R8P   !< \(\frac{1}{\rm{Re}}\) Inverse of Reynolds number (coefficient of viscous terms).
   real(R8P) :: Fr2_inv=1._R8P  !< \(\frac{1}{\rm{Fr}^2}\) Inverse of square of Froude number (coefficient of volume forces).
   real(R8P) :: PrRe_inv=1._R8P !< \(\frac{1}{\rm{Pr Re}}\) Inverse of Prandtl and Reynolds numbers (coef. of condution terms).
   contains
      ! public methods
      procedure, pass(self) :: compute_reference_values !< Compute derived reference values.
      procedure, pass(self) :: description              !< Return a pretty-formatted description of the numbers.
      procedure, pass(self) :: destroy                  !< Destroy numbers value.
      procedure, pass(self) :: initialize               !< Initialize numbers value.
      procedure, pass(self) :: load_from_file           !< Load from file.
      procedure, pass(self) :: save_into_file           !< Save into file.
      ! operators
      generic :: assignment(=) => adim_assign_adim !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: adim_assign_adim !< Operator `=`.
endtype non_dimensional_numbers_object

contains
   ! public methods
   elemental subroutine compute_reference_values(self)
   !< Compute derived reference values.
   class(non_dimensional_numbers_object), intent(inout) :: self !< Non dimensional numbers.

   self%mu0  = (self%r0 * self%v0 * self%L0) / self%Re
   self%f0   = (self%v0 * self%v0) / (self%L0 * self%Fr * self%Fr)
   self%k0   = (self%mu0 * self%c0) / self%Pr
   self%Dt0  = self%L0 / self%v0
   self%p0   = self%r0 * self%v0 * self%v0
   self%a0   = self%v0
   self%T0   = (self%v0 * self%v0) / self%c0
   self%E0   = self%v0 * self%v0
   self%q0   = (self%v0 * self%v0 * self%v0) / self%L0

   self%Re_inv   = 1._R8P / self%Re
   self%Fr2_inv  = 1._R8P / (self%Fr * self%Fr)
   self%PrRe_inv = 1._R8P / (self%Pr * self%Re)
   endsubroutine compute_reference_values

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the numbers.
   class(non_dimensional_numbers_object), intent(in)           :: self             !< Non dimensional numbers.
   character(*),                          intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                               :: desc             !< Description.
   character(len=:), allocatable                               :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter                                 :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//prefix_//'Re      : '//trim(str(n=self%Re      ))//NL
   desc = desc//prefix_//'Fr      : '//trim(str(n=self%Fr      ))//NL
   desc = desc//prefix_//'We      : '//trim(str(n=self%We      ))//NL
   desc = desc//prefix_//'Ma      : '//trim(str(n=self%Ma      ))//NL
   desc = desc//prefix_//'Pr      : '//trim(str(n=self%Pr      ))//NL
   desc = desc//prefix_//'L0      : '//trim(str(n=self%L0      ))//NL
   desc = desc//prefix_//'r0      : '//trim(str(n=self%r0      ))//NL
   desc = desc//prefix_//'v0      : '//trim(str(n=self%v0      ))//NL
   desc = desc//prefix_//'c0      : '//trim(str(n=self%c0      ))//NL
   desc = desc//prefix_//'mu0     : '//trim(str(n=self%mu0     ))//NL
   desc = desc//prefix_//'f0      : '//trim(str(n=self%f0      ))//NL
   desc = desc//prefix_//'k0      : '//trim(str(n=self%k0      ))//NL
   desc = desc//prefix_//'Dt0     : '//trim(str(n=self%Dt0     ))//NL
   desc = desc//prefix_//'p0      : '//trim(str(n=self%p0      ))//NL
   desc = desc//prefix_//'a0      : '//trim(str(n=self%a0      ))//NL
   desc = desc//prefix_//'T0      : '//trim(str(n=self%T0      ))//NL
   desc = desc//prefix_//'E0      : '//trim(str(n=self%E0      ))//NL
   desc = desc//prefix_//'q0      : '//trim(str(n=self%q0      ))//NL
   desc = desc//prefix_//'Re_inv  : '//trim(str(n=self%Re_inv  ))//NL
   desc = desc//prefix_//'Fr2_inv : '//trim(str(n=self%Fr2_inv ))//NL
   desc = desc//prefix_//'PrRe_inv: '//trim(str(n=self%PrRe_inv))
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy numbers value.
   class(non_dimensional_numbers_object), intent(inout) :: self  !< Non dimensional numbers.
   type(non_dimensional_numbers_object)                 :: fresh !< Fresh instance of non dimensional numbers.

   self = fresh
   endsubroutine destroy

   elemental subroutine initialize(self, adimensionals)
   !< Initialize numbers value.
   class(non_dimensional_numbers_object), intent(inout)        :: self          !< Non dimensional numbers.
   type(non_dimensional_numbers_object),  intent(in), optional :: adimensionals !< Non dimensional numbers values.

   call self%destroy
   if (present(adimensionals)) self = adimensionals
   endsubroutine initialize

   subroutine load_from_file(self, fini, go_on_fail)
   !< Load from file.
   class(non_dimensional_numbers_object), intent(inout)        :: self        !< Non dimensional numbers.
   type(file_ini),                        intent(in)           :: fini        !< Simulation parameters ini file handler.
   logical,                               intent(in), optional :: go_on_fail  !< Go on if load fails.
   logical                                                     :: go_on_fail_ !< Go on if load fails, local variable.

   go_on_fail_ = .true. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call fini%get(section_name=INI_SECTION_NAME, option_name='Re', val=self%Re, error=self%error%status)
   if (.not.go_on_fail_) call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(Re)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, option_name='Fr', val=self%Fr, error=self%error%status)
   if (.not.go_on_fail_) call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(Fr)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, option_name='We', val=self%We, error=self%error%status)
   if (.not.go_on_fail_) call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(We)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, option_name='Ma', val=self%Ma, error=self%error%status)
   if (.not.go_on_fail_) call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(Me)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, option_name='Pr', val=self%Pr, error=self%error%status)
   if (.not.go_on_fail_) call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(Pr)', is_severe=.not.go_on_fail_)

   call self%compute_reference_values
   endsubroutine load_from_file

   subroutine save_into_file(self, fini)
   !< Save into file.
   class(non_dimensional_numbers_object), intent(inout) :: self !< Simulation parameters.
   type(file_ini),                        intent(inout) :: fini !< Simulation parameters ini file handler.

   call fini%add(section_name=INI_SECTION_NAME, option_name='Re', val=self%Re, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='Fr', val=self%Fr, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='We', val=self%We, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='Ma', val=self%Ma, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='Pr', val=self%Pr, error=self%error%status)
   endsubroutine save_into_file

   ! private methods
   pure subroutine adim_assign_adim(lhs, rhs)
   !< Operator `=`.
   class(non_dimensional_numbers_object), intent(inout) :: lhs !< Left hand side.
   type(non_dimensional_numbers_object),  intent(in)    :: rhs !< Right hand side.

   lhs%error    = rhs%error
   lhs%Re       = rhs%Re
   lhs%Fr       = rhs%Fr
   lhs%We       = rhs%We
   lhs%Ma       = rhs%Ma
   lhs%Pr       = rhs%Pr
   lhs%L0       = rhs%L0
   lhs%r0       = rhs%r0
   lhs%v0       = rhs%v0
   lhs%c0       = rhs%c0
   lhs%mu0      = rhs%mu0
   lhs%f0       = rhs%f0
   lhs%k0       = rhs%k0
   lhs%Dt0      = rhs%Dt0
   lhs%p0       = rhs%p0
   lhs%a0       = rhs%a0
   lhs%T0       = rhs%T0
   lhs%E0       = rhs%E0
   lhs%q0       = rhs%q0
   lhs%Re_inv   = rhs%Re_inv
   lhs%Fr2_inv  = rhs%Fr2_inv
   lhs%PrRe_inv = rhs%PrRe_inv
   endsubroutine adim_assign_adim
endmodule off_non_dimensional_numbers_object
