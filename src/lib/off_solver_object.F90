#include "preprocessor_macros.h"
!< OFF solver object definition and implementation.

module off_solver_object
!< OFF solver object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_error_object, only : error_object
use finer, only : file_ini
use foodie, only : foodie_integrator_factory, integrand_object, integrator_object
use foreseer, only : riemann_solver_object, riemann_solver_compressible_hllc_id, foreseer_factory
use penf, only : I4P, R8P, str
use stringifor, only : string
use wenoof, only : interpolator_object, wenoof_create

implicit none
private
public :: solver_object

character(len=6), parameter :: INI_SECTION_NAME='solver' !< INI (config) file section name containing the solver parameters.

type :: solver_object
   !< Solver object class.
   !<
   !< Class designed to handle the general solver models parameters.
   type(error_object)                        :: error                         !< Errors handler.
   character(len=:),             allocatable :: time_integrator               !< Time integrator model: euler, rk2(3-10)...
   character(len=:),             allocatable :: convective_operator           !< Convective operator model: tvd, eno3, weno3...
   character(len=:),             allocatable :: diffusive_operator            !< Diffusive operator model: centered2...
   character(len=:),             allocatable :: turbulence_model              !< Turbulence_model: k-e, k-w, LES...
   type(string)                              :: riemann_solver_scheme         !< Riemann solver scheme: LLF, HLLC, Roe...
   real(R8P)                                 :: artificial_viscosity=0._R8P   !< Artifiical viscosity.
   real(R8P)                                 :: residuals_tolerance=0._R8P    !< Tolerance on residuals value.
   real(R8P)                                 :: pseudo_compressibility=0._R8P !< Pseudo compressibility.
   real(R8P)                                 :: chimera_forcing=0._R8P        !< Chimera forcing coefficient.
   class(integrator_object),     allocatable :: integrator                    !< Time integrator.
   class(interpolator_object),   allocatable :: interpolator                  !< WENO interpolator.
   class(riemann_solver_object), allocatable :: riemann_solver                !< Riemann solver.
   integer(I4P)                              :: gcu=0_I4P                     !< Number of ghost cells used (space accuracy).
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return a pretty-formatted description of solver parameters.
      procedure, pass(self) :: destroy        !< Destroy solver.
      procedure, pass(self) :: initialize     !< Initialize solver.
      procedure, pass(self) :: load_from_file !< Load from file.
      procedure, pass(self) :: save_into_file !< Save into file.
      ! operators
      generic :: assignment(=) => solver_assign_solver !< Overload `=`.
      ! private methods
      procedure, pass(self) :: set_convective_operator !< Set data of convective operator (interpolator).
      procedure, pass(lhs)  :: solver_assign_solver    !< Operator `=`.
endtype solver_object

contains
   ! public methods
   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of solver parameters.
   class(solver_object), intent(in)           :: self             !< Solver object.
   character(*),         intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable              :: desc             !< Description.
   character(len=:), allocatable              :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter                :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   if (allocated(self%integrator)) then
      desc=desc//prefix_//'Time integrator:'//NL
      desc=desc//self%integrator%description(prefix=prefix_//' ')//NL
   endif
   if (allocated(self%interpolator)) then
      desc=desc//prefix_//'Space interpolator:'//NL
      desc=desc//self%interpolator%description(prefix=prefix_//' ')//NL
   endif
   if(allocated(self%diffusive_operator   )) desc=desc//prefix_//'diffusive operator    : '//self%diffusive_operator//NL
   if(allocated(self%turbulence_model     )) desc=desc//prefix_//'turbulence model      : '//self%turbulence_model//NL
   if (allocated(self%riemann_solver)) then
      desc=desc//prefix_//'Riemann solver:'//NL
      desc=desc//self%riemann_solver%description(prefix=prefix_//' ')//NL
   endif

   desc=desc//prefix_//'artificial viscosity  : '//trim(str(self%artificial_viscosity))//NL
   desc=desc//prefix_//'residuals tolerance   : '//trim(str(self%residuals_tolerance))//NL
   desc=desc//prefix_//'pseudo compressibility: '//trim(str(self%pseudo_compressibility))//NL
   desc=desc//prefix_//'chimera forcing       : '//trim(str(self%chimera_forcing))
   endfunction description

   _ELEMENTAL_ subroutine destroy(self)
   !< Destroy solver.
   class(solver_object), intent(inout) :: self  !< Solver object.
   type(solver_object)                 :: fresh !< Fresh instance of solver object.

   self = fresh
   if (allocated(self%time_integrator)) deallocate(self%time_integrator)
   if (allocated(self%convective_operator)) deallocate(self%convective_operator)
   if (allocated(self%diffusive_operator)) deallocate(self%diffusive_operator)
   if (allocated(self%turbulence_model)) deallocate(self%turbulence_model)
   self%riemann_solver_scheme = ''
   if (allocated(self%interpolator)) then
      call self%interpolator%destroy
      deallocate(self%interpolator)
   endif
   if (allocated(self%riemann_solver)) then
      call self%riemann_solver%destroy
      deallocate(self%riemann_solver)
   endif
   endsubroutine destroy

   _ELEMENTAL_ subroutine initialize(self, solver)
   !< Initialize solver.
   class(solver_object), intent(inout)        :: self   !< Solver object.
   type(solver_object),  intent(in), optional :: solver !< Solver object.

   call self%destroy
   if (present(solver)) then
      self = solver
   else
      self%time_integrator       = 'runge_kutta_ssp_stages_1_order_1'
      self%convective_operator   = 'weno1'
      self%riemann_solver_scheme = 'llf'
   endif
   endsubroutine initialize

   subroutine load_from_file(self, fini, integrand_0, go_on_fail)
   !< Load from file.
   class(solver_object),    intent(inout)        :: self        !< Solver object.
   type(file_ini),          intent(in)           :: fini        !< Simulation parameters ini file handler.
   class(integrand_object), intent(in)           :: integrand_0 !< Initial conditions.
   logical,                 intent(in), optional :: go_on_fail  !< Go on if load fails.
   logical                                       :: go_on_fail_ !< Go on if load fails, local variable.
   character(999)                                :: buffer      !< Buffer string.

   go_on_fail_ = .true. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call fini%get(section_name=INI_SECTION_NAME, &
                 option_name='time_integrator', &
                 val=buffer,               &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(time_integrator)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%time_integrator = trim(adjustl(buffer))
   call foodie_integrator_factory(scheme=self%time_integrator, integrator=self%integrator, stages=1, &
                                  tolerance=1e2_R8P, iterations=1, autoupdate=.true., U=integrand_0)

   call fini%get(section_name=INI_SECTION_NAME,     &
                 option_name='convective_operator', &
                 val=buffer,                        &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(convective_operator)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%convective_operator = trim(adjustl(buffer))
   call self%set_convective_operator

   call fini%get(section_name=INI_SECTION_NAME,    &
                 option_name='diffusive_operator', &
                 val=buffer,                       &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(diffusive_operator)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%diffusive_operator = trim(adjustl(buffer))

   call fini%get(section_name=INI_SECTION_NAME,  &
                 option_name='turbulence_model', &
                 val=buffer,                     &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(turbulence_model)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%turbulence_model = trim(adjustl(buffer))

   call fini%get(section_name=INI_SECTION_NAME, &
                 option_name='riemann_solver',  &
                 val=buffer,                    &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(riemann_solver)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) then
      self%riemann_solver_scheme = trim(adjustl(buffer))
   else
      write(stderr, '(A)') 'warning: Riemann Solver not speciefied, set by default as LLF'
      self%riemann_solver_scheme = riemann_solver_compressible_hllc_id
   endif
   self%riemann_solver_scheme = self%riemann_solver_scheme%upper()
   call foreseer_factory(riemann_solver_scheme=self%riemann_solver_scheme%chars(), riemann_solver=self%riemann_solver)

   call fini%get(section_name=INI_SECTION_NAME,      &
                 option_name='artificial_viscosity', &
                 val=self%artificial_viscosity,      &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(artificial_viscosity)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME,     &
                 option_name='residuals_tolerance', &
                 val=self%residuals_tolerance,      &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(residuals_tolerance)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME,        &
                 option_name='pseudo_compressibility', &
                 val=self%pseudo_compressibility,      &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(pseudo_compressibility)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, &
                 option_name='chimera_forcing', &
                 val=self%chimera_forcing,      &
                 error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(chimera_forcing)', is_severe=.not.go_on_fail_)
   endsubroutine load_from_file

   subroutine save_into_file(self, fini)
   !< Save into file.
   class(solver_object), intent(inout) :: self !< Solver object.
   type(file_ini),       intent(inout) :: fini !< Simulation parameters ini file handler.

   call fini%add(section_name=INI_SECTION_NAME, option_name='time_integrator', val=self%time_integrator, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='convective_operator', val=self%convective_operator, &
                 error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='diffusive_operator', val=self%diffusive_operator, &
                 error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='turbulence_model', val=self%turbulence_model, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='artificial_viscosity', val=self%artificial_viscosity, &
                 error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='residuals_tolerance', val=self%residuals_tolerance, &
                 error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='pseudo_compressibility', val=self%pseudo_compressibility, &
                 error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='chimera_forcing', val=self%chimera_forcing, error=self%error%status)
   endsubroutine save_into_file

   ! private methods
   subroutine set_convective_operator(self)
   !< Set data of convective operator (interpolator).
   class(solver_object), intent(inout) :: self   !< Solver object.
   type(string)                        :: buffer !< Buffer string.

   if (allocated(self%convective_operator)) then
      buffer = trim(adjustl(self%convective_operator))
      buffer = buffer%upper()
      select case(buffer%chars())
      case('WENO1')
         self%gcu = 1
      case('WENO3')
         self%gcu = 2
      case('WENO5')
         self%gcu = 3
      case('WENO7')
         self%gcu = 4
      case('WENO9')
         self%gcu = 5
      case('WENO11')
         self%gcu = 6
      case('WENO13')
         self%gcu = 7
      case('WENO15')
         self%gcu = 8
      case('WENO17')
         self%gcu = 9
      case default
         error stop 'error: unknown convervective operator'
      endselect
   else
      self%gcu = 1
   endif
   call wenoof_create(interpolator_type='reconstructor-JS', S=self%gcu, interpolator=self%interpolator)
   endsubroutine set_convective_operator

   _PURE_ subroutine solver_assign_solver(lhs, rhs)
   !< Operator `=`.
   class(solver_object), intent(inout) :: lhs !< Left hand side.
   type(solver_object),  intent(in)    :: rhs !< Right hand side.

                                           lhs%error                  = rhs%error
   if (allocated(rhs%time_integrator))     lhs%time_integrator        = rhs%time_integrator
   if (allocated(rhs%convective_operator)) lhs%convective_operator    = rhs%convective_operator
   if (allocated(rhs%diffusive_operator))  lhs%diffusive_operator     = rhs%diffusive_operator
   if (allocated(rhs%turbulence_model))    lhs%turbulence_model       = rhs%turbulence_model
                                           lhs%riemann_solver_scheme  = rhs%riemann_solver_scheme
                                           lhs%artificial_viscosity   = rhs%artificial_viscosity
                                           lhs%residuals_tolerance    = rhs%residuals_tolerance
                                           lhs%pseudo_compressibility = rhs%pseudo_compressibility
                                           lhs%chimera_forcing        = rhs%chimera_forcing
   if (allocated(rhs%integrator)) then
      if (.not.allocated(lhs%integrator)) then
         allocate(lhs%integrator, mold=rhs%integrator)
      elseif (.not.same_type_as(lhs%integrator, rhs%integrator)) then
         deallocate(lhs%integrator)
         allocate(lhs%integrator, mold=rhs%integrator)
      endif
      lhs%integrator = rhs%integrator
   endif
   if (allocated(rhs%interpolator)) then
      if (.not.allocated(lhs%interpolator)) then
         allocate(lhs%interpolator, mold=rhs%interpolator)
      elseif (.not.same_type_as(lhs%interpolator, rhs%interpolator)) then
         deallocate(lhs%interpolator)
         allocate(lhs%interpolator, mold=rhs%interpolator)
      endif
      lhs%interpolator = rhs%interpolator
   endif
   if (allocated(rhs%riemann_solver)) then
      if (.not.allocated(lhs%riemann_solver)) then
         allocate(lhs%riemann_solver, mold=rhs%riemann_solver)
      elseif (.not.same_type_as(lhs%riemann_solver, rhs%riemann_solver)) then
         deallocate(lhs%riemann_solver)
         allocate(lhs%riemann_solver, mold=rhs%riemann_solver)
      endif
      lhs%riemann_solver = rhs%riemann_solver
   endif
   lhs%gcu = rhs%gcu
   endsubroutine solver_assign_solver
endmodule off_solver_object
