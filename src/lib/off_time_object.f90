!< OFF time object definition and implementation.

module off_time_object
!< OFF time object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
use off_error_object, only : error_object
use finer, only : file_ini
use penf, only : I8P, R8P, str

implicit none
private
public :: time_object

character(len=4), parameter :: INI_SECTION_NAME='time' !< INI (config) file section name containing the time parameters.

type :: time_object
   !< Time object class.
   !<
   !< Class designed to handle the timing data of the simulation.
   type(error_object) :: error              !< Errors handler.
   logical            :: is_unsteady=.true. !< Type of simulation: unsteady or steady.
   integer(I8P)       :: n=0                !< Time steps counter.
   real(R8P)          :: t=0._R8P           !< Time.
   integer(I8P)       :: n_max=0            !< Maximum number of time steps computed.
   real(R8P)          :: t_max=0._R8P       !< Maximum time of integration, ignored if `n_max>0`.
   real(R8P)          :: CFL=0.3_R8P        !< Courant-Friedrichs-Lewy stability coefficient.
   real(R8P)          :: dt=0._R8P          !< Global time step.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return a pretty-formatted description of time parameters.
      procedure, pass(self) :: destroy        !< Destroy time.
      procedure, pass(self) :: eta            !< Return ETA-related informations string.
      procedure, pass(self) :: initialize     !< Initialize time.
      procedure, pass(self) :: is_the_end     !< Return true if the end of simulation is reached.
      procedure, pass(self) :: print_eta      !< Print to stdout ETA-related informations.
      procedure, pass(self) :: progress       !< Return the progress of simulation.
      procedure, pass(self) :: load_from_file !< Load from file.
      procedure, pass(self) :: save_into_file !< Save into file.
      procedure, pass(self) :: set_stop       !< Set simulation stop condition.
      procedure, pass(self) :: update_dt      !< Update time step for the last iterate.
      procedure, pass(self) :: update_n       !< Update time step counter.
      procedure, pass(self) :: update_t       !< Update time.
      ! operators
      generic :: assignment(=) => time_assign_time !< Overload `=`.
      ! private methods
      procedure, pass(lhs) :: time_assign_time !< Operator `=`.
endtype time_object

contains
   ! public methods
   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of time parameters.
   class(time_object), intent(in)           :: self             !< Time object.
   character(*),       intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable            :: desc             !< Description.
   character(len=:), allocatable            :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter              :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   if (self%is_unsteady) then
      desc = desc//prefix_//'timing: unsteady'//NL
   else
      desc = desc//prefix_//'timing: steady'//NL
   endif
   desc = desc//prefix_//'n     : '//trim(str(n=self%n))//NL
   desc = desc//prefix_//'t     : '//trim(str(n=self%t))//NL
   desc = desc//prefix_//'n_max : '//trim(str(n=self%n_max))//NL
   desc = desc//prefix_//'t_max : '//trim(str(n=self%t_max))//NL
   desc = desc//prefix_//'CFL   : '//trim(str(n=self%CFL))//NL
   desc = desc//prefix_//'dt    : '//trim(str(n=self%dt))
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy time.
   class(time_object), intent(inout) :: self  !< Time object.
   type(time_object)                 :: fresh !< Fresh instance of time object.

   self = fresh
   endsubroutine destroy

   pure function eta(self, prefix)
   !< Return ETA-related informations string.
   class(time_object), intent(in)           :: self             !< Time object.
   character(*),       intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable            :: eta              !< ETA-realted informations string.
   character(len=:), allocatable            :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter              :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   eta = ''
   eta = eta//prefix_//'n         : '//trim(str(n=self%n))//NL
   eta = eta//prefix_//'t         : '//trim(str(n=self%t))//NL
   eta = eta//prefix_//'dt        : '//trim(str(n=self%dt))//NL
   eta = eta//prefix_//'progress %: '//trim(str(fm='(F7.2)', n=self%progress()))
   endfunction eta

   elemental subroutine initialize(self, is_unsteady, n, t, n_max, t_max, CFL, dt, time)
   !< Initialize time.
   class(time_object), intent(inout)        :: self        !< Time object.
   logical,            intent(in), optional :: is_unsteady !< Type of simulation: unsteady or steady.
   integer(I8P),       intent(in), optional :: n           !< Time steps counter.
   real(R8P),          intent(in), optional :: t           !< Time.
   integer(I8P),       intent(in), optional :: n_max       !< Maximum number of time steps computed.
   real(R8P),          intent(in), optional :: t_max       !< Maximum time of integration, ignored if `n_max>0`.
   real(R8P),          intent(in), optional :: CFL         !< Courant-Friedrichs-Lewy stability coefficient.
   real(R8P),          intent(in), optional :: dt          !< Global time step.
   type(time_object),  intent(in), optional :: time        !< Time object.

   call self%destroy
   if (present(is_unsteady)) self%is_unsteady = is_unsteady
   if (present(n)) self%n = n
   if (present(t)) self%t = t
   if (present(n_max)) self%n_max = n_max
   if (present(t_max)) self%t_max = t_max
   if (present(CFL)) self%CFL = CFL
   if (present(dt)) self%dt = dt
   if (present(time)) self = time
   endsubroutine initialize

   elemental function  is_the_end(self) result(yes)
   !< Return true if the end of simulation is reached.
   class(time_object), intent(in) :: self !< Time object.
   logical                        :: yes  !< Test result.

   if (self%n_max <= 0) then
      yes = (self%t>=self%t_max)
   else
      yes = (self%n==self%n_max)
   endif
   endfunction is_the_end

   subroutine  print_eta(self, prefix, frequency)
   !< Print to stdout ETA-related informations.
   class(time_object), intent(in)           :: self       !< Time object.
   character(*),       intent(in), optional :: prefix     !< Prefixing string.
   integer(I8P),       intent(in), optional :: frequency  !< Print frequency (on time steps).
   integer(I8P)                                frequency_ !< Print frequency (on time steps), local variable.

   frequency_ = 1 ; if (present(frequency)) frequency_ = frequency
   if (mod(self%n, frequency) == 0) then
      write(stdout, '(A)')
      write(stdout, '(A)') self%eta(prefix=prefix)
   endif
   endsubroutine  print_eta

   elemental function progress(self) result(prog)
   !< Return the progress of simulation.
   class(time_object), intent(in) :: self !< Time object.
   real(R8P)                      :: prog !< Actual progress value.

   if (self%n_max>0) then
      prog = self%n*100/(self%n_max*1._R8P)
   elseif (self%t_max>0._R8P) then
      prog = 100*self%t/self%t_max
   else
      prog = 0._R8P
   endif
   endfunction progress

   subroutine load_from_file(self, fini, go_on_fail)
   !< Load from file.
   class(time_object), intent(inout)        :: self        !< Time object.
   type(file_ini),     intent(in)           :: fini        !< Simulation parameters ini file handler.
   logical,            intent(in), optional :: go_on_fail  !< Go on if load fails.
   logical                                  :: go_on_fail_ !< Go on if load fails, local variable.

   go_on_fail_ = .true. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call fini%get(section_name=INI_SECTION_NAME, option_name='is_unsteady', val=self%is_unsteady, error=self%error%status)
   if (.not.go_on_fail_) &
      call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(is_unsteady)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, option_name='n_max', val=self%n_max, error=self%error%status)
   if (.not.go_on_fail_)call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(n_max)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, option_name='t_max', val=self%t_max, error=self%error%status)
   if (.not.go_on_fail_)call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(t_max)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, option_name='cfl', val=self%CFL, error=self%error%status)
   if (.not.go_on_fail_) call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(cfl)', is_severe=.not.go_on_fail_)

   call fini%get(section_name=INI_SECTION_NAME, option_name='dt', val=self%dt, error=self%error%status)
   if (.not.go_on_fail_) call self%error%check(message='failed to load ['//INI_SECTION_NAME//'].(dt)', is_severe=.not.go_on_fail_)
   endsubroutine load_from_file

   subroutine save_into_file(self, fini)
   !< Save into file.
   class(time_object), intent(inout) :: self !< Time object.
   type(file_ini),     intent(inout) :: fini !< Simulation parameters ini file handler.

   call fini%add(section_name=INI_SECTION_NAME, option_name='is_unsteady', val=self%is_unsteady, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='n_max', val=self%n_max, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='t_max', val=self%t_max, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='cfl', val=self%cfl, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='dt', val=self%cfl, error=self%error%status)
   endsubroutine save_into_file

   elemental subroutine set_stop(self)
   !< Set simulation stop condition.
   class(time_object), intent(inout) :: self !< Time object.

   if (self%is_unsteady) then
      if (self%n_max > 0_I8P) then
         self%t_max = -1._R8P ! the value of t_max is ignored because n_max>0
      else
         self%n_max = -1_I8P
      endif
   else
      self%t_max = -1._R8P ! the value of t_max is ignored because steady simulation
   endif
   endsubroutine set_stop

   elemental subroutine update_dt(self)
   !< Update time step for the last iterate.
   class(time_object), intent(inout) :: self !< Time object.

   if (self%is_unsteady) then
      ! for an unsteady accurate simulation each cell is updated by means of global minimum time step
      ! control for the last iterate
      if (self%n_max <= 0) then
         if ((self%t + self%dt) > self%t_max) then
            ! the global minimum time step is so high that the last iteration will go over t_max
            ! it is decreased in order to achieve exactly t_max
            self%dt = abs(self%t_max - self%t)
         endif
      endif
   endif
   endsubroutine update_dt

   elemental subroutine update_n(self)
   !< Update time step counter.
   class(time_object), intent(inout) :: self !< Time object.

   self%n = self%n + 1
   endsubroutine update_n

   elemental subroutine update_t(self)
   !< Update time.
   class(time_object), intent(inout) :: self !< Time object.

   self%t = self%t + self%dt
   endsubroutine update_t

   ! private methods
   pure subroutine time_assign_time(lhs, rhs)
   !< Operator `=`.
   class(time_object), intent(inout) :: lhs !< Left hand side.
   type(time_object),  intent(in)    :: rhs !< Right hand side.

   lhs%error       = rhs%error
   lhs%is_unsteady = rhs%is_unsteady
   lhs%n           = rhs%n
   lhs%t           = rhs%t
   lhs%n_max       = rhs%n_max
   lhs%t_max       = rhs%t_max
   lhs%CFL         = rhs%CFL
   lhs%dt          = rhs%dt
   endsubroutine time_assign_time
endmodule off_time_object
