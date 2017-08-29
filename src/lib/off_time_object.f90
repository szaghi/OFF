!< OFF time object definition and implementation.

module off_time_object
!< OFF time object definition and implementation.

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
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return a pretty-formatted description of time parameters.
      procedure, pass(self) :: destroy        !< Destroy time.
      procedure, pass(self) :: initialize     !< Initialize time.
      procedure, pass(self) :: is_the_end     !< Return true if the end of simulation is reached.
      procedure, pass(self) :: progress       !< Return the progress of simulation.
      procedure, pass(self) :: load_from_file !< Load from file.
      procedure, pass(self) :: save_into_file !< Save into file.
      procedure, pass(self) :: set_stop       !< Set simulation stop condition.
      procedure, pass(self) :: update         !< Update time.
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
   desc = desc//prefix_//'CFL   : '//trim(str(n=self%CFL))
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy time.
   class(time_object), intent(inout) :: self  !< Time object.
   type(time_object)                 :: fresh !< Fresh instance of time object.

   self = fresh
   endsubroutine destroy

   elemental subroutine initialize(self)
   !< Initialize time.
   class(time_object), intent(inout) :: self !< Time object.

   call self%destroy
   endsubroutine initialize

   elemental function  is_the_end(self) result(yes)
   !< Return true if the end of simulation is reached.
   class(time_object), intent(in) :: self !< Time object.
   logical                        :: yes  !< Test result.

   yes = ((self%t==self%t_max).or.(self%n==self%n_max))
   endfunction is_the_end

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
   endsubroutine load_from_file

   subroutine save_into_file(self, fini)
   !< Save into file.
   class(time_object), intent(inout) :: self !< Time object.
   type(file_ini),     intent(inout) :: fini !< Simulation parameters ini file handler.

   call fini%add(section_name=INI_SECTION_NAME, option_name='is_unsteady', val=self%is_unsteady, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='n_max', val=self%n_max, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='t_max', val=self%t_max, error=self%error%status)
   call fini%add(section_name=INI_SECTION_NAME, option_name='cfl', val=self%cfl, error=self%error%status)
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

   elemental subroutine update(self, global_min_dt)
   !< Update time.
   class(time_object), intent(inout) :: self          !< Time object.
   real(R8P),          intent(inout) :: global_min_dt !< Global (all processes/images, all blocks) minimum time step.

   if (self%is_unsteady) then
      ! for an unsteady accurate simulation each cell is updated by means of global minimum time step
      ! control for the last iterate
      if (self%n_max <= 0) then
         if ((self%t + global_min_dt) > self%t_max) then
            ! the global minimum time step is so high that the last iteration will go over t_max
            ! it is decreased in order to achieve exactly t_max
            global_min_dt = abs(self%t_max - self%t)
         endif
      endif
      self%t = self%t + global_min_dt
   endif
   endsubroutine update

   ! private methods
   pure subroutine time_assign_time(lhs, rhs)
   !< Operator `=`.
   class(time_object), intent(inout) :: lhs !< Left hand side.
   type(time_object),  intent(in)    :: rhs !< Right hand side.

   lhs%error       = rhs%error
   lhs%n           = rhs%n
   lhs%t           = rhs%t
   lhs%n_max       = rhs%n_max
   lhs%t_max       = rhs%t_max
   lhs%CFL         = rhs%CFL
   lhs%is_unsteady = rhs%is_unsteady
   endsubroutine time_assign_time
endmodule off_time_object
