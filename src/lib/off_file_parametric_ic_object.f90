!< OFF file parametric IC object definition and implementation.

module off_file_parametric_ic_object
!< OFF file parametric IC object definition and implementation.

!< The file parametric IC (Initial Conditions) is an ini-formatted file containing the definition of the initial conditions grid by
!< means of simple parametric shapes.
!< Its skeleton is the following
!<
!<```ini
!< [dimensions]
!< blocks_number = n
!<
!< [block_1]
!< id = id_number
!< level = ref_level
!< pressure = pressure
!< density = density
!< velocity = velocity_x velocity_y velocity_z
!<
!< [block_2]
!< id = id_number
!< level = ref_level
!< pressure = pressure
!< density = density
!< velocity = velocity_x velocity_y velocity_z
!< ...
!<```
!< where
!<+ `blocks_number = n` is integer;
!<+ `id = id_number` is the unique block ID;
!<+ `level = ref_level` is the refinement level;
!<+ `pressure = pressure` is real;
!<+ `density = density` is real;
!<+ `velocity = velocity_x velocity_y velocity_z` are 3 reals;
!<
!< [[file_parametric_ic_object]] provides standard API for loading and saving this file.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_object, only : block_object
use off_grid_dimensions_object, only : grid_dimensions_object
use off_file_object, only : file_object
use finer, only : file_ini
use flow, only : eos_compressible, primitive_compressible, primitive_to_conservative_compressible
use penf, only : I4P, I8P, R8P, str
use vecfor, only : ex, ey, ez

implicit none
private
public :: file_parametric_ic_object

type, extends(file_object) :: file_parametric_ic_object
   !< File parametric grid object class.
   type(file_ini) :: fini !< Parametric grid ini file handler.
   contains
      ! public methods
      procedure, pass(self) :: get_initial_conditions !< Return initial conditions from loaded parameters.
      procedure, pass(self) :: load_file              !< Load file.
      procedure, pass(self) :: save_file              !< Save file.
endtype file_parametric_ic_object

contains
   ! public methods
   subroutine get_initial_conditions(self, eos, blocks, file_name)
   !< Return initial conditions from loaded parameters.
   class(file_parametric_ic_object), intent(inout)        :: self            !< File object.
   type(eos_compressible),           intent(in)           :: eos             !< Equation of state.
   type(block_object),               intent(inout)        :: blocks(1:)      !< Blocks data.
   character(*),                     intent(in), optional :: file_name       !< File name.
   integer(I4P)                                           :: blocks_number   !< Blocks number.
   type(primitive_compressible)                           :: P               !< Primitive variables.
   real(R8P)                                              :: velocity_(3)    !< Velocity temporary array.
   integer(I4P)                                           :: b               !< Counter.

   if (.not.self%is_loaded) call self%load_file(file_name=file_name)
   call self%fini%get(section_name='dimensions', option_name='blocks_number', val=blocks_number, error=self%error%status)
   call self%error%check(message='failed to load [dimensions].(blocks_number)', is_severe=.true.)
   if (blocks_number>0.and.size(blocks, dim=1)>=blocks_number) then
      do b=1, blocks_number
         call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='pressure', &
                            val=P%pressure, error=self%error%status)
         call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(pressure)', is_severe=.true.)
         call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='density', &
                            val=P%density, error=self%error%status)
         call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(density)', is_severe=.true.)
         call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='velocity', &
                            val=velocity_, error=self%error%status)
         call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(velocity)', is_severe=.true.)
         P%velocity%x = velocity_(1)
         P%velocity%y = velocity_(2)
         P%velocity%z = velocity_(3)
         blocks(b)%cell%U = primitive_to_conservative_compressible(primitive=P, eos=eos)
      enddo
   endif
   endsubroutine get_initial_conditions

   subroutine load_file(self, file_name)
   !< Load file.
   class(file_parametric_ic_object), intent(inout)        :: self      !< File object.
   character(*),                     intent(in), optional :: file_name !< File name.

   call self%initialize(file_name=file_name)
   call self%fini%load(filename=self%file_name, error=self%error%status)
   call self%error%check(message='failed to load "'//trim(adjustl(self%file_name))//'"', is_severe=.true.)
   self%is_loaded = .true.
   endsubroutine load_file

   subroutine save_file(self, file_name)
   !< Save file.
   class(file_parametric_ic_object), intent(inout)        :: self       !< File object.
   character(*),                       intent(in), optional :: file_name  !< File name.

   error stop 'error: file_paremetric_ic_object%save_file to be implemented'
   endsubroutine save_file
endmodule off_file_parametric_ic_object
