!< OFF file solution object definition and implementation.

module off_file_solution_object
!< OFF file solution object definition and implementation.

!< The file solution is an unformatted, stream file containing the cells solution variable of the whole grid.
!< Its skeleton is the following
!<
!<```
!< # header
!< blocks_number
!< # for each block
!< id, level, gc, ni, nj, nk
!< # core
!< # for each block (for all nodes of block)
!< node%vertex%x, node%vertex%y, node%vertex%z
!<```
!<
!< Secondary, it can be used also in the *parametric* ascii mode, essentially for loading parametric initial conditions. Its
!< skeleton is the following:
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
!< [[file_solution_object]] provides standard API for loading and saving this file.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_object, only : block_object
use off_grid_dimensions_object, only : grid_dimensions_object
use off_file_object, only : file_object
use finer, only : file_ini
use flow, only : eos_compressible, primitive_compressible, primitive_to_conservative_compressible
use penf, only : I4P, I8P, R8P, str

implicit none
private
public :: file_solution_object

type, extends(file_object) :: file_solution_object
   !< File solution object class.
   contains
      ! public methods
      procedure, pass(self) :: load_conservatives_from_file !< Load conservative variables from file.
      procedure, pass(self) :: save_conservatives_into_file !< Save conservative variables into file.
endtype file_solution_object

contains
   ! public methods
   subroutine load_conservatives_from_file(self, grid_dimensions, blocks, file_name)
   !< Load conservative variables from file.
   class(file_solution_object),  intent(inout)        :: self            !< File object.
   type(grid_dimensions_object), intent(in)           :: grid_dimensions !< Grid dimensions off all blocks into file.
   type(block_object),           intent(inout)        :: blocks(1:)      !< Blocks storage.
   character(*),                 intent(in), optional :: file_name       !< File name.
   integer(I4P)                                       :: blocks_number   !< Blocks number.
   type(primitive_compressible)                       :: P               !< Primitive variables.
   real(R8P)                                          :: velocity_(3)    !< Velocity temporary array.
   character(len=:), allocatable                      :: emsg_suffix     !< Error message.
   integer(I4P)                                       :: b               !< Counter.

   if (present(file_name)) self%file_name = trim(adjustl(file_name))
   emsg_suffix = ' from file "'//self%file_name//'" in procedure "file_solution_object%load_conservative_from_file"'
   if (self%is_parametric) then
      call self%fini%load(filename=self%file_name, error=self%error%status)
      call self%fini%get(section_name='dimensions', option_name='blocks_number', val=blocks_number, error=self%error%status)
      call self%error%check(message='failed to load [dimensions].(blocks_number)'//emsg_suffix, is_severe=.true.)
      if (blocks_number>0.and.size(blocks, dim=1)>=blocks_number) then
         do b=1, blocks_number
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='pressure', &
                               val=P%pressure, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(pressure)'//emsg_suffix, &
                                  is_severe=.true.)
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='density', &
                               val=P%density, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(density)'//emsg_suffix, &
                                  is_severe=.true.)
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='velocity', &
                               val=velocity_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(velocity)'//emsg_suffix, &
                                  is_severe=.true.)
            P%velocity%x = velocity_(1)
            P%velocity%y = velocity_(2)
            P%velocity%z = velocity_(3)
            blocks(b)%cell%P = P
            call blocks(b)%primitive_to_conservative
         enddo
      endif
   else
      call self%open_file(action='read')
      do b=1, size(blocks, dim=1)
         ! call blocks(b)%load_conservatives_from_file(file_unit=self%file_unit, pos=grid_dimensions%iopos_block_nodes(b=b))
      enddo
      call self%close_file
   endif
   endsubroutine load_conservatives_from_file

   subroutine save_conservatives_into_file(self, grid_dimensions, blocks, file_name)
   !< Save conservative variables into file.
   class(file_solution_object),  intent(inout)        :: self            !< File object.
   type(grid_dimensions_object), intent(in)           :: grid_dimensions !< Grid dimensions off all blocks into file.
   type(block_object),           intent(inout)        :: blocks(1:)      !< Blocks storage.
   character(*),                 intent(in), optional :: file_name       !< File name.
   integer(I4P)                                       :: b               !< Counter.

   call self%open_file(file_name=file_name, action='write')
   do b=1, size(blocks, dim=1)
      ! call blocks(b)%save_conservatives_into_file(file_unit=self%file_unit, pos=grid_dimensions%iopos_block_nodes(b=b))
   enddo
   call self%close_file
   endsubroutine save_conservatives_into_file
endmodule off_file_solution_object
