!< OFF file parametric grid object definition and implementation.

module off_file_parametric_grid_object
!< OFF file parametric grid object definition and implementation.

!< The file parametric grid is an ini-formatted file containing the definition of the grid by means of simple parametric shapes.
!< Its skeleton is the following
!<
!<```ini
!< [dimensions]
!< blocks_number = n
!< ghost_cells_number = imin imax jmin jmax kmin kmax
!<
!< [block_1]
!< id = id_number
!< level = ref_level
!< dimensions = ni nj nk
!< extents = xmin ymin zmin xmax ymax zmax
!< boundary_conditions = xmin ymin zmin xmax ymax zmax
!<
!< [block_2]
!< id = id_number
!< level = ref_level
!< dimensions = ni nj nk
!< extents = xmin ymin zmin xmax ymax zmax
!< boundary_conditions = xmin ymin zmin xmax ymax zmax
!< ...
!<```
!< where
!<+ `blocks_number = n` is integer;
!<+ `ghost_cells_number = imin imax jmin jmax kmin kmax` are 6 integers;
!<+ `id = id_number` is the unique block ID;
!<+ `level = ref_level` is the refinement level;
!<+ `dimensions = ni nj nk` are 3 integers;
!<+ `extents = xmin ymin zmin xmax ymax zmax` are 6 reals;
!<+ `boundary_conditions = xmin ymin zmin xmax ymax zmax` are 6 strings with BC code as defined in [[boundary_conditions_object]].
!<
!< [[file_parametric_grid_object]] provides standard API for loading and saving this file.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_object, only : block_object
use off_grid_dimensions_object, only : grid_dimensions_object
use off_file_object, only : file_object
use finer, only : file_ini
use penf, only : I4P, I8P, R8P, str
use vecfor, only : ex, ey, ez

implicit none
private
public :: file_parametric_grid_object

type, extends(file_object) :: file_parametric_grid_object
   !< File parametric grid object class.
   type(file_ini) :: fini !< Parametric grid ini file handler.
   contains
      ! public methods
      procedure, pass(self) :: get_grid_dimensions !< Return grid dimensions from loaded paramters.
      procedure, pass(self) :: load_file           !< Load file.
      procedure, pass(self) :: save_file           !< Save file.
endtype file_parametric_grid_object

contains
   ! public methods
   subroutine get_grid_dimensions(self, grid_dimensions)
   !< Return grid dimensions from loaded paramters.
   class(file_parametric_grid_object), intent(inout) :: self            !< File object.
   type(grid_dimensions_object),       intent(out)   :: grid_dimensions !< Grid dimensions off all blocks.
   integer(I4P)                                      :: dimensions_(3)  !< Block dimensions buffer.
   real(R8P)                                         :: extents_(6)     !< Block extents buffer.
   integer(I8P)                                      :: id_             !< Block ID.
   integer(I4P)                                      :: level_          !< Block refinement level.
   integer(I4P)                                      :: b               !< Counter.

   if (self%is_loaded) then
      call self%fini%get(section_name='dimensions', option_name='blocks_number', val=grid_dimensions%blocks_number, &
                         error=self%error%status)
      call self%error%check(message='failed to load [dimensions].(blocks_number)', is_severe=.true.)
      if (grid_dimensions%blocks_number>0) then
         call grid_dimensions%alloc
         do b=1, grid_dimensions%blocks_number
            call self%fini%get(section_name='dimensions', option_name='ghost_cells_number', &
                               val=grid_dimensions%block_signature(b)%gc, error=self%error%status)
            call self%error%check(message='failed to load [dimensions].(ghost_cells_number)', is_severe=.true.)
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='id', &
                               val=id_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(id)', is_severe=.true.)
            grid_dimensions%block_signature(b)%id = id_
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='level', &
                               val=level_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(level)', is_severe=.true.)
            grid_dimensions%block_signature(b)%level = level_
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='dimensions', &
                               val=dimensions_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(dimensions)', is_severe=.true.)
            grid_dimensions%block_signature(b)%ni = dimensions_(1)
            grid_dimensions%block_signature(b)%nj = dimensions_(2)
            grid_dimensions%block_signature(b)%nk = dimensions_(3)
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='extents', &
                               val=extents_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b, no_sign=.true.))//'].(extents)', is_severe=.true.)
            grid_dimensions%block_signature(b)%emin = extents_(1) * ex +  extents_(2) * ey + extents_(3) * ez
            grid_dimensions%block_signature(b)%emax = extents_(4) * ex +  extents_(5) * ey + extents_(6) * ez
         enddo
      endif
   endif
   endsubroutine get_grid_dimensions

   subroutine load_file(self, file_name)
   !< Load file.
   class(file_parametric_grid_object), intent(inout)        :: self      !< File object.
   character(*),                       intent(in), optional :: file_name !< File name.

   call self%initialize(file_name=file_name)
   call self%fini%load(filename=self%file_name, error=self%error%status)
   call self%error%check(message='failed to load "'//trim(adjustl(self%file_name))//'"', is_severe=.true.)
   self%is_loaded = .true.
   endsubroutine load_file

   subroutine save_file(self, file_name)
   !< Save file.
   class(file_parametric_grid_object), intent(inout)        :: self       !< File object.
   character(*),                       intent(in), optional :: file_name  !< File name.

   error stop 'error: file_paremetric_grid_object%save_file to be implemented'
   endsubroutine save_file
endmodule off_file_parametric_grid_object
