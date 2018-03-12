!< OFF file grid object definition and implementation.

module off_file_grid_object
!< OFF file grid object definition and implementation.

!< The file grid is an unformatted, stream file containing the nodes coordinates of the whole grid.
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
!< Secondary, it can be used also in the *parametric* ascii mode, essentially for loading parametric grid definitions. Its skeleton
!< is the following:
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
!< [[file_grid_object]] provides standard API for loading and saving this file.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_object, only : block_object
use off_grid_dimensions_object, only : grid_dimensions_object
use off_file_object, only : file_object
use finer, only : file_ini
use penf, only : I4P, I8P, R8P, str
use vecfor, only : ex, ey, ez

implicit none
private
public :: file_grid_object

type, extends(file_object) :: file_grid_object
   !< File grid object class.
   contains
      ! public methods
      procedure, pass(self) :: load_grid_dimensions_from_file !< Load the grid dimensions of all blocks from file.
      procedure, pass(self) :: load_nodes_from_file           !< Load nodes coordinates from file.
      procedure, pass(self) :: save_grid_dimensions_into_file !< Save the grid dimensions of all blocks into file.
      procedure, pass(self) :: save_nodes_into_file           !< Save nodes coordinates into file.
endtype file_grid_object

contains
   ! public methods
   subroutine load_grid_dimensions_from_file(self, grid_dimensions, file_name)
   !< Load the grid dimensions of all blocks from file.
   class(file_grid_object),      intent(inout)        :: self            !< File object.
   type(grid_dimensions_object), intent(inout)        :: grid_dimensions !< Grid dimensions off all blocks into file.
   character(*),                 intent(in), optional :: file_name       !< File name.
   integer(I4P)                                       :: dimensions_(3)  !< Block dimensions buffer.
   real(R8P)                                          :: extents_(6)     !< Block extents buffer.
   character(32)                                      :: faces_bc_(6)    !< Faces boundary conditions.
   integer(I8P)                                       :: id_             !< Block ID.
   integer(I4P)                                       :: level_          !< Block refinement level.
   logical                                            :: is_cartesian_   !< Cartesian block sentinel.
   logical                                            :: is_null_xyz_    !< Nullify direction(s) sentinel.
   character(len=:), allocatable                      :: emsg_suffix     !< Error message.
   integer(I4P)                                       :: b               !< Counter.

   if (present(file_name)) self%file_name = trim(adjustl(file_name))
   emsg_suffix = ' from file "'//self%file_name//'" in procedure "file_grid_object%load_grid_dimensions_from_file"'
   if (self%is_parametric) then
      call self%fini%load(filename=self%file_name, error=self%error%status)
      call self%fini%get(section_name='dimensions', option_name='blocks_number', val=grid_dimensions%blocks_number, &
                         error=self%error%status)
      call self%error%check(message='failed to load [dimensions].(blocks_number)'//emsg_suffix, is_severe=.true.)
      if (grid_dimensions%blocks_number>0) then
         call grid_dimensions%alloc
         do b=1, grid_dimensions%blocks_number
            call self%fini%get(section_name='dimensions', option_name='ghost_cells_number', &
                               val=grid_dimensions%block_signature(b)%gc, error=self%error%status)
            call self%error%check(message='failed to load [dimensions].(ghost_cells_number)'//emsg_suffix, is_severe=.true.)

            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='id', &
                               val=id_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(id)'//emsg_suffix, &
                                  is_severe=.true.)
            grid_dimensions%block_signature(b)%id = id_

            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='level', &
                               val=level_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(level)'//emsg_suffix, &
                                  is_severe=.true.)
            grid_dimensions%block_signature(b)%level = level_

            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='dimensions', &
                               val=dimensions_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(dimensions)'//emsg_suffix, &
                                  is_severe=.true.)
            grid_dimensions%block_signature(b)%ni = dimensions_(1)
            grid_dimensions%block_signature(b)%nj = dimensions_(2)
            grid_dimensions%block_signature(b)%nk = dimensions_(3)

            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='extents', &
                               val=extents_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b, no_sign=.true.))//'].(extents)'//emsg_suffix, &
                                  is_severe=.true.)
            grid_dimensions%block_signature(b)%emin = extents_(1) * ex +  extents_(2) * ey + extents_(3) * ez
            grid_dimensions%block_signature(b)%emax = extents_(4) * ex +  extents_(5) * ey + extents_(6) * ez

            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='boundary_conditions', &
                               val=faces_bc_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b, no_sign=.true.))//'].(boundary_conditions)'// &
                                  emsg_suffix, is_severe=.false.)
            if (self%error%status==0) grid_dimensions%block_signature(b)%faces_bc = faces_bc_

            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='is_cartesian', &
                               val=is_cartesian_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b, no_sign=.true.))//'].(is_cartesian)'//emsg_suffix, &
                                  is_severe=.false.)
            if (self%error%status==0) grid_dimensions%block_signature(b)%is_cartesian = is_cartesian_
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='is_null_x', &
                               val=is_null_xyz_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b, no_sign=.true.))//'].(is_null_x)'//emsg_suffix, &
                                  is_severe=.false.)
            if (self%error%status==0) grid_dimensions%block_signature(b)%is_null_x = is_null_xyz_
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='is_null_y', &
                               val=is_null_xyz_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b, no_sign=.true.))//'].(is_null_y)'//emsg_suffix, &
                                  is_severe=.false.)
            if (self%error%status==0) grid_dimensions%block_signature(b)%is_null_y = is_null_xyz_
            call self%fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='is_null_z', &
                               val=is_null_xyz_, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b, no_sign=.true.))//'].(is_null_z)'//emsg_suffix, &
                                  is_severe=.false.)
            if (self%error%status==0) grid_dimensions%block_signature(b)%is_null_z = is_null_xyz_
         enddo
      endif
   else
      call self%open_file(action='read')
      call grid_dimensions%load_from_file(file_unit=self%file_unit)
      call self%close_file
   endif
   endsubroutine load_grid_dimensions_from_file

   subroutine load_nodes_from_file(self, grid_dimensions, blocks, file_name)
   !< Load nodes coordinates from file.
   class(file_grid_object),      intent(inout)        :: self            !< File object.
   type(grid_dimensions_object), intent(in)           :: grid_dimensions !< Grid dimensions off all blocks into file.
   type(block_object),           intent(inout)        :: blocks(1:)      !< Blocks storage.
   character(*),                 intent(in), optional :: file_name       !< File name.
   integer(I4P)                                       :: b               !< Counter.

   call self%open_file(file_name=file_name, action='read')
   do b=1, size(blocks, dim=1)
      call blocks(b)%load_nodes_from_file(file_unit=self%file_unit, pos=grid_dimensions%iopos_block_nodes(b=b))
   enddo
   call self%close_file
   endsubroutine load_nodes_from_file

   subroutine save_grid_dimensions_into_file(self, grid_dimensions, file_name)
   !< Load the grid dimensions of all blocks into file.
   class(file_grid_object),      intent(inout)        :: self            !< File object.
   type(grid_dimensions_object), intent(in)           :: grid_dimensions !< Grid dimensions off all blocks into file.
   character(*),                 intent(in), optional :: file_name       !< File name.

   call self%open_file(file_name=file_name, action='write')
   call grid_dimensions%save_into_file(file_unit=self%file_unit)
   call self%close_file
   endsubroutine save_grid_dimensions_into_file

   subroutine save_nodes_into_file(self, grid_dimensions, blocks, file_name)
   !< Save nodes coordinates into file.
   class(file_grid_object),      intent(inout)        :: self            !< File object.
   type(grid_dimensions_object), intent(in)           :: grid_dimensions !< Grid dimensions off all blocks into file.
   type(block_object),           intent(inout)        :: blocks(1:)      !< Blocks storage.
   character(*),                 intent(in), optional :: file_name       !< File name.
   integer(I4P)                                       :: b               !< Counter.

   call self%open_file(file_name=file_name, action='write')
   do b=1, size(blocks, dim=1)
      call blocks(b)%save_nodes_into_file(file_unit=self%file_unit, pos=grid_dimensions%iopos_block_nodes(b=b))
   enddo
   call self%close_file
   endsubroutine save_nodes_into_file
endmodule off_file_grid_object
