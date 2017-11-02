!< OFF file grid object definition and implementation.

module off_file_grid_object
!< OFF file grid object definition and implementation.

!< The file grid is an unformatted, stream file containing the nodes coordinates of the whole grid.
!< It skeleton is the following
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
!< [[file_grid_object]] provides standard API for loading and saving this file.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_object, only : block_object
use off_grid_dimensions_object, only : grid_dimensions_object
use off_file_object, only : file_object
use penf, only : I4P

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

   call self%open_file(file_name=file_name, action='read')
   call grid_dimensions%load_from_file(file_unit=self%file_unit)
   call self%close_file
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
