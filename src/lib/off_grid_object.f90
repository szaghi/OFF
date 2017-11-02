!< OFF grid object definition and implementation.

module off_grid_object
!< OFF grid object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_object, only : block_object
use off_file_grid_object, only : file_grid_object
use off_file_parametric_grid_object, only : file_parametric_grid_object
use off_grid_dimensions_object, only : grid_dimensions_object
use penf, only : I4P, str
use stringifor, only : string

implicit none
private
public :: grid_object

type :: grid_object
   !< grid object class.
   !<
   !< [[grid_object]] is a container for all grid data.
   type(file_grid_object)            :: file_grid            !< Grid file handler.
   type(file_parametric_grid_object) :: file_parametric_grid !< Parametric grid file handler.
   type(grid_dimensions_object)      :: grid_dimensions      !< Grid dimensions.
   type(block_object), allocatable   :: blocks(:)            !< Blocks list.
   contains
      ! public methods
      procedure, pass(self) :: allocate_blocks !< Allocate blocks accordingly to grid dimensions.
      procedure, pass(self) :: description     !< Return a pretty-formatted description of the grid.
      procedure, pass(self) :: destroy         !< Destroy grid.
      procedure, pass(self) :: initialize      !< Initialize grid.
      procedure, pass(self) :: load_from_file  !< Load from grid file.
      procedure, pass(self) :: save_into_file  !< Save into grid file.
endtype grid_object

contains
   ! public methods
   subroutine allocate_blocks(self)
   !< Allocate blocks accordingly to grid dimensions.
   class(grid_object), intent(inout) :: self !< Grid.
   integer(I4P)                      :: b    !< Counter.

   if (self%grid_dimensions%blocks_number > 0) then
      if (allocated(self%blocks)) then
         call self%blocks%destroy
         deallocate(self%blocks)
      endif
      allocate(self%blocks(1:self%grid_dimensions%blocks_number))
      do b=1, self%grid_dimensions%blocks_number
         call self%blocks(b)%initialize(signature=self%grid_dimensions%block_signature(b))
      enddo
   endif
   endsubroutine allocate_blocks

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the grid.
   class(grid_object), intent(in)           :: self             !< Grid.
   character(*),       intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable            :: desc             !< Description.
   character(len=:), allocatable            :: prefix_          !< Prefixing string, local variable.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//self%grid_dimensions%description(prefix=prefix_//'  ')
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy grid.
   class(grid_object), intent(inout) :: self !< Grid.

   call self%file_grid%destroy
   call self%grid_dimensions%destroy
   if (allocated(self%blocks)) then
      call self%blocks%destroy
      deallocate(self%blocks)
   endif
   endsubroutine destroy

   subroutine initialize(self, load_from_file, is_parametric, file_name, grid)
   !< Initialize grid.
   class(grid_object), intent(inout)        :: self            !< Grid.
   logical,            intent(in), optional :: load_from_file  !< Sentinel to load grid from file.
   logical,            intent(in), optional :: is_parametric   !< Sentinel to load grid from parametric grid file.
   character(*),       intent(in), optional :: file_name       !< Grid file name.
   type(grid_object),  intent(in), optional :: grid            !< Grid data.
   logical                                  :: load_from_file_ !< Sentinel to load grid from file, local variable.

   load_from_file_ = .false. ;  if (present(load_from_file)) load_from_file_ = load_from_file
   call self%destroy
   call self%grid_dimensions%initialize
   if (load_from_file_) then
      call self%load_from_file(is_parametric=is_parametric, file_name=file_name)
   else if (present(grid)) then
      call self%grid_dimensions%initialize(block_signature=grid%blocks%signature)
      allocate(self%blocks(1:size(grid%blocks, dim=1)), source=grid%blocks)
   endif
   endsubroutine initialize

   subroutine load_from_file(self, is_parametric, file_name)
   !< Load from grid file.
   class(grid_object), intent(inout)        :: self           !< Grid.
   logical,            intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   character(*),       intent(in), optional :: file_name      !< File name.
   logical                                  :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.

   is_parametric_ = .false. ;  if (present(is_parametric)) is_parametric_ = is_parametric
   if (is_parametric_) then
      call self%file_parametric_grid%initialize(file_name=file_name)
   else
      call self%file_grid%initialize(file_name=file_name)
      call self%file_grid%load_grid_dimensions_from_file(grid_dimensions=self%grid_dimensions)
      call self%allocate_blocks
      call self%file_grid%load_nodes_from_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
   endif
   endsubroutine load_from_file

   subroutine save_into_file(self, is_parametric, file_name, ascii, metrics, off, tecplot, vtk)
   !< Save grid file.
   class(grid_object), intent(inout)        :: self           !< Grid.
   logical,            intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   character(*),       intent(in), optional :: file_name      !< File name.
   logical,            intent(in), optional :: ascii          !< Ascii/binary output.
   logical,            intent(in), optional :: metrics        !< Save also metrics data.
   logical,            intent(in), optional :: off            !< Save in OFF format sentinel.
   logical,            intent(in), optional :: tecplot        !< Tecplot output format sentinel.
   logical,            intent(in), optional :: vtk            !< VTK output format sentinel.
   logical                                  :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.
   logical                                  :: off_           !< OFF format sentinel, local variable.
   logical                                  :: vtk_           !< VTK format sentinel, local variable.
   type(string)                             :: file_name_     !< File name buffer.
   integer(I4P)                             :: b              !< Counter.

   is_parametric_ = .false. ;  if (present(is_parametric)) is_parametric_ = is_parametric
   if (is_parametric_) then
      error stop 'error: grid_object%save_into_file(is_parametric=.true., ...) to be implemented'
   else
      off_ = .true.  ; if (present(off)) off_ = off
      vtk_    = .false. ; if (present(vtk   )) vtk_    = vtk

      if (off_) then
         if (present(file_name)) call self%file_grid%initialize(file_name=file_name)
         call self%file_grid%save_grid_dimensions_into_file(grid_dimensions=self%grid_dimensions)
         call self%file_grid%save_nodes_into_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
      endif

      if (vtk_) then
         do b=1, self%grid_dimensions%blocks_number
            file_name_ = trim(adjustl(file_name))
            file_name_ = file_name_%basename(strip_last_extension=.true.)
            file_name_ = file_name_//'-block'//                                             &
                         '-id_'//trim(str(n=self%blocks(b)%signature%id, no_sign=.true.))// &
                         '-lv_'//trim(str(n=self%blocks(b)%signature%level, no_sign=.true.))//'.vts'
            call self%blocks(b)%save_file_grid(file_name=file_name_%chars(), ascii=ascii, metrics=metrics, vtk=vtk)
         enddo
      endif
      endif
   endsubroutine save_into_file
endmodule off_grid_object
