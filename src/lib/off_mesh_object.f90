!< OFF mesh object definition and implementation.

module off_mesh_object
!< OFF mesh object definition and implementation.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use off_block_object, only : block_object
use off_file_grid_object, only : file_grid_object
use off_file_parametric_grid_object, only : file_parametric_grid_object
use off_file_parametric_ic_object, only : file_parametric_ic_object
use off_grid_dimensions_object, only : grid_dimensions_object
use flow, only : eos_compressible
use penf, only : I4P, str
use stringifor, only : string

implicit none
private
public :: mesh_object

type :: mesh_object
   !< mesh object class.
   !<
   !< [[mesh_object]] is a container for all mesh data.
   type(file_grid_object)            :: file_grid            !< Grid file handler.
   type(file_parametric_grid_object) :: file_parametric_grid !< Parametric grid file handler.
   type(file_parametric_ic_object)   :: file_parametric_ic   !< Parametric initial conditions file handler.
   type(grid_dimensions_object)      :: grid_dimensions      !< Grid dimensions.
   type(block_object), allocatable   :: blocks(:)            !< Blocks list.
   contains
      ! public methods
      procedure, pass(self) :: allocate_blocks     !< Allocate blocks accordingly to grid dimensions.
      procedure, pass(self) :: description         !< Return a pretty-formatted description of the mesh.
      procedure, pass(self) :: destroy             !< Destroy mesh.
      procedure, pass(self) :: initialize          !< Initialize mesh.
      procedure, pass(self) :: load_grid_from_file !< Load grid from file.
      procedure, pass(self) :: load_ic_from_file   !< Load initial conditions from file.
      procedure, pass(self) :: save_grid_into_file !< Save grid into file.
endtype mesh_object

contains
   ! public methods
   subroutine allocate_blocks(self)
   !< Allocate blocks accordingly to grid dimensions.
   class(mesh_object), intent(inout) :: self !< Mesh.
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
   !< Return a pretty-formatted description of the mesh.
   class(mesh_object), intent(in)           :: self    !< Mesh.
   character(*),       intent(in), optional :: prefix  !< Prefixing string.
   character(len=:), allocatable            :: desc    !< Description.
   character(len=:), allocatable            :: prefix_ !< Prefixing string, local variable.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//self%grid_dimensions%description(prefix=prefix_//'  ')
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy mesh.
   class(mesh_object), intent(inout) :: self !< Mesh.

   call self%file_grid%destroy
   call self%file_parametric_grid%destroy
   call self%file_parametric_ic%destroy
   call self%grid_dimensions%destroy
   if (allocated(self%blocks)) then
      call self%blocks%destroy
      deallocate(self%blocks)
   endif
   endsubroutine destroy

   subroutine initialize(self, grid_file_name, is_grid_file_parametric, ic_file_name, is_ic_file_parametric, mesh)
   !< Initialize mesh.
   !<
   !< The initialization can be done also loading from files (standard or paramatric ones).
   class(mesh_object), intent(inout)        :: self                    !< Mesh.
   character(*),       intent(in), optional :: grid_file_name          !< Grid file name.
   logical,            intent(in), optional :: is_grid_file_parametric !< Sentinel for parametric grid file.
   character(*),       intent(in), optional :: ic_file_name            !< Initial conditions file name.
   logical,            intent(in), optional :: is_ic_file_parametric   !< Sentinel for initial condition file.
   type(mesh_object),  intent(in), optional :: mesh                    !< mesh data.

   call self%destroy
   call self%grid_dimensions%initialize
   if (present(mesh)) then
      call self%grid_dimensions%initialize(block_signature=mesh%blocks%signature)
      allocate(self%blocks(1:size(mesh%blocks, dim=1)), source=mesh%blocks)
   else
      if (present(grid_file_name)) call self%load_grid_from_file(is_parametric=is_grid_file_parametric, file_name=grid_file_name)
      ! if (present(ic_file_name)) call self%load_ic_from_file(is_parametric=is_ic_file_parametric, file_name=ic_file_name)
   endif
   endsubroutine initialize

   subroutine load_grid_from_file(self, is_parametric, file_name)
   !< Load grid from file.
   class(mesh_object), intent(inout)        :: self           !< Mesh.
   logical,            intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   character(*),       intent(in), optional :: file_name      !< File name.
   logical                                  :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.
   integer(I4P)                             :: b              !< Counter.

   is_parametric_ = .false. ;  if (present(is_parametric)) is_parametric_ = is_parametric
   if (is_parametric_) then
      call self%file_parametric_grid%get_grid_dimensions(file_name=file_name, grid_dimensions=self%grid_dimensions)
      if (self%grid_dimensions%blocks_number>0) then
         call self%allocate_blocks
         do b=1, self%grid_dimensions%blocks_number
            call self%blocks(b)%create_linspace
         enddo
      endif
   else
      call self%file_grid%initialize(file_name=file_name)
      call self%file_grid%load_grid_dimensions_from_file(grid_dimensions=self%grid_dimensions)
      call self%allocate_blocks
      call self%file_grid%load_nodes_from_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
   endif
   endsubroutine load_grid_from_file

   subroutine load_ic_from_file(self, is_parametric, eos, file_name)
   !< Load initial conditions from file.
   class(mesh_object),     intent(inout)        :: self           !< Mesh.
   logical,                intent(in), optional :: is_parametric  !< Sentinel to load grid parametric grid file.
   type(eos_compressible), intent(in), optional :: eos            !< Equation of state.
   character(*),           intent(in), optional :: file_name      !< File name.
   logical                                      :: is_parametric_ !< Sentinel to load grid parametric grid file, local variable.

   is_parametric_ = .false. ;  if (present(is_parametric)) is_parametric_ = is_parametric
   if (is_parametric_.and.present(eos)) then
      call self%file_parametric_ic%load_file(file_name=file_name)
      if (self%grid_dimensions%blocks_number>0) then
         call self%file_parametric_ic%get_initial_conditions(eos=eos, blocks=self%blocks, file_name=file_name)
      else
         error stop 'error: mesh blocks have not been allocated!'
      endif
   else
      error stop 'error: loading initial conditions from restart files to be implemented!'
   endif
   endsubroutine load_ic_from_file

   subroutine save_grid_into_file(self, is_parametric, file_name, ascii, metrics, off, tecplot, vtk)
   !< Save grid into file.
   class(mesh_object), intent(inout)        :: self           !< Mesh.
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
      error stop 'error: mesh_object%save_grid_into_file(is_parametric=.true., ...) to be implemented'
   else
      off_ = .true.  ; if (present(off)) off_ = off
      vtk_ = .false. ; if (present(vtk)) vtk_ = vtk

      if (off_) then
         if (present(file_name)) call self%file_grid%initialize(file_name=file_name)
         call self%file_grid%save_grid_dimensions_into_file(grid_dimensions=self%grid_dimensions)
         call self%file_grid%save_nodes_into_file(grid_dimensions=self%grid_dimensions, blocks=self%blocks)
      endif

      if (vtk_) then
         do b=1, self%grid_dimensions%blocks_number
            file_name_ = trim(adjustl(file_name))
            if (file_name_%basename(strip_last_extension=.true.)/='') file_name_ = file_name_%basename(strip_last_extension=.true.)
            file_name_ = file_name_//'-block'//                                             &
                         '-id_'//trim(str(n=self%blocks(b)%signature%id, no_sign=.true.))// &
                         '-lv_'//trim(str(n=self%blocks(b)%signature%level, no_sign=.true.))//'.vts'
            call self%blocks(b)%save_file_grid(file_name=file_name_%chars(), ascii=ascii, metrics=metrics, vtk=vtk)
         enddo
      endif
      endif
   endsubroutine save_grid_into_file
endmodule off_mesh_object
