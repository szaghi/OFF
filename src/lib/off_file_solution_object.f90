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
   integer(I8P) :: save_frequency=1         !< Solution save frequency (on time steps).
   logical      :: save_metrics=.false.     !< Save metrics sentinel.
   logical      :: save_ghost_cells=.false. !< Save ghost cells sentinel.
   logical      :: ascii_format=.false.     !< ASCII file format sentinel.
   logical      :: off_format=.false.       !< OFF file format sentinel.
   logical      :: tecplot_format=.false.   !< Tecplot file format sentinel.
   logical      :: vtk_format=.false.       !< VTK file format sentinel.
   contains
      ! public methods
      procedure, pass(self) :: description                  !< Return a pretty-formatted description of the file.
      procedure, pass(self) :: destroy                      !< Destroy file.
      procedure, pass(self) :: load_parameters_from_file    !< Load file parameters from file.
      procedure, pass(self) :: load_conservatives_from_file !< Load conservative variables from file.
      procedure, pass(self) :: save_conservatives_into_file !< Save conservative variables into file.
      ! operators
      procedure, pass(lhs) :: file_assign_file !< Operator `=`.
endtype file_solution_object

contains
   ! public methods
   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted description of the file.
   class(file_solution_object), intent(in)           :: self             !< File object.
   character(*),                intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                     :: desc             !< Description.
   character(len=:), allocatable                     :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter                       :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = self%file_object%description(prefix=prefix_)//NL
   desc = desc//prefix_//'save frequency: '//trim(str(no_sign=.true., n=self%save_frequency))//NL
   desc = desc//prefix_//'save metrics: '//trim(str(self%save_metrics))//NL
   desc = desc//prefix_//'save ghost cells: '//trim(str(self%save_ghost_cells))//NL
   desc = desc//prefix_//'ascii format: '//trim(str(self%off_format))//NL
   desc = desc//prefix_//'off format: '//trim(str(self%off_format))//NL
   desc = desc//prefix_//'tecplot format: '//trim(str(self%tecplot_format))//NL
   desc = desc//prefix_//'vtk format: '//trim(str(self%vtk_format))//NL
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy file.
   class(file_solution_object), intent(inout) :: self  !< File object.
   type(file_object)                          :: fresh !< Fresh instance of file object.

   self = fresh
   endsubroutine destroy

   subroutine load_parameters_from_file(self, fini, options_prefix, go_on_fail)
   !< Load file parameters from file.
   class(file_solution_object),  intent(inout)        :: self           !< File object.
   type(file_ini),               intent(in)           :: fini           !< Solution parameters ini file handler.
   character(len=*),             intent(in)           :: options_prefix !< Prefix string of file options names.
   logical,                      intent(in), optional :: go_on_fail     !< Go on if load fails.
   logical                                            :: go_on_fail_    !< Go on if load fails, local variable.
   integer(I8P)                                       :: buffer_i       !< Buffer integer.
   logical                                            :: buffer_l       !< Buffer logical.
   character(999)                                     :: buffer_s       !< Buffer string.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call self%file_object%load_parameters_from_file(fini=fini, options_prefix=options_prefix, go_on_fail=go_on_fail_)

   call fini%get(section_name='files', option_name=options_prefix//'_save_frequency', val=buffer_i, error=self%error%status)
   call self%error%check(message='failed to load [files].('//options_prefix//'_save_frequency)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%save_frequency = buffer_i

   call fini%get(section_name='files', option_name=options_prefix//'_save_metrics', val=buffer_l, error=self%error%status)
   call self%error%check(message='failed to load [files].('//options_prefix//'_save_metrics)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%save_metrics = buffer_l

   call fini%get(section_name='files', option_name=options_prefix//'_save_ghost_cells', val=buffer_l, error=self%error%status)
   call self%error%check(message='failed to load [files].('//options_prefix//'_save_ghost_cells)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%save_ghost_cells = buffer_l

   call fini%get(section_name='files', option_name=options_prefix//'_ascii_format', val=buffer_l, error=self%error%status)
   call self%error%check(message='failed to load [files].('//options_prefix//'_ascii_format)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%ascii_format = buffer_l

   call fini%get(section_name='files', option_name=options_prefix//'_off_format', val=buffer_l, error=self%error%status)
   call self%error%check(message='failed to load [files].('//options_prefix//'_off_format)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%off_format = buffer_l

   call fini%get(section_name='files', option_name=options_prefix//'_tecplot_format', val=buffer_l, error=self%error%status)
   call self%error%check(message='failed to load [files].('//options_prefix//'_tecplot_format)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%tecplot_format = buffer_l

   call fini%get(section_name='files', option_name=options_prefix//'_vtk_format', val=buffer_l, error=self%error%status)
   call self%error%check(message='failed to load [files].('//options_prefix//'_vtk_format)', is_severe=.not.go_on_fail_)
   if (self%error%status <= 0) self%vtk_format = buffer_l
   endsubroutine load_parameters_from_file

   subroutine load_conservatives_from_file(self, grid_dimensions, blocks, file_name)
   !< Load conservative variables from file.
   class(file_solution_object),  intent(inout)        :: self            !< File object.
   type(grid_dimensions_object), intent(in)           :: grid_dimensions !< Grid dimensions off all blocks into file.
   type(block_object),           intent(inout)        :: blocks(1:)      !< Blocks storage.
   character(*),                 intent(in), optional :: file_name       !< File name.
   type(file_ini)                                     :: fini            !< Solution parameters ini file handler.
   integer(I4P)                                       :: blocks_number   !< Blocks number.
   type(primitive_compressible)                       :: P               !< Primitive variables.
   real(R8P)                                          :: velocity_(3)    !< Velocity temporary array.
   character(len=:), allocatable                      :: emsg_suffix     !< Error message.
   integer(I4P)                                       :: b               !< Counter.

   if (present(file_name)) self%file_name = trim(adjustl(file_name))
   emsg_suffix = ' from file "'//self%file_name//'" in procedure "file_solution_object%load_conservative_from_file"'
   if (self%is_parametric) then
      call fini%load(filename=self%file_name, error=self%error%status)
      call fini%get(section_name='dimensions', option_name='blocks_number', val=blocks_number, error=self%error%status)
      call self%error%check(message='failed to load [dimensions].(blocks_number)'//emsg_suffix, is_severe=.true.)
      if (blocks_number>0.and.size(blocks, dim=1)>=blocks_number) then
         do b=1, blocks_number
            call fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='pressure', &
                          val=P%pressure, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(pressure)'//emsg_suffix, &
                                  is_severe=.true.)
            call fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='density', &
                          val=P%density, error=self%error%status)
            call self%error%check(message='failed to load [block_'//trim(str(b,no_sign=.true.))//'].(density)'//emsg_suffix, &
                                  is_severe=.true.)
            call fini%get(section_name='block_'//trim(str(b, no_sign=.true.)), option_name='velocity', &
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

   ! operators
   pure subroutine file_assign_file(lhs, rhs)
   !< Operator `=`.
   class(file_solution_object), intent(inout) :: lhs !< Left hand side.
   class(file_object),          intent(in)    :: rhs !< Right hand side.

   call lhs%file_object%file_assign_file(rhs=rhs)
   select type(rhs)
   type is(file_solution_object)
      lhs%save_frequency = rhs%save_frequency
      lhs%save_metrics = rhs%save_metrics
      lhs%save_ghost_cells = rhs%save_ghost_cells
      lhs%ascii_format = rhs%ascii_format
      lhs%off_format = rhs%off_format
      lhs%tecplot_format = rhs%tecplot_format
      lhs%vtk_format = rhs%vtk_format
   endselect
   endsubroutine file_assign_file
endmodule off_file_solution_object
