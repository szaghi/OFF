!< OFF test: save (and re-load for checking) file grid.

program off_test_save_load_file_grid
!< OFF test: save (and re-load for checking) file grid.

use off_objects, only : block_object, simulation_object
use flap, only : command_line_interface
use penf, only : I4P, I8P
use vecfor, only : ex, ey, ez

implicit none
character(999)          :: file_basename       !< Basename of grid file.
type(simulation_object) :: simulation          !< Simulation data.
type(block_object)      :: blocks(2)           !< A block.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call cli_parse

call blocks(1)%initialize(id=1_I8P, level=1, gc=[2, 2, 2, 2, 2, 2], ni=8, nj=16, nk=32)
call blocks(1)%create_linspace(emin=(0*ex), emax=(ex+ey+ez))
call blocks(2)%initialize(id=2_I8P, level=1, gc=[2, 2, 2, 2, 2, 2], ni=16, nj=32, nk=64)
call blocks(2)%create_linspace(emin=(1*ex), emax=(2*ex+ey+ez))

call simulation%initialize(blocks=blocks)

call simulation%save_file_grid(file_basename=trim(adjustl(file_basename)), metrics=.true., off=.true., vtk=.true.)

! re-load for checking
call simulation%initialize

call simulation%load_file_grid(file_basename=trim(adjustl(file_basename)))

are_tests_passed(1) = simulation%blocks(2)%cells_number(with_ghosts=.false.) == 16_I4P * 32_I4P * 64_I4P

print '(A)', simulation%description()

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

! remove temporary files grid
call simulation%os%rm(file_name=trim(adjustl(file_basename))//'.grd')
call simulation%os%rm(file_name=trim(adjustl(file_basename))//'*.vts')

contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli   !< Test command line interface.
  integer(I4P)                 :: error !< Error trapping flag.

  call cli%init(progname='off_test_save_load_file_grid',                        &
                authors='S. Zaghi',                                             &
                help='Usage: ',                                                 &
                examples=["off_test_save_load_file_grid --grid-basename grid"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--grid-basename',     &
               switch_ab='-g',               &
               help='basename of grid file', &
               required=.false.,             &
               def='grid',                   &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--grid-basename', val=file_basename)
  endsubroutine cli_parse
endprogram off_test_save_load_file_grid
