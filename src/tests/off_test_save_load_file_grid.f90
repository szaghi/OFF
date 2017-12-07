!< OFF test: save (and re-load for checking) file grid.

program off_test_save_load_file_grid
!< OFF test: save (and re-load for checking) file grid.

use off_objects, only : mesh_object, simulation_object
use flap, only : command_line_interface
use penf, only : I4P, I8P
use vecfor, only : ex, ey, ez

implicit none
character(999)          :: file_name           !< File name of grid file.
type(simulation_object) :: simulation          !< Simulation data.
type(mesh_object)       :: mesh                !< A mesh.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call cli_parse

allocate(mesh%blocks(2))
call mesh%blocks(1)%initialize(id=1_I8P, level=1, gc=[2, 2, 2, 2, 2, 2], ni=8, nj=16, nk=32)
call mesh%blocks(1)%create_linspace(emin=(0*ex), emax=(ex+ey+ez))
call mesh%blocks(2)%initialize(id=2_I8P, level=1, gc=[2, 2, 2, 2, 2, 2], ni=16, nj=32, nk=64)
call mesh%blocks(2)%create_linspace(emin=(1*ex), emax=(2*ex+ey+ez))

call simulation%initialize(mesh=mesh)

call simulation%mesh%save_grid_into_file(file_name=trim(adjustl(file_name)), metrics=.true., off=.true., vtk=.true.)

! re-load for checking
call simulation%initialize

call simulation%mesh%load_grid_from_file(file_name=trim(adjustl(file_name)))

are_tests_passed(1) = simulation%mesh%blocks(2)%cells_number(with_ghosts=.false.) == 16_I4P * 32_I4P * 64_I4P

print '(A)', simulation%description()

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

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

  call cli%add(switch='--grid',          &
               switch_ab='-g',           &
               help='name of grid file', &
               required=.false.,         &
               def='grid.grd',           &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--grid', val=file_name)
  endsubroutine cli_parse
endprogram off_test_save_load_file_grid
