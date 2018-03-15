!< OFF test: load file parametric grid.

program off_test_load_file_parametric_grid
!< OFF test: load file parametric grid.

use off_objects, only : file_grid_object, mesh_object, simulation_object
use flap, only : command_line_interface
use penf, only : I4P, I8P
use vecfor, only : ex, ey, ez

implicit none
character(999)          :: file_name           !< File name of grid file.
type(file_grid_object)  :: file_grid           !< A grid file.
type(simulation_object) :: simulation          !< Simulation data.
type(mesh_object)       :: mesh                !< A mesh.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call cli_parse

call file_grid%initialize(file_name=trim(adjustl(file_name)), is_parametric=.true.)
call simulation%mesh%load_grid_from_file(file_grid=file_grid)

are_tests_passed(1) = simulation%mesh%blocks(2)%cells_number(with_ghosts=.false.) == 16_I4P * 32_I4P * 64_I4P

print '(A)', simulation%description()

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

contains
   subroutine cli_parse()
   !< Build and parse test cli.
   type(command_line_interface) :: cli   !< Test command line interface.
   integer(I4P)                 :: error !< Error trapping flag.

   call cli%init(progname='off_test_load_file_parametric_grid',                              &
                 authors='S. Zaghi',                                                         &
                 help='Usage: ',                                                             &
                 examples=["off_test_load_file_parametric_grid --grid parametric_grid.ini"], &
                 epilog=new_line('a')//"all done")

   call cli%add(switch='--grid',                          &
                switch_ab='-g',                           &
                help='name of parametric grid file',      &
                required=.false.,                         &
                def='file-skeletons/parametric_grid.ini', &
                act='store')

   call cli%parse(error=error) ; if (error/=0) stop

   call cli%get(switch='--grid', val=file_name)
   endsubroutine cli_parse
endprogram off_test_load_file_parametric_grid
