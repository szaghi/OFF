!< OFF test: save (and re-load for checking) file grid.

program off_test_save_load_file_grid
!< OFF test: save (and re-load for checking) file grid.

use off_objects, only : simulation_object
use penf, only : I4P, I8P
use vecfor, only : ex, ey, ez

implicit none
type(simulation_object) :: simulation          !< Simulation data.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call simulation%initialize(file_parameters='src/tests/file_grid_input/load_grid_parametric.ini')

call simulation%save_file_grid(file_grid=simulation%file_grid_output, metrics=.true., off=.true., vtk=.true.)

call simulation%initialize(file_parameters='src/tests/file_grid_input/load_grid.ini')

are_tests_passed(1) = simulation%mesh%blocks(3)%cells_number(with_ghosts=.false.) == 33 * 4 * 10

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
endprogram off_test_save_load_file_grid
