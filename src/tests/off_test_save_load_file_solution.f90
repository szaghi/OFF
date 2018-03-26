!< OFF test: save (and re-load for checking) file solution.

program off_test_save_load_file_solution
!< OFF test: save (and re-load for checking) file solution.

use off_objects, only : simulation_object
use penf, only : R8P

implicit none
type(simulation_object) :: simulation          !< Simulation data.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call simulation%initialize(file_parameters='src/tests/file_solution_input/test.ini')

call simulation%save_file_grid(file_grid=simulation%file_grid_output, metrics=.true., off=.true., vtk=.true., force=.true.)
call simulation%save_file_solution(file_solution=simulation%file_solution, metrics=.true., off=.true., vtk=.true., force=.true.)

call simulation%initialize(file_parameters='src/tests/file_solution_input/reload.ini')

are_tests_passed(1) = simulation%mesh%blocks(3)%cell(1,6,8)%U%density == 0.1_R8P

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
endprogram off_test_save_load_file_solution
