!< OFF test: shock tube tester.

program off_test_shock_tube
!< OFF test: shock tube tester.

use off_objects, only : simulation_object
use penf, only : I4P

implicit none
type(simulation_object) :: simulation          !< Simulation data.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call simulation%load_file_parameters(file_name='file-skeletons/simulation_parameters.ini')

are_tests_passed(1) = simulation%solver%time_integrator == 'rk5'

print '(A)', simulation%description()

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
endprogram off_test_shock_tube
