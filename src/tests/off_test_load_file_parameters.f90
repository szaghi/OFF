!< OFF test: load file of simulation parameters.

program off_test_load_file_parameters
!< OFF test: load file of simulation parameters.
!<
!< @note Execute test with `file-skeletons/simulation_parameters.ini` input:
!<```bash
!< off_test_load_file_parameters -par file-skeletons/simulation_parameters.ini
!<```

use off_objects, only : simulation_object
use penf, only : I4P

implicit none
logical                 :: go_on_fail          !< Go on if load fails.
type(simulation_object) :: simulation          !< Simulation data.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call simulation%load_file_parameters(file_name='file-skeletons/simulation_parameters.ini')

are_tests_passed(1) = simulation%solver%time_integrator == 'rk5'

print '(A)', simulation%description()

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
endprogram off_test_load_file_parameters
