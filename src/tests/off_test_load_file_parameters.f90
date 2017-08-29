!< OFF test: load file of simulation parameters.

program off_test_load_file_parameters
!< OFF test: load file of simulation parameters.

use off_objects, only : simulation_object
use flap, only : command_line_interface
use penf, only : I4P

implicit none
character(999)          :: file_parameters     !< Name of simulation parameters.
logical                 :: go_on_fail          !< Go on if load fails.
type(simulation_object) :: simulation          !< Simulation data.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call cli_parse

call simulation%load_file_parameters(file_name=file_parameters, go_on_fail=go_on_fail)

are_tests_passed(1) = simulation%solver%time_integrator == 'rk5'

print '(A)', simulation%description()

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli   !< Test command line interface.
  integer(I4P)                 :: error !< Error trapping flag.

  call cli%init(progname='off_test_load_file_parameters',                                   &
                authors='S. Zaghi',                                                         &
                help='Usage: ',                                                             &
                examples=["off_test_load_file_parameters --parameters sim_parameters.ini"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--parameters-file',                &
               switch_ab='-par',                          &
               help='name of simulation parameters file', &
               required=.true.,                           &
               act='store')

  call cli%add(switch='--go-on-fail',                &
               switch_ab='-gof',                     &
               help='go-on if load fails somewhere', &
               required=.false.,                     &
               def='.false.',                        &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--parameters-file', val=file_parameters)
  call cli%get(switch='--go-on-fail', val=go_on_fail)
  endsubroutine cli_parse
endprogram off_test_load_file_parameters
