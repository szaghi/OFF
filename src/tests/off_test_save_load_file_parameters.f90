!< OFF test: save (and re-load for checking) file of simulation parameters.

program off_test_save_load_file_parameters
!< OFF test: save (and re-load for checking) file of simulation parameters.

use off_objects, only : simulation_object
use flap, only : command_line_interface
use penf, only : I4P, R8P

implicit none
character(999)          :: file_parameters     !< Name of simulation parameters.
type(simulation_object) :: simulation          !< Simulation data.
logical                 :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

simulation%adimensionals%Ma = 3._R8P

call cli_parse

call simulation%save_file_parameters(file_name=file_parameters)

call simulation%initialize ! re-initialize

call simulation%load_file_parameters(file_name=file_parameters)

are_tests_passed(1) = simulation%adimensionals%Ma == 3._R8P

print '(A)', simulation%description()

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

call simulation%os%rm(file_name=trim(file_parameters)) ! remove temporary file parameters

contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli   !< Test command line interface.
  integer(I4P)                 :: error !< Error trapping flag.

  call cli%init(progname='off_test_save_load_file_parameters',                                   &
                authors='S. Zaghi',                                                              &
                help='Usage: ',                                                                  &
                examples=["off_test_save_load_file_parameters --parameters sim_parameters.ini"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--parameters-file',                &
               switch_ab='-par',                          &
               help='name of simulation parameters file', &
               required=.true.,                           &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--parameters-file', val=file_parameters)
  endsubroutine cli_parse
endprogram off_test_save_load_file_parameters
