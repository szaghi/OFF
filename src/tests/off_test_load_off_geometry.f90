!< OFF test: load file containing a geometry to be immersed into the grid.

program off_test_load_off_geometry
!< OFF test: load file containing a geometry to be immersed into the grid.
!<
!< @note Execute test with `file-skeletons/simulation_parameters.ini` input:
!<```bash
!< off_test_load_file_parameters -par file-skeletons/simulation_parameters.ini
!<```

use off_objects, only : mesh_object, simulation_object
use flap, only : command_line_interface
use penf, only : I4P, I8P, R8P
use vecfor, only : ex, ey, ez, vector

implicit none
character(999), allocatable :: file_name(:)        !< File name of OFF geometry file.
character(999)              :: outf_name           !< Output file name.
type(simulation_object)     :: simulation          !< Simulation data.
type(mesh_object)           :: mesh                !< A mesh.
integer(I4P)                :: f                   !< File counter.
logical                     :: are_tests_passed(1) !< Result of tests check.

are_tests_passed = .false.

call cli_parse
allocate(mesh%blocks(1))
call mesh%blocks(1)%initialize(id=1_I8P, level=1, gc=[0, 0, 0, 0, 0, 0], ni=128, nj=128, nk=64, &
                               interfaces_number=size(file_name, dim=1))
call mesh%blocks(1)%create_linspace(emin=(-180*ex-82*ey-53*ez), emax=(0*ex+100*ey+87*ez))
do f=1, size(file_name, dim=1)
   ! call mesh%blocks(1)%immerge_off_geometry(file_name=trim(adjustl(file_name(f))), n=f, &
                                            ! outside_reference=(1000*ex+1000*ey+100*ez), distance_sign_inverse=.false.)
enddo
call mesh%blocks(1)%update_level_set_distance
! call simulation%initialize(mesh=mesh)
call simulation%mesh%save_file_grid(file_name=trim(adjustl(outf_name)), metrics=.true., off=.false., vtk=.true.)

are_tests_passed(1) = simulation%mesh%blocks(1)%cell(1,1,1)%level_set%distance > 0._R8P

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

contains
  subroutine cli_parse()
  !< Build and parse test cli.
  type(command_line_interface) :: cli   !< Test command line interface.
  integer(I4P)                 :: error !< Error trapping flag.

  call cli%init(progname='off_test_load_off_geometry',                  &
                authors='S. Zaghi',                                     &
                help='Usage: ',                                         &
                examples=["off_test_load_off_geometry --off file.off"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--off',                         &
               help='name(s) of OFF geometry file(s)', &
               nargs='+',                              &
               required=.true.,                        &
               act='store')

  call cli%add(switch='--out',                   &
               help='output file name',          &
               required=.false.,                 &
               def='off_test_load_off_geometry', &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get_varying(switch='--off', val=file_name, error=error) ; if (error/=0) stop
  call cli%get(switch='--out', val=outf_name, error=error) ; if (error/=0) stop
  endsubroutine cli_parse
endprogram off_test_load_off_geometry
