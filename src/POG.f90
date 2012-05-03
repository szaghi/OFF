program POG
!-----------------------------------------------------------------------------------------------------------------------------------
!(doc)CNAME{POG, Post-processor Output Generator for OFF (Open Finite volume Fluid dynamics code)}
!(doc)AUTHORS{Stefano Zaghi}
!(doc)VERSION{v0.0.5}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                            ! Integers and reals precision definition.
USE Data_Type_Globals                       ! Definition of Type_Global and Type_Block.
USE Data_Type_OS, init_os=>init             ! Definition of Type_OS.
USE Data_Type_Time                          ! Definition of Type_Time.
USE Data_Type_Vector, set_vec=>set          ! Definition of Type_Vector.
USE Lib_IO_Misc                             ! Procedures for IO and strings operations.
USE Lib_Mesh,           only: node2center   ! Subroutine for computing cell center cell nodes ones.
USE Lib_PostProcessing, only: pp_format,  & ! Post-processing data format.
                              tec_output, & ! Function for writing Tecplot file.
                              vtk_output    ! Function for writing VTK file.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Global)::             global            ! Global-level data.
type(Type_Block), allocatable:: block(:)          ! Block-level data.
logical::                       meshonly = .true. ! Flag for post-process only mesh.
integer(I_P)::                  b,l               ! Counters.
integer(I_P)::                  err               ! Error trapping flag: 0 no errors, >0 error occurs.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! parsing command line for getting options
call parsing_command_line

! the number of grid levels/blocks are artificially set to 1: the post processor convert one level/block at time
global%mesh%Nl = 1 ; global%mesh%Nb = 1
! allocating block-level data
if (allocated(block)) then
  call free_block(block=block(1))
  deallocate(block)
endif
allocate(block(1))
! getting dimensions of block
err = load_bmesh_dims(myrank = 0, filename = global%file%File_Mesh, block = block(1))
! allocating block
call alloc_block(global=global,block=block(1))

write(stdout,'(A)',IOSTAT=err)'----------------------------------------------------------------------'
write(stdout,'(A)',IOSTAT=err)' Loading '//trim(global%file%File_Mesh)

err = load_bmesh(myrank = 0, filename = global%file%File_Mesh, block = block(1))
if (.not.meshonly) then
  write(stdout,'(A)',IOSTAT=err)' Loading '//trim(global%file%File_Sol)
  err = load_gfluid_Ns(binary = .true., myrank = 0, filename = global%file%File_Sol, global = global)
  ! allocating global fluid dynamic data
  call alloc_global_fluid(global=global)
  err = load_bfluid(myrank = 0, filename = global%file%File_Sol, global = global, block = block(1))
endif
write(stdout,*)

write(stdout,'(A)',IOSTAT=err)' Computing the mesh variables that are not loaded from input files'
! computing the cell-center coordinates
write(stdout,'(A)',IOSTAT=err)'   Computing cells centers coordinates'
call node2center(block = block(1))

write(stdout,'(A)',IOSTAT=err)'----------------------------------------------------------------------'
write(stdout,'(A)',IOSTAT=err)' Processing data'
write(stdout,'(A)',IOSTAT=err)'  Writing Tecplot and/or VTK output data'
! writing output
if (pp_format%tec) err = tec_output(meshonly=meshonly, global=global, block=block(1), filename=trim(global%file%File_Pout))
if (pp_format%vtk) err = vtk_output(meshonly=meshonly, global=global, block=block(1), filename=trim(global%file%File_Pout))
write(stdout,'(A)',IOSTAT=err)' Processing data Complete'
write(stdout,'(A)',IOSTAT=err)'----------------------------------------------------------------------'
write(stdout,*)
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine print_usage()
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for printing the correct use of the program.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stdout,'(A)')' POG'
  write(stdout,'(A)')' Post processing code for OFF code'
  write(stdout,'(A)')' Usage:'
  write(stdout,'(A)')'   POG -m mesh_file_name'
  write(stdout,'(A)') '     [-o output_file_name'
  write(stdout,'(A)') '      -s solution_file_name'
  write(stdout,'(A)') '      -bc'
  write(stdout,'(A)') '      -cell'
  write(stdout,'(A)') '      -ascii'
  write(stdout,'(A)') '      -tec yes/no'
  write(stdout,'(A)') '      -vtk yes/no'
  write(stdout,'(A)') '      -os UIX/WIN]'
  write(stdout,*)
  write(stdout,'(A)') ' Optional arguments and default values:'
  write(stdout,'(A)') '  -o output_file_name   => output file name; default is basename of mesh file with the proper extension'
  write(stdout,'(A)') '  -s solution_file_name => solution file name; if passed the solution variables are saved'
  write(stdout,'(A)') '  -bc                   => saving boundary conditions ghost cells (default no)'
  write(stdout,'(A)') '  -cell                 => all vars other than nodes coord. are cell centered (default no, node centered)'
  write(stdout,'(A)') '  -ascii                => write ascii output file (default no, write binary one)'
  write(stdout,'(A)') '  -tec yes/no           => write (or not) Tecplot file format (default yes)'
  write(stdout,'(A)') '  -vtk yes/no           => write (or not) VTK file format (default no)'
  write(stdout,'(A)') '  -os UIX/WIN           => type of Operating System write (default *UIX OS type)'
  write(stdout,*)
  write(stdout,'(A)')' Examples: '
  write(stdout,'(A)')'  Solution file'
  write(stdout,'(A)')'    Teplot) POG -m example.geo -s example.sol -o example'
  write(stdout,'(A)')'    VTK)    POG -m example.geo -s example.sol -o example'
  write(stdout,'(A)')'  Only mesh file'
  write(stdout,'(A)')'    Teplot) POG -m example.geo                -o example'
  write(stdout,'(A)')'    VTK)    POG -m example.geo                -o example'
  write(stdout,*)
  write(stdout,'(A)') ' Notes:'
  write(stdout,'(A)') '   1) the output file name extension is not necessary because it assigned according to the type of output:'
  write(stdout,'(A)') '      binary       Tecplot => .plt'
  write(stdout,'(A)') '      ascii        Tecplot => .dat'
  write(stdout,'(A)') '      binary/ascii VTK     => .vtm'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_usage

  subroutine parsing_command_line()
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for printing the correct use of the program.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P):: Nca = 0         ! Number of command line arguments.
  integer(I_P):: c               ! Counter for command line arguments.
  character(6):: ca_switch       ! Switch identifier.
  character(3):: yes             ! Yes (no) flag.
  character(3):: os_type = 'UIX' ! OS type 'UIX' or 'WIN'.
  character(20):: date           ! Actual date.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing solution file name to "undefined" in order to capture when only mesh is to be post-processed
  global%file%File_Sol = "undefined"
  Nca = command_argument_count()
  if (Nca==0) then
    write(stderr,'(A)')' Error: no argument has been passed to command line'
    call print_usage
    stop
  else
    ! processing switch
    c = 0
    do while (c<Nca)
      c = c + 1
      call get_command_argument(c, ca_switch)
      select case(adjustl(trim(ca_switch)))
      case('-m') ! file name of mesh file
        call get_command_argument(c+1, global%file%File_Mesh) ; c = c + 1
      case('-o') ! file name of output file
        call get_command_argument(c+1, global%file%File_Pout) ; c = c + 1
      case('-s') ! file name of solution file
        call get_command_argument(c+1, global%file%File_Sol) ; c = c + 1
      case('-bc') ! saving boundary conditions ghost cells
        pp_format%bc = .true.
      case('-cell') ! dependent variables are saved cell centered instead of interpolated at nodes
        pp_format%node = .false.
      case('-ascii') ! ascii output
        pp_format%binary = .false.
      case('-tec') ! Tecplot file format
        call get_command_argument(c+1,yes) ; c = c + 1 ; yes = Upper_Case(adjustl(trim(yes)))
        pp_format%tec = (adjustl(trim(yes))=='YES')
      case('-vtk') ! VKT file format
        call get_command_argument(c+1,yes) ; c = c + 1 ; yes = Upper_Case(adjustl(trim(yes)))
        pp_format%vtk = (adjustl(trim(yes))=='YES')
      case('-os') ! OS type
        call get_command_argument(c+1,os_type) ; c = c + 1 ; os_type = Upper_Case(os_type)
      case default
        write(stderr,'(A)') ' Error: switch "'//adjustl(trim(ca_switch))//'" is unknown'
        call print_usage
        stop
      endselect
    enddo
  endif
  if (pp_format%node.and.pp_format%bc) then
    write(stderr,'(A)') ' It is not possible to save bc ghost cells and node-centered interpolated variables!'
    write(stderr,*)
    stop
  endif
  ! set OS type
  OS=init_os(c_id=os_type)
  if (adjustl(trim(global%file%File_Mesh))=='unset') then
    write(stderr,'(A)') ' File name of mesh file must be provided!'
    write(stderr,*)
    call print_usage
    stop
  endif
  if (adjustl(trim(global%file%File_Pout))=='unset') then
    ! the name of mesh file is used as output file base name
    global%file%File_Pout = adjustl(trim(global%file%File_Mesh))
  endif
  global%file%File_Mesh = adjustl(trim(global%file%File_Mesh)) ; global%file%File_Mesh = string_OS_sep(global%file%File_Mesh)
  global%file%File_Sol  = adjustl(trim(global%file%File_Sol )) ; global%file%File_Sol  = string_OS_sep(global%file%File_Sol )
  global%file%File_Pout = adjustl(trim(global%file%File_Pout)) ; global%file%File_Pout = string_OS_sep(global%file%File_Pout)
  date = Get_Date_String()
  write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)' POG started on'
  write(stdout,'(A)',iostat=err)' '//date
  write(stdout,'(A)',iostat=err)' Input Files'
  write(stdout,'(A)',iostat=err)'  Mesh file:                '//trim(global%file%File_Mesh)
  write(stdout,'(A)',iostat=err)'  Boundary conditions file: '//trim(global%file%File_BC)
  write(stdout,'(A)',iostat=err)'  Solution file:            '//trim(global%file%File_Sol)
  write(stdout,'(A)',iostat=err)'  Output file:              '//trim(global%file%File_Pout)
  write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
  write(stdout,*)
  if (trim(global%file%File_Sol)/="unset") then
    ! a solution file name has been passed
    meshonly = .false.
  else
    ! no solution file name has been passed: post-processing of only mesh is assumed
    meshonly = .true.
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parsing_command_line
endprogram POG
