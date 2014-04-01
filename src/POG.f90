!> @ingroup Program
!> @{
!> @defgroup POGProgram POG
!> @}

!> @brief POG, Post-processor Output Generator for @off (Open Finite volume Fluid dynamics code).
!> This is an auxiliary tool useful for post-processing @off simulations outputs. It can manipulate @off outputs and it can produce
!> files ready to be visualized. Two different visualization standards are supported:
!> - Tecplot, Inc.: Tecplot is a wide-used visualization tool (http://www.tecplot.com/). POG can produce both ascii and binary files
!>   in Tecplot standard.
!> - VTK: The Visualization Toolkit (VTK) is an open-source, freely available software system for 3D computer graphics, image
!>   processing and visualization (http://www.vtk.org/). A lot of visualization tools support VTK standard. Among those tools
!>   Paraview (http://www.paraview.org/) seems to be one of the most complete. POG can produce both ascii and binary files in VTK
!>   standard. To this aim the Lib_VTK_IO is used.
!> @note In order to produce binary Tecplot file the Tecplot, Inc.'s library tecio.a (or tecio64.a) must be available during the
!>       compilation of POG. If you are using the makefile shipped with @off code there is a dedicated option: TECIO=yes/no. If
!>       during the compilation this option is set to yes (e.g. make POG TECIO=yes) the make search the correct library into ./lib/.
!>       Edit the makefile to point to the correct path where tecio.a (or tecio64.a) is placed. For more details see \ref Compiling
!>       "Compiling Instructions".
!> @todo \b DocImprove: Improve the documentation
!> @ingroup POGProgram
program POG
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                  ! Integers and reals precision definition.
USE Data_Type_Files,       only: Type_Files       ! Definition of Type_Files.
USE Data_Type_Global,      only: Type_Global      ! Definition of Type_Global.
USE Data_Type_OS,          only: Type_OS          ! Definition of Type_OS.
USE Data_Type_PostProcess, only: Type_PostProcess ! Definition of Type_PostProcess.
USE Data_Type_SBlock,      only: Type_SBlock      ! Definition of Type_SBlock.
USE Lib_IO_Misc                                   ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Global)::      global !< Global-level data.
type(Type_Files)::       IOFile !< Input/Output files.
type(Type_PostProcess):: pp     !< Post-process data.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! initializing the post-processing
call pog_init
! saving output files
if (pp%tec) then
  write(stdout,'(A)')'+--> Saving '//IOFile%tec%name
  call IOFile%tec%save(global=global)
  if (IOFile%tec%iostat/=0) then
    write(stderr,'(A)')'+--> '//IOFile%tec%iomsg
    stop
  endif
endif
if (pp%vtk) then
  write(stdout,'(A)')'+--> Saving '//IOFile%vtk%name
  call IOFile%vtk%save(global=global)
  if (IOFile%vtk%iostat/=0) then
    write(stderr,'(A)')'+--> '//IOFile%vtk%iomsg
    stop
  endif
endif
if (pp%gnu) then
  write(stdout,'(A)')'+--> Saving '//IOFile%gnu%name
  call IOFile%gnu%save(global=global)
  if (IOFile%gnu%iostat/=0) then
    write(stderr,'(A)')'+--> '//IOFile%gnu%iomsg
    stop
  endif
endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! Subroutine for printing the correct use of the program.
  subroutine print_usage()
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
  write(stdout,'(A)') '      -schl'
  write(stdout,'(A)') '      -mirrorX'
  write(stdout,'(A)') '      -mirrorY'
  write(stdout,'(A)') '      -mirrorZ'
  write(stdout,'(A)') '      -tec yes/no'
  write(stdout,'(A)') '      -vtk yes/no'
  write(stdout,'(A)') '      -gnu yes/no]'
  write(stdout,*)
  write(stdout,'(A)') ' Optional arguments and default values:'
  write(stdout,'(A)') '  -o output_file_name   => output file name; default is basename of mesh file with the proper extension'
  write(stdout,'(A)') '  -s solution_file_name => solution file name; if passed the solution variables are saved'
  write(stdout,'(A)') '  -bc                   => saving boundary conditions ghost cells (default no)'
  write(stdout,'(A)') '  -cell                 => all vars other than nodes coord. are cell centered (default no, node centered)'
  write(stdout,'(A)') '  -ascii                => write ascii output file (default no, write binary one)'
  write(stdout,'(A)') '  -schl                 => saving (pseudo) Schlieren field (default no)'
  write(stdout,'(A)') '  -mirrorX              => saving solution togheter a X-axis mirrored copy of flow fieds (default no)'
  write(stdout,'(A)') '  -mirrorY              => saving solution togheter a Y-axis mirrored copy of flow fieds (default no)'
  write(stdout,'(A)') '  -mirrorZ              => saving solution togheter a Z-axis mirrored copy of flow fieds (default no)'
  write(stdout,'(A)') '  -tec yes/no           => write (or not) Tecplot file format (default yes)'
  write(stdout,'(A)') '  -vtk yes/no           => write (or not) VTK file format (default no)'
  write(stdout,'(A)') '  -gnu yes/no           => write (or not) Gnuplot file format (default no)'
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
  write(stdout,'(A)') '      ascii        Gnuplot => .gnu'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_usage

  ! Subroutine for printing the correct use of the program.
  subroutine parsing_command_line()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P)::  Nca = 0   !< Number of command line arguments.
  integer(I4P)::  c         !< Counter for command line arguments.
  character(8)::  ca_switch !< Switch identifier.
  character(3)::  yes       !< Yes (no) flag.
  character(99):: filename  !< File names dummy string.
  character(99):: outfname  !< Output file name dummy string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing solution file name to "unset" in order to capture when only mesh is to be post-processed
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
      select case(trim(adjustl(ca_switch)))
      case('-m') ! file name of mesh file
        call get_command_argument(c+1, filename) ; c = c + 1
        call IOFile%mesh%set(name=global%OS%string_separator_fix(string=trim(adjustl(filename))))
      case('-o') ! file name of output file
        call get_command_argument(c+1, outfname) ; c = c + 1
        outfname=global%OS%string_separator_fix(string=trim(adjustl(outfname)))
      case('-s') ! file name of solution file
        call get_command_argument(c+1, filename) ; c = c + 1
        call IOFile%sol%set(name=global%OS%string_separator_fix(string=trim(adjustl(filename))))
      case('-proc') ! processes/blocks map
        call get_command_argument(c+1, filename) ; c = c + 1
        call IOFile%proc%set(name=global%OS%string_separator_fix(string=trim(adjustl(filename))))
      case('-bc') ! saving boundary conditions ghost cells
        pp%bc = .true.
      case('-cell') ! dependent variables are saved cell centered instead of interpolated at nodes
        pp%node = .false.
      case('-ascii') ! ascii output
        pp%binary = .false.
      case('-schl') ! Schlieren field
        pp%schlieren = .true.
      case('-mirrorX') ! mirroring solution using X versor
        pp%mirrorX = .true.
      case('-mirrorY') ! mirroring solution using Y versor
        pp%mirrorY = .true.
      case('-mirrorZ') ! mirroring solution using Z versor
        pp%mirrorZ = .true.
      case('-tec') ! Tecplot file format
        call get_command_argument(c+1,yes) ; c = c + 1 ; yes = Upper_Case(trim(adjustl(yes)))
        pp%tec = (trim(adjustl(yes))=='YES')
      case('-vtk') ! VKT file format
        call get_command_argument(c+1,yes) ; c = c + 1 ; yes = Upper_Case(trim(adjustl(yes)))
        pp%vtk = (trim(adjustl(yes))=='YES')
      case('-gnu') ! Gnuplot file format
        call get_command_argument(c+1,yes) ; c = c + 1 ; yes = Upper_Case(trim(adjustl(yes)))
        pp%gnu = (trim(adjustl(yes))=='YES')
      case default
        write(stderr,'(A)') ' Error: switch "'//trim(adjustl(ca_switch))//'" is unknown'
        call print_usage
        stop
      endselect
    enddo
  endif
  if (pp%tec) call IOFile%tec%set(name=trim(adjustl(outfname)))
  if (pp%vtk) call IOFile%vtk%set(name=trim(adjustl(outfname)))
  if (pp%gnu) call IOFile%gnu%set(name=trim(adjustl(outfname)))
  if (pp%node.and.pp%bc) then
    write(stderr,'(A)') ' It is not possible to save bc ghost cells and node-centered interpolated variables!'
    stop
  endif
  if (.not.allocated(IOFile%mesh%name)) then
    write(stderr,'(A)') ' File name of mesh file must be provided!'
    call print_usage
    stop
  endif
  ! the name of mesh file is used as output file base name
  if (.not.allocated(IOFile%tec%name)) call IOFile%tec%set(name=IOFile%mesh%name)
  if (.not.allocated(IOFile%vtk%name)) call IOFile%vtk%set(name=IOFile%mesh%name)
  if (.not.allocated(IOFile%gnu%name)) call IOFile%gnu%set(name=IOFile%mesh%name)
  pp%meshonly=.true. ; if (allocated(IOFile%sol%name)) pp%meshonly = .false. ! a solution file name has been passed
  IOFile%tec%pp = pp
  IOFile%vtk%pp = pp
  IOFile%gnu%pp = pp
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parsing_command_line

  !> @brief Procedure for initializing the post-processing.
  subroutine pog_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P):: b,l !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing IR_Precision module constants
  call IR_init
  ! initializing compiled code options collection
  call global%cco%init
  write(stdout,'(A)')'+--> Compiled code used options'
  call global%cco%print(unit=stdout,pref='|-->')
  ! parsing command line for getting global option file name
  call parsing_command_line
 ! setting files paths: the command line arguments have full paths
  if (pp%tec) call IOFile%tec%set( path_in='',path_out='')
  if (pp%vtk) call IOFile%vtk%set( path_in='',path_out='')
  if (pp%gnu) call IOFile%gnu%set( path_in='',path_out='')
              call IOFile%mesh%set(path_in='',path_out='')
              call IOFile%sol%set( path_in='',path_out='')
              call IOFile%proc%set(path_in='',path_out='')
  ! loading input files
  write(stdout,'(A)')'+--> Loading input files'
  if (allocated(IOFile%proc%name)) then
    call IOFile%proc%load(mesh_dims=global%mesh_dims,parallel=global%parallel)
    if (IOFile%proc%iostat/=0) then
      write(stderr,'(A)')'+--> '//IOFile%proc%iomsg
      stop
    endif
    call IOFile%mesh%load_header(mesh_dims=global%mesh_dims)
    if (IOFile%mesh%iostat/=0) then
      write(stderr,'(A)')'+--> '//IOFile%mesh%iomsg
      stop
    endif
  else
    call IOFile%mesh%load_header(mesh_dims=global%mesh_dims)
    if (IOFile%mesh%iostat/=0) then
      write(stderr,'(A)')'+--> '//IOFile%mesh%iomsg
      stop
    endif
    global%mesh_dims%Nb = global%mesh_dims%Nb_tot
    call global%parallel%set_serial(Nb_tot=global%mesh_dims%Nb_tot)
  endif
  call global%parallel%print(pref='|-->    ',unit=stdout)
  ! loading mesh file
  write(stdout,'(A)')'+-->   Loading '//IOFile%mesh%name
  call IOFile%mesh%load(global=global)
  if (IOFile%mesh%iostat/=0) then
    write(stderr,'(A)')'+--> '//IOFile%mesh%iomsg
    stop
  endif
  ! loading solution file
  if (.not.pp%meshonly) then
    write(stdout,'(A)')'+-->   Loading '//IOFile%sol%name
    call IOFile%sol%load(global=global)
    if (IOFile%sol%iostat/=0) then
      write(stderr,'(A)')'+--> '//IOFile%init%iomsg
      stop
    endif
  endif
  ! computing the mesh variables that are not loaded from input files
  do l=1,global%mesh_dims%Nl ; do b=1,global%mesh_dims%Nb
      call global%block(b,l)%metrics
      call global%block(b,l)%metrics_correction
  enddo ; enddo
  ! printing block infos
  write(stdout,'(A)')'+-->   Blocks infos'
  do l=1,global%mesh_dims%Nl ; do b=1,global%mesh_dims%Nb
    write(stdout,'(A)')'+-->     Block b='//trim(str(n=b))//' level l='//trim(str(n=l))
    call global%block(b,l)%print(unit=stdout,pref='|-->      ')
  enddo ; enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pog_init
endprogram POG
