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
USE IR_Precision                                                        ! Integers and reals precision definition.
USE Data_Type_Command_Line_Interface, only: Type_Command_Line_Interface ! Definition of Type_Command_Line_Interface.
USE Data_Type_Files,                  only: Type_Files                  ! Definition of Type_Files.
USE Data_Type_Global,                 only: Type_Global                 ! Definition of Type_Global.
USE Data_Type_OS,                     only: Type_OS                     ! Definition of Type_OS.
USE Data_Type_PostProcess,            only: Type_PostProcess            ! Definition of Type_PostProcess.
USE Data_Type_SBlock,                 only: Type_SBlock                 ! Definition of Type_SBlock.
USE Data_Type_Varying_String                                            ! Definition of Type_Varying_String.
USE Lib_IO_Misc                                                         ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Global)          :: global !< Global-level data.
type(Type_Files)           :: IOFile !< Input/Output files.
type(Type_PostProcess)     :: pp     !< Post-process data.
type(Type_SBlock), pointer :: block  !< Pointer for scanning global%block tree.
integer(I8P)               :: ID     !< Counter.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! initializing the post-processing
call pog_init
associate(rks=>global%parallel%rks)
    write(stdout,'(A)')'+-'//rks//'-> Loading input files'
    write(stdout,'(A)')'+-'//rks//'->   Loading '//IOFile%mesh%name
    ! call IOFile%mesh%load_header(global=global)
    call IOFile%mesh%load(global=global)
    if (.not.pp%meshonly) then
      write(stdout,'(A)')'+-'//rks//'->   Loading '//IOFile%sol%name
      call IOFile%sol%load(global=global)
    endif
    ! computing the mesh variables that are not loaded from input files
    print*,'Cazzo'
    do while(global%block%loopID(ID=ID))
      block => global%block%dat(ID=ID)
      call block%metrics
      call block%metrics_correction
      write(stdout,'(A)')'+-'//rks//'->     Block b='//trim(str(n=ID))
      call block%print(unit=stdout,pref='|-'//rks//'->      ')
    enddo
    print*,'Cazzo finale'
    ! saving output files
    if (pp%tec) then
      write(stdout,'(A)')'+-'//rks//'-> Saving '//IOFile%tec%name
      call IOFile%tec%save(global=global)
    endif
    if (pp%vtk) then
      write(stdout,'(A)')'+-'//rks//'-> Saving '//IOFile%vtk%name
      call IOFile%vtk%save(global=global)
    endif
    if (pp%gnu) then
      write(stdout,'(A)')'+-'//rks//'-> Saving '//IOFile%gnu%name
      call IOFile%gnu%save(global=global)
    endif
endassociate
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Procedure for parsing Command Line Arguments (CLA) implementing POG Command Line Interface (CLI).
  subroutine parsing_command_line()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(3)::                      yes      !< Yes (no) flag.
  character(99)::                     filename !< File names dummy string.
  character(len=:), allocatable::     outfname !< Output file name dummy string.
  type(Type_Command_Line_Interface):: cli      !< Command Line Interface (CLI).
  integer(I4P)::                      error    !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing CLI
  call cli%init(progname='POG',                                                            &
                examples=['POG -m example.geo -s example.sol -o example                 ', &
                          'POG -m example.geo -s example.sol -o example -ascii -cell -bc', &
                          'POG -m example.geo -ascii                                    '])
  ! setting CLAs
  call cli%add(switch='-m',      help='Mesh file name',               required=.true., act='store'                   ,error=error)
  call cli%add(switch='-o',      help='Output file name',             required=.false.,act='store',     def=''       ,error=error)
  call cli%add(switch='-s',      help='Solution file name',           required=.false.,act='store',     def=''       ,error=error)
  call cli%add(switch='-proc',   help='Processes/block map file name',required=.false.,act='store',     def=''       ,error=error)
  call cli%add(switch='-bc',     help='Save BC cells',                required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(switch='-cell',   help='Save data at cells center',    required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(switch='-ascii',  help='Save ascii output',            required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(switch='-schl',   help='Save (pseudo) Schlieren field',required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(switch='-mirrorX',help='Save also a X-mirrored copy',  required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(switch='-mirrorY',help='Save also a Y-mirrored copy',  required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(switch='-mirrorZ',help='Save also a Z-mirrored copy',  required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(switch='-tec',    help='Save output in Tecplot format',required=.false.,act='store',     def='yes'    ,error=error)
  call cli%add(switch='-vtk',    help='Save output in VTK format',    required=.false.,act='store',     def='no'     ,error=error)
  call cli%add(switch='-gnu',    help='Save output in VTK format',    required=.false.,act='store',     def='no'     ,error=error)
  associate(rks=>global%parallel%rks)
    ! parsing CLI
    write(stdout,'(A)')'+-'//rks//'-> Parsing Command Line Arguments'
    call cli%parse(error=error,pref='|-'//rks//'->') ; if (error/=0) stop
    ! using CLI data to set POG behaviour
    call cli%get(switch='-m',val=filename,pref='|-'//rks//'->',error=error)
    call IOFile%mesh%set(name=global%OS%string_separator_fix(string=trim(adjustl(filename))),errpref='+-'//rks//'->')
    if (cli%passed(switch='-o')) then
      call cli%get(switch='-o',val=filename,pref='|-'//rks//'->',error=error)
      outfname=global%OS%string_separator_fix(string=trim(adjustl(filename)))
    endif
    if (cli%passed(switch='-s')) then
      call cli%get(switch='-s',val=filename,pref='|-'//rks//'->',error=error)
      call IOFile%sol%set(name=global%OS%string_separator_fix(string=trim(adjustl(filename))),errpref='+-'//rks//'->')
    endif
    if (cli%passed(switch='-proc')) then
      call cli%get(switch='-proc',val=filename,pref='|-'//rks//'->',error=error)
      ! call IOFile%proc%set(name=global%OS%string_separator_fix(string=trim(adjustl(filename))),errpref='+-'//rks//'->')
    endif
    call cli%get(switch='-bc'     ,val=pp%bc       ,pref='|-'//rks//'->',error=error)
    call cli%get(switch='-cell'   ,val=pp%node     ,pref='|-'//rks//'->',error=error);pp%node = .not.pp%node
    call cli%get(switch='-ascii'  ,val=pp%binary   ,pref='|-'//rks//'->',error=error);pp%binary = .not.pp%binary
    call cli%get(switch='-schl'   ,val=pp%schlieren,pref='|-'//rks//'->',error=error)
    call cli%get(switch='-mirrorX',val=pp%mirrorX  ,pref='|-'//rks//'->',error=error)
    call cli%get(switch='-mirrorY',val=pp%mirrorY  ,pref='|-'//rks//'->',error=error)
    call cli%get(switch='-mirrorZ',val=pp%mirrorZ  ,pref='|-'//rks//'->',error=error)
    call cli%get(switch='-tec'    ,val=yes         ,pref='|-'//rks//'->',error=error);pp%tec=(Upper_Case(trim(adjustl(yes)))=='YES')
    call cli%get(switch='-vtk'    ,val=yes         ,pref='|-'//rks//'->',error=error);pp%vtk=(Upper_Case(trim(adjustl(yes)))=='YES')
    call cli%get(switch='-gnu'    ,val=yes         ,pref='|-'//rks//'->',error=error);pp%gnu=(Upper_Case(trim(adjustl(yes)))=='YES')
    if (allocated(outfname)) then
      if (pp%tec) call IOFile%tec%set(name=trim(adjustl(outfname)),errpref='+-'//rks//'->')
      if (pp%vtk) call IOFile%vtk%set(name=trim(adjustl(outfname)),errpref='+-'//rks//'->')
      if (pp%gnu) call IOFile%gnu%set(name=trim(adjustl(outfname)),errpref='+-'//rks//'->')
    endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parsing_command_line

  !> @brief Procedure for initializing the post-processing.
  subroutine pog_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing IR_Precision module constants
  call IR_init
  ! initializing parallel environments
  call global%parallel%init
  ! initializing compiled code options collection
  call global%cco%init
  ! parsing command line
  call parsing_command_line
  associate(rks=>global%parallel%rks)
    write(stdout,'(A)')'+-'//rks//'-> Compiled code used options'
    call global%cco%print(unit=stdout,pref='|-'//rks//'->')
    if (pp%node.and.pp%bc) then
      write(stderr,'(A)')'+-'//rks//'-> It is not possible to save bc ghost cells and node-centered interpolated variables!'
      stop
    endif
    ! the name of mesh file is used as output file base name if output file name has been exeplicitely declared
    if (.not.allocated(IOFile%tec%name)) call IOFile%tec%set(name=IOFile%mesh%name,errpref='+-'//rks//'->')
    if (.not.allocated(IOFile%vtk%name)) call IOFile%vtk%set(name=IOFile%mesh%name,errpref='+-'//rks//'->')
    if (.not.allocated(IOFile%gnu%name)) call IOFile%gnu%set(name=IOFile%mesh%name,errpref='+-'//rks//'->')
    pp%meshonly=.true. ; if (allocated(IOFile%sol%name)) pp%meshonly = .false. ! a solution file name has been passed
    IOFile%tec%pp = pp
    IOFile%vtk%pp = pp
    IOFile%gnu%pp = pp
    ! setting files paths: the command line arguments have full paths
    if (pp%tec) call IOFile%tec%set( path_in='',path_out='',errpref='+-'//rks//'->')
    if (pp%vtk) call IOFile%vtk%set( path_in='',path_out='',errpref='+-'//rks//'->')
    if (pp%gnu) call IOFile%gnu%set( path_in='',path_out='',errpref='+-'//rks//'->')
                call IOFile%mesh%set(path_in='',path_out='',errpref='+-'//rks//'->')
                call IOFile%sol%set( path_in='',path_out='',errpref='+-'//rks//'->')
                ! call IOFile%proc%set(path_in='',path_out='',errpref='+-'//rks//'->')
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine pog_init
endprogram POG
