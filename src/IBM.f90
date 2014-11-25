!> @ingroup Program
!> @{
!> @defgroup IBMProgram IBM
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup IBMPrivateVarPar IBM
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup IBMPrivateProcedure IBM
!> @}

!> @brief IBM, Initial and Boundary conditions, Mesh generator for @off (Open Finite volume Fluid dynamics code).
!> This is an auxiliary tool useful for building proper inputs for @off code. ICG can build Initial Conditions files, Mesh files and
!> Boundary Conditions files. It accepts two kinds of inputs: \n
!> - Direct Blocks Description: this is the simplest available input. The initial and boundary descriptions as well as the geometry
!>   are directly described by means of simple ascii files. This kind of inputs can describe only simple Cartesian grids.
!> - Ansys (http://www.ansys.com) IcemCFD Multiblock INFO importer: this is a more complex (but more flexible) input. The initial
!>   conditions are described by means of simple ascii files similar to the Direct Block Description input, but the boundary
!>   conditions and the geometry are loaded by Ansys IcemCFD Multiblock INFO files. These files can describe more complex scenario
!>   with general curvilinear grids.
!> @todo \b DocImprove: Improve the documentation
!> @ingroup IBMProgram
program IBM
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                        ! Integers and reals precision definition.
USE Data_Type_Command_Line_Interface, only: Type_Command_Line_Interface ! Definition of Type_Command_Line_Interface.
USE Data_Type_Files,                  only: Type_Files                  ! Definition of files structures.
USE Data_Type_Global,                 only: Type_Global                 ! Definition of Type_Global.
USE Data_Type_SBlock,                 only: Type_SBlock                 ! Definition of Type_SBlock.
USE Lib_IO_Misc                                                         ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
!> @ingroup IBMPrivateVarPar
!> @{
type(Type_Files)::         IOFile !< Input/Output files.
type(Type_Global),target:: global !< Global data.
integer(I4P)::             b,r    !< Counters.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! initializing IBM
call ibm_init
! loading the blocks description
associate(ibm_opts=>IOFile%ibm_opts,cart_blks=>IOFile%cart_blks,mesh=>IOFile%mesh,bc=>IOFile%bc,init=>IOFile%init,&
          species0=>global%species0,time_step=>global%time_step,rks=>global%parallel%rks)
  select case(ibm_opts%desc_type)
  case('BLOCKS')
    call cart_blks%load
    if (cart_blks%Nr>0) then
      do r=1,cart_blks%Nr
        write(stdout,'(A)')'+-'//rks//'->'//' Region r('//trim(str(.true.,r))//')'
        call cart_blks%region(r)%print(unit=stdout,pref='|-'//rks//'->  ')
      enddo
    endif
    call cart_blks%compute_mesh_dims(global=global,ref_ratio=ibm_opts%ref_ratio)
    call mesh%save_header(global=global)
    call mesh%save_block(mesh_nodes_tag_open=.true.)
    call   bc%save_block(        bc_tag_open=.true.)
    call init%save_header(species0=species0,time_step=time_step)
    call init%save_block(fluid_cells_tag_open=.true.)
    do b=1,cart_blks%Nb
      associate(block=>cart_blks%block(b))
        write(stdout,'(A)')'+-'//rks//'->'//' Block b('//trim(str(.true.,b))//')'
        call block%create_uniform_grid
        call block%alloc(members=.true.)
        call block%set_cells_bc(block_dims=global%block_dims)
        call block%metrics
        call block%metrics_correction
             block%C%P = block%IC
        call block%update_primitive(species0=species0)
        if (cart_blks%Nr>0) then
          do r=1,cart_blks%Nr
            associate(region=>cart_blks%region(r))
              if (region%center>=block%exts%emin.and.region%center<=block%exts%emax) then
                write(stdout,'(A)') '+-'//rks//'->'//' Region r('//trim(str(.true.,r))//&
                  ') overlapps block b('//trim(str(.true.,b))//')'
                call block%set_region_ic(region=region,species0=species0)
              endif
            endassociate
          enddo
        endif
        call mesh%save_block(block=block,ID=int(b,kind=I8P))
        call bc%save_block(block=block,ID=int(b,kind=I8P))
        call init%save_block(block=block,ID=int(b,kind=I8P))
        call block%print(unit=stdout,pref='|-'//rks//'->  ')
        call block%free(cells_data=.true.)
      endassociate
    enddo
    call mesh%save_block( mesh_nodes_tag_close=.true.)
    call   bc%save_block(         bc_tag_close=.true.)
    call init%save_block(fluid_cells_tag_close=.true.)
  case('ICEMCFD')
  case('GRIDGEN')
  case default
    write(stderr,'(A)')'+-'//rks//'->'//' Attentions: unknown type of blocks description! Valid values are:'
    write(stderr,'(A)')'|-'//rks//'->'//' "BLOCKS", "ICEMCFD", "GRIDGEN"'
  endselect
endassociate
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup IBMPrivateProcedure
  !> @{
  !> @brief Procedure for parsing Command Line Arguments (CLA) implementing IBM Command Line Interface (CLI).
  function parsing_command_line(File_Option) result(error)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(OUT)::         File_Option !< Options file name.
  type(Type_Command_Line_Interface):: cli         !< Command Line Interface (CLI).
  integer(I4P)::                      error       !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(rks=>global%parallel%rks)
    ! initializing CLI
    call cli%init(progname='IBM',examples=['IBM IBM_options_file'])
    ! setting CLAs
    call cli%add(pref='|-'//rks//'-> ',positional=.true.,position=1,help='IBM options file',required=.true.,error=error)
    if (error/=0) return
    ! parsing CLI
    write(stdout,'(A)')'+-'//rks//'-> Parsing Command Line Arguments'
    call cli%parse(error=error,pref='|-'//rks//'-> ') ; if (error/=0) return
    ! using CLI data to set POG behaviour
    call cli%get(position=1_I4P, val=File_Option, error=error,pref='|-'//rks//'-> ')
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction parsing_command_line

  !> Subroutine for initializing the initial conditions according to the input options.
  subroutine ibm_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(99):: File_Option !< Options file name.
  integer(I4P)::  error       !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing IR_Precision module constants
  call IR_init
  ! initializing parallel environments
  call global%parallel%init
  ! initializing compiled code options collection
  call global%cco%init
  ! parsing command line for getting global option file name
  error = parsing_command_line(File_Option=File_Option) ; if (error/=0) stop
  ! loading IBM options file
  associate(rks=>global%parallel%rks)
    write(stdout,'(A)')'+-'//rks//'-> Compiled code used options'
    call global%cco%print(unit=stdout,pref='|-'//rks//'->')
    ! setting IBM options file structures
    call IOFile%ibm_opts%set(name=trim(File_Option),path_in='',errpref='+-'//rks//'->')
    associate(OS=>global%OS,ibm_opts=>IOFile%ibm_opts)
      call ibm_opts%load(OS=OS)
      write(stdout,'(A)')'+-'//rks//'-> IBM options'
      call ibm_opts%print(pref='|-'//rks//'->  ',unit=stdout)
      ibm_opts%iostat=OS%make_dir(directory=ibm_opts%path_out)
      call ibm_opts%save
    endassociate
    ! setting files strucutures
    associate(ibm_opts=>IOFile%ibm_opts,cart_blks=>IOFile%cart_blks,spec=>IOFile%spec,mesh=>IOFile%mesh,bc=>IOFile%bc,&
              init=>IOFile%init)
      select case(ibm_opts%desc_type)
      case('BLOCKS')
        call cart_blks%set(name=ibm_opts%fn_blk_desc(1)%vs,path_in=ibm_opts%path_in,errpref='+-'//rks//'->')
      case('ICEMCFD')
      case('GRIDGEN')
      case default
        write(stderr,'(A)')'+-'//rks//'-> Attentions: unknown type of blocks description! Valid values are:'
        write(stderr,'(A)')'|-'//rks//'   1) "BLOCKS"'
        write(stderr,'(A)')'|-'//rks//'   2) "ICEMCFD"'
        write(stderr,'(A)')'|-'//rks//'   3) "GRIDGEN"'
        stop
      endselect
      call spec%set(name=ibm_opts%fn_spec,path_in =ibm_opts%path_in ,errpref='+-'//rks//'->')
      call mesh%set(name=ibm_opts%fn_mesh,path_out=ibm_opts%path_out,errpref='+-'//rks//'->')
      call bc%set(  name=ibm_opts%fn_bc  ,path_out=ibm_opts%path_out,errpref='+-'//rks//'->')
      call init%set(name=ibm_opts%fn_init,path_out=ibm_opts%path_out,errpref='+-'//rks//'->')
    endassociate
    ! loading initial species
    call IOFile%spec%load(species=global%species0)
    write(stdout,'(A)')'+-'//rks//'-> Initial species'
    call global%species0%print(pref='|-'//rks//'->  ',unit=stdout)
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine ibm_init
  !> @}
endprogram IBM
