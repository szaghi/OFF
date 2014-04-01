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
USE IR_Precision                        ! Integers and reals precision definition.
USE Data_Type_Files,  only: Type_Files  ! Definition of files structures.
USE Data_Type_Global, only: Type_Global ! Definition of Type_Global.
USE Data_Type_SBlock, only: Type_SBlock ! Definition of Type_SBlock.
USE Lib_IO_Misc                         ! Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
!> @ingroup IBMPrivateVarPar
!> @{
type(Type_Files)::         IOFile !< Input/Output files.
type(Type_Global),target:: global !< Global data.
integer(I4P)::             b,l,r  !< Counters.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! initializing IBM
call ibm_init
! loading the blocks description
associate(ibm_opts=>IOFile%ibm_opts,cart_blks=>IOFile%cart_blks,mesh=>IOFile%mesh,bc=>IOFile%bc,init=>IOFile%init,&
          mesh_dims=>IOFile%cart_blks%mesh_dims,Nl=>global%mesh_dims%Nl,species0=>global%species0,time_step=>global%time_step)
  select case(ibm_opts%desc_type)
  case('BLOCKS')
    call cart_blks%load ; if (cart_blks%iostat/=0) write(stderr,'(A)')trim(cart_blks%iomsg)
    if (cart_blks%Nr>0) then
      do r=1,cart_blks%Nr
        write(stdout,'(A)')'+-> Region r('//trim(str(.true.,r))//')'
        call cart_blks%region(r)%print(unit=stdout,pref='|->  ')
      enddo
    endif
    call cart_blks%compute_mesh_dims(Nl=Nl)
    call global%parallel%set_serial(Nb_tot=mesh_dims%Nb)
    call mesh%save_header(mesh_dims=cart_blks%mesh_dims)         ; if (mesh%iostat/=0) write(stderr,'(A)')trim(mesh%iomsg)
    call mesh%save_block(mesh_nodes_tag_open=.true.)             ; if (mesh%iostat/=0) write(stderr,'(A)')trim(mesh%iomsg)
    call   bc%save_block(        bc_tag_open=.true.)             ; if (  bc%iostat/=0) write(stderr,'(A)')trim(  bc%iomsg)
    call init%save_header(species0=species0,time_step=time_step) ; if (init%iostat/=0) write(stderr,'(A)')trim(init%iomsg)
    call init%save_block(fluid_cells_tag_open=.true.)            ; if (init%iostat/=0) write(stderr,'(A)')trim(init%iomsg)
    do l=1,mesh_dims%Nl ; do b=1,mesh_dims%Nb_tot
      associate(block=>IOFile%cart_blks%block(b))
        block%dims = mesh_dims%block_dims(b,l)
        call block%create_uniform_grid
        call block%alloc(members=.true.)
        call block%set_cells_bc(mesh_dims=mesh_dims,l=l)
        call block%metrics
        call block%metrics_correction
             block%C%P = block%IC
        call block%update_primitive(species0=species0)
        if (cart_blks%Nr>0) then
          do r=1,cart_blks%Nr
            associate(region=>cart_blks%region(r))
              if (region%center>=block%exts%emin.and.region%center<=block%exts%emax) then
                write(stdout,'(A)') '+-> Region r('//trim(str(.true.,r))//') overlapps block b('//trim(str(.true.,b))//') l('//&
                  trim(str(.true.,l))//')'
                call block%set_region_ic(region=region,species0=species0)
              endif
            endassociate
          enddo
        endif
        call mesh%save_block(block=block,b=b,l=l,parallel=global%parallel) ; if (mesh%iostat/=0) write(stderr,'(A)')trim(mesh%iomsg)
        call bc%save_block(  block=block,b=b,l=l,parallel=global%parallel) ; if (  bc%iostat/=0) write(stderr,'(A)')trim(  bc%iomsg)
        call init%save_block(block=block,b=b,l=l,parallel=global%parallel) ; if (init%iostat/=0) write(stderr,'(A)')trim(init%iomsg)
        write(stdout,'(A)')'+-> Block b('//trim(str(.true.,b))//') l('//trim(str(.true.,l))//')'
        call block%print(unit=stdout,pref='|->  ')
        call block%free(cells_data=.true.)
      endassociate
    enddo ; enddo
    call mesh%save_block( mesh_nodes_tag_close=.true.) ; if (mesh%iostat/=0) write(stderr,'(A)')trim(mesh%iomsg)
    call   bc%save_block(         bc_tag_close=.true.) ; if (  bc%iostat/=0) write(stderr,'(A)')trim(  bc%iomsg)
    call init%save_block(fluid_cells_tag_close=.true.) ; if (init%iostat/=0) write(stderr,'(A)')trim(init%iomsg)
    ! saving bc file
  case('ICEMCFD')
  case('GRIDGEN')
  case default
    write(stderr,'(A)')' Attentions: unknown type of blocks description! Valid values are "BLOCKS", "ICEMCFD", "GRIDGEN"'
  endselect
endassociate
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup IBMPrivateProcedure
  !> @{
  !> Subroutine for initializing the initial conditions according to the input options.
  subroutine ibm_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P)::  Nca = 0_I4P !< Number of command line arguments.
  character(99):: File_Option !< Options file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing IR_Precision module constants
  call IR_init
  ! parsing command line for getting global option file name
  Nca = command_argument_count()
  if (Nca==0) then
    write(stderr,'(A)')' A valid file name of the options file must be provided as command line argument'
    write(stderr,'(A)')' No argument has been passed to command line'
    write(stderr,'(A)')' Correct use is:'
    write(stderr,'(A)')'   IBM "valid_option_file_name"'
    stop
  else
    call get_command_argument(1, File_Option)
    File_Option = global%OS%string_separator_fix(string=trim(adjustl(File_Option)))
  endif
  ! setting IBM options file structures
  call IOFile%ibm_opts%set(name=trim(File_Option),path_in='')
  ! loading IBM options file
  associate(OS=>global%OS,Nl=>global%mesh_dims%Nl,ibm_opts=>IOFile%ibm_opts)
    call ibm_opts%load(OS=OS) ; if (ibm_opts%iostat/=0) write(stderr,'(A)')ibm_opts%iomsg ; Nl = ibm_opts%Nl
    write(stdout,'(A)')'+-> IBM options'
    call ibm_opts%print(pref='|->  ',unit=stdout)
    ibm_opts%iostat=OS%make_dir(directory=ibm_opts%path_out)
  endassociate
  ! setting files strucutures
  associate(ibm_opts=>IOFile%ibm_opts,cart_blks=>IOFile%cart_blks,spec=>IOFile%spec,mesh=>IOFile%mesh,bc=>IOFile%bc,&
            init=>IOFile%init)
    select case(ibm_opts%desc_type)
    case('BLOCKS')
      call cart_blks%set(name=ibm_opts%fn_blk_desc(1)%vs,path_in=ibm_opts%path_in)
    case('ICEMCFD')
    case('GRIDGEN')
    case default
      write(stderr,'(A)')' Attentions: unknown type of blocks description! Valid values are "BLOCKS", "ICEMCFD", "GRIDGEN"'
    endselect
    call spec%set(name=ibm_opts%fn_spec,path_in =ibm_opts%path_in )
    call mesh%set(name=ibm_opts%fn_mesh,path_out=ibm_opts%path_out)
    call bc%set(  name=ibm_opts%fn_bc  ,path_out=ibm_opts%path_out)
    call init%set(name=ibm_opts%fn_init,path_out=ibm_opts%path_out)
  endassociate
  ! loading initial species
  call IOFile%spec%load(species=global%species0) ; if (IOFile%spec%iostat/=0) write(stderr,'(A)')IOFile%spec%iomsg
  write(stdout,'(A)')'+-> Initial species'
  call global%species0%print(pref='|->  ',unit=stdout)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine ibm_init
  !> @}
endprogram IBM
