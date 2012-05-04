!> @brief @off is an Open source Finite volume Fluid dynamics code.
!> It is written in in standard (compliant) Fortran 2003 with highly modularity as design target. \n
!>
!> The aim of @off is to solve, numerically, the Navier-Stokes equations of fluid dynamics by means of Finite Volume technique.
!> The main features of @off code are the following:
!> - Finite Volume, Godunov-like scheme based on Euler conservation Laws written in fully conservative formulation:
!>   - the extension to viscous Navier-Stokes equations is under developing;
!> - Underling Riemann Problem solver for convective fluxes:
!>   - Approximate Riemann solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm;
!>   - Approximate Riemann solver based on Primitive Variables Linearization algorithm;
!>   - Approximate Riemann solver based on Two Rarefactions algorithm;
!>   - Approximate Riemann solver based on Two Shocks algorithm;
!>   - Approximate Riemann solver based on Adaptive (non iterative) PVL-TR-TS algorithm;
!>   - Approximate Riemann solver based on Adaptive (non iterative) LF-TR algorithm;
!>   - Approximate Riemann solver based on HLLC algorithm;
!>   - Approximate Riemann solver based on Roe linearization.
!>   - Exact Riemann solver based on iterative solution of u-function;
!> - Multi-Species fluids models:
!>   - Partial Densities species conservation (Standard Thermodynamic Model);
!>   - New multi-dimensional \f$\gamma,\eta,\chi\f$ conservation models of Favini, B. et al (under developing);
!> - Multi-Phases fluids models:
!>   - Fully-coupled Lagrangian particles transport model (under developing);
!> - Space numerical integration models:
!>   - \f$1^{st}\f$ order piece-wise constant reconstruction;
!>   - \f$2^{nd}\f$ order TVD linear-wise reconstruction;
!>   - \f$3^{rd},5^{th},7^{th}\f$ orders WENO non-linear reconstruction;
!> - Time approximation models:
!>   - \f$1^{st}\f$ order forward Euler integration;
!>   - \f$2^{nd},3^{rd},4^{th}\f$ orders Strong-Stability-Preserving explicit Runge-Kutta integration;
!> - Local pseudo-time convergence acceleration for steady simulations;
!> - Multi-grid time convergence acceleration:
!>   - Multi-grid model has been already developed, but it is affected by some not still recognized bugs. Testing and bugs fixing
!>     are in progress.
!> - Underling numerical grid models:
!>   - 3D, general curvilinear, body-fitted, structured multi-blocks mesh;
!>   - Adaptive Mesh Refinement, AMR model (under developing);
!>   - Blocks overlapping, overset (Chimera) model (to be developed in future);
!> - Computational parallelism ability:
!>   - Domain decomposition by means of Message Passing Interface (MPI) paradigm providing the ability to use distributed-memory
!>     cluster facilities;
!>   - Fine, local parallelism by means of OpenMP paradigm providing the ability to use shared-memory cluster facilities;
!>   - Fine, local parallelism by means of GPU programming (e.g. CUDA framework) providing the ability to use GPUs cluster
!>     facilities (to be developed in future).
!> @author    Stefano Zaghi
!> @version   0.0.5
!> @date      2012-04-24
!> @copyright GNU Public License version 3.
!> @bug <b>Multi-grid Models</b>: \n Multi-grid time convergence acceleration has been developed, but it is affected by some bugs
!>                                that <em>blow up</em> steady simulations.
!> @todo \b NavierStokesEq: Extension to viscous Navier-Stokes equations
!> @todo \b MultiSpeciesModel: Introducing new multi-dimensional \f$\gamma,\eta,\chi\f$ conservation models of Favini, B. et al
!> @todo \b MultiPhaseModel: Introducing fully-coupled Lagrangian particles transport model
!> @todo \b AMR: Introducing AMR (Adaptive Mesh Refinement) model
!> @todo \b Chimera: Introducing blocks overlapping, overset (Chimera) model
!> @todo \b GPU: Introducing fine, local parallelism by means of GPU programming (e.g. CUDA framework)
!> @todo \b DocImprove: Improve the documentation
!> @todo \b DocMakeFile: Create the documentation of makefile
program OFF
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                        ! Integers and reals precision definition.
USE Data_Type_BC, init_bc=>init                         ! Definition of Type_BC.
USE Data_Type_Globals                                   ! Definition of Type_Global and Type_Block.
USE Data_Type_OS                                        ! Definition of Type_OS.
USE Data_Type_Primitive, init_prim=>init, set_prim=>set ! Definition of Type_Primitive.
USE Data_Type_Probe                                     ! Definition of Type_Probe.
USE Data_Type_Tensor                                    ! Definition of Type_Tensor.
USE Data_Type_Time                                      ! Definition of Type_Time.
USE Lib_Fluidynamic, only: primitive2conservative, &    ! Function for converting primitive variables to conservative ones.
                           conservative2primitive, &    ! Function for converting conservative variables to primitive ones.
                           boundary_conditions,    &    ! Subroutine for imposing the boundary conditions.
                           solve_grl,              &    ! Subroutine for solving conservation eq. for a given grid level.
                           rk_init                      ! Subroutine for initializing Runge-Kutta coefficients.
USE Lib_Math,        only: digit                        ! Function for computing the number of digits of an integer.
USE Lib_Mesh,        only: mesh_metrics, &              ! Subroutine for computing metrics of the mesh.
                           mesh_metrics_correction      ! Subroutine for correcting metrics of the mesh according to bc.
USE Lib_IO_Misc                                         ! Procedures for IO and strings operations.
USE Lib_WENO,        only: weno_init                    ! Subroutine for initializing WENO coefficients.
#ifdef OPENMP
USE OMP_LIB                                             ! OpenMP runtime library.
#endif
USE Lib_Parallel,    only: Nthreads, &                  ! Number of threads.
                           Nproc,    &                  ! Number of processes.
                           blockmap, &                  ! Local/global blocks map.
                           procmap_load                 ! Function for loading the proc/blocks and local/global blocks maps.
#ifdef MPI2
USE MPI                                                 ! MPI runtime library.
USE Lib_Parallel,    only: Init_sendrecv                ! Subroutine for initialize send/receive communications.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_Global)::             global         !< Global-level data.
type(Type_Block), allocatable:: block(:,:)     !< Block-level data [1:Nb,1:Nl].
integer(I_P)::                  b              !< Blocks counter.
integer(I_P)::                  l              !< Grid levels counter.
integer(I_P)::                  err            !< Error trapping flag: 0 no errors, >0 error occurs.
integer(I_P)::                  lockfile       !< Locking unit file.
real::                          partial(1:10)  !< Partial time counters for code profiling
character(20)::                 date           !< Actual date.
integer(I_P)::                  myrank         !< Actual rank process.
integer(I_P)::                  Nprb = 0_I_P   !< Number of probes.
type(Type_Probe), allocatable:: probes(:)      !< Probes [1:Nprb].
integer(I_P)::                  unitprobe      !< Probes unit file.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
call off_init ! initializing the simulation

l = 1  ! grid level initializing: using only the finest grid for unsteady simulation

Temporal_Loop: do
  ! computing the solution for the actual time step
  partial(1) = Crono()
  call solve_grl(myrank = myrank, l = l, global = global, block= block(:,l))
  partial(2) = Crono(instant1=partial(1))
  ! saving probes
  if (Nprb>0) then
    if ((mod(global%fluid%n,global%file%probe_out)==0).OR. &
        (global%fluid%t==global%fluid%Tmax).OR.            &
        (global%fluid%n==global%fluid%Nmax)) then
      do b=1,Nprb
        unitprobe = Get_Unit()
        open(unit=unitprobe,&
             file=trim(global%file%Path_OutPut)//'probe'//trim(strz(4,b))//'-N_'//trim(strz(10,global%fluid%n))//'.dat')
        write(unitprobe,FR_P,iostat=err)global%fluid%t
        err = write(unitprobe,FR_P,block(blockmap(probes(b)%b),l)%fluid%P(probes(b)%i,probes(b)%j,probes(b)%k))
        close(unitprobe)
      enddo
    endif
  endif
  ! control sentinel for the temporal Loop
  if ((global%fluid%t==global%fluid%Tmax).OR. &
      (global%fluid%n==global%fluid%Nmax).OR. &
      (global%fluid%residual_stop)) exit Temporal_Loop
  !endif
enddo Temporal_Loop

! saving the final time step solution
do l=1,global%mesh%Nl
  ! converting conservative variables to primitive ones
  do b=1,global%mesh%Nb
    call conservative2primitive(global = global, block = block(b,l))
  enddo
  ! imposing the boundary conditions
  call boundary_conditions(myrank = myrank, l = l, global = global, block= block(:,l))
  ! saving the output file
  do b=1,global%mesh%Nb
    err = save_bfluid(filename=file_name(basename=trim(global%file%Path_OutPut)//global%file%File_Sol,&
                                         suffix='.sol',blk=blockmap(b),grl=l,n=global%fluid%n),       &
                      global=global,block=block(b,l))
  enddo
enddo

! the simulation is done: safe finalizing the simulation
if (myrank==0) then
  err = remove_file('lockfile') ! remove lockfile
  close(global%file%unit_res)   ! close log residuals file
endif
! finalizing parallel environments
#ifdef MPI2
call MPI_FINALIZE(err)
#endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Subroutine for initializing the simulation according to the input options.
  subroutine off_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P)::    b              !< Blocks counter.
  integer(I_P)::    l              !< Grid levels counter.
  integer(I_P)::    c              !< Cell and variable counter.
  integer(I_P)::    Nca = 0        !< Number of command line arguments.
  character(60)::   File_Option    !< Global option file name.
  character(500)::  varname_res    !< Variables name for the gnuplot residuals file.
  character(DI_P):: Ncstr          !< String containing current number id of conservative variables.
  real::            instant0 = 0.0 !< The Crono starting instant (used for timing the code).
  integer(I_P)::    UnitFree       !< Free logic unit.
  logical::         is_file        !< Flag for inquiring the presence of file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing parallel environments
#ifdef MPI2
  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, Nproc,err)
#else
  myrank = 0
#endif
#ifdef OPENMP
  !$OMP PARALLEL      &
  !$OMP DEFAULT(none) &
  !$OMP SHARED(Nthreads)
  Nthreads = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
#endif

  date = Get_Date_String()
  if (myrank==0) then
    ! inquiring the presence of a lockfile
    inquire(file='lockfile',exist=is_file,iostat=err)
    if (is_file) then
      ! a lockfile is present into the working directory; the simulation must be stopped
      write(stderr,'(A)')' Lockfile has been found into the working directory'
      write(stderr,'(A)')' Maybe a previous simulation has been aborted'
      write(stderr,'(A)')' Remove the lockfile before start the new simulation'
#ifdef MPI2
      call MPI_FINALIZE(err)
#endif
      stop
    else
      ! creating a lockfile
      lockfile = Get_Unit()
      open(unit=lockfile,file='lockfile')
      write(lockfile,'(A)')' Simulation started on'
      write(lockfile,'(A)')' '//date
      close(lockfile)
    endif
    ! printing machine precision information
    write(stdout,'(A)',   iostat=err)'----------------------------------------------------------------------'
    write(stdout,'(A)',   iostat=err)' Some information about the precision of runnig machine'
    call IR_Print()
    write(stdout,'(A)',   iostat=err)'----------------------------------------------------------------------'
    write(stdout,*)
    write(stdout,'(A)',   iostat=err)'----------------------------------------------------------------------'
    write(stdout,'(A,I3)',iostat=err)' Number of MPI processes:  ', Nproc
    write(stdout,'(A,I3)',iostat=err)' Number of OpenMP threads: ', Nthreads
    write(stdout,'(A)',   iostat=err)'----------------------------------------------------------------------'
    write(stdout,*)
  endif

  ! parsing command line for getting global option file name
  Nca = command_argument_count ()
  if (Nca==0) then
    write(stderr,'(A,I3)')' My RANK is: ',myrank
    write(stderr,'(A)')   ' A valid file name of the options file must be provided as command line argument'
    write(stderr,'(A)')   ' No argument has been passed to command line'
    write(stderr,'(A)')   ' Correct use is:'
    write(stderr,*)
    write(stderr,'(A)')   ' OFF "valid_option_file_name"'
#ifdef MPI2
    call MPI_FINALIZE(err)
#endif
    stop
  else
    call get_command_argument (1, File_Option)
    File_Option = string_OS_sep(File_Option) ; File_Option = adjustl(trim(File_Option))
  endif
  if (myrank==0) then
    write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
    write(stdout,'(A)',iostat=err)' Simulation started on'
    write(stdout,'(A)',iostat=err)' '//date
    write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
    write(stdout,*)
    write(stdout,'(A)',iostat=err)'----------------------------------------------------------------------'
    write(stdout,'(A)',iostat=err)' Loading input files'
    write(stdout,'(A)',iostat=err)'  Loading '//trim(File_Option)
  endif

  ! loading input options
  err = load_option_file(myrank = myrank, filename = File_Option, global = global)
  if (myrank==0) then
    write(stdout,'(A)',iostat=err)'  Loading '//trim(global%file%Path_InPut)//trim(global%file%File_Solver)
  endif
  ! loading solver options
  err = load_gfluid_soption(myrank=myrank, filename=trim(global%file%Path_InPut)//trim(global%file%File_Solver), global=global)
  ! loading processes/blocks map and computing the number global/local blocks
  err = procmap_load(myrank = myrank, filename = trim(global%file%Path_InPut)//'procmap.dat', global = global)
  ! loading the number of initial species from the first block, finest grid level, fluid file
  err = load_gfluid_Ns(binary   = .true.,                                                                         &
                       myrank   = myrank,                                                                         &
                       filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_Init), &
                                            suffix   = '.itc',                                                    &
                                            blk      = 1,                                                         &
                                            grl      = 1),                                                        &
                       global   = global)
  if (myrank==0) then
    write(stdout,*)
  endif

  ! allocating global fluidynamic data
  call alloc_global_fluid(global=global)
  ! allocating blocks data
  if (allocated(block)) then
    do l=lbound(block,dim=2),ubound(block,dim=2)
      do b=lbound(block,dim=1),ubound(block,dim=1)
        call free_block(block=block(b,l))
      enddo
    enddo
    deallocate(block)
  endif
  allocate(block(1:global%mesh%Nb,1:global%mesh%Nl))
  do l=1,global%mesh%Nl
    do b=1,global%mesh%Nb
      ! getting dimensions of block
      err = load_bmesh_dims(myrank   = myrank,                                                                         &
                            filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_Mesh), &
                                                 suffix   = '.geo',                                                    &
                                                 blk      = blockmap(b),                                               &
                                                 grl      = l),                                                        &
                            block    = block(b,l))
      ! allocating block
      call alloc_block(global=global,block=block(b,l))
    enddo
  enddo

  ! loading blocks data
  do l=1,global%mesh%Nl
    do b=1,global%mesh%Nb
      ! loading mesh
      err = load_bmesh(myrank   = myrank,                                                                         &
                       filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_Mesh), &
                                            suffix   = '.geo',                                                    &
                                            blk      = blockmap(b),                                               &
                                            grl      = l),                                                        &
                       block    = block(b,l))
      ! loading boundary conditions
      err = load_bbc(myrank   = myrank,                                                                       &
                     filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_BC), &
                                          suffix   = '.bco',                                                  &
                                          blk      = blockmap(b),                                             &
                                          grl      = l),                                                      &
                     block    = block(b,l))
      ! loading initial conditions
      err = load_bfluid(myrank   = myrank,                                                                         &
                        filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_Init), &
                                             suffix   = '.itc',                                                    &
                                             blk      = blockmap(b),                                               &
                                             grl      = l),                                                        &
                        global   = global,block=block(b,l))
      ! converting primitive variables to conservative ones
      call primitive2conservative(block = block(b,l))
      ! print some informations of the fluid data loaded
      err = print_info_bfluid(myrank=myrank,blk=b,grl=l,global=global,block=block(b,l))
    enddo
  enddo

  ! initializing Runge-Kutta coefficients
  call rk_init(global%fluid%rk_ord)

  select case(global%fluid%sp_ord)
  case(1_I1P) ! 1st order piecewise constant reconstruction
    global%mesh%gco = 1_I1P
  case(3_I1P) ! 3rd order WENO reconstruction
    global%mesh%gco = 2_I1P
  case(5_I1P) ! 5th order WENO reconstruction
    global%mesh%gco = 3_I1P
  case(7_I1P) ! 7th order WENO reconstruction
    global%mesh%gco = 4_I1P
  endselect

  ! coping input files in output path
  err = make_dir(myrank,global%file%Path_OutPut) ! creating the output directory
  if (myrank==0) then
    err = copy_file(source_file = trim(global%file%Path_InPut)//File_Option, &
                    target_file = trim(global%file%Path_OutPut)//File_Option)
    err = copy_file(source_file = trim(global%file%Path_InPut)//global%file%File_Solver, &
                    target_file = trim(global%file%Path_OutPut)//global%file%File_Solver)
#ifdef MPI2
    err = copy_file(source_file = trim(global%file%Path_InPut)//'procmap.dat', &
                    target_file = trim(global%file%Path_OutPut)//'procmap.dat')
#endif
  endif
  do l=1,global%mesh%Nl
    do b=1,global%mesh%Nb
      err = copy_file(source_file = file_name(basename=trim(global%file%Path_InPut)//trim(global%file%File_Mesh),suffix='.geo', &
                                              blk=blockmap(b),grl=l), &
                      target_file = file_name(basename=trim(global%file%Path_OutPut)//trim(global%file%File_Mesh),suffix='.geo', &
                                              blk=blockmap(b),grl=l))
      err = copy_file(source_file = file_name(basename=trim(global%file%Path_InPut)//trim(global%file%File_BC),suffix='.bco', &
                                              blk=blockmap(b),grl=l), &
                      target_file = file_name(basename=trim(global%file%Path_OutPut)//trim(global%file%File_BC),suffix='.bco', &
                                              blk=blockmap(b),grl=l))
      err = copy_file(source_file = file_name(basename=trim(global%file%Path_InPut)//trim(global%file%File_Init),suffix='.itc', &
                                              blk=blockmap(b),grl=l), &
                      target_file = file_name(basename=trim(global%file%Path_OutPut)//trim(global%file%File_Init),suffix='.itc', &
                                              blk=blockmap(b),grl=l))
    enddo
  enddo

  ! print some informations of the initial data loaded
  !err = fluid_print_info(myrank)

  ! computing the mesh variables that are not loaded from input files
  do l=1,global%mesh%Nl
    do b=1,global%mesh%Nb
      call mesh_metrics(           block= block(b,l))
      call mesh_metrics_correction(block= block(b,l))
      ! print some informations of the mesh data loaded
      err = print_info_bmesh(myrank=myrank,blk=b,grl=l,global=global,block=block(b,l))
    enddo
  enddo

  ! initializing WENO coefficients
  call weno_init(myrank=myrank,global=global,block=block,S=global%mesh%gco)

  ! initializing the log file of residuals
  if (myrank==0) then
    Ncstr = adjustl(trim(str(.true.,global%fluid%Nc)))
    ! creating the gnuplot script file for visualizing the residuals log file
    global%file%unit_res = Get_Unit()
    open(unit=global%file%unit_res,file=trim(global%file%Path_OutPut)//'gplot_res')
    write(global%file%unit_res,'(A)')'set xlabel "Iteration"'
    write(global%file%unit_res,'(A)')'set ylabel "Residuals"'
    write(global%file%unit_res,'(A)')'set log y'
    write(global%file%unit_res,'(A)')"p 'residuals.log' u 1:3 w l title 'R1', "//char(92)
    do c=2,global%fluid%Nc-1
      global%file%varform_res = '(A,I'//trim(str(.true.,digit(2+c)))//',A,I'//trim(str(.true.,digit(c)))//',A)'
      write(global%file%unit_res,trim(global%file%varform_res)) "  'residuals.log' u 1:",2+c," w l title 'R",c,"', "//char(92)
    enddo
    global%file%varform_res = '(A,I'//trim(str(.true.,digit(2+global%fluid%Nc)))// &
                              ',A,I'//trim(str(.true.,digit(global%fluid%Nc)))//',A)'
    write(global%file%unit_res,trim(global%file%varform_res)) "  'residuals.log' u 1:", &
                                                              2+global%fluid%Nc,        &
                                                              " w l title 'R",global%fluid%Nc,"'"
    close(global%file%unit_res)
    ! initialize gnuplot log file of residuals
    if (global%fluid%n>0) then
      global%file%unit_res = Get_Unit()
      open(unit=global%file%unit_res,file=trim(global%file%Path_OutPut)//'residuals.log',position='APPEND')
    else
      global%file%unit_res = Get_Unit()
      open(unit=global%file%unit_res,file=trim(global%file%Path_OutPut)//'residuals.log')
    endif
    ! initialize header
    global%file%varform_res = '('//trim(str(.true.,global%fluid%Nc))//&
      '(A,I'//trim(str(.true.,digit(global%fluid%Nc)))//'.'//trim(str(.true.,digit(global%fluid%Nc)))//',A))'
    write(varname_res,trim(global%file%varform_res))('"R',c,'",',c=1,global%fluid%Nc)
    varname_res = varname_res(1:len_trim(varname_res)-1)
    write(global%file%unit_res,'(A)')'# L2 norm of residuals'
    write(global%file%unit_res,'(A)')'# Simulation started on '//date
    write(global%file%unit_res,'(A)')'# "n","t",'//adjustl(trim(varname_res))
    ! initialize output format
    global%file%varform_res ='('//FI8P//',1X,'//adjustl(trim(str(.true.,global%fluid%Nc+1)))//'('//FR_P//',1X))'
  endif

  ! initialize the multi-processes send/recive comunications and doing the first comunication if necessary
#ifdef MPI2
  call Init_sendrecv(myrank=myrank,global=global,block=block)
#endif

  ! allocate multigrid variables if necessary
  !if (Nl>1) then
    !call alloc_multigrid(Nc = Nc, bb1 = bb(1,1,2), bb2 = bb(2,Nb,Nl))
  !endif

  ! checking the presence of probes file
  inquire(file=trim(global%file%Path_InPut)//'probes.dat',exist=is_file,iostat=err)
  if (is_file) then
    UnitFree = Get_Unit() ; open(unit = UnitFree, file = trim(global%file%Path_InPut)//'probes.dat', action = 'READ')
    read(UnitFree,*)Nprb
    if (allocated(probes)) deallocate(probes) ; allocate(probes(1:Nprb))
    do b=1,Nprb
      read(UnitFree,*)probes(b)%b,probes(b)%i,probes(b)%j,probes(b)%k
    enddo
    close(UnitFree)
  endif

  instant0 = Crono(start=.true.) ! code timing start
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine off_init

  !> @brief Function for loading global file options.
  function load_option_file(myrank,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where option variables are saved.
  type(Type_Global), intent(INOUT):: global   !< Global-level data.
  integer(I_P)::                     err      !< Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of file.
  character(3)::                     os_type  !< Type operating system.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) then
      call File_Not_Found(myrank,filename,'load_option_file')
  endif
  UnitFree = Get_Unit()
  open(unit = UnitFree, file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err) !               OS
  read(UnitFree,*,iostat=err)os_type
  read(UnitFree,*,iostat=err) !               INPUT OPTIONS
  read(UnitFree,*,iostat=err)global%file%Path_InPut
  read(UnitFree,*,iostat=err)global%file%File_Solver
  read(UnitFree,*,iostat=err)global%mesh%Nl
  read(UnitFree,*,iostat=err)global%file%File_Mesh
  read(UnitFree,*,iostat=err)global%file%File_BC
  read(UnitFree,*,iostat=err)global%file%File_Init
  read(UnitFree,*,iostat=err) !               OUTPUT OPTIONS
  read(UnitFree,*,iostat=err)global%file%Path_OutPut
  read(UnitFree,*,iostat=err)global%file%File_Sol
  read(UnitFree,*,iostat=err)global%file%Screen_out
  read(UnitFree,*,iostat=err)global%file%sol_out
  read(UnitFree,*,iostat=err)global%file%restart_out
  read(UnitFree,*,iostat=err)global%file%probe_out
  close(UnitFree)
  os_type = Upper_Case(os_type) ; OS = init(c_id=os_type)

  global%file%Path_InPut =string_OS_sep(global%file%Path_InPut ) ; global%file%Path_InPut  = adjustl(trim(global%file%Path_InPut ))
  global%file%File_Solver=string_OS_sep(global%file%File_Solver) ; global%file%File_Solver = adjustl(trim(global%file%File_Solver))
  global%file%File_Mesh  =string_OS_sep(global%file%File_Mesh  ) ; global%file%File_Mesh   = adjustl(trim(global%file%File_Mesh  ))
  global%file%File_BC    =string_OS_sep(global%file%File_BC    ) ; global%file%File_BC     = adjustl(trim(global%file%File_BC    ))
  global%file%File_Init  =string_OS_sep(global%file%File_Init  ) ; global%file%File_Init   = adjustl(trim(global%file%File_Init  ))
  global%file%Path_OutPut=string_OS_sep(global%file%Path_OutPut) ; global%file%Path_OutPut = adjustl(trim(global%file%Path_OutPut))
  global%file%File_Sol   =string_OS_sep(global%file%File_Sol   ) ; global%file%File_Sol    = adjustl(trim(global%file%File_Sol   ))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_option_file
endprogram OFF
