!> @addtogroup Program Programs
!> List of excutable programs.
!> @addtogroup DerivedType Derived Types
!> List of derived data types.
!> @addtogroup GlobalVarPar Global Variables and Parameters
!> List of global variables and parameters.
!> @addtogroup PrivateVarPar Private Variables and Parameters
!> List of private variables and parameters.
!> @addtogroup Interface Interfaces
!> List of explicitly defined interface.
!> @addtogroup Library Modules Libraries
!> List of modules containing libraries of procedures.
!> @addtogroup PublicProcedure Public Procedures
!> List of public procedures.
!> @addtogroup PrivateProcedure Private Procedures
!> List of private procedures.

!> @ingroup Program
!> @{
!> @defgroup OFFProgram OFF
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup OFFPrivateVarPar OFF
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup OFFPrivateProcedure OFF
!> @}

!> @brief @off is an Open source Finite volumes Fluid dynamics code.
!> It is written in in standard (compliant) Fortran 2003 with highly modularity as design target. \n
!>
!> The aim of @off is to solve, numerically, the compressible Navier-Stokes equations of fluid dynamics
!> (\ref Equations "see governing equations") by means of Finite Volumes techniques.
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
!>   - \f$3^{rd}\f$,\f$5^{th}\f$,\f$7^{th}\f$ order WENO non-linear reconstruction;
!> - Time approximation models:
!>   - \f$1^{st}\f$ order forward Euler integration;
!>   - \f$2^{nd}\f$,\f$3^{rd}\f$,\f$4^{th}\f$ order Strong-Stability-Preserving explicit Runge-Kutta integration;
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
!> @note
!> <b>Compiling Instructions</b> \n
!> @off is shipped with a makefile for compiling the codes on Unix/GNU Linux architectures. Other OS are not supported.
!> For more details see \ref Compiling "Compiling Instructions".
!> @bug <b>Multi-grid Models</b>: \n Multi-grid time convergence acceleration has been developed, but it is affected by some bugs
!>                                that <em>blow up</em> steady simulations.
!> @todo \b MultiSpeciesModel: Introducing new multi-dimensional \f$\gamma,\eta,\chi\f$ conservation models of Favini, B. et al
!> @todo \b MultiPhaseModel: Introducing fully-coupled Lagrangian particles transport model
!> @todo \b AMR: Introducing AMR (Adaptive Mesh Refinement) model
!> @todo \b Chimera: Introducing blocks overlapping, overset (Chimera) model
!> @todo \b GPU: Introducing fine, local parallelism by means of GPU programming (e.g. CUDA framework)
!> @todo \b DocImprove: Improve the documentation
!> @ingroup OFFProgram
program OFF
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                        ! Integers and reals precision definition.
USE Data_Type_AMRBlock                                  ! Definition of Type_AMRBlock.
USE Data_Type_BC                                        ! Definition of Type_BC.
USE Data_Type_Global                                    ! Definition of Type_Global.
USE Data_Type_OS                                        ! Definition of Type_OS.
USE Data_Type_Primitive                                 ! Definition of Type_Primitive.
USE Data_Type_Probe                                     ! Definition of Type_Probe.
USE Data_Type_SBlock                                    ! Definition of Type_SBlock.
USE Data_Type_Tensor                                    ! Definition of Type_Tensor.
USE Data_Type_Time                                      ! Definition of Type_Time.
USE Lib_Fluidynamic, only: primitive2conservative, &    ! Function for converting primitive variables to conservative ones.
                           conservative2primitive, &    ! Function for converting conservative variables to primitive ones.
                           boundary_conditions,    &    ! Subroutine for imposing the boundary conditions.
                           solve_grl                    ! Subroutine for solving conservation eq. for a given grid level.
USE Lib_IO_Misc                                         ! Procedures for IO and strings operations.
USE Lib_Math,        only: digit                        ! Function for computing the number of digits of an integer.
#ifdef PROFILING
USE Lib_Profiling                                       ! Procedures for profiling the code.
#endif
USE Lib_Runge_Kutta, only: rk_init                      ! Subroutine for initializing Runge-Kutta coefficients.
USE Lib_WENO,        only: weno_init                    ! Subroutine for initializing WENO coefficients.
USE Lib_Parallel,    only: Nthreads, &                  ! Number of threads.
                           Nproc,    &                  ! Number of processes.
                           blockmap, &                  ! Local/global blocks map.
                           procmap_load                 ! Function for loading the proc/blocks and local/global blocks maps.
#ifdef OPENMP
USE OMP_LIB                                             ! OpenMP runtime library.
#endif
#ifdef MPI2
USE MPI                                                 ! MPI runtime library.
USE Lib_Parallel,    only: Init_sendrecv                ! Subroutine for initialize send/receive communications.
#endif
USE Lib_Morton   !  Procedure for Morton's encoding.
USE Data_Type_HashID   !  Procedure for Morton's encoding.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
!> @ingroup OFFPrivateVarPar
!> @{
type(Type_Global), target::      global        !< Global-level data.
type(Type_SBlock), allocatable:: block(:,:)    !< Block-level data [1:Nb,1:Nl].
integer(I_P)::                   b             !< Blocks counter.
integer(I_P)::                   l             !< Grid levels counter.
integer(I_P)::                   err           !< Error trapping flag: 0 no errors, >0 error occurs.
integer(I_P)::                   lockfile      !< Locking unit file.
character(20)::                  date          !< Actual date.
integer(I_P)::                   Nprb = 0_I_P  !< Number of probes.
type(Type_Probe), allocatable::  probes(:)     !< Probes [1:Nprb].
integer(I_P)::                   unitprobe     !< Probes unit file.
type(Type_AMRBlock), allocatable:: amrblock(:) !< AMR grids [1:Nb].
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
call off_init ! initializing the simulation

l = 1  ! grid level initializing: using only the finest grid for unsteady simulation

#ifdef PROFILING
call profile(p=7,pstart=.true.,myrank=global%myrank)
#endif

Temporal_Loop: do
  ! computing the solution for the actual time step
#ifdef PROFILING
  call profile(p=1,pstart=.true.,myrank=global%myrank)
#endif
  call solve_grl(l = l, block= block(:,l))
#ifdef PROFILING
  call profile(p=1,pstop=.true.,myrank=global%myrank)
#endif
  ! saving probes
  if (Nprb>0) then
    if ((mod(global%n,global%file%probe_out)==0).OR. &
        (global%t==global%Tmax).OR.            &
        (global%n==global%Nmax)) then
      do b=1,Nprb
        open(unit=Get_Unit(unitprobe),&
             file=trim(global%file%Path_OutPut)//'probe'//trim(strz(4,b))//'-N_'//trim(strz(10,global%n))//'.dat')
        write(unitprobe,FR_P,iostat=err)global%t
        err = write_primitive(scalar=block(blockmap(probes(b)%b),l)%C(probes(b)%i,probes(b)%j,probes(b)%k)%P, &
                              unit=unitprobe,format=FR_P)
        close(unitprobe)
      enddo
    endif
  endif
  ! control sentinel for the temporal Loop
  if ((global%t==global%Tmax).OR.(global%n==global%Nmax).OR.(global%residual_stop)) exit Temporal_Loop
  !endif
enddo Temporal_Loop

! saving the final time step solution
do l=1,global%Nl
  ! converting conservative variables to primitive ones
  do b=1,global%Nb
    call conservative2primitive(block = block(b,l))
  enddo
  ! imposing the boundary conditions
  call boundary_conditions(l = l, block = block(:,l))
  ! saving the output file
  do b=1,global%Nb
    err = block(b,l)%save_fluid(filename=file_name(basename=trim(global%file%Path_OutPut)//global%file%File_Sol,&
                                                   suffix='.sol',blk=blockmap(b),grl=l,n=global%n))
  enddo
enddo

! the simulation is done: safe finalizing the simulation
if (global%myrank==0) then
  err = remove_file('lockfile') ! remove lockfile
  close(global%file%unit_res)   ! close log residuals file
endif
! finalizing parallel environments
#ifdef MPI2
call MPI_FINALIZE(err)
#endif
#ifdef PROFILING
! finalizing profiling
call profile(p=7,pstop=.true.,myrank=global%myrank)
call profile(finalize=.true.,myrank=global%myrank)
#endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup OFFPrivateProcedure
  !> @{
  !> @brief Subroutine for initializing the simulation according to the input options.
  subroutine off_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P)::    b                 !< Blocks counter.
  integer(I_P)::    l                 !< Grid levels counter.
  integer(I_P)::    c                 !< Cell and variable counter.
  integer(I_P)::    i,j,k             !< Space counters.
  integer(I_P)::    Nca = 0           !< Number of command line arguments.
  character(60)::   File_Option       !< Global option file name.
  character(500)::  varname_res       !< Variables name for the gnuplot residuals file.
  character(DI_P):: Ncstr             !< String containing current number id of conservative variables.
  integer(I_P)::    UnitFree          !< Free logic unit.
  logical::         is_file           !< Flag for inquiring the presence of file.
#ifndef PROFILING
  real(R8P)::       instant0 = 0._R8P !< The Crono starting instant used for profing the code.
#endif
  ! cazzo
  integer(I8P):: i64 , j64, k64
  !integer(I4P):: i4  , j4 , k4
  integer(I2P):: i2  , j2 , k2
  !integer(I1P):: i1  , j1 , k1
  type(Type_HashID)::                   ID    !< ID value.
  ! cazzo
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing parallel environments
#ifdef MPI2
  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,global%myrank,err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,err)
#else
  global%myrank = 0
#endif
#ifdef OPENMP
  !$OMP PARALLEL      &
  !$OMP DEFAULT(none) &
  !$OMP SHARED(Nthreads)
  Nthreads = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
#endif

  date = Get_Date_String()
  if (global%myrank==0) then
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
      open(unit=Get_Unit(lockfile),file='lockfile')
      write(lockfile,'(A)')' Simulation started on'
      write(lockfile,'(A)')' '//date
      close(lockfile)
    endif
    ! printing machine precision information
    write(stdout,'(A)',   iostat=err)'----------------------------------------------------------------------'
    write(stdout,'(A)',   iostat=err)' Some information about the precision of runnig machine'
    call IR_Print()
    ! cazzo
    !allocate(amrblock(1))
    !amrblock(1)%gc = 0
    !amrblock(1)%Ni = 4
    !amrblock(1)%Nj = 4
    !amrblock(1)%Nk = 4
    !call amrblock(1)%init(1_I2P)
    !do k=1,amrblock(1)%Nk
    !  do j=1,amrblock(1)%Nj
    !    do i=1,amrblock(1)%Ni
    !      call ID%build(b=1_I2P,i=i,j=j,k=k,l=0_I1P,p=0_I8P)
    !      call demorton3(ID%bcl,i2,j2,k2)
    !      write(stdout,'(4I10)')ID%bcl,i2,j2,k2
    !    enddo
    !  enddo
    !enddo
    !read(stdout,'(B8.8)') i4
    !write(stdout,'(I10)') i4
    !!read(stdout,'(B64.64)') i64
    !!write(stdout,'(Z64.64)') i64
    !!write(stdout,'(I10)') int(log10(4._R_P)/log10(2._R_P),I1P)
    !i1 = 6 ; j1 = 11 ; k1 = 30
    !i64 = morton2(i=i1,j=j1)
    !write(stdout,'(4I10,1X)') i1,j1,k1,i64
    !call demorton2(key=i64,i=i1,j=j1)
    !write(stdout,'(4I10,1X)') i1,j1,k1
    !i2 = 16 ; j2 = 101 ; k2 = 300
    !i64 = morton2(i=i2,j=j2)
    !write(stdout,'(4I10,1X)') i2,j2,k2,i64
    !call demorton2(key=i64,i=i2,j=j2)
    !write(stdout,'(4I10,1X)') i2,j2,k2
    !i4 = 600 ; j4 = 1001 ; k4 = 3000
    !i64 = morton2(i=i4,j=j4)
    !write(stdout,'(4I10,1X)') i4,j4,k4,i64
    !call demorton2(key=i64,i=i4,j=j4)
    !write(stdout,'(4I10,1X)') i4,j4,k4
    !i64 = treedim(3,5)*100/10**6
    !write(stdout,'(I10)') i64
    !stop 'cazzo'
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
    write(stderr,'(A,I3)')' My RANK is: ',global%myrank
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
  if (global%myrank==0) then
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
  err = load_off_option_file(filename = File_Option, global = global)
  if (global%myrank==0) then
    write(stdout,'(A)',iostat=err)'  Loading '//trim(global%file%Path_InPut)//trim(global%file%File_Solver)
  endif
  ! loading solver options
  err = global%load_fluid_soption(filename=trim(global%file%Path_InPut)//trim(global%file%File_Solver))
  ! loading processes/blocks map and computing the number global/local blocks
  err = procmap_load(filename = trim(global%file%Path_InPut)//'procmap.dat', global = global)
  ! loading the number of initial species from the first block, finest grid level, fluid file
  err = global%load_fluid_Ns(binary   = .true.,                                                                         &
                             filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_Init), &
                                                  suffix   = '.itc',                                                    &
                                                  blk      = 1,                                                         &
                                                  grl      = 1))
  if (global%myrank==0) then
    write(stdout,*)
  endif

  ! allocating global fluidynamic data
  call global%alloc_fluid
  ! allocating blocks data
  if (allocated(block)) then
    do l=lbound(block,dim=2),ubound(block,dim=2)
      do b=lbound(block,dim=1),ubound(block,dim=1)
        call block(b,l)%free
      enddo
    enddo
    deallocate(block)
  endif
  allocate(block(1:global%Nb,1:global%Nl))
  do l=1,global%Nl
    do b=1,global%Nb
      ! updating block pointer to global-level data
      block(b,l)%global => global
      ! getting dimensions of block
      err = block(b,l)%load_mesh_dims(filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_Mesh), &
                                                           suffix   = '.geo',                                                    &
                                                           blk      = blockmap(b),                                               &
                                                           grl      = l))
      ! allocating block
      call block(b,l)%alloc
    enddo
  enddo

  ! loading blocks data
  do l=1,global%Nl
    do b=1,global%Nb
      ! loading mesh
      err = block(b,l)%load_mesh(filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_Mesh), &
                                                      suffix   = '.geo',                                                    &
                                                      blk      = blockmap(b),                                               &
                                                      grl      = l))
      ! loading boundary conditions
      err = block(b,l)%load_bc(filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_BC), &
                                                    suffix   = '.bco',                                                  &
                                                    blk      = blockmap(b),                                             &
                                                    grl      = l))
      ! loading initial conditions
      err = block(b,l)%load_fluid(filename = file_name(basename = trim(global%file%Path_InPut)//trim(global%file%File_Init), &
                                                       suffix   = '.itc',                                                    &
                                                       blk      = blockmap(b),                                               &
                                                       grl      = l))
      ! converting primitive variables to conservative ones
      call primitive2conservative(block=block(b,l))
      ! print some informations of the fluid data loaded
      err = block(b,l)%print_info_fluid(blk=b,grl=l)
    enddo
  enddo

  ! loading inflow boundary conditions if necessary
  ! inflow 1
  do l=1,global%Nl ; do b=1,global%Nb
    do k=1-block(b,l)%gc(5),block(b,l)%Nk+block(b,l)%gc(6)
      do j=1-block(b,l)%gc(3),block(b,l)%Nj+block(b,l)%gc(4)
        do i=0-block(b,l)%gc(1),block(b,l)%Ni+block(b,l)%gc(2)
          if (block(b,l)%Fi(i,j,k)%BC%tp==bc_in1) then
            global%Nin1 = max(global%Nin1,block(b,l)%Fi(i,j,k)%BC%inf)
          endif
        enddo
      enddo
    enddo
    do k=1-block(b,l)%gc(5),block(b,l)%Nk+block(b,l)%gc(6)
      do j=0-block(b,l)%gc(3),block(b,l)%Nj+block(b,l)%gc(4)
        do i=1-block(b,l)%gc(1),block(b,l)%Ni+block(b,l)%gc(2)
          if (block(b,l)%Fj(i,j,k)%BC%tp==bc_in1) then
            global%Nin1 = max(global%Nin1,block(b,l)%Fj(i,j,k)%BC%inf)
          endif
        enddo
      enddo
    enddo
    do k=0-block(b,l)%gc(5),block(b,l)%Nk+block(b,l)%gc(6)
      do j=1-block(b,l)%gc(3),block(b,l)%Nj+block(b,l)%gc(4)
        do i=1-block(b,l)%gc(1),block(b,l)%Ni+block(b,l)%gc(2)
          if (block(b,l)%Fk(i,j,k)%BC%tp==bc_in1) then
            global%Nin1 = max(global%Nin1,block(b,l)%Fk(i,j,k)%BC%inf)
          endif
        enddo
      enddo
    enddo
  enddo ; enddo
  if (global%Nin1>0) then
    if (global%myrank==0) then
      write(stdout,'(A)')'rank0----------------------------------------------------------------------'
      write(stdout,'(A)')'rank0 There are Nin1='//trim(str(.true.,global%Nin1))//' "inflow 1"-type boundary conditions'
      write(stdout,'(A)')'rank0----------------------------------------------------------------------'
      write(stdout,*)
    endif
    call global%alloc_bc
    do b=1,global%Nin1
      if (global%myrank==0) then
        write(stdout,'(A)')'rank0----------------------------------------------------------------------'
        write(stdout,'(A)')'rank0 Loading file "'//trim(global%file%Path_InPut)//'in1.'//trim(strz(3,b))//'.bco"'
        write(stdout,'(A)')'rank0----------------------------------------------------------------------'
        write(stdout,*)
      endif
      err = global%load_bc_in1(filename=trim(global%file%Path_InPut)//'in1.'//trim(strz(3,b))//'.bco',in1=b)
    enddo
  endif

  ! initializing Runge-Kutta coefficients
  call rk_init(global%rk_ord)

  select case(global%sp_ord)
  case(1_I1P) ! 1st order piecewise constant reconstruction
    global%gco = 1_I1P
  case(3_I1P) ! 3rd order WENO reconstruction
    global%gco = 2_I1P
  case(5_I1P) ! 5th order WENO reconstruction
    global%gco = 3_I1P
  case(7_I1P) ! 7th order WENO reconstruction
    global%gco = 4_I1P
  endselect

  ! coping input files in output path
  err = make_dir(global%myrank,global%file%Path_OutPut) ! creating the output directory
  if (global%myrank==0) then
    err = copy_file(source_file = trim(global%file%Path_InPut)//File_Option, &
                    target_file = trim(global%file%Path_OutPut)//File_Option)
    err = copy_file(source_file = trim(global%file%Path_InPut)//global%file%File_Solver, &
                    target_file = trim(global%file%Path_OutPut)//global%file%File_Solver)
    err = copy_file(source_file = trim(global%file%Path_InPut)//'procmap.dat', &
                    target_file = trim(global%file%Path_OutPut)//'procmap.dat')
  endif
  do l=1,global%Nl ; do b=1,global%Nb
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
  enddo ; enddo

  ! print some informations of the initial data loaded
  !err = fluid_print_info(myrank)

  ! computing the mesh variables that are not loaded from input files
  do l=1,global%Nl
    do b=1,global%Nb
      call block(b,l)%metrics
      call block(b,l)%metrics_correction
      ! print some informations of the mesh data loaded
      err = block(b,l)%print_info_mesh(blk=b,grl=l)
    enddo
  enddo

  ! initializing WENO coefficients
  call weno_init(block=block(:,1),S=global%gco)

  ! initializing the log file of residuals
  if (global%myrank==0) then
    Ncstr = adjustl(trim(str(.true.,global%Nc)))
    ! creating the gnuplot script file for visualizing the residuals log file
    open(unit=Get_Unit(global%file%unit_res),file=trim(global%file%Path_OutPut)//'gplot_res')
    write(global%file%unit_res,'(A)')'set xlabel "Iteration"'
    write(global%file%unit_res,'(A)')'set ylabel "Residuals"'
    write(global%file%unit_res,'(A)')'set log y'
    write(global%file%unit_res,'(A)')"p 'residuals.log' u 1:3 w l title 'R1', "//char(92)
    do c=2,global%Nc-1
      global%file%varform_res = '(A,I'//trim(str(.true.,digit(2+c)))//',A,I'//trim(str(.true.,digit(c)))//',A)'
      write(global%file%unit_res,trim(global%file%varform_res)) "  'residuals.log' u 1:",2+c," w l title 'R",c,"', "//char(92)
    enddo
    global%file%varform_res = '(A,I'//trim(str(.true.,digit(2+global%Nc)))//',A,I'//trim(str(.true.,digit(global%Nc)))//',A)'
    write(global%file%unit_res,trim(global%file%varform_res)) "  'residuals.log' u 1:",2+global%Nc," w l title 'R",global%Nc,"'"
    close(global%file%unit_res)
    ! initialize gnuplot log file of residuals
    if (global%n>0) then
      open(unit=Get_Unit(global%file%unit_res),file=trim(global%file%Path_OutPut)//'residuals.log',position='APPEND')
    else
      open(unit=Get_Unit(global%file%unit_res),file=trim(global%file%Path_OutPut)//'residuals.log')
    endif
    ! initialize header
    global%file%varform_res = '('//trim(str(.true.,global%Nc))//&
      '(A,I'//trim(str(.true.,digit(global%Nc)))//'.'//trim(str(.true.,digit(global%Nc)))//',A))'
    write(varname_res,trim(global%file%varform_res))('"R',c,'",',c=1,global%Nc)
    varname_res = varname_res(1:len_trim(varname_res)-1)
    write(global%file%unit_res,'(A)')'# L2 norm of residuals'
    write(global%file%unit_res,'(A)')'# Simulation started on '//date
    write(global%file%unit_res,'(A)')'# "n","t",'//adjustl(trim(varname_res))
    ! initialize output format
    global%file%varform_res ='('//FI8P//',1X,'//adjustl(trim(str(.true.,global%Nc+1)))//'('//FR_P//',1X))'
  endif

  ! initialize the multi-processes send/recive comunications and doing the first comunication if necessary
#ifdef MPI2
  call Init_sendrecv(block=block)
#endif

  ! allocate multigrid variables if necessary
  !if (Nl>1) then
    !call alloc_multigrid(Nc = Nc, bb1 = bb(1,1,2), bb2 = bb(2,Nb,Nl))
  !endif

  ! checking the presence of probes file
  inquire(file=trim(global%file%Path_InPut)//'probes.dat',exist=is_file,iostat=err)
  if (is_file) then
     open(unit = Get_Unit(UnitFree), file = trim(global%file%Path_InPut)//'probes.dat', action = 'READ')
    read(UnitFree,*)Nprb
    if (allocated(probes)) deallocate(probes) ; allocate(probes(1:Nprb))
    do b=1,Nprb
      read(UnitFree,*)probes(b)%b,probes(b)%i,probes(b)%j,probes(b)%k
    enddo
    close(UnitFree)
  endif

#ifdef PROFILING
  ! code profiling initialization
  call profile(Np=7,header=['ZONE T="solve_grl"             ', &
                            'ZONE T="conservative2primitive"', &
                            'ZONE T="boundary_conditions"   ', &
                            'ZONE T="compute_time"          ', &
                            'ZONE T="residuals"             ', &
                            'ZONE T="rk_time_integration"   ', &
                            'ZONE T="OFF"                   '],myrank=global%myrank)
#else
  instant0 = Crono(start=.true.)
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine off_init

  !> @brief Function for loading global file options.
  function load_off_option_file(filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),      intent(IN)::    filename !< Name of file where option variables are saved.
  type(Type_Global), intent(INOUT):: global   !< Global-level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of file.
  character(3)::                     os_type  !< Type operating system.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) then
      call File_Not_Found(global%myrank,filename,'load_off_option_file')
  endif
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err) !               OS
  read(UnitFree,*,iostat=err)os_type
  read(UnitFree,*,iostat=err) !               INPUT OPTIONS
  read(UnitFree,*,iostat=err)global%file%Path_InPut
  read(UnitFree,*,iostat=err)global%file%File_Solver
  read(UnitFree,*,iostat=err)global%Nl
  read(UnitFree,*,iostat=err)global%file%File_Mesh
  read(UnitFree,*,iostat=err)global%file%File_BC
  read(UnitFree,*,iostat=err)global%file%File_Init
  read(UnitFree,*,iostat=err) !               OUTPUT OPTIONS
  read(UnitFree,*,iostat=err)global%file%Path_OutPut
  read(UnitFree,*,iostat=err)global%file%File_Sol
  read(UnitFree,*,iostat=err)global%file%screen_out
  read(UnitFree,*,iostat=err)global%file%sol_out
  read(UnitFree,*,iostat=err)global%file%restart_out
  read(UnitFree,*,iostat=err)global%file%probe_out
  close(UnitFree)
  os_type = Upper_Case(os_type) ; call OS%init(c_id=os_type,myrank=global%myrank)

  global%file%Path_InPut =string_OS_sep(global%file%Path_InPut ) ; global%file%Path_InPut  = adjustl(trim(global%file%Path_InPut ))
  global%file%File_Solver=string_OS_sep(global%file%File_Solver) ; global%file%File_Solver = adjustl(trim(global%file%File_Solver))
  global%file%File_Mesh  =string_OS_sep(global%file%File_Mesh  ) ; global%file%File_Mesh   = adjustl(trim(global%file%File_Mesh  ))
  global%file%File_BC    =string_OS_sep(global%file%File_BC    ) ; global%file%File_BC     = adjustl(trim(global%file%File_BC    ))
  global%file%File_Init  =string_OS_sep(global%file%File_Init  ) ; global%file%File_Init   = adjustl(trim(global%file%File_Init  ))
  global%file%Path_OutPut=string_OS_sep(global%file%Path_OutPut) ; global%file%Path_OutPut = adjustl(trim(global%file%Path_OutPut))
  global%file%File_Sol   =string_OS_sep(global%file%File_Sol   ) ; global%file%File_Sol    = adjustl(trim(global%file%File_Sol   ))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_off_option_file
  !> @}
endprogram OFF
