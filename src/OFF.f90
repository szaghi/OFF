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
USE IR_Precision                                                         ! Integers and reals precision definition.
USE Data_Type_Command_Line_Interface, only: Type_Command_Line_Interface  ! Definition of Type_Command_Line_Interface.
USE Data_Type_Files,                  only: Type_Files                   ! Definition of Type_Files.
USE Data_Type_Global,                 only: Type_Global                  ! Definition of Type_Global.
USE Data_Type_Hash_Table,             only: Type_Hash_Table              ! Definition of Type_Hash_Table.
USE Data_Type_Time,                   only: Type_Time                    ! Definition of Type_Time.
USE Data_Type_Tree,                   only: Type_Tree                    ! Definition of Type_Tree.
USE Data_Type_SBlock,                 only: Type_SBlock                  ! Definition of Type_SBlock.
USE Lib_IO_Misc                                                          ! Procedures for IO and strings operations.
USE Lib_Fluxes_Convective,            only: set_interface_reconstruction ! Procedure for initializing reconstruction algorithm.
USE Lib_Riemann_Solvers,              only: set_riemann_solver           ! Procedure for initializing Riemann solver.
USE Lib_Runge_Kutta,                  only: rk_init                      ! Procedure for initializing Runge-Kutta coefficients.
USE Lib_WENO,                         only: weno_init,weno_print         ! Procedure for initializing WENO coefficients.
#ifdef OPENMP
USE OMP_LIB                                                              ! OpenMP runtime library.
#endif
#ifdef _MPI
USE MPI                                                                  ! MPI runtime library.
USE Lib_Parallel,                     only: init_MPI_maps,print_MPI_maps ! Library for send/receive data for parallel MPI.
#endif
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
!> @ingroup OFFPrivateVarPar
!> @{
type(Type_Files)::  IOFile        !< Input/Output files.
type(Type_Global):: global        !< Global data.
type(Type_Time)::   time,time_res !< Code timing: current-elapsed and residual time.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! initializing the simulation
call off_init
#ifdef PROFILING
call IOFile%prof%profile(p=7,pstart=.true.,&
                         myrank=global%parallel%myrank,Nthreads=global%parallel%Nthreads,Nproc=global%parallel%Nproc)
#endif

Temporal_Loop: do
  ! computing the solution for the actual time step
#ifdef PROFILING
  call IOFile%prof%profile(p=1,pstart=.true.,&
                           myrank=global%parallel%myrank,Nthreads=global%parallel%Nthreads,Nproc=global%parallel%Nproc)
#endif
  call global%solve(prof=IOFile%prof)
#ifdef PROFILING
  call IOFile%prof%profile(p=1,pstop=.true.,&
                           myrank=global%parallel%myrank,Nthreads=global%parallel%Nthreads,Nproc=global%parallel%Nproc)
#endif
  ! saving the actual solution
  if (IOFile%sol%is_to_save(global%time_step%n).or.global%time_step%is_to_save()) then
    call IOFile%sol%save(global=global,n=global%time_step%n)
    if (IOFile%sol%iostat/=0) then
      write(stderr,'(A)')'+-'//global%parallel%rks//'-> '//IOFile%sol%iomsg
      call off_stop
    endif
  endif
  ! updating shell output
  if (global%parallel%myrank==0) then
    if (mod(global%time_step%n,IOFile%off_opts%shl_out)==0.or.global%time_step%is_the_end()) then
      call time%chronos
      time_res%Seconds = 100._R8P/(global%time_step%progress()) - time%Seconds
      call time%sec2dhms(seconds=time%Seconds)
      call time_res%sec2dhms(seconds=time_res%Seconds)
      associate(rks=>global%parallel%rks)
        write(stdout,'(A)')'+-'//rks//'-> Simulation progress   p:'//trim(str('(F6.2)',global%time_step%progress()))//'%'
        write(stdout,'(A)')'|-'//rks//'-> Time step number      n:'//trim(str(.true.,global%time_step%n))
        write(stdout,'(A)')'|-'//rks//'-> Time                  t:'//trim(str(.true.,global%time_step%t))
        call time%print(    pref='|-'//rks//'-> Elapsed  time',unit=stdout)
        call time_res%print(pref='|-'//rks//'-> Residual time',unit=stdout)
      endassociate
    endif
  endif
  ! control sentinel for the temporal Loop
  if (global%time_step%is_the_end()) exit Temporal_Loop
enddo Temporal_Loop
#ifdef PROFILING
call IOFile%prof%profile(p=7,pstop=.true.,&
                         myrank=global%parallel%myrank,Nthreads=global%parallel%Nthreads,Nproc=global%parallel%Nproc)
#endif
! finalizing the simulation
call off_finalize
! stopping the code
call off_stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup OFFPrivateProcedure
  !> @{
  !> @brief Procedure for parsing Command Line Arguments (CLA) implementing OFF Command Line Interface (CLI).
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
    call cli%init(progname='OFF',examples=['OFF OFF_options_file'])
    ! setting CLAs
    call cli%add(pref='|-'//rks//'-> ',positional=.true.,position=1,help='OFF options file',required=.true.,error=error)
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

  !> @brief Procedure for initializing the simulation.
  subroutine off_init()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P)::               Nca = 0_I4P    !< Number of command line arguments.
  character(99)::              File_Option    !< Options file name.
  integer(I4P)::               err            !< Error traping flag.
  real(R8P)::                  min_space_step !< Minimum space step.
  type(Type_SBlock), pointer:: block          !< Pointer for scanning global%block tree.
  integer(I8P)::               ID             !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! initializing IR_Precision module constants
  call IR_init
  ! initialzing timing
  call time%chronos(start=.true.)
  ! initializing parallel environments
  call global%parallel%init
  associate(rks=>global%parallel%rks)
    ! checking/creating lockfile
    if (global%parallel%is_master()) then
      call IOFile%lockfile%lock(errpref='+-'//rks//'->',outpref='+-'//rks//'->')
      if (IOFile%lockfile%iostat/=0) then
        write(stderr,'(A)')'+-'//rks//'->'//IOFile%lockfile%iomsg
        call off_stop
      endif
    endif
    ! parsing command line for getting global option file name
    err = parsing_command_line(File_Option=File_Option) ; if (err/=0) call off_stop
    write(stdout,'(A)')'+-'//rks//'-> Compiled code used options'
    call global%cco%print(unit=stdout,pref='|-'//rks//'->')
    ! printing architecture informations
    if (global%parallel%is_master()) then
      write(stdout,'(A)')'+-'//rks//'-> Running architecture general informations'
      call global%OS%print(pref='|-'//rks//'->  ',unit=stdout)
      write(stdout,'(A)')'+-'//rks//'-> Precision of running architecture'
      call IR_Print(pref='|-'//rks//'->  ',unit=stdout)
      write(stdout,'(A)')'+-'//rks//'-> Number of MPI processes:  '//trim(str(.true.,global%parallel%Nproc))
      write(stdout,'(A)')'+-'//rks//'-> Number of OpenMP threads: '//trim(str(.true.,global%parallel%Nthreads))
    endif
    ! setting OFF options file structures
    call IOFile%off_opts%set(name=trim(File_Option),path_in='',errpref='+-'//rks//'->',outpref='+-'//rks//'->')
    ! loading OFF options file
    associate(OS=>global%OS,off_opts=>IOFile%off_opts)
      call off_opts%load(OS=OS)
      if (global%parallel%is_master()) then
        write(stdout,'(A)')'+-'//rks//'-> OFF options'
        call off_opts%print(pref='|-'//rks//'->  ',unit=stdout)
      endif
      err=OS%make_dir(directory=off_opts%path_out)
    endassociate
    ! setting files strucutures
    associate(off_opts=>IOFile%off_opts,solv_opts=>IOFile%solv_opts,mesh=>IOFile%mesh,bc=>IOFile%bc,init=>IOFile%init,&
              !sol=>IOFile%sol,proc=>IOFile%proc,prof=>IOFile%prof)
              sol=>IOFile%sol,prof=>IOFile%prof)
      call solv_opts%set(name=off_opts%fn_solv,path_in=off_opts%path_in,path_out=off_opts%path_out,errpref='+-'//rks//'->',&
                         final_call=off_stop)
      call mesh%set(     name=off_opts%fn_mesh,path_in=off_opts%path_in,path_out=off_opts%path_out,errpref='+-'//rks//'->',&
                         final_call=off_stop)
      call bc%set(       name=off_opts%fn_bc  ,path_in=off_opts%path_in,path_out=off_opts%path_out,errpref='+-'//rks//'->',&
                         final_call=off_stop)
      call init%set(     name=off_opts%fn_init,path_in=off_opts%path_in,path_out=off_opts%path_out,errpref='+-'//rks//'->',&
                         final_call=off_stop)
      call sol%set(      name=off_opts%fn_sol ,path_in=off_opts%path_in,path_out=off_opts%path_out,errpref='+-'//rks//'->',&
                         final_call=off_stop,fout=off_opts%sol_out)
      !call proc%set(     name=off_opts%fn_proc,path_in=off_opts%path_in,path_out=off_opts%path_out)
      call prof%set(     name=off_opts%fn_prof,path_in=off_opts%path_in,path_out=off_opts%path_out,errpref='+-'//rks//'->',&
                         final_call=off_stop)
    endassociate
    ! loading input files
    if (global%parallel%is_master()) then
      write(stdout,'(A)')'+-'//rks//'-> Loading input files'
    endif
    call IOFile%solv_opts%load
    if (global%parallel%is_master()) then
      write(stdout,'(A)')'+-'//rks//'->   Loading '//IOFile%solv_opts%name
      call IOFile%solv_opts%print(pref='|-'//rks//'->    ',unit=stdout)
    endif
    call set_riemann_solver(solver=IOFile%solv_opts%RSU)
    call set_interface_reconstruction(reconstruction_type=IOFile%solv_opts%recon_tp)
    global%time_step  = IOFile%solv_opts%time_step
    global%space_step = IOFile%solv_opts%space_step
    global%adim       = IOFile%solv_opts%adim
    ! initializing Runge-Kutta coefficients
    call rk_init(global%time_step%rk_ord)
   !! loading processes/blocks map and computing the number global/local blocks
   !if (global%parallel%myrank==0) then
   !  write(stdout,'(A)')'+-'//rks//'->   Loading '//IOFile%proc%name
   !endif
   !call IOFile%proc%load(mesh_dims=global%mesh_dims,parallel=global%parallel)
   !if (IOFile%proc%iostat/=0) then
   !  write(stderr,'(A)')'+-'//rks//'-> '//IOFile%proc%iomsg
   !  call off_stop
   !endif
   !call global%parallel%print(pref='|-'//rks//'->    ',unit=stdout)
    if (global%parallel%is_master()) write(stdout,'(A)')'+-'//rks//'->   Loading '//IOFile%mesh%name
    call IOFile%mesh%load(global=global)
    if (global%parallel%is_master()) write(stdout,'(A)')'+-'//rks//'->   Loading '//IOFile%bc%name
    call IOFile%bc%load(global=global)
    if (global%parallel%is_master()) write(stdout,'(A)')'+-'//rks//'->   Loading '//IOFile%init%name
    call IOFile%init%load(global=global)
    ! computing the mesh variables that are not loaded from input files
    min_space_step = MaxR8P
    do while(global%block%loopID(ID=ID))
      block => global%block%dat(ID=ID)
      call block%metrics
      call block%metrics_correction
      min_space_step = min(min_space_step,block%min_space_step())
      call block%print(unit=stdout,pref='|-'//rks//'->  ')
    enddo
    ! initializing WENO coefficients
    call weno_init(S=global%space_step%gco,min_space_step=min_space_step)
    if (global%parallel%is_master()) then
      write(stdout,'(A)')'+-'//rks//'->   WENO settings'
      call weno_print(unit=stdout,pref='|-'//rks//'->     ')
    endif
    ! printing block infos
    if (global%parallel%is_master()) then
      write(stdout,'(A)')'+-'//rks//'->   Blocks infos'
      do while(global%block%loopID(ID=ID))
        block => global%block%dat(ID=ID)
        write(stdout,'(A)')'+-'//rks//'->     Block ID='//trim(str(n=ID))
        call block%print(unit=stdout,pref='|-'//rks//'->      ')
      enddo
    endif
    ! coping input files in output path
#ifdef _MPI
    ! syncronizing MPI processes
    call MPI_BARRIER(MPI_COMM_WORLD,err)
#endif
    if (global%parallel%is_master()) then
      call IOFile%off_opts%backup(OS=global%OS)
      call IOFile%solv_opts%backup(OS=global%OS)
      !call IOFile%proc%backup(OS=global%OS)
      call IOFile%mesh%backup(OS=global%OS)
      call IOFile%bc%backup(OS=global%OS)
      call IOFile%init%backup(OS=global%OS)
    endif
    ! initialize the multi-processes send/recive comunications and doing the first comunication if necessary
#ifdef _MPI
    call init_MPI_maps(parallel=global%parallel,mesh_dims=global%mesh_dims,block=global%block)
    if (global%parallel%myrank==0) write(stdout,'(A)')'+-'//rks//'->   MPI send/recv maps'
    call print_MPI_maps(parallel=global%parallel,pref='|-'//rks//'->    ',unit=stdout)
#endif
#ifdef PROFILING
    ! code profiling initialization
    call IOFile%prof%profile(Np=7,                                      &
                             header=['ZONE T="solve_grl"             ', &
                                     'ZONE T="conservative2primitive"', &
                                     'ZONE T="boundary_conditions"   ', &
                                     'ZONE T="compute_time"          ', &
                                     'ZONE T="residuals"             ', &
                                     'ZONE T="rk_time_integration"   ', &
                                     'ZONE T="OFF"                   '],&
                             myrank=global%parallel%myrank,Nthreads=global%parallel%Nthreads,Nproc=global%parallel%Nproc)
#endif

 !! loading inflow boundary conditions if necessary
 !! inflow 1
 !do l=1,global%mesh_dims%Nl ; do b=1,global%Nb
 !  do k=1-block(b,l)%gc(5),block(b,l)%Nk+block(b,l)%gc(6)
 !    do j=1-block(b,l)%gc(3),block(b,l)%Nj+block(b,l)%gc(4)
 !      do i=0-block(b,l)%gc(1),block(b,l)%Ni+block(b,l)%gc(2)
 !        if (block(b,l)%Fi(i,j,k)%BC%tp==bc_in1) then
 !          global%Nin1 = max(global%Nin1,block(b,l)%Fi(i,j,k)%BC%inf)
 !        endif
 !      enddo
 !    enddo
 !  enddo
 !  do k=1-block(b,l)%gc(5),block(b,l)%Nk+block(b,l)%gc(6)
 !    do j=0-block(b,l)%gc(3),block(b,l)%Nj+block(b,l)%gc(4)
 !      do i=1-block(b,l)%gc(1),block(b,l)%Ni+block(b,l)%gc(2)
 !        if (block(b,l)%Fj(i,j,k)%BC%tp==bc_in1) then
 !          global%Nin1 = max(global%Nin1,block(b,l)%Fj(i,j,k)%BC%inf)
 !        endif
 !      enddo
 !    enddo
 !  enddo
 !  do k=0-block(b,l)%gc(5),block(b,l)%Nk+block(b,l)%gc(6)
 !    do j=1-block(b,l)%gc(3),block(b,l)%Nj+block(b,l)%gc(4)
 !      do i=1-block(b,l)%gc(1),block(b,l)%Ni+block(b,l)%gc(2)
 !        if (block(b,l)%Fk(i,j,k)%BC%tp==bc_in1) then
 !          global%Nin1 = max(global%Nin1,block(b,l)%Fk(i,j,k)%BC%inf)
 !        endif
 !      enddo
 !    enddo
 !  enddo
 !enddo ; enddo
 !if (global%Nin1>0) then
 !  if (global%myrank==0) then
 !    write(stdout,'(A)',iostat=err)trim(rks)//' There are Nin1='//trim(str(.true.,global%Nin1))//' "inflow 1"-type BC'
 !  endif
 !  call global%alloc_bc
 !  do b=1,global%Nin1
 !    if (global%myrank==0) then
 !    write(stdout,'(A)',iostat=err)trim(rks)//'   Loading file "'//trim(global%dfile%Path_InPut)//'in1.'//trim(strz(3,b))//'.bco"'
 !    endif
 !    err = global%load_bc_in1(filename=trim(global%dfile%Path_InPut)//'in1.'//trim(strz(3,b))//'.bco',in1=b)
 !    if (err/=0) call off_stop
 !  enddo
 !endif

 !! initializing the log file of residuals
 !if (global%myrank==0) then
 !  err = global%dfile%init_file_res(Nc=global%Nc,n=global%n,date=date)
 !endif

  ! checking the presence of probes file
 !inquire(file=trim(global%dfile%Path_InPut)//'probes.dat',exist=is_file,iostat=err)
 !if (is_file) then
 !   open(unit = Get_Unit(UnitFree), file = trim(global%dfile%Path_InPut)//'probes.dat', action = 'READ')
 !  read(UnitFree,*)Nprb
 !  if (allocated(probes)) deallocate(probes) ; allocate(probes(1:Nprb))
 !  do b=1,Nprb
 !    read(UnitFree,*)probes(b)%b,probes(b)%i,probes(b)%j,probes(b)%k
 !  enddo
 !  close(UnitFree)
 !endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine off_init

  !> @brief Procedure for initializing the simulation according to the input options.
  subroutine off_finalize()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
 !integer(I4P)::    b                 !< Blocks counter.
 !integer(I4P)::    l                 !< Grid levels counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! saving the final time step solution
 !do l=1,global%mesh_dims%Nl
 !  ! converting conservative variables to primitive ones
 !  do b=1,global%Nb
 !    call conservative2primitive(block = block(b,l))
 !  enddo
 !  ! imposing the boundary conditions
 !  call boundary_conditions(l = l, block = block(:,l))
 !  ! saving the output file
 !  do b=1,global%Nb
 !    err = block(b,l)%save_fluid(filename=file_name(basename=trim(global%dfile%Path_OutPut)//global%dfile%File_Sol,&
 !                                                   suffix=trim(global%dfile%Sol_Ext),blk=global%blockmap(b),grl=l,n=global%n))
 !  enddo
 !enddo
  ! the simulation is done: safe finalizing the simulation
  if (global%parallel%is_master()) then
    call IOFile%lockfile%unlock
 !  close(global%dfile%unit_res)     ! close log residuals file
  endif
#ifdef PROFILING
  call IOFile%prof%profile(finalize=.true.,&
                           myrank=global%parallel%myrank,Nthreads=global%parallel%Nthreads,Nproc=global%parallel%Nproc)
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine off_finalize

  !> @brief Procedure for stopping the code wiht safe calling to MPI finalize.
  subroutine off_stop
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
#ifdef _MPI
  integer(I4P):: err !< Error traping flag.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
#ifdef _MPI
  call MPI_FINALIZE(err)
#endif
  stop
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine off_stop
  !> @}
endprogram OFF
