!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_GlobalDerivedType Data_Type_Global
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_GlobalInterface Data_Type_Global
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_GlobalPublicProcedure Data_Type_Global
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_GlobalPrivateProcedure Data_Type_Global
!> @}

!> @brief Module Data_Type_Global contains the definition of Type_Global and useful procedures for its handling.
!> Global-level data are referred to those informations of global interest.
module Data_Type_Global
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision           !< Integers and reals precision definition.
USE Data_Type_Primitive    !< Definition of Type_Primitive.
USE Lib_IO_Misc            !< Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: file_name
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Data_Type_GlobalDerivedType
!> @{
!> @brief Derived type containing the global-level file data.
!> Global-level file data are referred to those informations concerning with files of global interest.
type, public:: Type_File
  character(60)::  Path_InPut  = ''     !< Path of input files.
  character(60)::  Path_OutPut = ''     !< Path of output files.
  character(60)::  File_BC     ='unset' !< Base name of boundary conditions file.
  character(60)::  File_Init   ='unset' !< Base name of initial conditions file.
  character(60)::  File_Mesh   ='unset' !< Base name of mesh file.
  character(60)::  File_Spec   ='unset' !< Name of initial species file.
  character(60)::  File_Sol    ='unset' !< Base name of solution file.
  character(4)::   Sol_Ext     ='null'  !< Solution file extension.
  character(60)::  File_Solver ='unset' !< Name of solver options file.
  character(60)::  File_Pout   ='unset' !< Post-processed output file name.
  character(60)::  File_Prof   ='unset' !< Prefix name of profiling file names.
  character(500):: varform_res = ''     !< Gnuplot residuals writing format.
  integer(I8P)::   screen_out  = 1      !< Console refresh frequency.
  integer(I8P)::   sol_out     = 0      !< Actual solution writing frequency (if 0 only restart solution is saved).
  integer(I8P)::   restart_out = 1      !< Restart Solution writing frequency.
  integer(I8P)::   probe_out   = 1      !< Probes writing frequency.
  integer(I_P)::   unit_res    = 10     !< Logical unit of gnuplot log file of residuals.
endtype Type_File
!> Derived type containing global-level non-dimensional numbers and reference values.
type, public:: Type_Adimensional
  ! non dimensional numbers loaded from input file
  real(R_P):: Re = 1._R_P !< \f$\rm{Re}=\frac{\rho_0 v_0 L_0}{\mu_0}\f$ Reynolds number.
  real(R_P):: Fr = 1._R_P !< \f$\rm{Fr}=\sqrt{\frac{v_0^2}{f_0 L_0}}\f$ Froude number.
  real(R_P):: Pr = 1._R_P !< \f$\rm{Pr}=\frac{\mu_0 c_p}{k_0}\f$ Prandtl number.
  ! reference values loaded from input file
  real(R_P):: L0 = 1._R_P !< Reference length.
  real(R_P):: r0 = 1._R_P !< Reference density.
  real(R_P):: v0 = 1._R_P !< Reference velocity.
  real(R_P):: c0 = 1._R_P !< Reference specific heats (\f$cp_0 = cv_0 = R_0 = c_0\f$).
  ! reference values computed by means of previous values
  real(R_P):: mu0  = 1._R_P !< \f$\mu_0= \frac{\rho_0 v_0 L_0}{\rm{Re}}\f$ Reference dynamic viscosity.
  real(R_P):: f0   = 1._R_P !< \f$f_0= \frac{v_0^2}{L_0 \rm{Fr}^2}\f$ Reference specific force.
  real(R_P):: k0   = 1._R_P !< \f$k_0= \frac{\mu_0 c_0}{\rm{Pr}}\f$ Reference thermal conductivity coefficient.
  real(R_P):: Dt0  = 1._R_P !< \f$Dt_0=\frac{L_0}{v_0}\f$ Reference time interval.
  real(R_P):: p0   = 1._R_P !< \f$p_0=\rho_0 v_0^2\f$ Reference pressure.
  real(R_P):: a0   = 1._R_P !< \f$a_0=v_0\f$ Reference speed of sound.
  real(R_P):: T0   = 1._R_P !< \f$T_0=\frac{v_0^2}{c_0}\f$ Reference temperature.
  real(R_P):: E0   = 1._R_P !< \f$E_0=v_0^2\f$ Reference specific energy.
  real(R_P):: q0   = 1._R_P !< \f$q_0=\frac{v_0^3}{L_0}\f$ Reference specific heat.
  ! equations coefficients computed by means of previous values
  real(R_P):: Re_inv   = 1._R_P !< \f$\frac{1}{\rm{Re}}\f$ Inverse of Reynolds number (coefficient of viscous terms).
  real(R_P):: Fr2_inv  = 1._R_P !< \f$\frac{1}{\rm{Fr}^2}\f$ Inverse of square of Froude number (coefficient of volume forces).
  real(R_P):: PrRe_inv = 1._R_P !< \f$\frac{1}{\rm{Pr Re}}\f$ Inverse of Prandtl and Reynolds numbers (coef. of condution terms).
endtype Type_Adimensional
!> @brief Derived type containing the global-level data.
type, public:: Type_Global
  integer(I_P):: myrank = 0_I_P !< Rank of the process which data belongs to.
  type(Type_File):: file !< File data.
  ! mesh data
  integer(I_P):: Nl     = 1_I_P !< Number of grid levels.
  integer(I_P):: Nb     = 0_I_P !< Number of blocks.
  integer(I_P):: Nb_tot = 0_I_P !< Number of total blocks (sum over each process).
  integer(I1P):: gco    = 1_I1P !< Number of ghost cells necessary to achieve the space reconstruction order.
  ! boundary conditions data
  integer(I_P)::                      Nin1 = 0_I_P !< Number of inflow 1 boundary conditions.
  type(Type_Primitive), allocatable:: in1(:)       !< Inflow 1 boundary conditions primitive variables [1:Nin1].
  ! fluid dynamic data
  integer(I8P)::            n             = 0_I8P    !< Time steps counter.
  real(R_P)::               t             = 0._R_P   !< Time.
  integer(I_P)::            Ns            = 1_I_P    !< Number of species.
  integer(I_P)::            Np            = 7_I_P    !< Number of primitive variables    (Np = Ns + 6).
  integer(I_P)::            Nc            = 5_I_P    !< Number of conservative variables (Nc = Ns + 4).
  logical::                 inviscid      = .true.   !< Type of simulation: inviscid (Euler's eq.) or viscous (Navier-Stokes eq.).
  logical::                 unsteady      = .true.   !< Type of simulation: unsteady or not.
  integer(I8P)::            Nmax          = 0_I8P    !< Max number of iterations.
  real(R_P)::               Tmax          = 0._R_P   !< Max time, ignored if Nmax>0.
  integer(I1P)::            sp_ord        = 1_I1P    !< Order of space convergence (number of ghost cells).
  integer(I1P)::            rk_ord        = 1_I1P    !< Order of time convergence (number of Runge-Kutta stages).
  real(R_P)::               CFL           = 0.3_R_P  !< Value of stability coefficient.
  real(R_P)::               residual_toll = 0.01_R_P !< Tolerance for residuals vanishing evaluation.
  logical::                 residual_stop = .false.  !< Sentinel for stopping steady simulation when residuals vanish.
  real(R_P), allocatable::  cp0(:)                   !< Initial specific heat cp for each specie [1:Ns].
  real(R_P), allocatable::  cv0(:)                   !< Initial specific heat cv for each specie [1:Ns].
  type(Type_Adimensional):: adim                     !< Non-dimensionalization data.
  contains
    procedure, non_overridable:: alloc_bc => alloc_gbc                       ! Procedure for allocating bc memory.
    procedure, non_overridable:: load_bc_in1 => load_gbc_in1                 ! Procedure for loading the inflow1 bc data.
    procedure, non_overridable:: load_fluid_soption => load_gfluid_soption   ! Procedure for loading the fluidynamic solver options.
    procedure, non_overridable:: load_fluid_Ns => load_gfluid_Ns             ! Procedure for loading the number of species.
    procedure, non_overridable:: load_fluid_0species => load_gfluid_0species ! Procedure for loading the initial species.
    procedure, non_overridable:: alloc_fluid => alloc_gfluid                 ! Procedure for allocating the fluidynamic data.
endtype Type_Global
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Generic procedure for building file names.
!> Three different types of output file name can be build by means of this generic interface: \n
!> - File name of \b block-clean data; this type of file name is referred to those files containing block-level data; calling
!>   signature is of the type:
!>   @code fname = file_name(basename,suffix,blk,grl) @endcode
!> - File name of \b block-time_step data; this type of file name is referred to those files containing block-level data that
!>   varying with time; calling signature is of the type:
!>   @code fname = file_name(basename,suffix,blk,grl,n) @endcode
!> - File name of \b block-flip_flop data; this type of file name is referred to those files containing block-level data stored as
!>   backup flip/flop file; calling signature is of the type:
!>   @code fname = file_name(basename,suffix,blk,grl,flip) @endcode
!> @ingroup Data_Type_GlobalInterface
interface file_name
  module procedure Block_File_Name,Block_Step_File_Name,Block_Flip_File_Name
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_GlobalPrivateProcedure
  !> @{
  !> Function for building block-clean file name.
  !> @return \b filename character(len_trim(basename)+4+5+len_trim(suffix)) variable.
  function Block_File_Name(basename,suffix,blk,grl) result(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN)::                           basename !< Base name of file name.
  character(*), intent(IN)::                           suffix   !< Suffix    of file name.
  integer(I_P), intent(IN)::                           blk      !< Block number.
  integer(I_P), intent(IN)::                           grl      !< Grid refinement level.
  character(len_trim(basename)+4+5+len_trim(suffix)):: filename !< Output file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  filename = trim(basename)//'.g'//trim(strz(2,grl))//'.b'//trim(strz(3,blk))//trim(suffix)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Block_File_Name

  !> Function for building block-time_step file name.
  !> @return \b filename character(len_trim(basename)+4+5+13+len_trim(suffix)) variable.
  function Block_Step_File_Name(basename,suffix,blk,grl,n) result(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN)::                              basename !< Base name of file name.
  character(*), intent(IN)::                              suffix   !< Suffix   of file name.
  integer(I_P), intent(IN)::                              blk      !< Block number.
  integer(I_P), intent(IN)::                              grl      !< Grid refinement level.
  integer(I8P), intent(IN)::                              n        !< Time step number.
  character(len_trim(basename)+4+5+13+len_trim(suffix)):: filename !< Output file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  filename = trim(basename)//'.g'//trim(strz(2,grl))//'.b'//trim(strz(3,blk))//'-N_'//trim(strz(10,n))//trim(suffix)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Block_Step_File_Name

  !> Function for building block-flip_flop file name.
  !> @return \b filename character(len_trim(basename)+3+4+5+len_trim(suffix)) variable.
  function Block_Flip_File_Name(basename,suffix,blk,grl,flip) result(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN)::                             basename !< Base name of file name.
  character(*), intent(IN)::                             suffix   !< Suffix   of file name.
  integer(I_P), intent(IN)::                             blk      !< Block number.
  integer(I_P), intent(IN)::                             grl      !< Grid refinement level.
  integer(I1P), intent(IN)::                             flip     !< Flip-flop number.
  character(len_trim(basename)+3+4+5+len_trim(suffix)):: filename !< Output file name.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  filename = trim(basename)//'-f'//trim(strz(1,flip))//'.g'//trim(strz(2,grl))//'.b'//trim(strz(3,blk))//trim(suffix)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Block_Flip_File_Name

  !> Subroutine for computing the reference values for non-dimensional quantities.
  subroutine compute_values0(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Global), intent(INOUT):: global !< Global level data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! reference values
  global%adim%mu0  = (global%adim%r0*global%adim%v0*global%adim%L0)/global%adim%Re
  global%adim%f0   = (global%adim%v0*global%adim%v0)/(global%adim%L0*global%adim%Fr*global%adim%Fr)
  global%adim%k0   = (global%adim%mu0*global%adim%c0)/global%adim%Pr
  global%adim%Dt0  = global%adim%L0/global%adim%v0
  global%adim%p0   = global%adim%r0*global%adim%v0*global%adim%v0
  global%adim%a0   = global%adim%v0
  global%adim%T0   = (global%adim%v0*global%adim%v0)/global%adim%c0
  global%adim%E0   = global%adim%v0*global%adim%v0
  global%adim%q0   = (global%adim%v0*global%adim%v0*global%adim%v0)/global%adim%L0
  ! equations coefficients
  global%adim%Re_inv   = 1._R_P/global%adim%Re
  global%adim%Fr2_inv  = 1._R_P/(global%adim%Fr*global%adim%Fr)
  global%adim%PrRe_inv = 1._R_P/(global%adim%Pr*global%adim%Re)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_values0

  !> Subroutine for allocating dynamic data of Type_Global boundary conditions variables.
  subroutine alloc_gbc(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(global%in1)) then
    !call global%in1%free ; deallocate(global%in1)
  endif
  if (global%Nin1>0) allocate(global%in1(1:global%Nin1)) ; call global%in1%init(Ns=global%Ns)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_gbc

  !> Subroutine for allocating dynamic data of Type_Global fluid dynamic variables.
  subroutine alloc_gfluid(global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global !< Global data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(global%cp0)) deallocate(global%cp0) ; allocate(global%cp0(1:global%Ns))
  if (allocated(global%cv0)) deallocate(global%cv0) ; allocate(global%cv0(1:global%Ns))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_gfluid

  !> Function for loading inflow 1 boundary conditions.
  !> @return \b err integer(I4P) variable.
  function load_gbc_in1(global,filename,in1) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global   !< Global level data.
  character(*),       intent(IN)::    filename !< Name of file where boundary conditions variables are saved.
  integer(I_P),       intent(IN)::    in1      !< Actual in1 index.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of boundary conditions file.
  character(DI_P)::                   rks      !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(global%myrank,filename,'load_gbc_in1')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  err = read_primitive(scalar=global%in1(in1),format='*',unit=UnitFree)
  close(UnitFree)
  rks = 'rank'//trim(str(.true.,global%myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Inflow 1 primitive variables loaded:'
  err = global%in1(in1)%pprint(myrank=global%myrank,unit=stdout)
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,*,iostat=err)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gbc_in1

  !> Function for loading solver option variables.
  !> @return \b err integer(I4P) variable.
  function load_gfluid_soption(global,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global   !< Global level data.
  character(*),       intent(IN)::    filename !< Name of file where option variables are saved.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of option file.
  character(8)::                      timing   !< timing = 'unsteady' simulation otherwise steady one.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 1
  inquire(file=adjustl(trim(filename)),exist=is_file)
  if (.NOT.is_file) call File_Not_Found(global%myrank,filename,'load_gfluid_soption')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err)timing ; timing = Upper_Case(timing) ; global%unsteady = (adjustl(trim(timing))=='UNSTEADY')
  read(UnitFree,*,iostat=err)global%Nmax
  read(UnitFree,*,iostat=err)global%Tmax
  read(UnitFree,*,iostat=err)global%sp_ord
  read(UnitFree,*,iostat=err)global%rk_ord
  read(UnitFree,*,iostat=err)global%CFL
  read(UnitFree,*,iostat=err)global%residual_toll
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err)global%adim%Re
  read(UnitFree,*,iostat=err)global%adim%Fr
  read(UnitFree,*,iostat=err)global%adim%Pr
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err)global%adim%L0
  read(UnitFree,*,iostat=err)global%adim%r0
  read(UnitFree,*,iostat=err)global%adim%v0
  read(UnitFree,*,iostat=err)global%adim%c0
  close(UnitFree)
  if (global%unsteady) then
    if (global%Nmax>0_I8P) then
      global%Tmax=-1._R_P ! the value of Tmax is ignored because Nmax>0
    endif
  else
    global%Tmax=-1._R_P ! the value of Tmax is ignored because steady simulation
  endif
  call compute_values0(global) ! computing the reference values for non-dimensional quantities
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_soption

  !> Function for loading the number of species from fluid dynamic file "filename".
  !> @return \b err integer(I4P) variable.
  function load_gfluid_Ns(global,binary,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT)::        global   !< Global level data.
  logical,            intent(IN), optional:: binary   !< Flag for binary of ascii input file.
  character(*),       intent(IN)::           filename !< Name of file where mesh variables are saved.
  integer(I_P)::                             err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                             UnitFree !< Free logic unit.
  logical::                                  is_file  !< Flag for inquiring the presence of fluid dynamic file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(global%myrank,filename,'load_gfluid_Ns')
  if (present(binary)) then
    open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)global%Ns
  else
    open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err)global%Ns
  endif
  close(UnitFree)
  global%Np = global%Ns + 6
  global%Nc = global%Ns + 4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_Ns

  !> Function for loading initial species from ascii file.
  !> @return \b err integer(I4P) variable.
  function load_gfluid_0species(global,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global   !< Global level data.
  character(*),       intent(IN)::    filename !< Name of file where initial species are saved.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of initial species file.
  integer(I_P)::                      s        !< Species counter.
  real(R_P)::                         g        !< Specific heats ratio (cp/cv).
  real(R_P)::                         R        !< Specific heats difference (cp-cv).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(global%myrank,filename,'load_gfluid_0species')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) ! record skipped
  read(UnitFree,*,iostat=err) ! record skipped
  read(UnitFree,*,iostat=err)global%Ns ; call global%alloc_fluid
  do s=1,global%Ns
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err)g
    read(UnitFree,*,iostat=err)R
    read(UnitFree,*,iostat=err)global%cp0(s)
    read(UnitFree,*,iostat=err)global%cv0(s)
    if ((g>0.0_R_P).AND.(R>0.0_R_P)) then
      global%cv0(s)=R/(g-1._R_P)
      global%cp0(s)=g*global%cv0(s)
    elseif ((g>0.0_R_P).AND.(global%cp0(s)>0.0_R_P)) then
      global%cv0(s)=global%cp0(s)/g
    elseif ((g>0.0_R_P).AND.(global%cv0(s)>0.0_R_P)) then
      global%cp0(s)=global%cv0(s)*g
    elseif ((R>0.0_R_P).AND.(global%cp0(s)>0.0_R_P)) then
      global%cv0(s)=global%cp0(s)-R
    elseif ((R>0.0_R_P).AND.(global%cv0(s)>0.0_R_P)) then
      global%cp0(s)=global%cv0(s)+R
    endif
  enddo
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_0species
  !> @}
endmodule Data_Type_Global
