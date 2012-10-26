!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_GlobalsPublicProcedure Data_Type_Globals
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_GlobalsPrivateProcedure Data_Type_Globals
!> @}

!> @brief Module Data_Type_Globals contains the definition of global types variables and useful procedures for their handling.
!> Two main derived type are defined: \n
!> - Type_Global: derived type containing global-level data.
!> - Type_Block: derived type containing block-level data.
module Data_Type_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision           !< Integers and reals precision definition.
USE Data_Type_BC           !< Definition of Type_BC.
USE Data_Type_Cell         !< Definition of Type_Cell.
USE Data_Type_Conservative !< Definition of Type_Conservative.
USE Data_Type_Primitive    !< Definition of Type_Primitive.
USE Data_Type_Vector       !< Definition of Type_Vector.
USE Lib_IO_Misc            !< Procedures for IO and strings operations.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: file_name
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the global-level file data.
!> Global-level file data are referred to those informations concerning with files of global interest.
!> @ingroup DerivedType
type, public:: Type_File_Global
  character(60)::  Path_InPut           !< Path of input files.
  character(60)::  Path_OutPut          !< Path of output files.
  character(60)::  File_BC     ='unset' !< Base name of boundary conditions file.
  character(60)::  File_Init   ='unset' !< Base name of initial conditions file.
  character(60)::  File_Mesh   ='unset' !< Base name of mesh file.
  character(60)::  File_Spec   ='unset' !< Name of initial species file.
  character(60)::  File_Sol    ='unset' !< Base name of solution file.
  character(60)::  File_Solver ='unset' !< Name of solver options file.
  character(60)::  File_Pout   ='unset' !< Post-processed output file name.
  character(500):: varform_res          !< Gnuplot residuals writing format.
  integer(I8P)::   screen_out  = 1      !< Console refresh frequency.
  integer(I8P)::   sol_out     = 0      !< Actual solution writing frequency (if 0 only restart solution is saved).
  integer(I8P)::   restart_out = 1      !< Restart Solution writing frequency.
  integer(I8P)::   probe_out   = 1      !< Probes writing frequency.
  integer(I_P)::   unit_res             !< Logical unit of gnuplot log file of residuals.
endtype Type_File_Global

!> @brief Derived type containing the global-level mesh data.
!> Global-level mesh data are referred to those informations concerning with mesh's details of global interest.
!> @ingroup DerivedType
type, public:: Type_Mesh_Global
  integer(I_P):: Nl     = 1_I_P !< Number of grid levels.
  integer(I_P):: Nb     = 0_I_P !< Number of blocks.
  integer(I_P):: Nb_tot = 0_I_P !< Number of total blocks (sum over each process).
  integer(I1P):: gco    = 1_I_P !< Number of ghost cells necessary to achieve the space reconstruction order.
endtype Type_Mesh_Global
!> @brief Derived type containing the block-level mesh data.
!> Block-level mesh data are referred to those informations concerning with mesh's details of blocks interest.
type, public:: Type_Mesh_Block
  integer(I_P)::                   Ni     = 0_I_P    !< Number of cells in i direction.
  integer(I_P)::                   Nj     = 0_I_P    !< Number of cells in j direction.
  integer(I_P)::                   Nk     = 0_I_P    !< Number of cells in k direction.
  integer(I1P)::                   gc(1:6)=&
                                          (/1_I1P, & ! gc(1) => left  i.
                                            1_I1P, & ! gc(2) => right i.
                                            1_I1P, & ! gc(3) => left  j.
                                            1_I1P, & ! gc(4) => right j.
                                            1_I1P, & ! gc(5) => left  k.
                                            1_I1P  & ! gc(6) => right k.
                                            /)       !< Number of ghost cells for the 6 faces of the block.
  type(Type_Vector), allocatable:: node(:,:,:)       !< Cell nodes coordinates [0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  type(Type_Vector), allocatable:: NFi(:,:,:)        !< Face i normals, versor [0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Vector), allocatable:: NFj(:,:,:)        !< Face j normals, versor [1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Vector), allocatable:: NFk(:,:,:)        !< Face k normals, versor [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  real(R_P),         allocatable:: Si(:,:,:)         !< Face i area            [0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  real(R_P),         allocatable:: Sj(:,:,:)         !< Face j area            [1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  real(R_P),         allocatable:: Sk(:,:,:)         !< Face k area            [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)].
  real(R_P),         allocatable:: V(:,:,:)          !< Cell volume            [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Cell),   allocatable:: cell(:,:,:)       !< Cell data informations [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  type(Type_Vector), allocatable:: cent(:,:,:)       !< Cell center coordinates[1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
  ! Adaptive Mesh Refinement
  integer(I1P)::                   octree = 0_I_P    !< Octree refinement level.
  type(Type_Mesh_Block), pointer:: amr(:,:,:)=>null()!< Mesh data of refined cells.
endtype Type_Mesh_Block

!> @brief Derived type containing the global-level boundary conditions data.
!> Global-level boundary conditions data are referred to those informations concerning with bc's details of global interest.
!> @ingroup DerivedType
type, public:: Type_BC_Global
  integer(I_P)::                      Nin1 = 0_I_P !< Number of inflow 1 boundary conditions.
  type(Type_Primitive), allocatable:: in1(:)       !< Inflow 1 boundary conditions primitive variables [1:Nin1].
endtype Type_BC_Global
!> @brief Derived type containing the block-level boundary conditions data.
!> Block-level boundary conditions data are referred to those informations concerning with bc's details of blocks interest.
!> @ingroup DerivedType
type, public:: Type_BC_Block
  type(Type_BC), allocatable:: BCi(:,:,:) !< Boundary conditions of i faces [0-gc(1):1+gc(2),1-gc(3):1+gc(4),1-gc(5):1+gc(6)].
  type(Type_BC), allocatable:: BCj(:,:,:) !< Boundary conditions of j faces [1-gc(1):1+gc(2),0-gc(3):1+gc(4),1-gc(5):1+gc(6)].
  type(Type_BC), allocatable:: BCk(:,:,:) !< Boundary conditions of k faces [1-gc(1):1+gc(2),1-gc(3):1+gc(4),0-gc(5):1+gc(6)].
endtype Type_BC_Block

!> Derived type containing the global-level fluid dynamic data.
!> Global-level fluid dynamic data are referred to those informations concerning with fluids details of global interest.
!> @ingroup DerivedType
type, public:: Type_Fluid_Global
  integer(I8P)::           n             = 0_I8P    !< Time steps counter.
  real(R_P)::              t             = 0._R_P   !< Time.
  integer(I_P)::           Ns            = 1_I_P    !< Number of species.
  integer(I_P)::           Np            = 7_I_P    !< Number of primitive variables    (Np = Ns + 6).
  integer(I_P)::           Nc            = 5_I_P    !< Number of conservative variables (Nc = Ns + 4).
  logical::                inviscid      = .true.   !< Type of simulation: inviscid (Euler's eq.) or viscous (Navier-Stokes eq.).
  logical::                unsteady      = .true.   !< Type of simulation: unsteady or not.
  integer(I8P)::           Nmax          = 0_I8P    !< Max number of iterations.
  real(R_P)::              Tmax          = 0._R_P   !< Max time, ignored if Nmax>0.
  integer(I1P)::           sp_ord        = 1_I_P    !< Order of space convergence (number of ghost cells).
  integer(I1P)::           rk_ord        = 1_I_P    !< Order of time convergence (number of Runge-Kutta stages).
  real(R_P)::              CFL           = 0.3_R_P  !< Value of stability coefficient.
  real(R_P)::              residual_toll = 0.01_R_P !< Tolerance for residuals vanishing evaluation.
  logical::                residual_stop = .false.  !< Sentinel for stopping steady simulation when residuals vanish.
  real(R_P), allocatable:: cp0(:)                   !< Initial specific heat cp for each specie [1:Ns].
  real(R_P), allocatable:: cv0(:)                   !< Initial specific heat cv for each specie [1:Ns].
endtype Type_Fluid_Global
!> Derived type containing the block-level fluid dynamic data.
!> Block-level fluid dynamic data are referred to those informations concerning with fluids details of blocks interest.
!> @note Dimensions of 3D array are [1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)].
!> @ingroup DerivedType
type, public:: Type_Fluid_Block
  real(R_P),               allocatable:: Dt(:,:,:)   !< Local time step.
  type(Type_Primitive),    allocatable:: P (:,:,:)   !< Primitive variables.
  type(Type_Conservative), allocatable:: U (:,:,:)   !< Conservative variables.
  type(Type_Conservative), allocatable:: KS(:,:,:,:) !< Runge-Kutta stages of conservative variables [1:rk_ord].
endtype Type_Fluid_Block

!> Derived type containing global-level non-dimensional numbers and reference values.
!> @ingroup DerivedType
type, public:: Type_Adimensional_Global
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
endtype Type_Adimensional_Global

!> @brief Derived type containing the global-level data.
!> @ingroup DerivedType
type, public:: Type_Global
  type(Type_File_Global)::         file  !< File data.
  type(Type_Mesh_Global)::         mesh  !< Mesh data.
  type(Type_BC_Global)::           bc    !< Boundary conditions data.
  type(Type_Fluid_Global)::        fluid !< Fluid dynamic data.
  type(Type_Adimensional_Global):: adim  !< Non-dimensionalization data.
  contains
    procedure, non_overridable:: alloc_bc => alloc_gbc                       ! Procedure for allocating bc memory.
    procedure, non_overridable:: load_bc_in1 => load_gbc_in1                 ! Procedure for loading the inflow1 bc data.
    procedure, non_overridable:: load_fluid_soption => load_gfluid_soption   ! Procedure for loading the fuidynamic solver options.
    procedure, non_overridable:: load_fluid_Ns => load_gfluid_Ns             ! Procedure for loading the number of species.
    procedure, non_overridable:: load_fluid_0species => load_gfluid_0species ! Procedure for loading the initial species.
    procedure, non_overridable:: alloc_fluid => alloc_gfluid                 ! Procedure for allocating the fluidynamic data.
endtype Type_Global
!> Derived type containing the block-level data.
!> @ingroup DerivedType
type, public:: Type_Block
  type(Type_Mesh_Block)::  mesh  !< Mesh data.
  type(Type_BC_Block)::    bc    !< Boundary conditions data.
  type(Type_Fluid_Block):: fluid !< Fluid dynamic data.
  contains
    procedure, non_overridable:: free => free_block                    ! Procedure for freeing memory.
    procedure, non_overridable:: alloc => alloc_block                  ! Procedure for allocating memory.
    procedure, non_overridable:: load_mesh_dims => load_bmesh_dims     ! Procedure for loading the mesh data dimensions.
    procedure, non_overridable:: load_mesh => load_bmesh               ! Procedure for loading the mesh data.
    procedure, non_overridable:: save_mesh => save_bmesh               ! Procedure for saving the mesh data.
    procedure, non_overridable:: print_info_mesh => print_info_bmesh   ! Procedure for printing of the mesh data.
    procedure, non_overridable:: load_bc => load_bbc                   ! Procedure for loading the bc data.
    procedure, non_overridable:: save_bc => save_bbc                   ! Procedure for saving the bc data.
    procedure, non_overridable:: load_fluid => load_bfluid             ! Procedure for loading the fluidynamic data.
    procedure, non_overridable:: save_fluid => save_bfluid             ! Procedure for saving the fluidynamic data.
    procedure, non_overridable:: print_info_fluid => print_info_bfluid ! Procedure for printing of the mesh data.
endtype Type_Block
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
!> @ingroup Interface,Data_Type_GlobalsPublicProcedure
interface file_name
  module procedure Block_File_Name,Block_Step_File_Name,Block_Flip_File_Name
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_GlobalsPrivateProcedure
  !> @{
  !> Function for building block-clean file name.
  !> @return \b filename character(len_trim(basename)+4+5+len_trim(suffix)) variable.
  function Block_File_Name(basename,suffix,blk,grl) result(filename)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN)::                           basename !< Base name of file name.
  character(*), intent(IN)::                           suffix   !< Suffix   of file name.
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
  if (allocated(global%bc%in1)) then
    call global%bc%in1%free ; deallocate(global%bc%in1)
  endif
  if (global%bc%Nin1>0) allocate(global%bc%in1(1:global%bc%Nin1)) ; call global%bc%in1%init(Ns=global%fluid%Ns)
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
  if (allocated(global%fluid%cp0)) deallocate(global%fluid%cp0) ; allocate(global%fluid%cp0(1:global%fluid%Ns))
  if (allocated(global%fluid%cv0)) deallocate(global%fluid%cv0) ; allocate(global%fluid%cv0(1:global%fluid%Ns))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_gfluid

  !> Subroutine for freeing dynamic data of Type_Block variables.
  subroutine free_block(block)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_block), intent(INOUT):: block !< Block data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Mesh data
  if (allocated(block%mesh%node)) deallocate(block%mesh%node)
  if (allocated(block%mesh%NFi )) deallocate(block%mesh%NFi )
  if (allocated(block%mesh%NFj )) deallocate(block%mesh%NFj )
  if (allocated(block%mesh%NFk )) deallocate(block%mesh%NFk )
  if (allocated(block%mesh%Si  )) deallocate(block%mesh%Si  )
  if (allocated(block%mesh%Sj  )) deallocate(block%mesh%Sj  )
  if (allocated(block%mesh%Sk  )) deallocate(block%mesh%Sk  )
  if (allocated(block%mesh%V   )) deallocate(block%mesh%V   )
  if (allocated(block%mesh%cent)) deallocate(block%mesh%cent)
  if (allocated(block%mesh%cell)) deallocate(block%mesh%cell)
  ! Boundary conditions data
  if (allocated(block%bc%BCi)) then
    call block%bc%BCi%free ; deallocate(block%bc%BCi)
  endif
  if (allocated(block%bc%BCj)) then
    call block%bc%BCj%free ; deallocate(block%bc%BCj)
  endif
  if (allocated(block%bc%BCk)) then
    call block%bc%BCk%free ; deallocate(block%bc%BCk)
  endif
  ! Fluid dynamic data
  if (allocated(block%fluid%Dt)) deallocate(block%fluid%Dt)
  if (allocated(block%fluid%P)) then
    call block%fluid%P%free ; deallocate(block%fluid%P)
  endif
  if (allocated(block%fluid%U)) then
    call block%fluid%U%free ; deallocate(block%fluid%U)
  endif
  if (allocated(block%fluid%KS)) then
    call block%fluid%KS%free ; deallocate(block%fluid%KS)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_block

  !> Subroutine for allocating dynamic data of Type_Block variables.
  subroutine alloc_block(block,global)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_block), intent(INOUT):: block    !< Block data.
  type(Type_Global), intent(IN)::    global   !< Global data.
  integer(I_P)::                     Ni,Nj,Nk !< Temporary variables for storing blocks dimensions.
  integer(I_P)::                     gc(1:6)  !< Temporary variable  for storing blocks ghost cells number.
  integer(I_P)::                     Ns       !< Temporary variable  for storing number of species.
  integer(I_P)::                     rk_ord   !< Temporary variable  for storing rk_ord.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! Free dynamic data block if previously allocated
  call block%free
  ! Storing block dimensions into temporary variables to simplify the code
  Ni      = block%mesh%Ni
  Nj      = block%mesh%Nj
  Nk      = block%mesh%Nk
  gc(1:6) = block%mesh%gc(1:6)
  Ns      = global%fluid%Ns
  rk_ord  = global%fluid%rk_ord
  ! Mesh data
  allocate(block%mesh%node(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  allocate(block%mesh%NFi (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%mesh%NFj (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%mesh%NFk (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  allocate(block%mesh%Si  (0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%Si  =0._R_P
  allocate(block%mesh%Sj  (1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%Sj  =0._R_P
  allocate(block%mesh%Sk  (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; block%mesh%Sk  =0._R_P
  allocate(block%mesh%V   (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; block%mesh%V   =0._R_P
  allocate(block%mesh%cell(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%mesh%cent(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  ! Boundary conditions data
  allocate(block%bc%BCi(0-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%bc%BCj(1-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  allocate(block%bc%BCk(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  ! Fluid dynamic data
  allocate(block%fluid%Dt(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; block%fluid%Dt = 0._R_P
  allocate(block%fluid%P (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; call block%fluid%P%init(Ns=Ns)
  allocate(block%fluid%U (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))          ; call block%fluid%U%init(Ns=Ns)
  allocate(block%fluid%KS(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6),1:rk_ord)) ; call block%fluid%KS%init(Ns=Ns)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_block

  !> Function for loading the mesh data dimensions of block from the mesh file "filename".
  !> @return \b err integer(I4P) variable.
  function load_bmesh_dims(block,ascii,myrank,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(INOUT):: block    !< Block level data.
  logical, optional, intent(IN)::    ascii    !< Flag for ascii file.
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where mesh variables are saved.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of mesh file.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_bmesh_dims')
  if (present(ascii)) then
    open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err)block%mesh%gc
    read(UnitFree,*,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree ), file = trim(filename), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)block%mesh%gc
    read(UnitFree,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bmesh_dims

  !> Function for loading block mesh file.
  !> @return \b err integer(I4P) variable.
  function load_bmesh(block,ascii,myrank,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(INOUT):: block    !< Block level data.
  logical, optional, intent(IN)::    ascii    !< Flag for ascii file.
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where mesh variables are saved.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of mesh file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_bmesh')
  if (present(ascii)) then
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err)block%mesh%gc
    read(UnitFree,*,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    err = read_vector(array3D=block%mesh%node,format='*',unit=UnitFree)
    err = read_cell(array3D=block%mesh%cell,format='*',unit=UnitFree)
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)block%mesh%gc
    read(UnitFree,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    err = read_vector(array3D=block%mesh%node,unit=UnitFree)
    err = read_cell(array3D=block%mesh%cell,unit=UnitFree)
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bmesh

  !> Function for saving block mesh file.
  !> @return \b err integer(I4P) variable.
  function save_bmesh(block,ascii,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(IN):: block    !< Block level data.
  logical, optional, intent(IN):: ascii    !< Flag for ascii file.
  character(*),      intent(IN):: filename !< Name of file where mesh variables are saved.
  integer(I_P)::                  err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                  UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (present(ascii)) then
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), action = 'WRITE', form = 'FORMATTED')
    write(UnitFree,'(6('//FI1P//',1X))',iostat=err)block%mesh%gc
    write(UnitFree,'(3('//FI_P//',1X))',iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    err = write_vector(array3D=block%mesh%node,format='*',unit=UnitFree)
    err = write_cell(array3D=block%mesh%cell,format='*',unit=UnitFree)
    close(UnitFree)
  else
    open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), action = 'WRITE', form = 'UNFORMATTED')
    write(UnitFree,iostat=err)block%mesh%gc
    write(UnitFree,iostat=err)block%mesh%Ni,block%mesh%Nj,block%mesh%Nk
    err = write_vector(array3D=block%mesh%node,unit=UnitFree)
    err = write_cell(array3D=block%mesh%cell,unit=UnitFree)
    close(UnitFree)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bmesh

  !> Function for printing to standard output info of block mesh data.
  !> @return \b err integer(I4P) variable.
  function print_info_bmesh(block,myrank,blk,grl,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(IN)::  block        !< Block level data.
  integer(I_P),      intent(IN)::  myrank       !< Rank process identification for MPI communications.
  integer(I_P),      intent(IN)::  blk          !< Actual block number.
  integer(I_P),      intent(IN)::  grl          !< Actual grid level number.
  type(Type_Global), intent(IN)::  global       !< Global level data.
  integer(I_P)::                   err          !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DR_P)::                vmax,vmin    !< String for printing max and min of variables.
  character(DR_P)::                xmm,ymm,zmm  !< String for printing max and min of variables.
  type(Type_Vector), allocatable:: NFi(:,:,:)   !< |
  type(Type_Vector), allocatable:: NFj(:,:,:)   !< |
  type(Type_Vector), allocatable:: NFk(:,:,:)   !< |
  real(R_P),         allocatable:: Si (:,:,:)   !< | Dummy variables for printing only internal cells info.
  real(R_P),         allocatable:: Sj (:,:,:)   !< |
  real(R_P),         allocatable:: Sk (:,:,:)   !< |
  real(R_P),         allocatable:: V  (:,:,:)   !< |
  integer(I_P)::                   Ni,Nj,Nk     !< Temporary variables for storing block dimensions.
  character(DI_P)::                rks          !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Mesh infos of myrank '//trim(rks)
  write(stdout,'(A)',iostat=err)trim(rks)//'  Global number of blocks '//trim(str(.true.,global%mesh%Nb_tot))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Block '//trim(str(.true.,blk))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Grid level '//trim(str(.true.,grl))
  Ni = block%mesh%Ni
  Nj = block%mesh%Nj
  Nk = block%mesh%Nk
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  i '//trim(str(.true.,block%mesh%gc(1)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right i '//trim(str(.true.,block%mesh%gc(2)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  j '//trim(str(.true.,block%mesh%gc(3)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right j '//trim(str(.true.,block%mesh%gc(4)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: left  k '//trim(str(.true.,block%mesh%gc(5)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Ghost cells: right k '//trim(str(.true.,block%mesh%gc(6)))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Ni '//trim(str(.true.,Ni))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Nj '//trim(str(.true.,Nj))
  write(stdout,'(A)',iostat=err)trim(rks)//'   Number of cells Nk '//trim(str(.true.,Nk))
  if (allocated(NFi)) deallocate(NFi) ; allocate(NFi(0:Ni,1:Nj,1:Nk)) ; NFi(0:Ni,1:Nj,1:Nk) = block%mesh%NFi(0:Ni,1:Nj,1:Nk)
  if (allocated(NFj)) deallocate(NFj) ; allocate(NFj(1:Ni,0:Nj,1:Nk)) ; NFj(1:Ni,0:Nj,1:Nk) = block%mesh%NFj(1:Ni,0:Nj,1:Nk)
  if (allocated(NFk)) deallocate(NFk) ; allocate(NFk(1:Ni,1:Nj,0:Nk)) ; NFk(1:Ni,1:Nj,0:Nk) = block%mesh%NFk(1:Ni,1:Nj,0:Nk)
  if (allocated(Si )) deallocate(Si ) ; allocate(Si (0:Ni,1:Nj,1:Nk)) ; Si (0:Ni,1:Nj,1:Nk) = block%mesh%Si (0:Ni,1:Nj,1:Nk)
  if (allocated(Sj )) deallocate(Sj ) ; allocate(Sj (1:Ni,0:Nj,1:Nk)) ; Sj (1:Ni,0:Nj,1:Nk) = block%mesh%Sj (1:Ni,0:Nj,1:Nk)
  if (allocated(Sk )) deallocate(Sk ) ; allocate(Sk (1:Ni,1:Nj,0:Nk)) ; Sk (1:Ni,1:Nj,0:Nk) = block%mesh%Sk (1:Ni,1:Nj,0:Nk)
  if (allocated(V  )) deallocate(V  ) ; allocate(V  (1:Ni,1:Nj,1:Nk)) ; V  (1:Ni,1:Nj,1:Nk) = block%mesh%V  (1:Ni,1:Nj,1:Nk)
  xmm=trim(str(n=maxval(NFi(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFi(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFi(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFi-x-max '//trim(xmm)//' NFi-y-max '//trim(ymm)//' NFi-z-max '//trim(zmm)
  xmm=trim(str(n=minval(NFi(:,:,:)%x))) ; ymm=trim(str(n=minval(NFi(:,:,:)%y))) ; zmm=trim(str(n=minval(NFi(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFi-x-min '//trim(xmm)//' NFi-y-min '//trim(ymm)//' NFi-z-min '//trim(zmm)
  xmm=trim(str(n=maxval(NFj(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFj(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFj(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFj-x-max '//trim(xmm)//' NFj-y-max '//trim(ymm)//' NFj-z-max '//trim(zmm)
  xmm=trim(str(n=minval(NFj(:,:,:)%x))) ; ymm=trim(str(n=minval(NFj(:,:,:)%y))) ; zmm=trim(str(n=minval(NFj(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFj-x-min '//trim(xmm)//' NFj-y-min '//trim(ymm)//' NFj-z-min '//trim(zmm)
  xmm=trim(str(n=maxval(NFk(:,:,:)%x))) ; ymm=trim(str(n=maxval(NFk(:,:,:)%y))) ; zmm=trim(str(n=maxval(NFk(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFk-x-max '//trim(xmm)//' NFk-y-max '//trim(ymm)//' NFk-z-max '//trim(zmm)
  xmm=trim(str(n=minval(NFk(:,:,:)%x))) ; ymm=trim(str(n=minval(NFk(:,:,:)%y))) ; zmm=trim(str(n=minval(NFk(:,:,:)%z)))
  write(stdout,'(A)',iostat=err)trim(rks)//'  NFk-x-min '//trim(xmm)//' NFk-y-min '//trim(ymm)//' NFk-z-min '//trim(zmm)
  vmax=trim(str(n=maxval(Si(:,:,:)))) ; vmin=trim(str(n=minval(Si(:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Si-max '//trim(vmax)//' Si-min '//trim(vmin)
  vmax=trim(str(n=maxval(Sj(:,:,:)))) ; vmin=trim(str(n=minval(Sj(:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Sj-max '//trim(vmax)//' Sj-min '//trim(vmin)
  vmax=trim(str(n=maxval(Sk(:,:,:)))) ; vmin=trim(str(n=minval(Sk(:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Sk-max '//trim(vmax)//' Sk-min '//trim(vmin)
  vmax=trim(str(n=maxval(V(:,:,:)))) ; vmin=trim(str(n=minval(V(:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  V-max '//trim(vmax)//' V-min '//trim(vmin)
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction print_info_bmesh

  !> Function for loading inflow 1 boundary conditions.
  !> @return \b err integer(I4P) variable.
  function load_gbc_in1(global,myrank,filename,in1) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global   !< Global level data.
  integer(I_P),       intent(IN)::    myrank   !< Actual rank process.
  character(*),       intent(IN)::    filename !< Name of file where boundary conditions variables are saved.
  integer(I_P),       intent(IN)::    in1      !< Actual in1 index.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of boundary conditions file.
  character(DI_P)::                   rks      !< String containing myrank.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_gbc_in1')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  err = read_primitive(scalar=global%bc%in1(in1),format='*',unit=UnitFree)
  close(UnitFree)
  rks = 'rank'//trim(str(.true.,myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Inflow 1 primitive variables loaded:'
  err = global%bc%in1(in1)%pprint(myrank=myrank,unit=stdout)
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,*,iostat=err)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gbc_in1

  !> Function for loading block boundary conditions file.
  !> @return \b err integer(I4P) variable.
  function load_bbc(block,myrank,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(INOUT):: block    !< Block level data.
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where mesh variables are saved.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of boundary conditions file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_bbc')
  open(unit=Get_Unit(UnitFree), file=adjustl(trim(filename)), status='OLD', action='READ', form='UNFORMATTED')
  err = read_bc(array3D=block%bc%BCi,unit=UnitFree)
  err = read_bc(array3D=block%bc%BCj,unit=UnitFree)
  err = read_bc(array3D=block%bc%BCk,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bbc

  !> Function for saving block boundary conditions file.
  !> @return \b err integer(I4P) variable.
  function save_bbc(block,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(IN):: block    !< Block level data.
  character(*),      intent(IN):: filename !< Name of file where mesh variables are saved.
  integer(I_P)::                  err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                  UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  open(unit=Get_Unit(UnitFree), file=adjustl(trim(filename)), action='WRITE', form='UNFORMATTED')
  err = write_bc(array3D=block%bc%BCi,unit=UnitFree)
  err = write_bc(array3D=block%bc%BCj,unit=UnitFree)
  err = write_bc(array3D=block%bc%BCk,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bbc

  !> Function for loading solver option variables.
  !> @return \b err integer(I4P) variable.
  function load_gfluid_soption(global,myrank,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global   !< Global level data.
  integer(I_P),       intent(IN)::    myrank   !< Rank process identification for MPI communications.
  character(*),       intent(IN)::    filename !< Name of file where option variables are saved.
  integer(I_P)::                      err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                      UnitFree !< Free logic unit.
  logical::                           is_file  !< Flag for inquiring the presence of option file.
  character(8)::                      timing   !< timing = 'unsteady' simulation otherwise steady one.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  err = 1
  inquire(file=adjustl(trim(filename)),exist=is_file)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_gfluid_soption')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err) ! record skipped because unnecessary
  read(UnitFree,*,iostat=err)timing ; timing = Upper_Case(timing) ; global%fluid%unsteady = (adjustl(trim(timing))=='UNSTEADY')
  read(UnitFree,*,iostat=err)global%fluid%Nmax
  read(UnitFree,*,iostat=err)global%fluid%Tmax
  read(UnitFree,*,iostat=err)global%fluid%sp_ord
  read(UnitFree,*,iostat=err)global%fluid%rk_ord
  read(UnitFree,*,iostat=err)global%fluid%CFL
  read(UnitFree,*,iostat=err)global%fluid%residual_toll
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
  if (global%fluid%unsteady) then
    if (global%fluid%Nmax>0_I8P) then
      global%fluid%Tmax=-1._R_P ! the value of Tmax is ignored because Nmax>0
    endif
  else
    global%fluid%Tmax=-1._R_P ! the value of Tmax is ignored because steady simulation
  endif
  call compute_values0(global) ! computing the reference values for non-dimensional quantities
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_soption

  !> Function for loading the number of species from fluid dynamic file "filename".
  !> @return \b err integer(I4P) variable.
  function load_gfluid_Ns(global,binary,myrank,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT)::        global   !< Global level data.
  logical,            intent(IN), optional:: binary   !< Flag for binary of ascii input file.
  integer(I_P),       intent(IN)::           myrank   !< Actual rank process.
  character(*),       intent(IN)::           filename !< Name of file where mesh variables are saved.
  integer(I_P)::                             err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                             UnitFree !< Free logic unit.
  logical::                                  is_file  !< Flag for inquiring the presence of fluid dynamic file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=trim(filename),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_gfluid_Ns')
  if (present(binary)) then
    open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
    read(UnitFree,iostat=err)global%fluid%Ns
  else
    open(unit = Get_Unit(UnitFree), file = trim(filename), status = 'OLD', action = 'READ', form = 'FORMATTED')
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err)global%fluid%Ns
  endif
  close(UnitFree)
  global%fluid%Np = global%fluid%Ns + 6
  global%fluid%Nc = global%fluid%Ns + 4
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_Ns

  !> Function for loading initial species from ascii file.
  !> @return \b err integer(I4P) variable.
  function load_gfluid_0species(global,myrank,filename) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Global), intent(INOUT):: global   !< Global level data.
  integer(I_P),       intent(IN)::    myrank   !< Actual rank process.
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
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_gfluid_0species')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'FORMATTED')
  read(UnitFree,*,iostat=err) ! record skipped
  read(UnitFree,*,iostat=err) ! record skipped
  read(UnitFree,*,iostat=err)global%fluid%Ns ; call global%alloc_fluid
  do s=1,global%fluid%Ns
    read(UnitFree,*,iostat=err) ! record skipped
    read(UnitFree,*,iostat=err)g
    read(UnitFree,*,iostat=err)R
    read(UnitFree,*,iostat=err)global%fluid%cp0(s)
    read(UnitFree,*,iostat=err)global%fluid%cv0(s)
    if ((g>0.0_R_P).AND.(R>0.0_R_P)) then
      global%fluid%cv0(s)=R/(g-1._R_P)
      global%fluid%cp0(s)=g*global%fluid%cv0(s)
    elseif ((g>0.0_R_P).AND.(global%fluid%cp0(s)>0.0_R_P)) then
      global%fluid%cv0(s)=global%fluid%cp0(s)/g
    elseif ((g>0.0_R_P).AND.(global%fluid%cv0(s)>0.0_R_P)) then
      global%fluid%cp0(s)=global%fluid%cv0(s)*g
    elseif ((R>0.0_R_P).AND.(global%fluid%cp0(s)>0.0_R_P)) then
      global%fluid%cv0(s)=global%fluid%cp0(s)-R
    elseif ((R>0.0_R_P).AND.(global%fluid%cv0(s)>0.0_R_P)) then
      global%fluid%cp0(s)=global%fluid%cv0(s)+R
    endif
  enddo
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_gfluid_0species

  !> Function for loading fluid dynamic file.
  !> @return \b err integer(I4P) variable.
  function load_bfluid(block,myrank,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(INOUT):: block    !< Block level data.
  integer(I_P),      intent(IN)::    myrank   !< Actual rank process.
  character(*),      intent(IN)::    filename !< Name of file where mesh variables are saved.
  type(Type_Global), intent(INOUT):: global   !< Global level data.
  integer(I_P)::                     err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                     UnitFree !< Free logic unit.
  logical::                          is_file  !< Flag for inquiring the presence of fluid dynamic file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(filename)),exist=is_file,iostat=err)
  if (.NOT.is_file) call File_Not_Found(myrank,filename,'load_bfluid')
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), status = 'OLD', action = 'READ', form = 'UNFORMATTED')
  read(UnitFree,iostat=err)global%fluid%Ns
  read(UnitFree,iostat=err)global%fluid%cp0,global%fluid%cv0
  read(UnitFree,iostat=err)global%fluid%n
  read(UnitFree,iostat=err)global%fluid%t
  read(UnitFree,iostat=err)block%fluid%Dt
  err = read_primitive(array3D=block%fluid%P,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction load_bfluid

  !> Function for saving fluid dynamic file.
  !> @return \b err integer(I4P) variable.
  function save_bfluid(block,filename,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(IN):: block    !< Block level data.
  character(*),      intent(IN):: filename !< Name of file where mesh variables are saved.
  type(Type_Global), intent(IN):: global   !< Global level data.
  integer(I_P)::                  err      !< Error trapping flag: 0 no errors, >0 error occurs.
  integer(I_P)::                  UnitFree !< Free logic unit.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  open(unit = Get_Unit(UnitFree), file = adjustl(trim(filename)), action = 'WRITE', form = 'UNFORMATTED')
  write(UnitFree,iostat=err)global%fluid%Ns
  write(UnitFree,iostat=err)global%fluid%cp0,global%fluid%cv0
  write(UnitFree,iostat=err)global%fluid%n
  write(UnitFree,iostat=err)global%fluid%t
  write(UnitFree,iostat=err)block%fluid%Dt
  err = write_primitive(array3D=block%fluid%P,unit=UnitFree)
  close(UnitFree)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_bfluid

  !> Function for printing to standard output info of block fluid data.
  !> @return \b err integer(I4P) variable.
  function print_info_bfluid(block,myrank,blk,grl,global) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Block), intent(IN):: block          !< Block level data.
  integer(I_P),      intent(IN):: myrank         !< Rank process identification for MPI communications.
  integer(I_P),      intent(IN):: blk            !< Actual block number.
  integer(I_P),      intent(IN):: grl            !< Actual grid level number.
  type(Type_Global), intent(IN):: global         !< Global level data.
  integer(I_P)::                  err            !< Error trapping flag: 0 no errors, >0 error occurs.
  character(DR_P)::               vmax,vmin      !< String for printing max and min of variables.
  real(R_P), allocatable::        dummy(:,:,:,:) !<  Dummy variables for printing only internal cells info.
  integer(I_P)::                  Ni,Nj,Nk,Np,Nc !< Temporary variables for storing block dimensions.
  character(DI_P)::               rks            !< String containing myrank.
  integer(I_P)::                  i,j,k,s        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rks = 'rank'//trim(str(.true.,myrank))
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,'(A)',iostat=err)trim(rks)//' Fluid infos of myrank '//trim(rks)
  write(stdout,'(A)',iostat=err)trim(rks)//'  Global number of blocks '//trim(str(.true.,global%mesh%Nb_tot))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Number of initial species '//trim(str(.true.,global%fluid%Ns))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Block '//trim(str(.true.,blk))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Grid level '//trim(str(.true.,grl))
  Ni = block%mesh%Ni
  Nj = block%mesh%Nj
  Nk = block%mesh%Nk
  Np = global%fluid%Np
  Nc = global%fluid%Nc
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Nc,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(:,i,j,k) = block%fluid%U(i,j,k)%cons2array()
      enddo
    enddo
  enddo
  write(stdout,'(A)',iostat=err)trim(rks)//'  Conservative variables'
  do s=1,global%fluid%Nc
    vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
    write(stdout,'(A)',iostat=err)trim(rks)//'  U('//trim(str(.true.,s))//'): max '//trim(vmax)//' min '//trim(vmin)
  enddo
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Np,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(:,i,j,k) = block%fluid%P(i,j,k)%prim2array()
      enddo
    enddo
  enddo
  write(stdout,'(A)',iostat=err)trim(rks)//'  Primitive variables'
  do s=1,global%fluid%Np
    vmax=trim(str(n=maxval(dummy(s,:,:,:)))) ; vmin=trim(str(n=minval(dummy(s,:,:,:))))
    write(stdout,'(A)',iostat=err)trim(rks)//'  P('//trim(str(.true.,s))//'): max '//trim(vmax)//' min '//trim(vmin)
  enddo
  if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1,1:Ni,1:Nj,1:Nk)) ; dummy = 0._R_P
  do k=1,Nk
    do j=1,Nj
      do i=1,Ni
        dummy(1,i,j,k) = block%fluid%Dt(i,j,k)
      enddo
    enddo
  enddo
  vmax=trim(str(n=maxval(dummy(1,:,:,:)))) ; vmin=trim(str(n=minval(dummy(1,:,:,:))))
  write(stdout,'(A)',iostat=err)trim(rks)//'  Time steps Dt: max '//trim(vmax)//' min '//trim(vmin)
  write(stdout,'(A)',iostat=err)trim(rks)//'----------------------------------------------------------------------'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction print_info_bfluid
  !> @}
endmodule Data_Type_Globals
