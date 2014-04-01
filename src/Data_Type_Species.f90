!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_SpeciesDerivedType Data_Type_Species
!> Module definition of fluid specie, Type_Species
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Data_Type_SpeciesInterface Data_Type_Species
!> Module definition of Type_Species
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_SpeciesPrivateProcedure Data_Type_Species
!> Module definition of fluid specie, Type_Species
!> @}

!> @brief Module Data_Type_Species contains the definition of Type_Species, that defines the fluid specie properties.
module Data_Type_Species
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                        ! Integers and reals precision definition.
USE Data_Type_Specie, only: Type_Specie ! Definition of Type_Specie.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the definition of Type_Species.
!> @ingroup Data_Type_SpeciesDerivedType
type, public:: Type_Species
  integer(I4P)::                   Np = 7_I4P !< Number of primitive variables    (Np = Ns + 6).
  integer(I4P)::                   Nc = 5_I4P !< Number of conservative variables (Nc = Ns + 4).
  integer(I4P)::                   Ns = 1_I4P !< Number of species.
  type(Type_Specie), allocatable:: heats(:)   !< Specific heats for each specie [1:Ns].
  contains
    procedure:: free        => free_specie        ! Procedure for freeing dynamic memory.
    procedure:: alloc       => alloc_specie       ! Procedure for allocating dynamic memory.
    procedure:: compute_Npc => compute_Np_Nc      ! Procedure for computing number of primitive and conservative variables.
    procedure:: print       => print_species_self ! Procedure for printing species with a pretty format.
    final::     finalize                          ! Procedure for freeing dynamic memory when finalizing.
    ! operators overloading
    generic:: assignment(=) => assign_species
    ! private procedures
    procedure, pass(spec1), private:: assign_species
endtype Type_Species
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_SpeciesPrivateProcedure
  !> @{
  !> Subroutine for freeing dynamic data of Type_Species variables.
  elemental subroutine free_specie(species)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Species), intent(INOUT):: species !< Species data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(species%heats)) deallocate(species%heats)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine free_specie

  !> @brief Procedure for freeing dynamic memory when finalizing.
  elemental subroutine finalize(species)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_Species), intent(INOUT):: species !< Species data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call species%free
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize

  !> Subroutine for allocating dynamic data of Type_Species variables.
  elemental subroutine alloc_specie(species)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Species), intent(INOUT):: species !< Species data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call species%free ; allocate(species%heats(1:species%Ns))
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alloc_specie

  !> @brief Procedure for computing number of primitive and conservative variables.
  elemental subroutine compute_Np_Nc(species)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Species), intent(INOUT):: species !< Species data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  species%Np = species%Ns + 6_I4P ; species%Nc = species%Ns + 4_I4P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_Np_Nc

  !> @brief Procedure for printing species with a pretty format.
  subroutine print_species_self(spec,pref,iostat,iomsg,unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Species),    intent(IN)::  spec    !< Specie.
  character(*), optional, intent(IN)::  pref    !< Prefixing string.
  integer(I4P), optional, intent(OUT):: iostat  !< IO error.
  character(*), optional, intent(OUT):: iomsg   !< IO error message.
  integer(I4P),           intent(IN)::  unit    !< Logic unit.
  character(len=:), allocatable::       prefd   !< Prefixing string.
  integer(I4P)::                        iostatd !< IO error.
  character(500)::                      iomsgd  !< Temporary variable for IO error message.
  integer(I4P)::                        s       !< Species counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  prefd = '' ; if (present(pref)) prefd = pref
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Np='//str(n=spec%Np)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Nc='//str(n=spec%Nc)
  write(unit=unit,fmt='(A)',iostat=iostatd,iomsg=iomsgd)prefd//' Ns='//str(n=spec%Ns)
  do s=1,spec%Ns
    call spec%heats(s)%print(unit=unit,iostat=iostatd,iomsg=iomsgd,pref=prefd//'  s('//trim(str(.true.,s))//'):')
  enddo
  if (present(iostat)) iostat = iostatd
  if (present(iomsg))  iomsg  = iomsgd
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_species_self

  ! Assignment (=)
  !> @brief Procedure for assignment between two species variables.
  elemental subroutine assign_species(spec1,spec2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_Species), intent(INOUT):: spec1
  type(Type_Species),  intent(IN)::    spec2
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(spec1%heats)) deallocate(spec1%heats) ; allocate(spec1%heats(1:spec2%Ns))
  spec1%Np    = spec2%Np
  spec1%Nc    = spec2%Nc
  spec1%Ns    = spec2%Ns
  spec1%heats = spec2%heats
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_species
  !> @}
endmodule Data_Type_Species
