!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_PostProcessDerivedType Data_Type_PostProcess
!> @}

!> @ingroup GlobalVarPar
!> @{
!> @defgroup Data_Type_PostProcessGlobalVarPar Data_Type_PostProcess
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Data_Type_PostProcessPublicProcedure Data_Type_PostProcess
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Data_Type_PostProcessPrivateProcedure Data_Type_PostProcess
!> @}

!> This module contains the definition of procedures and variables useful for post-process the OFF data.
!> @ingroup Data_Type_PostProcessDerivedType
module Data_Type_PostProcess
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                        ! Integers and reals precision definition.
USE Data_Type_SBlock, only: Type_SBlock ! Definition of Type_SBlock.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing the post-processing options.
!> @ingroup Data_Type_PostProcessDerivedType
type, public:: Type_PostProcess
  logical:: binary    = .true.  !< Binary or ascii post-process file.
  logical:: node      = .true.  !< Node or cell data location.
  logical:: bc        = .false. !< Saving or not boundary conditions cells.
  logical:: meshonly  = .true.  !< Flag for post-process only mesh.
  logical:: schlieren = .false. !< Saving or not (pseudo) Schlieren field.
  logical:: mirrorX   = .false. !< saving solution togheter a X-axis mirrored copy of flow fieds.
  logical:: mirrorY   = .false. !< saving solution togheter a Y-axis mirrored copy of flow fieds.
  logical:: mirrorZ   = .false. !< saving solution togheter a Z-axis mirrored copy of flow fieds.
  logical:: tec       = .true.  !< Tecplot file.
  logical:: vtk       = .false. !< VTK file.
  logical:: gnu       = .false. !< Gnuplot file.
  contains
    procedure:: compute_dimensions ! Procedure for computing the dimensions of the domain to be post-processed.
endtype Type_PostProcess
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_PostProcessPrivateProcedure
  !> @{
  !> @brief Procedure for computing the dimensions of the domain to be post-processed.
  subroutine compute_dimensions(pp,block,ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess), intent(IN)::  pp    !< Post-processing options.
  type(Type_SBlock),       intent(IN)::  block !< Block-level data.
  integer(I4P),            intent(OUT):: ni1   !< Minimum I index of node-centered data.
  integer(I4P),            intent(OUT):: ni2   !< Maximum I index of node-centered data.
  integer(I4P),            intent(OUT):: nj1   !< Minimum J index of node-centered data.
  integer(I4P),            intent(OUT):: nj2   !< Maximum J index of node-centered data.
  integer(I4P),            intent(OUT):: nk1   !< Minimum K index of node-centered data.
  integer(I4P),            intent(OUT):: nk2   !< Maximum K index of node-centered data.
  integer(I4P),            intent(OUT):: ci1   !< Minimum I index of cell-centered data.
  integer(I4P),            intent(OUT):: ci2   !< Maximum I index of cell-centered data.
  integer(I4P),            intent(OUT):: cj1   !< Minimum J index of cell-centered data.
  integer(I4P),            intent(OUT):: cj2   !< Maximum J index of cell-centered data.
  integer(I4P),            intent(OUT):: ck1   !< Minimum K index of cell-centered data.
  integer(I4P),            intent(OUT):: ck2   !< Maximum K index of cell-centered data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(gc=>block%dims%gc,Ni=>block%dims%Ni,Nj=>block%dims%Nj,Nk=>block%dims%Nk)
    if (pp%node) then
      ni1 = 0 ; ni2 = Ni
      nj1 = 0 ; nj2 = Nj
      nk1 = 0 ; nk2 = Nk

      ci1 = 0 ; ci2 = Ni
      cj1 = 0 ; cj2 = Nj
      ck1 = 0 ; ck2 = Nk
    else
      if (pp%bc) then
        ni1 = 0 - gc(1) ; ni2 = Ni + gc(2)
        nj1 = 0 - gc(3) ; nj2 = Nj + gc(4)
        nk1 = 0 - gc(5) ; nk2 = Nk + gc(6)

        ci1 = 1 - gc(1) ; ci2 = Ni + gc(2)
        cj1 = 1 - gc(3) ; cj2 = Nj + gc(4)
        ck1 = 1 - gc(5) ; ck2 = Nk + gc(6)
      else
        ni1 = 0         ; ni2 = Ni
        nj1 = 0         ; nj2 = Nj
        nk1 = 0         ; nk2 = Nk

        ci1 = 1         ; ci2 = Ni
        cj1 = 1         ; cj2 = Nj
        ck1 = 1         ; ck2 = Nk
      endif
    endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_dimensions
  !> @}
endmodule Data_Type_PostProcess
