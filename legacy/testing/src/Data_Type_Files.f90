!> @ingroup DerivedType
!> @{
!> @defgroup Data_Type_FilesDerivedType Data_Type_Files
!> @}

!> @brief Module Data_Type_Files aggregate the definition of all files structures of OFF project.
!> @ingroup Data_Type_FilesDerivedType
module Data_Type_Files
!-----------------------------------------------------------------------------------------------------------------------------------
USE Data_Type_File_BC,               only: Type_File_BC               ! Definition of Type_File_BC.
USE Data_Type_File_Fluid,            only: Type_File_Fluid            ! Definition of Type_File_Fluid.
USE Data_Type_File_Blocks_Cartesian, only: Type_File_Blocks_Cartesian ! Definition of Type_File_Blocks_Cartesian.
USE Data_Type_File_GNU,              only: Type_File_GNU              ! Definition of Type_File_GNU.
USE Data_Type_File_IBM_Options,      only: Type_File_IBM_Options      ! Definition of Type_File_IBM_Options.
USE Data_Type_File_Lock,             only: Type_File_Lock             ! Definition of Type_File_Lock.
USE Data_Type_File_Mesh,             only: Type_File_Mesh             ! Definition of Type_File_Mesh.
USE Data_Type_File_OFF_Options,      only: Type_File_OFF_Options      ! Definition of Type_File_OFF_Options.
!USE Data_Type_File_Procmap,          only: Type_File_Procmap          ! Definition of Type_File_Procmap.
USE Data_Type_File_Profile,          only: Type_File_Profile          ! Definition of Type_File_Profile.
USE Data_Type_File_Species,          only: Type_File_Species          ! Definition of Type_File_Species.
USE Data_Type_File_Solver_Options,   only: Type_File_Solver_Options   ! Definition of Type_File_Solver_Options.
USE Data_Type_File_Tec,              only: Type_File_Tec              ! Definition of Type_File_Tec.
USE Data_Type_File_VTK,              only: Type_File_VTK              ! Definition of Type_File_VTK.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Derived type containing the all files structures of OFF project.
!> @ingroup Data_Type_FilesDerivedType
type, public:: Type_Files
  type(Type_File_IBM_Options)::      ibm_opts  !< IBM options file.
  type(Type_File_OFF_Options)::      off_opts  !< OFF options file.
  type(Type_File_Solver_Options)::   solv_opts !< Solver options file.
  type(Type_File_Blocks_Cartesian):: cart_blks !< IBM blocks description file.
  type(Type_File_Species)::          spec      !< Species file.
  !type(Type_File_Procmap)::          proc      !< Processes/blocks map file.
  type(Type_File_Mesh)::             mesh      !< Mesh file.
  type(Type_File_BC)::               bc        !< Boundary conditions file.
  type(Type_File_Fluid)::            init      !< Initial conditions file.
  type(Type_File_Fluid)::            sol       !< Solution file.
  type(Type_File_Profile)::          prof      !< Profiling file.
  type(Type_File_Lock)::             lockfile  !< Lockfile.
  type(Type_File_Tec)::              tec       !< Tecplot post-precessing file.
  type(Type_File_VTK)::              vtk       !< VTK post-precessing file.
  type(Type_File_GNU)::              gnu       !< Gnuplot post-precessing file.
endtype Type_Files
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule Data_Type_Files
