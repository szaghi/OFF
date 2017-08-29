!< OFF objects collection.

module off_objects
!< OFF objects collection.
!<
!< This is convenient entry point to access all OFF objects.

use off_block_object, only : block_object
use off_block_signature_object, only : block_signature_object
use off_cell_object, only : cell_object
use off_error_object, only : error_object
use off_face_object, only : face_object
use off_file_object, only : file_object
use off_file_grid_object, only : file_grid_object
use off_files_collection_object, only : files_collection_object
use off_free_conditions_object, only : free_conditions_object
use off_grid_dimensions_object, only : grid_dimensions_object
use off_node_object, only : node_object
use off_non_dimensional_numbers_object, only : non_dimensional_numbers_object
use off_os_object, only : os_object
use off_simulation_object, only : simulation_object
use off_solver_object, only : solver_object
use off_time_object, only : time_object

implicit none
private
public :: block_object
public :: block_signature_object
public :: cell_object
public :: error_object
public :: face_object
public :: file_object
public :: file_grid_object
public :: files_collection_object
public :: free_conditions_object
public :: grid_dimensions_object
public :: node_object
public :: non_dimensional_numbers_object
public :: os_object
public :: simulation_object
public :: solver_object
public :: time_object
endmodule off_objects
