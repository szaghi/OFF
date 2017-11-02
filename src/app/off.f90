!< OFF, Open source Finite volumes Fluid dynamics code.

program off
!< OFF, Open source Finite volumes Fluid dynamics code.
use off_objects, only : simulation_object

implicit none
type(simulation_object) :: simulations !< Simulation data.

call simulations%initialize()

endprogram off
