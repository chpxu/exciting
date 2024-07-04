!> Import this module for using xgrid codes.
module xgrid
  use regular_mesh, only: regular_mesh_type
  use regular_grid, only: regular_grid_type
  use regular_grid_constructors, only: setup_unitcell_grid, setup_fft_grid

  private
  public :: &
          regular_mesh_type, &
          regular_grid_type, &
          setup_unitcell_grid, &
          setup_fft_grid
end module xgrid