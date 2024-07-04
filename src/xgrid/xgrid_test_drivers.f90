!> Collect all unit test drivers for structure modules

module xgrid_test_drivers
   use modmpi, only: mpiinfo

   ! Load test drivers here
   use regular_mesh_test, only: regular_mesh_test_driver
   use regular_grid_test, only: regular_grid_test_driver
   use regular_grid_constructors_test, only: regular_grid_constructors_test_driver

   private
   public :: xgrid_test_driver

contains

   subroutine xgrid_test_driver(mpiglobal, kill_on_failure)
      !> mpi information
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails
      logical, optional :: kill_on_failure

      ! Call test drivers here
      call regular_mesh_test_driver(mpiglobal, kill_on_failure)
      call regular_grid_test_driver(mpiglobal, kill_on_failure)
      call regular_grid_constructors_test_driver(mpiglobal, kill_on_failure)
   end subroutine xgrid_test_driver

end module xgrid_test_drivers

