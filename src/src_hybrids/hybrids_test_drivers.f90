!> Module for collecting unit test drivers for modules in the hybrids directory.
module hybrids_test_drivers
  use modmpi, only: mpiinfo
  ! Load test drivers here
  use hse_singularity_test, only: hse_singularity_test_driver
  private
  public :: hybrids_test_driver

  contains

  subroutine hybrids_test_driver(mpiglobal, kill_on_failure)
    !> mpi information
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program before the test driver finishes
    !> if an assertion fails 
    logical, optional :: kill_on_failure 

    ! Call test drivers here
    call hse_singularity_test_driver(mpiglobal, kill_on_failure)
  end subroutine hybrids_test_driver

end module hybrids_test_drivers
