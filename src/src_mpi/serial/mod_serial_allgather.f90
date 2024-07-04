!> Serial routine overloads for exciting's MPI wrappers.
!> These typically do nothing other than ensure the code
!> will compile in serial, without needing to dress all MPI 
!> calls throughout the code in preprocessor variables.
module mod_serial_allgather
  use mod_mpi_env, only: mpiinfo
  use precision, only: sp, dp, i32

#ifdef MPI
  use mpi
#endif   

  implicit none 
  private

  public :: xmpi_allgather, xmpi_allgatherv

  !> Wrappers for mpi_allgather.
  interface xmpi_allgather
    module procedure :: &
            mpi_allgather_integer_i32_rank_0_inplace
  end interface

  !> Wrappers for mpi_allgatherv.
  interface xmpi_allgatherv
    module procedure :: &
            mpi_allgatherv_in_place_rank3_complex_dp
  end interface

contains

  !> Dummy routine for serial version for mpi_gather
  subroutine mpi_allgather_integer_i32_rank_0_inplace(mpi_env, send_buffer, receive_buffer)
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Integer to be send by the current rank.
    integer(i32), intent(in) :: send_buffer
    !> Array that contains the the integers, recieved by all ranks
    integer(i32), allocatable, intent(out) :: receive_buffer(:)
  end subroutine mpi_allgather_integer_i32_rank_0_inplace

  !> Dummy routine for serial version for mpi_gatherv.
  subroutine mpi_allgatherv_in_place_rank3_complex_dp(mpi_env, buffer, offset, chunk_shape)
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer. On input, contains the data of the current rank in the range `offset` till `offset + chunk_shape - 1`.
    !> On output, contains the data of all ranks.
    complex(dp), intent(inout) :: buffer(:, :, :)
    !> Offset of the data of the current rank in `[[buffer]]`.
    integer(i32), intent(in) :: offset(3)
    !> Shape of the data chunk, handled by the current rank.
    integer(i32), intent(in) :: chunk_shape(3)
  end subroutine mpi_allgatherv_in_place_rank3_complex_dp

end module mod_serial_allgather