!> Module for mpi_allgather wrappers.
module mod_mpi_allgather
  use mod_mpi_env, only: mpiinfo
  use precision, only: sp, dp, i32

#ifdef MPI
  use mpi
#endif   

  implicit none 
  private

  
  public :: xmpi_allgather, xmpi_allgatherv


  !> Wrappers for mpi_allgether.
  interface xmpi_allgather
    module procedure :: &
            mpi_allgather_integer_i32_rank_0_inplace
  end interface


  !> Wrappers for mpi_allgetherv.
  interface xmpi_allgatherv
    module procedure :: &
            mpi_allgatherv_in_place_rank3_complex_dp
  end interface

contains


  !> Wrapper for mpi_allgather for an `integer(i32)` scalar.
  !>
  !> Gather data from all tasks and send combined data to all tasks.
  subroutine mpi_allgather_integer_i32_rank_0_inplace(mpi_env, send_buffer, receive_buffer)
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Integer to be send by the current rank.
    integer(i32), intent(in) :: send_buffer
    !> Array that contains the the integers, recieved by all ranks
    integer(i32), allocatable, intent(out) :: receive_buffer(:)
#ifdef MPI    
    integer(i32) :: mpi_err
    allocate(receive_buffer(mpi_env%procs))
    receive_buffer(mpi_env%rank + 1) = send_buffer
    call mpi_allgather(send_buffer, 1, MPI_INTEGER, &
                       receive_buffer, 1, MPI_INTEGER, &
                       mpi_env%comm, mpi_err)
#endif
  end subroutine mpi_allgather_integer_i32_rank_0_inplace


  !> Wrapper for mpi_allgatherv for a 3-rank `complex(dp)` array.
  !>
  !> Gather data from all tasks and send combined data to all tasks.
  subroutine mpi_allgatherv_in_place_rank3_complex_dp(mpi_env, buffer, offset, chunk_shape)
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Buffer. On input, contains the data of the current rank in the range `offset` till `offset + chunk_shape - 1`.
    !> On output, contains the data of all ranks.
    complex(dp), intent(inout) :: buffer(:, :, :)
    !> Offset or start adress of the data chunk handled by the current rank with reference to `[[buffer]]`.
    integer(i32), intent(in) :: offset(3)
    !> Shape of the data chunk, handled by the current rank.
    integer(i32), intent(in) :: chunk_shape(3)
#ifdef MPI
    integer(i32) :: mpi_error
    integer(i32) :: send_count
    integer(i32), allocatable :: receive_counts(:), displacements(:)

    send_count = product(chunk_shape)

    call xmpi_allgather(mpi_env, send_count, receive_counts)
    call calculate_displacements(mpi_env, receive_counts, displacements)
    
    call mpi_allgatherv( &
            MPI_IN_PLACE, &
            0, &
            MPI_DATATYPE_NULL, &
            buffer, &
            receive_counts, &
            displacements, &
            MPI_DOUBLE_COMPLEX, &
            mpi_env%comm, &
            mpi_error &
          )
#endif    
  end subroutine mpi_allgatherv_in_place_rank3_complex_dp


  !> Calculate the number of elements, the data chunk of each rank is displaced by.
  subroutine calculate_displacements(mpi_env, receive_counts, displacements)
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env
    !> Array that holds the size of the data chunk for each rank.
    integer(i32), intent(in) :: receive_counts(:)
    !> Number of elements, the data chunk of each rank is displaced by.
    integer(i32), allocatable, intent(out) :: displacements(:)

    integer(i32) :: rank

    allocate(displacements(mpi_env%procs))

    do rank=0, mpi_env%procs-1
      displacements(rank + 1) = sum(receive_counts(1 : rank))
    end do
  end subroutine 
  

end module mod_mpi_allgather