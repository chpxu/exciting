! Copyright (C) 2015-2023 exciting team (Riga and Berlin)
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
subroutine GSortho(ld,m,n,x)
!> The subroutine orthogonalises the vectors provided in the input.  
!> Note, it is orthogonalisation with the inner product x^T * y and not x^T * S * y.
!> GSortho assumes that the first m (number provided in the input) vectors are already orthogonal.
      Use constants, Only : zzero,zone
      Use modmpi, Only: rank
      Use mod_gkvector
      implicit none
      integer, intent(in) :: ld             ! Leading dimension, size of the vectors
      integer, intent(in) :: m              ! The number of vectors that are already orthogonalised
      integer, intent(in) :: n              ! The total number of vectors
      complex(8), intent(inout) :: x(ld,n)  ! Matrix of vectors

! local variables
      complex(8) :: zdotc
      external :: zdotc
      integer :: i
      real(8) :: ta,tb
      complex(8), allocatable :: A(:,:), update(:,:),normlist(:),zbuf(:)
      integer :: lwork,lrwork,liwork,info
      complex(8), allocatable :: work(:)
      real(8), allocatable :: w(:),rwork(:)
      integer, allocatable :: iwork(:)

      call timesec(ta)
! The first m vectors are already orthogonal among themselves.
! Orthogonalize the n-m successors to the first m vectors.
      if (m.gt.0) then
        Allocate(A(m,n-m))
! inner products
        call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                   m, &          ! M ... rows of op( A ) = rows of C
                   n-m, &           ! N ... cols of op( B ) = cols of C
                   ld, &          ! K ... cols of op( A ) = rows of op( B )
                   zone, &          ! alpha
                   x(1,1), &           ! A
                   ld,&           ! LDA ... leading dimension of A
                   x(1,m+1), &           ! B
                   ld, &          ! LDB ... leading dimension of B
                   zzero, &          ! beta
                   A, &  ! C
                   m &      ! LDC ... leading dimension of C
                  )
        allocate(update(ld,n-m))
! calculate the projections
        call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                   ld, &          ! M ... rows of op( A ) = rows of C
                   n-m, &           ! N ... cols of op( B ) = cols of C
                   m, &          ! K ... cols of op( A ) = rows of op( B )
                   zone, &          ! alpha
                   x(1,1), & ! A
                   ld,           & ! LDA ... leading dimension of A
                   A, &           ! B
                   m,&           ! LDB ... leading dimension of B
                   zzero, &          ! beta
                   update, &  ! C
                   ld &      ! LDC ... leading dimension of C
                  )
! subtract the projections
        x(1:ld,m+1:n)=x(1:ld,m+1:n)-update(1:ld,1:n-m)

! It should be voila at this point, but it's numerics. Nothing is exact...
! Test the orthogonality, if you need to debug or just to verify
#ifdef TEST_GSORTHO
        deallocate(update)
        Allocate(update(m,n-m))

        call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                   m, &          ! M ... rows of op( A ) = rows of C
                   n-m, &           ! N ... cols of op( B ) = cols of C
                   ld, &          ! K ... cols of op( A ) = rows of op( B )
                   zone, &          ! alpha
                   x(1,1), &           ! A
                   ld,&           ! LDA ... leading dimension of A
                   x(1,m+1), &           ! B
                   ld, &          ! LDB ... leading dimension of B
                   zzero, &          ! beta
                   A, &  ! C
                   m &      ! LDC ... leading dimension of C
                  )
        A=abs(A)
        write(*,*) 'GSortho test1',sum(A)
#endif 

        deallocate(A,update)
      endif

! Now, we orthogonalize the (n-m) successors
      allocate(normlist(n-m))
      allocate(zbuf(n-m))

! Normalize the vectors first
      do i=1,n-m
        normlist(i)=zdotc(ld,x(1,m+i),1,x(1,m+i),1)
      enddo
      do i=1,n-m
        x(1:ld,m+i)=x(1:ld,m+i)/sqrt(dble(normlist(i)))
      enddo
      deallocate(normlist,zbuf)

! Calculate inner products (overlaps)
      Allocate(A(n-m,n-m))
      Allocate(update(n-m,n-m))
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 n-m, &          ! M ... rows of op( A ) = rows of C
                 n-m, &           ! N ... cols of op( B ) = cols of C
                 ld, &          ! K ... cols of op( A ) = rows of op( B )
                 zone, &          ! alpha
                 x(1,m+1), &           ! A
                 ld,&           ! LDA ... leading dimension of A
                 x(1,m+1), &           ! B
                 ld, &          ! LDB ... leading dimension of B
                 zzero, &          ! beta
                 A, &  ! C
                 n-m &      ! LDC ... leading dimension of C
                )
! Diagonalise the overlap matrix. The eigenvectors will be orthogonal.
      allocate(w(n-m))
! This is a dry run to determine the necessary array dimensions
      allocate(work(1))
      allocate(rwork(1))
      allocate(iwork(1))
      lwork=-1
      lrwork=-1
      liwork=-1
      call zheevd ('V', 'U', n-m, A, n-m, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
      lwork=int(work(1))
      lrwork=int(rwork(1))
      liwork=iwork(1)
      deallocate(work,rwork,iwork)
! Now the proper diagonalisation
      allocate(work(lwork))
      allocate(rwork(lrwork))
      allocate(iwork(liwork))
      call zheevd ('V', 'U', n-m, A, n-m, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
      deallocate(work,rwork,iwork)
    
      deallocate(update)
      allocate(update(ld,n-m))
! Transform the vectors
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 ld, &          ! M ... rows of op( A ) = rows of C
                 n-m, &           ! N ... cols of op( B ) = cols of C
                 n-m, &          ! K ... cols of op( A ) = rows of op( B )
                 zone, &          ! alpha
                 x(1,m+1), & ! A
                 ld,           & ! LDA ... leading dimension of A
                 A, &           ! B
                 n-m,&           ! LDB ... leading dimension of B
                 zzero, &          ! beta
                 update, &  ! C
                 ld &      ! LDC ... leading dimension of C
                )
#ifdef aggressive
      do i=1,n-m
        x(1:ld,m+i)=update(1:ld,i)/sqrt(w(i))
      enddo
#else
! Normalise the vectors
      x(1:ld,m+1:n)=update(1:ld,1:n-m)
      allocate(normlist(n-m))
      allocate(zbuf(n-m))

      do i=1,n-m
        normlist(i)=zdotc(ld,x(1,m+i),1,x(1,m+i),1)
      enddo
      do i=1,n-m
        x(1:ld,m+i)=x(1:ld,m+i)/sqrt(dble(normlist(i)))
      enddo
      deallocate(normlist,zbuf)
#endif

      deallocate(A,update)

! Again, test the orthogonality if necessary.
#ifdef TEST_GSORTHO
      Allocate(update(n,n))
      Allocate(A(n,n))
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  n, &          ! M ... rows of op( A ) = rows of C
                  n, &           ! N ... cols of op( B ) = cols of C
                  ld, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  x(1,1), &           ! A
                  ld,&           ! LDA ... leading dimension of A
                  x(1,1), &           ! B
                  ld, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  A, &  ! C
                  n &      ! LDC ... leading dimension of C
                 )
      A=abs(A)
      write(*,*) 'GSortho test2',n,sum(A)
      deallocate(A,update,w)
#endif

      call timesec(tb)
#ifdef TIMINGS
      if (rank.eq.0) write(*,*) 'GSortho', tb-ta
#endif
end subroutine GSortho



