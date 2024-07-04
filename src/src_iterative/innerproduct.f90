! Copyright (C) 2015-2023 exciting team (Berlin and Riga)
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> The subroutine calculates the matrix product C = A**H * B.  
!> 
subroutine innerproduct(l,m,n,A,B,C)
      Use constants, Only : zzero, zone
      use modmpi
      use mod_gkvector
      implicit none
      integer, intent(in) :: l            ! leading dimension of matrices A and B
      integer, intent(in) :: m            ! 
      integer, intent(in) :: n            ! 
      complex(8), intent(in) :: A(l,m)    ! matrix A
      complex(8), intent(in) :: B(l,n)    ! matrix B
      complex(8), intent(out) :: C(m,n)

      real(8) :: ta,tb


      call timesec(ta)
      C=zzero
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  m, &          ! M ... rows of op( A ) = rows of C
                  n, &           ! N ... cols of op( B ) = cols of C
                  l, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  A(1,1), &           ! A
                  l,&           ! LDA ... leading dimension of A
                  B(1,1), &           ! B
                  l, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  C, &  ! C
                  m &      ! LDC ... leading dimension of C
                  )

      call timesec(tb)
#ifdef TIMINGS
      write(*,*) 'innerproduct', tb-ta
#endif

end subroutine innerproduct



